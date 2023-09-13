## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 05: Strategies.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame

# 1. Load data ==========
par <- import(here("data","parameters_manual.Rdata"))
WPP <- import(here("data","pop","WPP.Rdata"))
bWPP <- import(here("data","basepop.Rdata"))

# 2. Strategies ==========
# 2.1 Baseline parameters
parms_baseline <- par %>% 
  filter(!parameter %in% c('omega_ini', 'iota_cln_ini', 'phi_cln_ini', 
                           'rho_ini', 'rho_fin','sigma_ini','sigma_fin')) %>% 
  mutate(parameter = case_when(
    parameter == 'omega_fin' ~ 'omega', parameter == 'iota_cln_fin' ~ 'iota_cln', 
    parameter == 'phi_cln_fin' ~ 'phi_cln', TRUE ~ parameter)) %>% 
  add_row(parameter = "iota_min", low = 0, med = 0, hig = 0) %>% # Diagnosis Minimal
  add_row(parameter = "iota_sub", low = 0, med = 0, hig = 0) %>% # Diagnosis Subclinical
  add_row(parameter = "phi_min", low = 0, med = 0, hig = 0) %>%  # Treatment failure Minimal
  add_row(parameter = "phi_sub", low = 0, med = 0, hig = 0) %>%  # Treatment failure Subclinical
  select(parameter, med)
parms <- parms_baseline

# 2.2 Parameter set function
parms_set <- function(parms_base, parms_acf, acf_year, end_time) {
  sets <- list()
  
  for (year in 2020:end_time) {
    params <- parms_base
    
    if (year %in% acf_year) {
      for (i in seq_along(parms_acf$parameter)) {
        parm_name <- parms_acf$parameter[i]
        parm_value <- parms_acf$med[i]
        
        parms$med[parms$parameter == parm_name] <- parm_value
      }
    }
    sets[[as.character(year)]] <- parms
  }
  return(sets)
}

# 2.3 Baseline Scenario 
parms_base <- parms_set(parms_baseline, NULL, NULL, end_time = 2050)
parms_base <- as.data.frame(do.call(rbind, parms_base)) %>% 
  rownames_to_column(var = "year") %>%
  mutate(year = str_extract(year, "^\\d+"))

# 2.4 ACF Scenario A - Xpert only, Xpert+ get treatment
parms_ACFA <- data.frame(
  parameter = c("iota_sub", "phi_sub"), 
  med = c(0.80, 0.05))

acfa_year <- c(2025, 2030) 

parms_ACFA <- parms_set(parms_baseline, parms_ACFA, acfa_year, end_time = 2050)
parms_ACFA <- as.data.frame(do.call(rbind, parms_ACFA)) %>% 
  rownames_to_column(var = "year") %>%
  mutate(year = str_extract(year, "^\\d+"))

# 2.5 ACF Scenario B - CXR only, CXR+ get Xpert, Xpert+ get treatment
parms_ACFB <- data.frame(
  parameter = c("iota_sub", "phi_sub", "iota_min", "phi_min"), 
  med = c(0.5, 0.05, 0.8, 0.05))

acfb_year <- c(2025, 2030) 

parms_ACFB <- parms_set(parms_baseline, parms_ACFB, acfb_year, end_time = 2050)
parms_ACFB <- as.data.frame(do.call(rbind, parms_ACFB)) %>% 
  rownames_to_column(var = "year") %>%
  mutate(year = str_extract(year, "^\\d+"))

# 2.6 ACF Scenario C - CXR and Xpert for all, Xpert+ get treatment
parms_ACFC <- data.frame(
  parameter = c("iota_sub", "phi_sub", "iota_min", "phi_min"), 
  med = c(0.8, 0.05, 0.8, 0.05))

acfc_year <- c(2025, 2030) 

parms_ACFC <- parms_set(parms_baseline, parms_ACFC, acfb_year, end_time = 2050)
parms_ACFC <- as.data.frame(do.call(rbind, parms_ACFC)) %>% 
  rownames_to_column(var = "year") %>%
  mutate(year = str_extract(year, "^\\d+"))

# 2.7 Demographic parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)

# 3. Models ==========
ode <- function(parms, end_time = 2050) {
  
  # Parameter functions
  beta <- function(times){
    beta <- parms[parms$parameter == 'beta' & parms$year == floor(times), 'med']
    return(beta)
  }
  kappa <- function(times){
    kappa <- parms[parms$parameter == 'kappa' & parms$year == floor(times), 'med']
    return(kappa)
  }
  gamma_infcle <- function(times){
    gamma_infcle <- parms[parms$parameter == 'gamma_infcle' & parms$year == floor(times), 'med']
    return(gamma_infcle)
  }
  lambda_infmin <- function(times){
    lambda_infmin <- parms[parms$parameter == 'lambda_infmin' & parms$year == floor(times), 'med']
    return(lambda_infmin)
  }
  gamma_mincle <- function(times){
    gamma_mincle <- parms[parms$parameter == 'gamma_mincle' & parms$year == floor(times), 'med']
    return(gamma_mincle)
  }
  theta_cleinf <- function(times){
    theta_cleinf <- parms[parms$parameter == 'theta_cleinf' & parms$year == floor(times), 'med']
    return(theta_cleinf)
  }
  iota_min <- function(times){
    iota_min <- parms[parms$parameter == 'iota_min' & parms$year == floor(times), 'med']
    return(iota_min)
  }
  phi_min <- function(times){
    phi_min <- parms[parms$parameter == 'phi_min' & parms$year == floor(times), 'med']
    return(phi_min)
  }
  lambda_minsub <- function(times){
    lambda_minsub <- parms[parms$parameter == 'lambda_minsub' & parms$year == floor(times), 'med']
    return(lambda_minsub)
  }
  lambda_infsub <- function(times){
    lambda_infsub <- parms[parms$parameter == 'lambda_infsub' & parms$year == floor(times), 'med']
    return(lambda_infsub)
  }
  gamma_submin <- function(times){
    gamma_submin <- parms[parms$parameter == 'gamma_submin' & parms$year == floor(times), 'med']
    return(gamma_submin)
  }
  iota_sub <- function(times){
    iota_sub <- parms[parms$parameter == 'iota_sub' & parms$year == floor(times), 'med']
    return(iota_sub)
  }
  phi_sub <- function(times){
    phi_sub <- parms[parms$parameter == 'phi_sub' & parms$year == floor(times), 'med']
    return(phi_sub)
  }
  lambda_subcln <- function(times){
    lambda_subcln <- parms[parms$parameter == 'lambda_subcln' & parms$year == floor(times), 'med']
    return(lambda_subcln)
  }
  gamma_clnsub <- function(times){
    gamma_clnsub <- parms[parms$parameter == 'gamma_clnsub' & parms$year == floor(times), 'med']
    return(gamma_clnsub)
  }
  omega <- function(times){
    omega <- parms[parms$parameter == 'omega' & parms$year == floor(times), 'med']
    return(omega)
  }
  iota_cln <- function(times){
    iota_cln <- parms[parms$parameter == 'iota_cln' & parms$year == floor(times), 'med']
    return(iota_cln)
  }
  phi_cln <- function(times){
    phi_cln <- parms[parms$parameter == 'phi_cln' & parms$year == floor(times), 'med']
    return(phi_cln)
  }
  tau_min <- function(times){
    tau_min <- parms[parms$parameter == 'tau_min' & parms$year == floor(times), 'med']
    return(tau_min)
  }
  tau_sub <- function(times){
    tau_sub <- parms[parms$parameter == 'tau_sub' & parms$year == floor(times), 'med']
    return(tau_sub)
  }
  
  # Static parameters  
  delta <- 1      # Treatment year
  theta_recinf <- 1   # REINF: Recovered -> Infected (No protection)
  
  # Time-dependent parameters
  mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
  nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
  chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)
  
  des <- function(times, state, parms) {
    with(as.list(c(times, state, parms)), {
      
      dN_RL  = nu(times)*(N_RL+I_RL+O_RL+M_RL+RM_RL+S_RL+RS_RL+C_RL+RC_RL+P_RL) - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - mu(times)*N_RL
      dI_RL  = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf(times))+(P_RL*theta_recinf)) - gamma_infcle(times)*I_RL - lambda_infmin(times)*I_RL - lambda_infsub(times)*I_RL - mu(times)*I_RL
      dO_RL  = gamma_infcle(times)*I_RL + gamma_mincle(times)*M_RL - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf(times) - mu(times)*O_RL
      dM_RL  = lambda_infmin(times)*I_RL + gamma_submin(times)*S_RL - gamma_mincle(times)*M_RL - lambda_minsub(times)*M_RL - iota_min(times)*M_RL + phi_min(times)*RM_RL + tau_min(times)*P_RL - mu(times)*M_RL
      dRM_RL = iota_min(times)*M_RL - phi_min(times)*RM_RL - delta*RM_RL - mu(times)*RM_RL
      dS_RL  = lambda_infsub(times)*I_RL + lambda_minsub(times)*M_RL + gamma_clnsub(times)*C_RL - gamma_submin(times)*S_RL - lambda_subcln(times)*S_RL - iota_sub(times)*S_RL + phi_sub(times)*RS_RL + tau_sub(times)*P_RL - mu(times)*S_RL
      dRS_RL = iota_sub(times)*S_RL - phi_sub(times)*RS_RL - delta*RS_RL - mu(times)*RS_RL
      dC_RL  = lambda_subcln(times)*S_RL - gamma_clnsub(times)*C_RL - omega(times)*C_RL - mu(times)*C_RL - iota_cln(times)*C_RL + phi_cln(times)*RC_RL
      dRC_RL = iota_cln(times)*C_RL - phi_cln(times)*RC_RL - delta*RC_RL - mu(times)*RC_RL
      dP_RL  = delta*(RM_RL+RS_RL+RC_RL) - tau_min(times)*P_RL - tau_sub(times)*P_RL - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu(times)*P_RL
      
      dN_RH  = nu(times)*(N_RH+I_RH+O_RH+M_RH+RM_RH+S_RH+RS_RH+C_RH+RC_RH+P_RH) - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu(times)*N_RH
      dI_RH  = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(O_RH*theta_cleinf(times))+(P_RH*theta_recinf)) - gamma_infcle(times)*I_RH - lambda_infmin(times)*I_RH - lambda_infsub(times)*I_RH - mu(times)*I_RH
      dO_RH  = gamma_infcle(times)*I_RH + gamma_mincle(times)*M_RH - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_cleinf(times) - mu(times)*O_RH
      dM_RH  = lambda_infmin(times)*I_RH + gamma_submin(times)*S_RH - gamma_mincle(times)*M_RH - lambda_minsub(times)*M_RH - iota_min(times)*M_RH + phi_min(times)*RM_RH + tau_min(times)*P_RH - mu(times)*M_RH
      dRM_RH = iota_min(times)*M_RH - phi_min(times)*RM_RH - delta*RM_RH - mu(times)*RM_RH
      dS_RH  = lambda_infsub(times)*I_RH + lambda_minsub(times)*M_RH + gamma_clnsub(times)*C_RH - gamma_submin(times)*S_RH - lambda_subcln(times)*S_RH - iota_sub(times)*S_RH + phi_sub(times)*RS_RH + tau_sub(times)*P_RH - mu(times)*S_RH
      dRS_RH = iota_sub(times)*S_RH - phi_sub(times)*RS_RH - delta*RS_RH - mu(times)*RS_RH
      dC_RH  = lambda_subcln(times)*S_RH - gamma_clnsub(times)*C_RH - omega(times)*C_RH - mu(times)*C_RH - iota_cln(times)*C_RH + phi_cln(times)*RC_RH
      dRC_RH = iota_cln(times)*C_RH - phi_cln(times)*RC_RH - delta*RC_RH - mu(times)*RC_RH
      dP_RH  = delta*(RM_RH+RS_RH+RC_RH) - tau_min(times)*P_RH - tau_sub(times)*P_RH - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_recinf - mu(times)*P_RH
      
      dN_UL  = nu(times)*(N_UL+I_UL+O_UL+M_UL+RM_UL+S_UL+RS_UL+C_UL+RC_UL+P_UL) - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu(times)*N_UL
      dI_UL  = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(O_UL*theta_cleinf(times))+(P_UL*theta_recinf)) - gamma_infcle(times)*I_UL - lambda_infmin(times)*I_UL - lambda_infsub(times)*I_UL - mu(times)*I_UL
      dO_UL  = gamma_infcle(times)*I_UL + gamma_mincle(times)*M_UL - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_cleinf(times) - mu(times)*O_UL
      dM_UL  = lambda_infmin(times)*I_UL + gamma_submin(times)*S_UL - gamma_mincle(times)*M_UL - lambda_minsub(times)*M_UL - iota_min(times)*M_UL + phi_min(times)*RM_UL + tau_min(times)*P_UL - mu(times)*M_UL
      dRM_UL = iota_min(times)*M_UL - phi_min(times)*RM_UL - delta*RM_UL - mu(times)*RM_UL
      dS_UL  = lambda_infsub(times)*I_UL + lambda_minsub(times)*M_UL + gamma_clnsub(times)*C_UL - gamma_submin(times)*S_UL - lambda_subcln(times)*S_UL - iota_sub(times)*S_UL + phi_sub(times)*RS_UL + tau_sub(times)*P_UL - mu(times)*S_UL
      dRS_UL = iota_sub(times)*S_UL - phi_sub(times)*RS_UL - delta*RS_UL - mu(times)*RS_UL
      dC_UL  = lambda_subcln(times)*S_UL - gamma_clnsub(times)*C_UL - omega(times)*C_UL - mu(times)*C_UL - iota_cln(times)*C_UL + phi_cln(times)*RC_UL
      dRC_UL = iota_cln(times)*C_UL - phi_cln(times)*RC_UL - delta*RC_UL - mu(times)*RC_UL
      dP_UL  = delta*(RM_UL+RS_UL+RC_UL) - tau_min(times)*P_UL - tau_sub(times)*P_UL - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_recinf - mu(times)*P_UL
      
      dN_UH  = nu(times)*(N_UH+I_UH+O_UH+M_UH+RM_UH+S_UH+RS_UH+C_UH+RC_UH+P_UH) - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu(times)*N_UH
      dI_UH  = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(O_UH*theta_cleinf(times))+(P_UH*theta_recinf)) - gamma_infcle(times)*I_UH - lambda_infmin(times)*I_UH - lambda_infsub(times)*I_UH - mu(times)*I_UH
      dO_UH  = gamma_infcle(times)*I_UH + gamma_mincle(times)*M_UH - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_cleinf(times) - mu(times)*O_UH
      dM_UH  = lambda_infmin(times)*I_UH + gamma_submin(times)*S_UH - gamma_mincle(times)*M_UH - lambda_minsub(times)*M_UH - iota_min(times)*M_UH + phi_min(times)*RM_UH + tau_min(times)*P_UH - mu(times)*M_UH
      dRM_UH = iota_min(times)*M_UH - phi_min(times)*RM_UH - delta*RM_UH - mu(times)*RM_UH
      dS_UH  = lambda_infsub(times)*I_UH + lambda_minsub(times)*M_UH + gamma_clnsub(times)*C_UH - gamma_submin(times)*S_UH - lambda_subcln(times)*S_UH - iota_sub(times)*S_UH + phi_sub(times)*RS_UH + tau_sub(times)*P_UH - mu(times)*S_UH
      dRS_UH = iota_sub(times)*S_UH - phi_sub(times)*RS_UH - delta*RS_UH - mu(times)*RS_UH
      dC_UH  = lambda_subcln(times)*S_UH - gamma_clnsub(times)*C_UH - omega(times)*C_UH - mu(times)*C_UH - iota_cln(times)*C_UH + phi_cln(times)*RC_UH
      dRC_UH = iota_cln(times)*C_UH - phi_cln(times)*RC_UH - delta*RC_UH - mu(times)*RC_UH
      dP_UH  = delta*(RM_UH+RS_UH+RC_UH) - tau_min(times)*P_UH - tau_sub(times)*P_UH - ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_recinf - mu(times)*P_UH
      
      return(list(c(
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH, dRM_RL, dRM_RH, dRM_UL, dRM_UH,
        dS_RL, dS_RH, dS_UL, dS_UH, dRS_RL, dRS_RH, dRS_UL, dRS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dP_RL, dP_RH, dP_UL, dP_UH),
        Pop   = (N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+O_RL+O_RH+O_UL+O_UH+M_RL+M_RH+M_UL+M_UH+dRM_RL+dRM_RH+dRM_UL+dRM_UH+S_RL+S_RH+S_UL+S_UH+RS_RL+RS_RH+RS_UL+RS_UH+C_RL+C_RH+C_UL+C_UH+RC_RL+RC_RH+RC_UL+RC_UH+P_RL+P_RH+P_UL+P_UH), # Total population
        PRur  = (N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH), # Population rural
        PUrb  = (N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH), # Population urban
        PHig  = (N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH), # Population high SES
        PLow  = (N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL), # Population low SES
        PRL   = (N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL), # Population Rural - Low SES
        PRH   = (N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH), # Population Rural - High SES
        PUL   = (N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL), # Population Urban - Low SES
        PUH   = (N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH), # Population Urban - High SES
        Sub   = ((S_RL+S_RH+S_UL+S_UH)/chi(times)*1e5), # Subclinical TB (per 100k)
        SubLo = ((S_RL+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Subclinical TB in low SES (per 100k)
        SubHi = ((S_RH+S_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Subclinical TB in high SES (per 100k)
        SubUr = ((S_UH+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Subclinical TB in urban (per 100k)
        SubRu = ((S_RH+S_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Subclinical TB in urban (per 100k)
        SubRL = (S_RL/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Subclinical TB Rural - Low SES (per 100k)
        SubRH = (S_RH/(N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Subclinical TB Rural - High SES (per 100k)
        SubUL = (S_UL/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL)*1e5), # Subclinical TB Urban - Low SES (per 100k)
        SubUH = (S_UH/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Subclinical TB Urban - High SES (per 100k)
        Cln   = ((C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5), # Clinical TB (per 100k)
        ClnLo = ((S_RL+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Clinical TB in low SES (per 100k)
        ClnHi = ((S_RH+S_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Clinical TB in high SES (per 100k)
        ClnUr = ((S_UH+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Clinical TB in urban (per 100k)
        ClnRu = ((S_RH+S_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Clinical TB in urban (per 100k)
        ClnRL = (C_RL/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Clinical TB Rural - Low SES (per 100k)
        ClnRH = (C_RH/(N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Clinical TB Rural - High SES (per 100k)
        ClnUL = (C_UL/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL)*1e5), # Clinical TB Urban - Low SES (per 100k)
        ClnUH = (C_UH/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Clinical TB Urban - High SES (per 100k)
        TBc   = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5, # Infectious TB (per 100k)
        TBcLo = ((S_RL+S_UL+C_RL+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Infectious TB in low SES (per 100k)
        TBcHi = ((S_RH+S_UH+C_RH+C_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Infectious TB in high SES (per 100k)
        TBcUr = ((S_UH+S_UL+C_UH+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Infectious TB in urban (per 100k)
        TBcRu = ((S_RH+S_RL+C_RH+C_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Infectious TB in urban (per 100k)
        TBcRL = ((S_RL+C_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+RS_RL+C_RL+RC_RL+P_RL)*1e5), # Infectious TB Rural - Low SES (per 100k)
        TBcRH = ((S_RH+C_RH)/(N_RH+I_RH+O_RH+M_RH+RM_RH+RS_RH+C_RH+RC_RH+P_RH)*1e5), # Infectious TB Rural - High SES (per 100k)
        TBcUL = ((S_UL+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+RS_UL+C_UL+RC_UL+P_UL)*1e5), # Infectious TB Urban - Low SES (per 100k)
        TBcUH = ((S_UH+C_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+RS_UH+C_UH+RC_UH+P_UH)*1e5), # Infectious TB Urban - High SES (per 100k)
        Mor   = (omega(times)*(C_RL+C_RH+C_UL+C_UH))/chi(times)*1e5, # Clinical TB mortality per time (per 100k)
        Dxs   = (iota_cln(times)*(C_RL+C_RH+C_UL+C_UH))/chi(times)*1e5, # Notifications cTB per time in adults (per 100k)
        Spr   = (S_RL+S_RH+S_UL+S_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB/infectious TB
        MinIf = (M_RL+M_RH+M_UL+M_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion minTB
        SubIf = (S_RL+S_RH+S_UL+S_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB
        ClnIf = (C_RL+C_RH+C_UL+C_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion clnTB
        URs   = (S_UL+S_UH)/(S_RL+S_RH), # Relative urban/rural in scTB
        URc   = (C_UL+C_UH)/(C_RL+C_RH), # Relative urban/rural in cTB
        URt   = (S_UL+S_UH+C_UL+C_UH)/(S_UL+S_UH+C_UL+C_UH+S_RL+S_RH+C_RL+C_RH), # Proportion infectious TB in urban
        HLs   = (S_RH+S_UH)/(S_RL+S_UL), # Relative high/low SES in scTB
        HLc   = (C_RH+C_UH)/(C_RL+C_UL), # Relative high/low SES in cTB
        HLt   = (S_RH+S_UH+C_RH+C_UH)/(S_RL+S_UL+C_RL+C_UL+S_RH+S_UH+C_RH+C_UH), # Proportion infectious TB in high SES
        ARIsi = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+N_RH+N_UL+N_UH), # ARI: Susceptible -> Infected (%) 
        ARIoi = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_RL+O_RH+O_UL+O_UH)*theta_cleinf(times), # ARI: Cleared -> Infected (%) 
        ARIpi = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_RL+P_RH+P_UL+P_UH)*theta_recinf, # ARI: Recovered -> Infected (%) 
        ARI   = ((beta(times)/chi(times))*((kappa(times)*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
    })
  }
  
  yini <- c(N_RL = bWPP[, "N_RL"], N_RH = bWPP[, "N_RH"], N_UL = bWPP[, "N_UL"], N_UH = bWPP[, "N_UH"], 
            I_RL = bWPP[, "I_RL"], I_RH = bWPP[, "I_RH"], I_UL = bWPP[, "I_UL"], I_UH = bWPP[, "I_UH"], 
            O_RL = bWPP[, "O_RL"], O_RH = bWPP[, "O_RH"], O_UL = bWPP[, "O_UL"], O_UH = bWPP[, "O_UH"],
            M_RL = bWPP[, "M_RL"], M_RH = bWPP[, "M_RH"], M_UL = bWPP[, "M_UL"], M_UH = bWPP[, "M_UH"], 
            RM_RL = bWPP[, "RM_RL"], RM_RH = bWPP[, "RM_RH"], RM_UL = bWPP[, "RM_UL"], RM_UH = bWPP[, "RM_UH"],
            S_RL = bWPP[, "S_RL"], S_RH = bWPP[, "S_RH"], S_UL = bWPP[, "S_UL"], S_UH = bWPP[, "S_UH"], 
            RS_RL = bWPP[, "RS_RL"], RS_RH = bWPP[, "RS_RH"], RS_UL = bWPP[, "RS_UL"], RS_UH = bWPP[, "RS_UH"],
            C_RL = bWPP[, "C_RL"], C_RH = bWPP[, "C_RH"], C_UL = bWPP[, "C_UL"], C_UH = bWPP[, "C_UH"], 
            RC_RL = bWPP[, "RC_RL"], RC_RH = bWPP[, "RC_RH"], RC_UL = bWPP[, "RC_UL"], RC_UH = bWPP[, "RC_UH"],
            P_RL = bWPP[, "P_RL"], P_RH = bWPP[, "P_RH"], P_UL = bWPP[, "P_UL"], P_UH = bWPP[, "P_UH"]) 
  
  times <- seq(2020, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}  

# 4. Outputs ==========
outbase <- as.data.frame(ode(parms_base)) %>% 
  mutate(type = "base")

outacfa <- as.data.frame(ode(parms_ACFA)) %>% 
  mutate(type = "acfa")

outacfb <- as.data.frame(ode(parms_ACFB)) %>% 
  mutate(type = "acfb")

outacfc <- as.data.frame(ode(parms_ACFC)) %>% 
  mutate(type = "acfc")

# 5. Data curation ==========
outs <- rbind(outbase, outacfa, outacfb, outacfc) %>% 
  pivot_longer(cols = -c(time, type), names_to = "var", values_to = "val")

# 6. Plots ==========
prev_targets <- c(150, 100, 50, 25)
acf_year <- unique(acfa_year, acfb_year, acfc_year)
dis_state <- factor(c("Minimal","Subclinical","Clinical"), levels = c("Minimal","Subclinical","Clinical"))
scenarios <- c("Scenario 1", "Scenario 2", "Scenario 3", "Baseline")
type_label <- labeller(type = c("acfa" = "Scenario 1", "acfb" = "Scenario 2", "acfc" = "Scenario 3", "base" = "Baseline"))
dem_urbrur <- c("Rural", "Urban")
dem_ses <- c("High","Low")

# TB prevalence per scenarios 
tiff(here("plots", "TBprev.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "TBc"), aes(x = time, y = val, colour = type)) +
  scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,200,25), limits = c(0,200), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_hline(yintercept = prev_targets, linetype = "dashed", color = "gray") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Prevalence rate (urban vs rural)
tiff(here("plots", "Prop_urbvrur.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "URt"), aes(x = time, y = val, colour = type)) +
  scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB urban vs rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Total prevalence (low vs high SES)
tiff(here("plots", "Prop_lovhi.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "HLt"), aes(x = time, y = (1-val), colour = type)) +
  scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB low vs high") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion subclinical
tiff(here("plots", "Prop_scTB.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "Spr"), aes(x = time, y = val, colour = type)) +
  scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion subclinical TB (scTB/infTB)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Disease states
tiff(here("plots", "Prop_disstates.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  facet_wrap(~type, labeller = type_label) +
  geom_area(data = filter(outs, var %in% c("MinIf","SubIf","ClnIf")), aes(x = time, y = val, fill = reorder(var, -val)), position = "fill") +
  scale_fill_manual(values = c("#4DC4CB","#FDC75D","#F58B65"), name = "State:", labels = dis_state) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion disease state") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()
