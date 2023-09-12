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

# 2.4 ACF Scenario A
parms_ACFA <- data.frame(
  parameter = c("iota_min", "phi_min"), 
  med = c(0.1, 0.05))

acf_year <- c(2025, 2027) 

parms_ACFA <- parms_set(parms_baseline, parms_ACFA, acf_year, end_time = 2050)
parms_ACFA <- as.data.frame(do.call(rbind, parms_ACFA)) %>% 
  rownames_to_column(var = "year") %>%
  mutate(year = str_extract(year, "^\\d+"))

parms <- parms_base
times <- 2020

mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)

chi(times)

# 3. Models ==========
ode <- function(times, state, parms) {
  
  # Parameters
  beta <- parms[parms$parameter == 'beta' & parms$year == floor(times), 'med']
  kappa <- parms[parms$parameter == 'kappa' & parms$year == floor(times), 'med']
  gamma_infcle <- parms[parms$parameter == 'gamma_infcle' & parms$year == floor(times), 'med']
  lambda_infmin <- parms[parms$parameter == 'lambda_infmin' & parms$year == floor(times), 'med']
  gamma_mincle <- parms[parms$parameter == 'gamma_mincle' & parms$year == floor(times), 'med']
  theta_cleinf <- parms[parms$parameter == 'theta_cleinf' & parms$year == floor(times), 'med']
  iota_min <- parms[parms$parameter == 'iota_min' & parms$year == floor(times), 'med']
  phi_min <- parms[parms$parameter == 'phi_min', 'med']
  lambda_minsub <- parms[parms$parameter == 'lambda_minsub' & parms$year == floor(times), 'med']
  lambda_infsub <- parms[parms$parameter == 'lambda_infsub' & parms$year == floor(times), 'med']
  gamma_submin <- parms[parms$parameter == 'gamma_submin' & parms$year == floor(times), 'med']
  iota_sub <- parms[parms$parameter == 'iota_sub' & parms$year == floor(times), 'med']
  phi_sub <- parms[parms$parameter == 'phi_sub' & parms$year == floor(times), 'med']
  lambda_subcln <- parms[parms$parameter == 'lambda_subcln' & parms$year == floor(times), 'med']
  gamma_clnsub <- parms[parms$parameter == 'gamma_clnsub' & parms$year == floor(times), 'med']
  omega <- parms[parms$parameter == 'omega' & parms$year == floor(times), 'med']
  iota_cln <- parms[parms$parameter == 'iota_cln' & parms$year == floor(times), 'med']
  phi_cln <- parms[parms$parameter == 'phi_cln' & parms$year == floor(times), 'med']
  tau_min <- parms[parms$parameter == 'tau_min' & parms$year == floor(times), 'med']
  tau_sub <- parms[parms$parameter == 'tau_sub' & parms$year == floor(times), 'med']
  
  # Static parameters  
  delta <- 1      # Treatment year
  theta_recinf <- 1   # REINF: Recovered -> Infected (No protection)
  
  # Time-dependent parameters
  mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
  nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
  chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)
  
  with(as.list(c(times, state, parms)), {
    
    dN_RL  = nu(times)*(N_RL+I_RL+O_RL+M_RL+RM_RL+S_RL+RS_RL+C_RL+RC_RL+P_RL) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - mu(times)*N_RL
    dI_RL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf)+(P_RL*theta_recinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - mu(times)*I_RL
    dO_RL  = gamma_infcle*I_RL + gamma_mincle*M_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf - mu(times)*O_RL
    dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - iota_min*M_RL + phi_min*RM_RL + tau_min*P_RL - mu(times)*M_RL
    dRM_RL = iota_min*M_RL - phi_min*RM_RL - delta*RM_RL - mu(times)*RM_RL
    dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - iota_sub*S_RL + phi_sub*RS_RL + tau_sub*P_RL - mu(times)*S_RL
    dRS_RL = iota_sub*S_RL - phi_sub*RS_RL - delta*RS_RL - mu(times)*RS_RL
    dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - omega*C_RL - mu(times)*C_RL - iota_cln*C_RL + phi_cln*RC_RL
    dRC_RL = iota_cln*C_RL - phi_cln*RC_RL - delta*RC_RL - mu(times)*RC_RL
    dP_RL  = delta*(RM_RL+RS_RL+RC_RL) - tau_min*P_RL - tau_sub*P_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu(times)*P_RL
    
    dN_RH  = nu(times)*(N_RH+I_RH+O_RH+M_RH+RM_RH+S_RH+RS_RH+C_RH+RC_RH+P_RH) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu(times)*N_RH
    dI_RH  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(O_RH*theta_cleinf)+(P_RH*theta_recinf)) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu(times)*I_RH
    dO_RH  = gamma_infcle*I_RH + gamma_mincle*M_RH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_cleinf - mu(times)*O_RH
    dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH - iota_min*M_RH + phi_min*RM_RH + tau_min*P_RH - mu(times)*M_RH
    dRM_RH = iota_min*M_RH - phi_min*RM_RH - delta*RM_RH - mu(times)*RM_RH
    dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH - iota_sub*S_RH + phi_sub*RS_RH + tau_sub*P_RH - mu(times)*S_RH
    dRS_RH = iota_sub*S_RH - phi_sub*RS_RH - delta*RS_RH - mu(times)*RS_RH
    dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - omega*C_RH - mu(times)*C_RH - iota_cln*C_RH + phi_cln*RC_RH
    dRC_RH = iota_cln*C_RH - phi_cln*RC_RH - delta*RC_RH - mu(times)*RC_RH
    dP_RH  = delta*(RM_RH+RS_RH+RC_RH) - tau_min*P_RH - tau_sub*P_RH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_recinf - mu(times)*P_RH
    
    dN_UL  = nu(times)*(N_UL+I_UL+O_UL+M_UL+RM_UL+S_UL+RS_UL+C_UL+RC_UL+P_UL) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu(times)*N_UL
    dI_UL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(O_UL*theta_cleinf)+(P_UL*theta_recinf)) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu(times)*I_UL
    dO_UL  = gamma_infcle*I_UL + gamma_mincle*M_UL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_cleinf - mu(times)*O_UL
    dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL - iota_min*M_UL + phi_min*RM_UL + tau_min*P_UL - mu(times)*M_UL
    dRM_UL = iota_min*M_UL - phi_min*RM_UL - delta*RM_UL - mu(times)*RM_UL
    dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL - iota_sub*S_UL + phi_sub*RS_UL + tau_sub*P_UL - mu(times)*S_UL
    dRS_UL = iota_sub*S_UL - phi_sub*RS_UL - delta*RS_UL - mu(times)*RS_UL
    dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - omega*C_UL - mu(times)*C_UL - iota_cln*C_UL + phi_cln*RC_UL
    dRC_UL = iota_cln*C_UL - phi_cln*RC_UL - delta*RC_UL - mu(times)*RC_UL
    dP_UL  = delta*(RM_UL+RS_UL+RC_UL) - tau_min*P_UL - tau_sub*P_UL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_recinf - mu(times)*P_UL
    
    dN_UH  = nu(times)*(N_UH+I_UH+O_UH+M_UH+RM_UH+S_UH+RS_UH+C_UH+RC_UH+P_UH) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu(times)*N_UH
    dI_UH  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(O_UH*theta_cleinf)+(P_UH*theta_recinf)) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu(times)*I_UH
    dO_UH  = gamma_infcle*I_UH + gamma_mincle*M_UH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_cleinf - mu(times)*O_UH
    dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH - iota_min*M_UH + phi_min*RM_UH + tau_min*P_UH - mu(times)*M_UH
    dRM_UH = iota_min*M_UH - phi_min*RM_UH - delta*RM_UH - mu(times)*RM_UH
    dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH - iota_sub*S_UH + phi_sub*RS_UH + tau_sub*P_UH - mu(times)*S_UH
    dRS_UH = iota_sub*S_UH - phi_sub*RS_UH - delta*RS_UH - mu(times)*RS_UH
    dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - omega*C_UH - mu(times)*C_UH - iota_cln*C_UH + phi_cln*RC_UH
    dRC_UH = iota_cln*C_UH - phi_cln*RC_UH - delta*RC_UH - mu(times)*RC_UH
    dP_UH  = delta*(RM_UH+RS_UH+RC_UH) - tau_min*P_UH - tau_sub*P_UH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_recinf - mu(times)*P_UH
    
    return(list(c(
      dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH, dRM_RL, dRM_RH, dRM_UL, dRM_UH,
      dS_RL, dS_RH, dS_UL, dS_UH, dRS_RL, dRS_RH, dRS_UL, dRS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dP_RL, dP_RH, dP_UL, dP_UH),
      Pop = (N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+O_RL+O_RH+O_UL+O_UH+M_RL+M_RH+M_UL+M_UH+dRM_RL+dRM_RH+dRM_UL+dRM_UH+S_RL+S_RH+S_UL+S_UH+RS_RL+RS_RH+RS_UL+RS_UH+C_RL+C_RH+C_UL+C_UH+RC_RL+RC_RH+RC_UL+RC_UH+P_RL+P_RH+P_UL+P_UH), # Total population
      Sub = (S_RL+S_RH+S_UL+S_UH)/chi(times)*100000, # Subclinical TB (per 100k)
      Cln = (C_RL+C_RH+C_UL+C_UH)/chi(times)*100000, # Clinical TB (per 100k)
      TBc = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/chi(times)*100000, # All TB (per 100k)
      Mor = (omega*(C_RL+C_RH+C_UL+C_UH))/chi(times)*100000, # Clinical TB mortality per time (per 100k)
      Dxs = (iota_cln*(C_RL+C_RH+C_UL+C_UH))/chi(times)*100000, # Notifications cTB per time in adults (per 100k)
      Spr = (S_RL+S_RH+S_UL+S_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB
      URs = (S_UL+S_UH)/(S_RL+S_RH), # Relative urban/rural in scTB
      URc = (C_UL+C_UH)/(C_RL+C_RH), # Relative urban/rural in cTB
      HLs = (S_RH+S_UH)/(S_RL+S_UL), # Relative high/low SES in scTB
      HLc = (C_RH+C_UH)/(C_RL+C_UL), # Relative high/low SES in cTB
      ARIsi = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+N_RH+N_UL+N_UH), # ARI: Susceptible -> Infected (%) 
      ARIoi = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_RL+O_RH+O_UL+O_UH)*theta_cleinf, # ARI: Cleared -> Infected (%) 
      ARIpi = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_RL+P_RH+P_UL+P_UH)*theta_recinf, # ARI: Recovered -> Infected (%) 
      ARI = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))) # ARI
    )
  })
}

# 4. Outputs ==========

end_time <- 2050

times <- seq(2020, end_time, by = 1)

state <- c(N_RL = bWPP[, "N_RL"], N_RH = bWPP[, "N_RH"], N_UL = bWPP[, "N_UL"], N_UH = bWPP[, "N_UH"], 
           I_RL = bWPP[, "I_RL"], I_RH = bWPP[, "I_RH"], I_UL = bWPP[, "I_UL"], I_UH = bWPP[, "I_UH"], 
           O_RL = bWPP[, "O_RL"], O_RH = bWPP[, "O_RH"], O_UL = bWPP[, "O_UL"], O_UH = bWPP[, "O_UH"],
           M_RL = bWPP[, "M_RL"], M_RH = bWPP[, "M_RH"], M_UL = bWPP[, "M_UL"], M_UH = bWPP[, "M_UH"], 
           RM_RL = bWPP[, "RM_RL"], RM_RH = bWPP[, "RM_RH"], RM_UL = bWPP[, "RM_UL"], RM_UH = bWPP[, "RM_UH"],
           S_RL = bWPP[, "S_RL"], S_RH = bWPP[, "S_RH"], S_UL = bWPP[, "S_UL"], S_UH = bWPP[, "S_UH"], 
           RS_RL = bWPP[, "RS_RL"], RS_RH = bWPP[, "RS_RH"], RS_UL = bWPP[, "RS_UL"], RS_UH = bWPP[, "RS_UH"],
           C_RL = bWPP[, "C_RL"], C_RH = bWPP[, "C_RH"], C_UL = bWPP[, "C_UL"], C_UH = bWPP[, "C_UH"], 
           RC_RL = bWPP[, "RC_RL"], RC_RH = bWPP[, "RC_RH"], RC_UL = bWPP[, "RC_UL"], RC_UH = bWPP[, "RC_UH"],
           P_RL = bWPP[, "P_RL"], P_RH = bWPP[, "P_RH"], P_UL = bWPP[, "P_UL"], P_UH = bWPP[, "P_UH"]) 

out <- deSolve::ode(state, times, ode, parms_base)

# 5. Plots ==========
ggplot() +
  geom_line(data = base, aes(x = time, y = TBc), colour = "#CE2931") +
  scale_y_continuous(limits = c(0,200)) +
  labs(x = "Year", y = "TB prevalence rate (per 100K)") +
  theme_minimal()
