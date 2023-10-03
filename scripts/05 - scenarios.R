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
parms <- import(here("data","parameters_manual.Rdata"))
WPP <- import(here("data","pop","WPP.Rdata"))
WUP <- import(here("data","pop","WUP.Rdata"))
bWPP <- import(here("data","basepop.Rdata"))

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  filter(!parameter %in% c('omega_ini', 'iota_cln_ini', 'phi_cln_ini', 
                           'rho_ini', 'rho_fin','sigma_ini','sigma_fin')) %>% 
  mutate(parameter = case_when(
    parameter == 'omega_fin' ~ 'omega', parameter == 'iota_cln_fin' ~ 'iota_cln', 
    parameter == 'phi_cln_fin' ~ 'phi_cln', TRUE ~ parameter)) %>% 
  add_row(parameter = "phi_min", low = 0, med = 0, hig = 0) %>%  # Treatment failure Minimal
  add_row(parameter = "phi_sub", low = 0, med = 0, hig = 0) %>%  # Treatment failure Subclinical
  mutate(med = ifelse(parameter == "tau_sub", 0.005, med)) %>%
  mutate(med = ifelse(parameter == "tau_min", 0.01, med)) %>%
  select(parameter, med)

# 2.2 Intervention parameters
acf_year <- seq(2025,2034,2)

prop_target <- 1 # Proportion of population targeted for ACF (NB: Check for annual rate)
prop_reached <- 0.8 # Proportion of population participating in ACF
prop_notreached <- 1-prop_reached # Proportion of population not reached (NB: See difference per sociodemographic status)

# 2.3 Xpert
xpert_sens_cln <- 0.91 # Xpert Ultra (95%CI: 0.86-0.95) [Zifodya et al 2021]
xpert_sens_sub <- 0.80 # 0.69
xpert_sens_min <- 0.01 
xpert_fp_nds <- 0.01

prov_sputum <- 1

# 2.4 Chest X-ray
cxr_sens_cln <- 0.95
cxr_sens_sub <- 0.95
cxr_sens_min <- 0.40
cxr_fp_nds <- 0.17

# 2.3 ACF Scenario A - Xpert only, Xpert+ get treatment
acfa <- data.frame(parameter = character(), med = numeric()) %>% 
  add_row(parameter = "alpha_min", med = prop_target*prop_reached*prov_sputum*xpert_sens_min) %>% 
  add_row(parameter = "alpha_sub", med = prop_target*prop_reached*prov_sputum*xpert_sens_sub) %>%
  add_row(parameter = "alpha_cln", med = prop_target*prop_reached*prov_sputum*xpert_sens_cln) %>% 
  add_row(parameter = "alpha_nds", med = prop_target*prop_reached*prov_sputum*xpert_fp_nds)

# 2.4 ACF Scenario B - CXR only, CXR+ get Xpert, Xpert+ get treatment
acfb <- data.frame(parameter = character(), med = numeric()) %>% 
  add_row(parameter = "alpha_min", med = prop_target*prop_reached*cxr_sens_min) %>% 
  add_row(parameter = "alpha_sub", med = prop_target*prop_reached*cxr_sens_sub) %>%
  add_row(parameter = "alpha_cln", med = prop_target*prop_reached*cxr_sens_cln) %>% 
  add_row(parameter = "alpha_nds", med = prop_target*prop_reached*cxr_fp_nds)

# 2.5 ACF Scenario C - CXR and Xpert for all, Xpert+ get treatment
acfc <- data.frame(parameter = character(), med = numeric()) %>% 
  add_row(parameter = "alpha_min", med = prop_target*prop_reached*cxr_sens_min) %>% 
  add_row(parameter = "alpha_sub", med = prop_target*prop_reached*cxr_sens_sub) %>%
  add_row(parameter = "alpha_cln", med = prop_target*prop_reached*cxr_sens_cln) %>% 
  add_row(parameter = "alpha_nds", med = prop_target*prop_reached*cxr_fp_nds)

# 2.6 Demographic parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)
rhoU <- approxfun(WUP$year, WUP$urbprop, method = 'linear', rule = 2)
rhoR <- approxfun(WUP$year, WUP$rurprop, method = 'linear', rule = 2)

# 3. Models ==========
ode <- function(parms, interv = NULL, acf_times = NULL, end_time = 2050) {
  
  times <- seq(2020, end_time, by = 1)
  
  # Parameter functions
  beta <- parms[parms$parameter == 'beta', 'med']
  kappa <- parms[parms$parameter == 'kappa', 'med']
  gamma_infcle <- parms[parms$parameter == 'gamma_infcle', 'med']
  lambda_infmin <- parms[parms$parameter == 'lambda_infmin', 'med']
  gamma_mincle <- parms[parms$parameter == 'gamma_mincle', 'med']
  theta_cleinf <- parms[parms$parameter == 'theta_cleinf', 'med']
  iota_min <- parms[parms$parameter == 'iota_min', 'med']
  phi_min <- parms[parms$parameter == 'phi_min', 'med']
  lambda_minsub <- parms[parms$parameter == 'lambda_minsub', 'med']
  lambda_infsub <- parms[parms$parameter == 'lambda_infsub', 'med']
  gamma_submin <- parms[parms$parameter == 'gamma_submin', 'med']
  iota_sub <- parms[parms$parameter == 'iota_sub', 'med']
  phi_sub <- parms[parms$parameter == 'phi_sub', 'med']
  lambda_subcln <- parms[parms$parameter == 'lambda_subcln', 'med']
  gamma_clnsub <- parms[parms$parameter == 'gamma_clnsub', 'med']
  omega <- parms[parms$parameter == 'omega', 'med']
  iota_cln <- parms[parms$parameter == 'iota_cln', 'med']
  phi_cln <- parms[parms$parameter == 'phi_cln', 'med']
  tau_min <- parms[parms$parameter == 'tau_min', 'med']
  tau_sub <- parms[parms$parameter == 'tau_sub', 'med']
  
  # Intervention parameters
  if(is.null(interv)){
    alpha_min <- 0
    alpha_sub <- 0
    alpha_cln <- 0
  } else { 
  alpha_min <- interv[interv$parameter == 'alpha_min', 'med']
  alpha_sub <- interv[interv$parameter == 'alpha_sub', 'med']
  alpha_cln <- interv[interv$parameter == 'alpha_cln', 'med']
  }

  # Active-case finding switch
  if(is.null(acf_times)){
    acf <- function(t) 0
  } else { 
    values <- ifelse(floor(times) %in% acf_times, 1, 0)
    acf <- approxfun(times, values, rule = 2)
  }
  
  # Static parameters  
  delta <- 1      # Treatment year
  theta_recinf <- 1   # REINF: Recovered -> Infected (No protection)
  
  des <- function(times, state, parms) {
    with(as.list(c(times, state, parms)), {
      
      dN_RL  = nu(times)*rhoR(times)*(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - mu(times)*N_RL
      dI_RL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf)+(P_RL*theta_recinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - mu(times)*I_RL
      dO_RL  = gamma_infcle*I_RL + gamma_mincle*M_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf - mu(times)*O_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - acf(floor(times))*alpha_min*M_RL + phi_min*SM_RL + tau_min*P_RL - mu(times)*M_RL
      dSM_RL = acf(floor(times))*alpha_min*M_RL - phi_min*SM_RL - delta*SM_RL - mu(times)*SM_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - acf(floor(times))*alpha_sub*S_RL + phi_sub*SS_RL + tau_sub*P_RL - mu(times)*S_RL
      dSS_RL = acf(floor(times))*alpha_sub*S_RL - phi_sub*SS_RL - delta*SS_RL - mu(times)*SS_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - iota_cln*C_RL + phi_cln*RC_RL - acf(floor(times))*alpha_cln*C_RL + phi_cln*SC_RL - omega*C_RL - mu(times)*C_RL
      dRC_RL = iota_cln*C_RL - phi_cln*RC_RL - delta*RC_RL - mu(times)*RC_RL
      dSC_RL = acf(floor(times))*alpha_cln*C_RL - phi_cln*SC_RL - delta*SC_RL - mu(times)*SC_RL
      dP_RL  = delta*(SM_RL+SS_RL+RC_RL+SC_RL) - tau_min*P_RL - tau_sub*P_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu(times)*P_RL
      
      return(list(c(
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dO_RL, dO_RH, dO_UL, dO_UH, 
        dM_RL, dM_RH, dM_UL, dM_UH, dRM_RL, dRM_RH, dRM_UL, dRM_UH, dSM_RL, dSM_RH, dSM_UL, dSM_UH, 
        dS_RL, dS_RH, dS_UL, dS_UH, dRS_RL, dRS_RH, dRS_UL, dRS_UH, dSS_RL, dSS_RH, dSS_UL, dSS_UH,
        dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dSC_RL, dSC_RH, dSC_UL, dSC_UH,
        dP_RL, dP_RH, dP_UL, dP_UH),
        Pop   = (N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+O_RL+O_RH+O_UL+O_UH+
                   M_RL+M_RH+M_UL+M_UH+RM_RL+RM_RH+RM_UL+RM_UH+SM_RL+SM_RH+SM_UL+SM_UH+
                   S_RL+S_RH+S_UL+S_UH+RS_RL+RS_RH+RS_UL+RS_UH+SS_RL+SS_RH+SS_UL+SS_UH+
                   C_RL+C_RH+C_UL+C_UH+RC_RL+RC_RH+RC_UL+RC_UH+SC_RL+SC_RH+SC_UL+SC_UH+
                   P_RL+P_RH+P_UL+P_UH), # Total population
        PRur  = (N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH), # Population rural
        PUrb  = (N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH), # Population urban
        PHig  = (N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH), # Population high SES
        PLow  = (N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL), # Population low SES
        PRL   = (N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL), # Population Rural - Low SES
        PRH   = (N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH), # Population Rural - High SES
        PUL   = (N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL), # Population Urban - Low SES
        PUH   = (N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH), # Population Urban - High SES
        Sub   = ((S_RL+S_RH+S_UL+S_UH)/chi(times)*1e5), # Subclinical TB (per 100k)
        SubLo = ((S_RL+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Subclinical TB in low SES (per 100k)
        SubHi = ((S_RH+S_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Subclinical TB in high SES (per 100k)
        SubUr = ((S_UH+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Subclinical TB in urban (per 100k)
        SubRu = ((S_RH+S_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Subclinical TB in urban (per 100k)
        SubRL = (S_RL/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Subclinical TB Rural - Low SES (per 100k)
        SubRH = (S_RH/(N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Subclinical TB Rural - High SES (per 100k)
        SubUL = (S_UL/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL)*1e5), # Subclinical TB Urban - Low SES (per 100k)
        SubUH = (S_UH/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Subclinical TB Urban - High SES (per 100k)
        Cln   = ((C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5), # Clinical TB (per 100k)
        ClnLo = ((S_RL+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Clinical TB in low SES (per 100k)
        ClnHi = ((S_RH+S_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Clinical TB in high SES (per 100k)
        ClnUr = ((S_UH+S_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Clinical TB in urban (per 100k)
        ClnRu = ((S_RH+S_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Clinical TB in urban (per 100k)
        ClnRL = (C_RL/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Clinical TB Rural - Low SES (per 100k)
        ClnRH = (C_RH/(N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Clinical TB Rural - High SES (per 100k)
        ClnUL = (C_UL/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL)*1e5), # Clinical TB Urban - Low SES (per 100k)
        ClnUH = (C_UH/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Clinical TB Urban - High SES (per 100k)
        TBc   = ((S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5), # Infectious TB (per 100k)
        TBcLo = ((S_RL+S_UL+C_RL+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Infectious TB in low SES (per 100k)
        TBcHi = ((S_RH+S_UH+C_RH+C_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Infectious TB in high SES (per 100k)
        TBcUr = ((S_UH+S_UL+C_UH+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Infectious TB in urban (per 100k)
        TBcRu = ((S_RH+S_RL+C_RH+C_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Infectious TB in urban (per 100k)
        TBcRL = ((S_RL+C_RL)/(N_RL+I_RL+O_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL)*1e5), # Infectious TB Rural - Low SES (per 100k)
        TBcRH = ((S_RH+C_RH)/(N_RH+I_RH+O_RH+M_RH+RM_RH+SM_RH+S_RH+RS_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH)*1e5), # Infectious TB Rural - High SES (per 100k)
        TBcUL = ((S_UL+C_UL)/(N_UL+I_UL+O_UL+M_UL+RM_UL+SM_UL+S_UL+RS_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL)*1e5), # Infectious TB Urban - Low SES (per 100k)
        TBcUH = ((S_UH+C_UH)/(N_UH+I_UH+O_UH+M_UH+RM_UH+SM_UH+S_UH+RS_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH)*1e5), # Infectious TB Urban - High SES (per 100k)
        rMor  = (omega*(C_RL+C_RH+C_UL+C_UH))/chi(times)*1e5, # Clinical TB mortality per time (per 100k)
        tMor  = (omega*(C_RL+C_RH+C_UL+C_UH)), # Clinical TB mortality per time
        Dxs   = (iota_cln*(C_RL+C_RH+C_UL+C_UH))/chi(times)*1e5, # Notifications cTB per time in adults (per 100k)
        Spr   = (S_RL+S_RH+S_UL+S_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB/infectious TB
        Cpr   = (C_RL+C_RH+C_UL+C_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion cTB/infectious TB
        MinIf = (M_RL+M_RH+M_UL+M_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion minTB
        SubIf = (S_RL+S_RH+S_UL+S_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB
        ClnIf = (C_RL+C_RH+C_UL+C_UH)/(M_RL+M_RH+M_UL+M_UH+S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion clnTB
        URs   = (S_UL+S_UH)/(S_RL+S_RH), # Relative urban/rural in scTB
        URc   = (C_UL+C_UH)/(C_RL+C_RH), # Relative urban/rural in cTB
        URt   = (S_UL+S_UH+C_UL+C_UH)/(S_UL+S_UH+C_UL+C_UH+S_RL+S_RH+C_RL+C_RH), # Proportion infectious TB in urban
        HLs   = (S_RH+S_UH)/(S_RL+S_UL), # Relative high/low SES in scTB
        HLc   = (C_RH+C_UH)/(C_RL+C_UL), # Relative high/low SES in cTB
        HLt   = (S_RH+S_UH+C_RH+C_UH)/(S_RL+S_UL+C_RL+C_UL+S_RH+S_UH+C_RH+C_UH), # Proportion infectious TB in high SES
        ARI   = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
    })
  }
  
  yini <- c(N_RL = bWPP[, "N_RL"], N_RH = bWPP[, "N_RH"], N_UL = bWPP[, "N_UL"], N_UH = bWPP[, "N_UH"], 
            I_RL = bWPP[, "I_RL"], I_RH = bWPP[, "I_RH"], I_UL = bWPP[, "I_UL"], I_UH = bWPP[, "I_UH"], 
            O_RL = bWPP[, "O_RL"], O_RH = bWPP[, "O_RH"], O_UL = bWPP[, "O_UL"], O_UH = bWPP[, "O_UH"],
            M_RL = bWPP[, "M_RL"], M_RH = bWPP[, "M_RH"], M_UL = bWPP[, "M_UL"], M_UH = bWPP[, "M_UH"], 
            RM_RL = bWPP[, "RM_RL"], RM_RH = bWPP[, "RM_RH"], RM_UL = bWPP[, "RM_UL"], RM_UH = bWPP[, "RM_UH"],
            SM_RL = 0, SM_RH = 0, SM_UL = 0, SM_UH = 0,
            S_RL = bWPP[, "S_RL"], S_RH = bWPP[, "S_RH"], S_UL = bWPP[, "S_UL"], S_UH = bWPP[, "S_UH"], 
            RS_RL = bWPP[, "RS_RL"], RS_RH = bWPP[, "RS_RH"], RS_UL = bWPP[, "RS_UL"], RS_UH = bWPP[, "RS_UH"],
            SS_RL = 0, SS_RH = 0, SS_UL = 0, SS_UH = 0,
            C_RL = bWPP[, "C_RL"], C_RH = bWPP[, "C_RH"], C_UL = bWPP[, "C_UL"], C_UH = bWPP[, "C_UH"], 
            RC_RL = bWPP[, "RC_RL"], RC_RH = bWPP[, "RC_RH"], RC_UL = bWPP[, "RC_UL"], RC_UH = bWPP[, "RC_UH"],
            SC_RL = 0, SC_RH = 0, SC_UL = 0, SC_UH = 0,
            P_RL = bWPP[, "P_RL"], P_RH = bWPP[, "P_RH"], P_UL = bWPP[, "P_UL"], P_UH = bWPP[, "P_UH"]) 
  
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}  

# 4. Outputs ==========
outbase <- as.data.frame(ode(parms)) %>% 
  mutate(type = "base") 

outacfa <- as.data.frame(ode(parms, acfa, acf_year)) %>% 
  mutate(type = "acfa") 

outacfb <- as.data.frame(ode(parms, acfb, acf_year)) %>% 
  mutate(type = "acfb") 

# outacfc <- as.data.frame(ode(parms, acfc, acf_year)) %>% 
#   mutate(type = "acfc") 

# 5. Data curation ==========
#outs <- rbind(outbase, outacfa, outacfb, outacfc) %>% 
outs <- rbind(outbase, outacfa, outacfb) %>% 
  arrange(time) %>% 
  group_by(type) %>% 
  mutate(cumMor = cumsum(tMor)) %>% 
  group_by(time) %>% 
  mutate(dMor = tMor - tMor[type == 'base'], dcumMor = cumMor - cumMor[type == 'base'],
         pTBc = TBc/TBc[type == 'base']) %>% 
  pivot_longer(cols = -c(time, type), names_to = "var", values_to = "val")

#outpop <- rbind(outbase, outacfa, outacfb, outacfc) %>% 
outpop <- rbind(outbase, outacfa, outacfb) %>%
  select(time, type, Pop, PRur, PUrb, PHig, PLow, PRL, PRH, PUL, PUH) %>% 
  mutate(pPRur = PRur/Pop, pPUrb = PUrb/Pop, pPHig = PHig/Pop, pPLow = PLow/Pop,
         pPRL = PRL/Pop, pPRH = PRH/Pop, pPUL = PUL/Pop, pPUH = PUH/Pop) %>% 
  pivot_longer(cols = -c(time, type), names_to = "var", values_to = "val")

DALY = TB*multip
  
# 6. Plots ==========
prev_targets <- c(100, 50, 25, 10)
dis_state <- factor(c("Clinical","Subclinical","Minimal"))
inf_dis <- c("Clinical", "Subclinical")
#scenarios <- c("Scenario 1", "Scenario 2", "Scenario 3", "Baseline")
scenarios <- c("ACF: Xpert", "ACF: CXR", "Baseline")
#type_label <- labeller(type = c("acfa" = "Scenario 1", "acfb" = "Scenario 2", "acfc" = "Scenario 3", "base" = "Baseline"))
type_label <- labeller(type = c("acfa" = "ACF: Xpert", "acfb" = "ACF: CXR", "acfc", "base" = "Baseline"))
dem_urbrur <- c("Rural", "Urban")
dem_ses <- c("High","Low")

# TB prevalence per scenarios 
tiff(here("plots", "TBprev.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "TBc"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
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

# Annual risk of infection 
tiff(here("plots", "ARI.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "ARI"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(limits = c(0,0.0125), labels = scales::percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Annual risk of infection") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Proportion reduction of TB prevalence
tiff(here("plots", "TBreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pTBc"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(min(acf_year),2050,5), limits = c(min(acf_year),2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(min(acf_year),2050)) + 
  labs(x = "Year", y = "Proportion reduction of TB prevalence") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_segment(aes(x = 2028 - 2, xend = 2028 + 2, y = 0.31, yend = 0.31), color = "gray", linetype = "dashed") + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Proportion prevalence (urban vs rural)
tiff(here("plots", "Prop_urbvrur.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "URt"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB urban vs rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion prevalence (low vs high SES)
tiff(here("plots", "Prop_lovhi.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "HLt"), aes(x = time, y = (1-val), colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB low vs high") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Urban
tiff(here("plots", "Prop_UH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outpop) +
  geom_line(data = filter(outpop, var == "pPUH"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Rural
tiff(here("plots", "Prop_RH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outpop) +
  geom_line(data = filter(outpop, var == "pPRH"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Urban
tiff(here("plots", "Prop_UL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outpop) +
  geom_line(data = filter(outpop, var == "pPUL"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931", "#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion Low SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Rural
tiff(here("plots", "Prop_RL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outpop) +
  geom_line(data = filter(outpop, var == "pPRL"), aes(x = time, y = val, colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion Low SES - Rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion subclinical
tiff(here("plots", "Prop_scTB.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  facet_wrap(~type, labeller = type_label) +
  geom_area(data = filter(outs, var %in% c("Cpr","Spr")), aes(x = time, y = val, fill = var), position = "fill") +
  scale_fill_manual(values = c("#F58B65","#FDC75D"), name = NULL, labels = inf_dis) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion infectious TB") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Disease states
tiff(here("plots", "Prop_disstates.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  facet_wrap(~type, labeller = type_label) +
  geom_area(data = filter(outs, var %in% c("ClnIf","SubIf","MinIf")), aes(x = time, y = val, fill = reorder(var, -val, decreasing = TRUE)), position = "fill") +
  scale_fill_manual(values = c("#F58B65","#FDC75D","#4DC4CB"), name = "State:", labels = dis_state) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion disease state") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB deaths averted
tiff(here("plots", "TBdeathsavert.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dcumMor"), aes(x = time, y = abs(val), colour = type)) +
  #scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_color_manual(values = c("#CE2931","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative TB deaths averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

TBdeaths <- outs %>% 
  filter(time %in% c(2025,2030,2040,2050)) %>% 
  filter(var %in% c("tMor","dMor"))

TBdeathsaverttb <- TBdeaths %>%
  pivot_wider(names_from = time, values_from = val) %>% 
  slice(1,4,6,8) %>% 
  mutate_if(is.numeric, ~ round(., digits = 0))
