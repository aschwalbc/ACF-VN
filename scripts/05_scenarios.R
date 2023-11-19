## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 05: Scenarios.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame
library(progress) # Displays progress bar

# 1. Load data ==========
parms <- import(here("outputs","pts","fitpts.Rdata"))
parms <- parms %>% sample_frac(0.02)
WPP <- import(here("data","pop","WPP.Rdata"))
WUP <- import(here("data","pop","WUP.Rdata"))
GDP <- import(here("data","pop","GDP.Rdata"))
base <- import(here("data","state_base.Rdata"))
DALY <- import(here("data","pop","DALYs.xlsx"))

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  select(-omega_ini, -iota_cln_ini, -phi_cln_ini, -rho_ini, -rho_fin, -sigma_ini, -sigma_fin) %>% 
  rename(omega = omega_fin, iota_cln = iota_cln_fin, phi_cln = phi_cln_fin, gamma_minrec = gamma_mincle) 

# 2.2 Intervention parameters
prop_target <- 1 # Proportion of population targeted for ACF
prop_reached <- 1 # Proportion of population participating in ACF

# 2.3 Xpert Ultra
xpert_fp_sic <- c(val = 0.04, lo = 0.03, hi = 0.07)
xpert_fp_rec <- c(val = 0.04, lo = 0.03, hi = 0.07)
xpert_sens_min <- c(val = 0.04, lo = 0.03, hi = 0.07)
xpert_sens_sub <- c(val = 0.78, lo = 0.68, hi = 0.86)
xpert_sens_cln <- c(val = 0.91, lo = 0.86, hi = 0.95)
xpert_fp_tre <- c(val = 0.12, lo = 0.03, hi = 0.30)

# 2.4 Chest X-ray
cxr_fp_sic <- c(val = 0.085, lo = 0.069, hi = 0.134)
cxr_fp_rec <- c(val = 0.16, lo = 0.14, hi = 0.27)
cxr_sens_min <- c(val = 0.85, lo = 0.80, hi = 0.90)
cxr_sens_sub <- c(val = 0.91, lo = 0.90, hi = 0.92)
cxr_sens_cln <- c(val = 0.91, lo = 0.90, hi = 0.92)
cxr_fp_tre <- c(val = 0.16, lo = 0.14, hi = 0.27)

# 2.3 ACF Scenario A - Xpert only, Xpert+ get treatment
acfa <- data.frame(parameter = character(), val = numeric(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", val = xpert_fp_sic['val'], lo = xpert_fp_sic['lo'], hi = xpert_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", val = xpert_fp_rec['val'], lo = xpert_fp_rec['lo'], hi = xpert_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", val = xpert_sens_min['val'], lo = xpert_sens_min['lo'], hi = xpert_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", val = xpert_sens_sub['val'], lo = xpert_sens_sub['lo'], hi = xpert_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", val = xpert_sens_cln['val'], lo = xpert_sens_cln['lo'], hi = xpert_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", val = xpert_fp_tre['val'], lo = xpert_fp_tre['lo'], hi = xpert_fp_tre['hi']) %>% 
  mutate_if(is.numeric, ~ . * (prop_target*prop_reached))

# 2.4 ACF Scenario B - CXR only, CXR+ get Xpert, Xpert+ get treatment
acfb <- data.frame(parameter = character(), val = numeric(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", val = xpert_fp_sic['val']*cxr_fp_sic['val'], lo = xpert_fp_sic['lo']*cxr_fp_sic['lo'], hi = xpert_fp_sic['hi']*cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", val = xpert_fp_rec['val']*cxr_fp_rec['val'], lo = xpert_fp_rec['lo']*cxr_fp_rec['lo'], hi = xpert_fp_rec['hi']*cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", val = xpert_sens_min['val']*cxr_sens_min['val'], lo = xpert_sens_min['lo']*cxr_sens_min['lo'], hi = xpert_sens_min['hi']*cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", val = xpert_sens_sub['val']*cxr_sens_sub['val'], lo = xpert_sens_sub['lo']*cxr_sens_sub['lo'], hi = xpert_sens_sub['hi']*cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", val = xpert_sens_cln['val']*cxr_sens_cln['val'], lo = xpert_sens_cln['lo']*cxr_sens_cln['lo'], hi = xpert_sens_cln['hi']*cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", val = xpert_fp_tre['val']*cxr_fp_tre['val'], lo = xpert_fp_tre['lo']*cxr_fp_tre['lo'], hi = xpert_fp_tre['hi']*cxr_fp_tre['hi']) %>% 
  mutate_if(is.numeric, ~ . * (prop_target*prop_reached))

# 2.5 ACF Scenario C - CXR only, CXR+ get treatment
acfc <- data.frame(parameter = character(), val = numeric(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", val = cxr_fp_sic['val'], lo = cxr_fp_sic['lo'], hi = cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", val = cxr_fp_rec['val'], lo = cxr_fp_rec['lo'], hi = cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", val = cxr_sens_min['val'], lo = cxr_sens_min['lo'], hi = cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", val = cxr_sens_sub['val'], lo = cxr_sens_sub['lo'], hi = cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", val = cxr_sens_cln['val'], lo = cxr_sens_cln['lo'], hi = cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", val = cxr_fp_tre['val'], lo = cxr_fp_tre['lo'], hi = cxr_fp_tre['hi']) %>% 
  mutate_if(is.numeric, ~ . * (prop_target*prop_reached))

rm(list = ls(pattern = "^(xpert|cxr|prop)"))

# 2.6 Sociodemographic parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
rhoU <- approxfun(WUP$year, WUP$urbprop, method = 'linear', rule = 2)
rhoR <- approxfun(WUP$year, WUP$rurprop, method = 'linear', rule = 2)
sigL <- approxfun(GDP$year, GDP$lowses, method = 'linear', rule = 2)
sigH <- approxfun(GDP$year, (1-GDP$lowses), method = 'linear', rule = 2)

rm(WPP, WUP, GDP)

# 2.7 Cost data
screen_pcf_DSTB <- 104 # Passive case-finding DSTB
screen_pcf_DRTB <- 500 # Passive case-finding DRTB

screen_acf_xpert <- 8 # ACFA: Xpert
screen_acf_cxrxpert <- 1.7 # ACFB: CXR -> Xpert
screen_acf_cxr <- 1 # ACFC: CXR

rx_DSTB <- 81 # Treatment DSTB
rx_DRTB <- 973 # Treatment DRTB

# 2.8 DALYs
daly <- approxfun(DALY$year, DALY$daly_pc, method = 'linear', rule = 2)

# 2.9 MDR proportion
mdr <- 0.05 # DR-TB incidence over total incidence (2015-2021)

# 3. Models ==========
ode <- function(parms, base, interv = NULL, acf_times = NULL, end_time = 2050) {
  
  times <- seq(2020, end_time, by = 1/12)
  
  # Intervention parameters
  if(is.null(interv)){
    alpha_sic <- 0
    alpha_rec <- 0
    alpha_min <- 0
    alpha_sub <- 0
    alpha_cln <- 0
    alpha_tre <- 0
  } else { 
    alpha_sic <- runif(1, min = interv[interv$parameter == 'alpha_sic', 'lo'], max = interv[interv$parameter == 'alpha_sic', 'hi'])
    alpha_rec <- runif(1, min = interv[interv$parameter == 'alpha_rec', 'lo'], max = interv[interv$parameter == 'alpha_rec', 'hi'])
    alpha_min <- runif(1, min = interv[interv$parameter == 'alpha_min', 'lo'], max = interv[interv$parameter == 'alpha_min', 'hi'])
    alpha_sub <- runif(1, min = interv[interv$parameter == 'alpha_sub', 'lo'], max = interv[interv$parameter == 'alpha_sub', 'hi'])
    alpha_cln <- runif(1, min = interv[interv$parameter == 'alpha_cln', 'lo'], max = interv[interv$parameter == 'alpha_cln', 'hi'])
    alpha_tre <- runif(1, min = interv[interv$parameter == 'alpha_tre', 'lo'], max = interv[interv$parameter == 'alpha_tre', 'hi'])
  }

  # Active-case finding switch
  if(is.null(acf_times)){
    acf <- function(times) 0
  } else { 
    values <- ifelse(floor(times) %in% acf_times, 1, 0)
    acf <- approxfun(times, values, rule = 2)
  }
  
  # Static parameters  
  delta <- 1         # Treatment year
  theta_recinf <- 1  # REINF: Recovered -> Infected (No protection)
  theta_treinf <- 1  # REINF: Treated -> Infected (No protection)
  phi_min <- 0.0     # Treatment failure Minimal
  phi_sub <- 0.0     # Treatment failure Subclinical
  
  des <- function(times, state, parms) {
    with(as.list(c(times, state, parms)), {
      
      PopT   = (N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+
                N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)
      
      dN_RL  = nu(times)*rhoR(times)*sigL(times)*(PopT) - ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - acf(floor(times))*alpha_sic*N_RL + delta*SN_RL - mu(times)*N_RL
      dSN_RL = acf(floor(times))*alpha_sic*N_RL - delta*SN_RL - mu(times)*SN_RL
      dI_RL  = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(W_RL*theta_cleinf)+(O_RL*theta_recinf)+(P_RL*theta_treinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - acf(floor(times))*alpha_sic*I_RL - mu(times)*I_RL
      dSI_RL = acf(floor(times))*alpha_sic*I_RL - delta*SI_RL - mu(times)*SI_RL
      dW_RL  = gamma_infcle*I_RL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(W_RL*theta_cleinf)) - acf(floor(times))*alpha_sic*W_RL + delta*(SI_RL+SW_RL) - mu(times)*SI_RL
      dSW_RL = acf(floor(times))*alpha_sic*W_RL - delta*SW_RL - mu(times)*SW_RL
      dO_RL  = gamma_minrec*M_RL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_RL*theta_recinf)) - acf(floor(times))*alpha_rec*O_RL + delta*SO_RL - mu(times)*O_RL
      dSO_RL = acf(floor(times))*alpha_rec*O_RL - delta*SO_RL - mu(times)*SO_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_minrec*M_RL - lambda_minsub*M_RL - acf(floor(times))*alpha_min*M_RL + phi_min*SM_RL + tau_min*P_RL - mu(times)*M_RL
      dSM_RL = acf(floor(times))*alpha_min*M_RL - phi_min*SM_RL - delta*SM_RL - mu(times)*SM_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - acf(floor(times))*alpha_sub*S_RL + phi_sub*SS_RL + tau_sub*P_RL - mu(times)*S_RL
      dSS_RL = acf(floor(times))*alpha_sub*S_RL - phi_sub*SS_RL - delta*SS_RL - mu(times)*SS_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - iota_cln*C_RL + phi_cln*(RC_RL+SC_RL) - acf(floor(times))*alpha_cln*C_RL - omega*C_RL - mu(times)*C_RL
      dRC_RL = iota_cln*C_RL - phi_cln*RC_RL - delta*RC_RL - mu(times)*RC_RL
      dSC_RL = acf(floor(times))*alpha_cln*C_RL - phi_cln*SC_RL - delta*SC_RL - mu(times)*SC_RL
      dP_RL  = delta*(SM_RL+SS_RL+RC_RL+SC_RL+SP_RL) - tau_min*P_RL - tau_sub*P_RL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_RL*theta_recinf)) - acf(floor(times))*alpha_tre*P_RL - mu(times)*P_RL
      dSP_RL = acf(floor(times))*alpha_tre*P_RL - delta*SP_RL - mu(times)*SP_RL

      dN_RH  = nu(times)*rhoR(times)*sigH(times)*(PopT) - ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - acf(floor(times))*alpha_sic*N_RH + delta*SN_RH - mu(times)*N_RH
      dSN_RH = acf(floor(times))*alpha_sic*N_RH - delta*SN_RH - mu(times)*SN_RH
      dI_RH  = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(W_RH*theta_cleinf)+(O_RH*theta_recinf)+(P_RH*theta_treinf)) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - acf(floor(times))*alpha_sic*I_RH - mu(times)*I_RH
      dSI_RH = acf(floor(times))*alpha_sic*I_RH - delta*SI_RH - mu(times)*SI_RH
      dW_RH  = gamma_infcle*I_RH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(W_RH*theta_cleinf)) - acf(floor(times))*alpha_sic*W_RH + delta*(SI_RH+SW_RH) - mu(times)*SI_RH
      dSW_RH = acf(floor(times))*alpha_sic*W_RH - delta*SW_RH - mu(times)*SW_RH
      dO_RH  = gamma_minrec*M_RH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_RH*theta_recinf)) - acf(floor(times))*alpha_rec*O_RH + delta*SO_RH - mu(times)*O_RH
      dSO_RH = acf(floor(times))*alpha_rec*O_RH - delta*SO_RH - mu(times)*SO_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_minrec*M_RH - lambda_minsub*M_RH - acf(floor(times))*alpha_min*M_RH + phi_min*SM_RH + tau_min*P_RH - mu(times)*M_RH
      dSM_RH = acf(floor(times))*alpha_min*M_RH - phi_min*SM_RH - delta*SM_RH - mu(times)*SM_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH - acf(floor(times))*alpha_sub*S_RH + phi_sub*SS_RH + tau_sub*P_RH - mu(times)*S_RH
      dSS_RH = acf(floor(times))*alpha_sub*S_RH - phi_sub*SS_RH - delta*SS_RH - mu(times)*SS_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - iota_cln*C_RH + phi_cln*(RC_RH+SC_RH) - acf(floor(times))*alpha_cln*C_RH - omega*C_RH - mu(times)*C_RH
      dRC_RH = iota_cln*C_RH - phi_cln*RC_RH - delta*RC_RH - mu(times)*RC_RH
      dSC_RH = acf(floor(times))*alpha_cln*C_RH - phi_cln*SC_RH - delta*SC_RH - mu(times)*SC_RH
      dP_RH  = delta*(SM_RH+SS_RH+RC_RH+SC_RH+SP_RH) - tau_min*P_RH - tau_sub*P_RH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_RH*theta_recinf)) - acf(floor(times))*alpha_tre*P_RH - mu(times)*P_RH
      dSP_RH = acf(floor(times))*alpha_tre*P_RH - delta*SP_RH - mu(times)*SP_RH
      
      dN_UL  = nu(times)*rhoU(times)*sigL(times)*(PopT) - ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - acf(floor(times))*alpha_sic*N_UL + delta*SN_UL - mu(times)*N_UL
      dSN_UL = acf(floor(times))*alpha_sic*N_UL - delta*SN_UL - mu(times)*SN_UL
      dI_UL  = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(W_UL*theta_cleinf)+(O_UL*theta_recinf)+(P_UL*theta_treinf)) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - acf(floor(times))*alpha_sic*I_UL - mu(times)*I_UL
      dSI_UL = acf(floor(times))*alpha_sic*I_UL - delta*SI_UL - mu(times)*SI_UL
      dW_UL  = gamma_infcle*I_UL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(W_UL*theta_cleinf)) - acf(floor(times))*alpha_sic*W_UL + delta*(SI_UL+SW_UL) - mu(times)*SI_UL
      dSW_UL = acf(floor(times))*alpha_sic*W_UL - delta*SW_UL - mu(times)*SW_UL
      dO_UL  = gamma_minrec*M_UL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_UL*theta_recinf)) - acf(floor(times))*alpha_rec*O_UL + delta*SO_UL - mu(times)*O_UL
      dSO_UL = acf(floor(times))*alpha_rec*O_UL - delta*SO_UL - mu(times)*SO_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_minrec*M_UL - lambda_minsub*M_UL - acf(floor(times))*alpha_min*M_UL + phi_min*SM_UL + tau_min*P_UL - mu(times)*M_UL
      dSM_UL = acf(floor(times))*alpha_min*M_UL - phi_min*SM_UL - delta*SM_UL - mu(times)*SM_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL - acf(floor(times))*alpha_sub*S_UL + phi_sub*SS_UL + tau_sub*P_UL - mu(times)*S_UL
      dSS_UL = acf(floor(times))*alpha_sub*S_UL - phi_sub*SS_UL - delta*SS_UL - mu(times)*SS_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - iota_cln*C_UL + phi_cln*(RC_UL+SC_UL) - acf(floor(times))*alpha_cln*C_UL - omega*C_UL - mu(times)*C_UL
      dRC_UL = iota_cln*C_UL - phi_cln*RC_UL - delta*RC_UL - mu(times)*RC_UL
      dSC_UL = acf(floor(times))*alpha_cln*C_UL - phi_cln*SC_UL - delta*SC_UL - mu(times)*SC_UL
      dP_UL  = delta*(SM_UL+SS_UL+RC_UL+SC_UL+SP_UL) - tau_min*P_UL - tau_sub*P_UL - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_UL*theta_recinf)) - acf(floor(times))*alpha_tre*P_UL - mu(times)*P_UL
      dSP_UL = acf(floor(times))*alpha_tre*P_UL - delta*SP_UL - mu(times)*SP_UL
      
      dN_UH  = nu(times)*rhoU(times)*sigH(times)*(PopT) - ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - acf(floor(times))*alpha_sic*N_UH + delta*SN_UH - mu(times)*N_UH
      dSN_UH = acf(floor(times))*alpha_sic*N_UH - delta*SN_UH - mu(times)*SN_UH
      dI_UH  = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(W_UH*theta_cleinf)+(O_UH*theta_recinf)+(P_UH*theta_treinf)) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - acf(floor(times))*alpha_sic*I_UH - mu(times)*I_UH
      dSI_UH = acf(floor(times))*alpha_sic*I_UH - delta*SI_UH - mu(times)*SI_UH
      dW_UH  = gamma_infcle*I_UH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(W_UH*theta_cleinf)) - acf(floor(times))*alpha_sic*W_UH + delta*(SI_UH+SW_UH) - mu(times)*SI_UH
      dSW_UH = acf(floor(times))*alpha_sic*W_UH - delta*SW_UH - mu(times)*SW_UH
      dO_UH  = gamma_minrec*M_UH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_UH*theta_recinf)) - acf(floor(times))*alpha_rec*O_UH + delta*SO_UH - mu(times)*O_UH
      dSO_UH = acf(floor(times))*alpha_rec*O_UH - delta*SO_UH - mu(times)*SO_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_minrec*M_UH - lambda_minsub*M_UH - acf(floor(times))*alpha_min*M_UH + phi_min*SM_UH + tau_min*P_UH - mu(times)*M_UH
      dSM_UH = acf(floor(times))*alpha_min*M_UH - phi_min*SM_UH - delta*SM_UH - mu(times)*SM_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH - acf(floor(times))*alpha_sub*S_UH + phi_sub*SS_UH + tau_sub*P_UH - mu(times)*S_UH
      dSS_UH = acf(floor(times))*alpha_sub*S_UH - phi_sub*SS_UH - delta*SS_UH - mu(times)*SS_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - iota_cln*C_UH + phi_cln*(RC_UH+SC_UH) - acf(floor(times))*alpha_cln*C_UH - omega*C_UH - mu(times)*C_UH
      dRC_UH = iota_cln*C_UH - phi_cln*RC_UH - delta*RC_UH - mu(times)*RC_UH
      dSC_UH = acf(floor(times))*alpha_cln*C_UH - phi_cln*SC_UH - delta*SC_UH - mu(times)*SC_UH
      dP_UH  = delta*(SM_UH+SS_UH+RC_UH+SC_UH+SP_UH) - tau_min*P_UH - tau_sub*P_UH - (((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_UH*theta_recinf)) - acf(floor(times))*alpha_tre*P_UH - mu(times)*P_UH
      dSP_UH = acf(floor(times))*alpha_tre*P_UH - delta*SP_UH - mu(times)*SP_UH
      
      return(list(c(
        dN_RL, dN_RH, dN_UL, dN_UH, dSN_RL, dSN_RH, dSN_UL, dSN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dSI_RL, dSI_RH, dSI_UL, dSI_UH,
        dW_RL, dW_RH, dW_UL, dW_UH, dSW_RL, dSW_RH, dSW_UL, dSW_UH, dO_RL, dO_RH, dO_UL, dO_UH, dSO_RL, dSO_RH, dSO_UL, dSO_UH, 
        dM_RL, dM_RH, dM_UL, dM_UH, dSM_RL, dSM_RH, dSM_UL, dSM_UH, dS_RL, dS_RH, dS_UL, dS_UH, dSS_RL, dSS_RH, dSS_UL, dSS_UH, 
        dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dSC_RL, dSC_RH, dSC_UL, dSC_UH, 
        dP_RL, dP_RH, dP_UL, dP_UH, dSP_RL, dSP_RH, dSP_UL, dSP_UH),
        rMin    = ((M_RL+M_RH+M_UL+M_UH)/PopT*1e5), # Minimal TB (per 100k)
        tMin    = (M_RL+M_RH+M_UL+M_UH), # Total minimal TB
        rSub    = ((S_RL+S_RH+S_UL+S_UH)/PopT*1e5), # Subclinical TB (per 100k)
        tSub    = (S_RL+S_RH+S_UL+S_UH), # Total subclinical TB
        rCln    = ((C_RL+C_RH+C_UL+C_UH)/PopT*1e5), # Clinical TB (per 100k)
        tCln    = (C_RL+C_RH+C_UL+C_UH), # Total clinical TB
        rTBc    = ((S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/PopT*1e5), # Infectious TB (per 100k)
        tTBc    = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Total infectious TB
        rMor    = ((omega*(C_RL+C_RH+C_UL+C_UH))/PopT*1e5), # Clinical TB mortality per time (per 100k)
        tMor    = (omega*(C_RL+C_RH+C_UL+C_UH)), # Clinical TB mortality per time
        rDxs    = ((iota_cln*(C_RL+C_RH+C_UL+C_UH))/PopT*1e5), # Notifications cTB per time in adults (per 100k)
        tDxs    = (iota_cln*(C_RL+C_RH+C_UL+C_UH)), # Total notifications cTB per time in adults
        tFPos   = acf(floor(times))*((alpha_sic*(N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+W_RL+W_RH+W_UL+W_UH))+(alpha_rec*(O_RL+O_RH+O_UL+O_UH))+(alpha_tre*(P_RL+P_RH+P_UL+P_UH))), # FP diagnosed
        tTPos   = acf(floor(times))*((alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH))), # TP diagnosed
        cPCFs   = (((1-mdr)*((iota_cln)*(C_RL+C_RH+C_UL+C_UH)))*screen_pcf_DSTB),
        cPCFr   = (((mdr)*((iota_cln)*(C_RL+C_RH+C_UL+C_UH)))*screen_pcf_DRTB),
        cACF1   = acf(floor(times))*((alpha_sic*(N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+W_RL+W_RH+W_UL+W_UH))+(alpha_rec*(O_RL+O_RH+O_UL+O_UH))(alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH)))*screen_acf_xpert,
        cACF2   = acf(floor(times))*((alpha_sic*(N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+W_RL+W_RH+W_UL+W_UH))+(alpha_rec*(O_RL+O_RH+O_UL+O_UH))(alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH)))*screen_acf_cxrxpert,
        cACF3   = acf(floor(times))*((alpha_sic*(N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+W_RL+W_RH+W_UL+W_UH))+(alpha_rec*(O_RL+O_RH+O_UL+O_UH))(alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH)))*screen_acf_cxr,
        cRxPCFs = (((1-mdr)*((delta)*(RC_RL+RC_RH+RC_UL+RC_UH)))*rx_DSTB),
        cRxPCFr = (((mdr)*((delta)*(RC_RL+RC_RH+RC_UL+RC_UH)))*rx_DRTB),
        cRxTPs  = (((1-mdr)*((delta)*(SM_RL+SM_RH+SM_UL+SM_UH+SS_RL+SS_RH+SS_UL+SS_UH+SC_RL+SC_RH+SC_UL+SC_UH)))*rx_DSTB),
        cRxTPr  = (((mdr)*((delta)*(SM_RL+SM_RH+SM_UL+SM_UH+SS_RL+SS_RH+SS_UL+SS_UH+SC_RL+SC_RH+SC_UL+SC_UH)))*rx_DRTB),
        cRxFPs  = (((1-mdr)*((delta)*(SN_RL+SN_RH+SN_UL+SN_UH+SI_RL+SI_RH+SI_UL+SI_UH+SW_RL+SW_RH+SW_UL+SW_UH+SO_RL+SO_RH+SO_UL+SO_UH+SP_RL+SP_RH+SP_UL+SP_UH)))*rx_DSTB),
        cRxFPr  = (((mdr)*((delta)*(SN_RL+SN_RH+SN_UL+SN_UH+SI_RL+SI_RH+SI_UL+SI_UH+SW_RL+SW_RH+SW_UL+SW_UH+SO_RL+SO_RH+SO_UL+SO_UH+SP_RL+SP_RH+SP_UL+SP_UH)))*rx_DRTB),
        DALYs   = daly(floor(times))*(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # DALY estimates
        ARI     = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
    })
  }
  
  yini <- c(N_RL = base[, "N_RL"], N_RH = base[, "N_RH"], N_UL = base[, "N_UL"], N_UH = base[, "N_UH"],
            SN_RL = 0, SN_RH = 0, SN_UL = 0, SN_UH = 0,
            I_RL = base[, "I_RL"], I_RH = base[, "I_RH"], I_UL = base[, "I_UL"], I_UH = base[, "I_UH"],
            SI_RL = 0, SI_RH = 0, SI_UL = 0, SI_UH = 0,
            W_RL = base[, "W_RL"], W_RH = base[, "W_RH"], W_UL = base[, "W_UL"], W_UH = base[, "W_UH"],
            SW_RL = 0, SW_RH = 0, SW_UL = 0, SW_UH = 0,
            O_RL = base[, "O_RL"], O_RH = base[, "O_RH"], O_UL = base[, "O_UL"], O_UH = base[, "O_UH"],
            SO_RL = 0, SO_RH = 0, SO_UL = 0, SO_UH = 0,
            M_RL = base[, "M_RL"], M_RH = base[, "M_RH"], M_UL = base[, "M_UL"], M_UH = base[, "M_UH"],
            SM_RL = 0, SM_RH = 0, SM_UL = 0, SM_UH = 0,
            S_RL = base[, "S_RL"], S_RH = base[, "S_RH"], S_UL = base[, "S_UL"], S_UH = base[, "S_UH"],
            SS_RL = 0, SS_RH = 0, SS_UL = 0, SS_UH = 0,
            C_RL = base[, "C_RL"], C_RH = base[, "C_RH"], C_UL = base[, "C_UL"], C_UH = base[, "C_UH"],
            RC_RL = base[, "RC_RL"], RC_RH = base[, "RC_RH"], RC_UL = base[, "RC_UL"], RC_UH = base[, "RC_UH"],
            SC_RL = 0, SC_RH = 0, SC_UL = 0, SC_UH = 0,
            P_RL = base[, "P_RL"], P_RH = base[, "P_RH"], P_UL = base[, "P_UL"], P_UH = base[, "P_UH"],
            SP_RL = 0, SP_RH = 0, SP_UL = 0, SP_UH = 0)

  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}  

# 4. Outputs ==========
# 4.1 Baseline
outbase <- list()
pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))
for (i in 1:nrow(parms)) {
  curr_parms <- as.data.frame(parms[i,])
  curr_base <- as.data.frame(base[i,-1])
  
  outbase[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base))
  outbase[[i]] <- outbase[[i]] %>% mutate(type = 'base', run = i)
  pb$tick()
}
r00_outbase_df <- do.call(rbind, outbase)
r00_outbase_df <- r00_outbase_df %>% mutate(round = "00")
export(r00_outbase_df, here("outputs", "results", "r00_base.Rdata"))

# 4.2 Scenarios
acf_year <- list(
  "01" = 2025,
  "02" = seq(2025, 2026, 1),
  "03" = seq(2025, 2027, 1),
  "05" = seq(2025, 2029, 1),
  "10" = seq(2025, 2034, 1),
  "12" = seq(2025, 2036, 1))

outacfa <- list()
outacfb <- list()
outacfc <- list()

for (j in names(acf_year)) {
  rounds <- acf_year[[j]]
  print(acf_year[j])
  
  outacfa[[j]] <- list()
  outacfb[[j]] <- list()
  outacfc[[j]] <- list()
  
  pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))
  
  for (i in 1:nrow(parms)) {
    curr_parms <- as.data.frame(parms[i,])
    curr_base <- as.data.frame(base[i,-1])
    
    outacfa[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = acf_year[[j]]))
    outacfa[[j]][[i]] <- outacfa[[j]][[i]] %>% mutate(type = 'acfa', run = i, round = j)
    
    outacfb[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfb, acf_times = acf_year[[j]]))
    outacfb[[j]][[i]] <- outacfb[[j]][[i]] %>% mutate(type = 'acfb', run = i, round = j)
    
    outacfc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfc, acf_times = acf_year[[j]]))
    outacfc[[j]][[i]] <- outacfc[[j]][[i]] %>% mutate(type = 'acfc', run = i, round = j)
    
    pb$tick()
  }
}

for (k in 1:length(acf_year)) { 
  acfa_name <- paste0("r", names(acf_year[k]), "_outacfa_df")
  assign(acfa_name, do.call(rbind, outacfa[[k]]))
  export(get(acfa_name), here("outputs", "results", paste0(acfa_name, ".Rdata")))
  
  acfb_name <- paste0("r", names(acf_year[k]), "_outacfb_df")
  assign(acfb_name, do.call(rbind, outacfb[[k]]))
  export(get(acfb_name), here("outputs", "results", paste0(acfb_name, ".Rdata")))
  
  acfc_name <- paste0("r", names(acf_year[k]), "_outacfc_df")
  assign(acfc_name, do.call(rbind, outacfc[[k]]))
  export(get(acfc_name), here("outputs", "results", paste0(acfc_name, ".Rdata")))
}

rm(list = ls())
