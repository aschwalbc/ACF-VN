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
parms <- parms %>% sample_n(nrow(parms)/40)
WPP <- import(here("data","pop","WPP.Rdata"))
WUP <- import(here("data","pop","WUP.Rdata"))
GDP <- import(here("data","pop","GDP.Rdata"))
base <- import(here("data","state_base.Rdata"))

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  select(-omega_ini, -iota_cln_ini, -phi_cln_ini, -rho_ini, -rho_fin, -sigma_ini, -sigma_fin) %>% 
  rename(omega = omega_fin, iota_cln = iota_cln_fin, phi_cln = phi_cln_fin, gamma_minrec = gamma_mincle) 

# 2.2 Intervention parameters
acf_year <- seq(2025,2028,1)

prop_target <- 1 # Proportion of population targeted for ACF
prop_reached <- 1 # Proportion of population participating in ACF

# 2.3 Xpert
# 2.3.1 MTB/RIF
# xpert_fp_sic <- c(val = 0.01, lo = 0.01, hi = 0.03)
# xpert_fp_rec <- c(val = 0.01, lo = 0.01, hi = 0.03)
# xpert_sens_min <- c(val = 0.01, lo = 0.01, hi = 0.03)
# xpert_sens_sub <- c(val = 0.69, lo = 0.48, hi = 0.86)
# xpert_sens_cln <- c(val = 0.85, lo = 0.79, hi = 0.90)
# xpert_fp_tre <- c(val = 0.03, lo = 0.01, hi = 0.08)

# 2.3.2 Ultra
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

# 2.6 Sociodemographic parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
rhoU <- approxfun(WUP$year, WUP$urbprop, method = 'linear', rule = 2)
rhoR <- approxfun(WUP$year, WUP$rurprop, method = 'linear', rule = 2)
sigL <- approxfun(GDP$year, GDP$lowses, method = 'linear', rule = 2)
sigH <- approxfun(GDP$year, (1-GDP$lowses), method = 'linear', rule = 2)

# 2.7 Cost data
screen_pcf_DSTB <- 104 # Passive case-finding DSTB
screen_pcf_DRTB <- 500 # Passive case-finding DRTB

screen_acf_xpert <- 8 # ACFA: Xpert
screen_acf_cxrxpert <- 1.7 # ACFB: CXR -> Xpert
screen_acf_cxr <- 1 # ACFC: CXR

rx_DSTB <- 81 # Treatment DSTB
rx_DRTB <- 973 # Treatment DRTB

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
        Pop   = PopT, # Total population
        PRur  = (N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH), # Population rural
        PUrb  = (N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population urban
        PHig  = (N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population high SES
        PLow  = (N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL), # Population low SES
        PRL   = (N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL), # Population Rural - Low SES
        PRH   = (N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH), # Population Rural - High SES
        PUL   = (N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL), # Population Urban - Low SES
        PUH   = (N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population Urban - High SES
        Sub   = ((S_RL+S_RH+S_UL+S_UH)/PopT*1e5), # Subclinical TB (per 100k)
        SubLo = ((S_RL+S_UL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Subclinical TB in low SES (per 100k)
        SubHi = ((S_RH+S_UH)/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB in high SES (per 100k)
        SubUr = ((S_UH+S_UL)/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB in urban (per 100k)
        SubRu = ((S_RH+S_RL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Subclinical TB in rural (per 100k)
        SubRL = (S_RL/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Subclinical TB Rural - Low SES (per 100k)
        SubRH = (S_RH/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Subclinical TB Rural - High SES (per 100k)
        SubUL = (S_UL/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Subclinical TB Urban - Low SES (per 100k)
        SubUH = (S_UH/(N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB Urban - High SES (per 100k)
        Cln   = ((C_RL+C_RH+C_UL+C_UH)/PopT*1e5), # Clinical TB (per 100k)
        ClnLo = ((C_RL+C_UL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Clinical TB in low SES (per 100k)
        ClnHi = ((C_RH+C_UH)/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB in high SES (per 100k)
        ClnUr = ((C_UH+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB in urban (per 100k)
        ClnRu = ((C_RH+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Clinical TB in rural (per 100k)
        ClnRL = (C_RL/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Clinical TB Rural - Low SES (per 100k)
        ClnRH = (C_RH/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Clinical TB Rural - High SES (per 100k)
        ClnUL = (C_UL/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Clinical TB Urban - Low SES (per 100k)
        ClnUH = (C_UH/(N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB Urban - High SES (per 100k)
        TBc   = ((S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/PopT*1e5), # Infectious TB (per 100k)
        TBcLo = ((S_RL+S_UL+C_RL+C_UL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Infectious TB in low SES (per 100k)
        TBcHi = ((S_RH+S_UH+C_RH+C_UH)/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB in high SES (per 100k)
        TBcUr = ((S_UH+S_UL+C_UH+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB in urban (per 100k)
        TBcRu = ((S_RH+S_RL+C_RH+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Infectious TB in rural (per 100k)
        TBcRL = ((S_RL+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+W_RL+SW_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Infectious TB Rural - Low SES (per 100k)
        TBcRH = ((S_RH+C_RH)/(N_RH+SN_RH+I_RH+SI_RH+W_RH+SW_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Infectious TB Rural - High SES (per 100k)
        TBcUL = ((S_UL+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+W_UL+SW_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Infectious TB Urban - Low SES (per 100k)
        TBcUH = ((S_UH+C_UH)/(N_UH+SN_UH+I_UH+SI_UH+W_UH+SW_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB Urban - High SES (per 100k)
        rMor  = (omega*(C_RL+C_RH+C_UL+C_UH))/PopT*1e5, # Clinical TB mortality per time (per 100k)
        tMor  = (omega*(C_RL+C_RH+C_UL+C_UH)), # Clinical TB mortality per time
        Dxs   = (iota_cln*(C_RL+C_RH+C_UL+C_UH))/PopT*1e5, # Notifications cTB per time in adults (per 100k)
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
        NumSC = acf(floor(times))*((alpha_sic*(SN_RL+SN_RH+SN_UL+SN_UH+SI_RL+SI_RH+SI_UL+SI_UH+SW_RL+SW_RH+SW_UL+SW_UH))+(alpha_rec*(SO_RL+SO_RH+SO_UL+SO_UH))+(alpha_min*(SM_RL+SM_RH+SM_UL+SM_UH))+(alpha_sub*(SS_RL+SS_RH+SS_UL+SS_UH))+(alpha_cln*(SC_RL+SC_RH+SC_UL+SC_UH))+(alpha_tre*(SP_RL+SP_RH+SP_UL+SP_UH))),
        FPnds = acf(floor(times))*((alpha_sic*(SN_RL+SN_RH+SN_UL+SN_UH+SI_RL+SI_RH+SI_UL+SI_UH+SW_RL+SW_RH+SW_UL+SW_UH))+(alpha_rec*(SO_RL+SO_RH+SO_UL+SO_UH))+(alpha_tre*(SP_RL+SP_RH+SP_UL+SP_UH))),
        TPdis = acf(floor(times))*((alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH))),
        SCmin = acf(floor(times))*alpha_min*(M_RL+M_RH+M_UL+M_UH),
        SCsub = acf(floor(times))*alpha_sub*(S_RL+S_RH+S_UL+S_UH),
        SCcln = acf(floor(times))*alpha_cln*(C_RL+C_RH+C_UL+C_UH),
        DXcln = iota_cln*(C_RL+C_RH+C_UL+C_UH),
        FNdis = acf(floor(times))*(((1-alpha_min)*(M_RL+M_RH+M_UL+M_UH))+((1-alpha_sub)*(S_RL+S_RH+S_UL+S_UH))+((1-alpha_cln)*(C_RL+C_RH+C_UL+C_UH))),
        FNmin = acf(floor(times))*((1-alpha_min)*(M_RL+M_RH+M_UL+M_UH)),
        FNsub = acf(floor(times))*((1-alpha_sub)*(S_RL+S_RH+S_UL+S_UH)),
        FNcln = acf(floor(times))*((1-alpha_cln)*(C_RL+C_RH+C_UL+C_UH)),
        ARI   = ((beta/PopT)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
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
outbase <- list()
outacfa <- list()
outacfb <- list()
outacfc <- list()

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))

for (i in 1:nrow(parms)) {
  curr_parms <- as.data.frame(parms[i,])
  curr_base <- as.data.frame(base[i,-1])
  
  outbase[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base))
  outbase[[i]] <- outbase[[i]] %>% mutate(type = 'base', run = i)
  
  outacfa[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = acf_year))
  outacfa[[i]] <- outacfa[[i]] %>% mutate(type = 'acfa', run = i)

  outacfb[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfb, acf_times = acf_year))
  outacfb[[i]] <- outacfb[[i]] %>% mutate(type = 'acfb', run = i)

  outacfc[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfc, acf_times = acf_year))
  outacfc[[i]] <- outacfc[[i]] %>% mutate(type = 'acfc', run = i)
  
  pb$tick()
}

outbase_df <- do.call(rbind, outbase)
# export(outbase_df, here("outputs","outbase_df.Rdata"))
outacfa_df <- do.call(rbind, outacfa)
# export(outacfa_df, here("outputs","outacfa_df.Rdata"))
outacfb_df <- do.call(rbind, outacfb)
# export(outacfb_df, here("outputs","outacfb_df.Rdata"))
outacfc_df <- do.call(rbind, outacfc)
# export(outacfc_df, here("outputs","outacfc_df.Rdata"))

# 5. Data curation ==========
outs <- rbind(outbase_df, outacfa_df, outacfb_df, outacfc_df) %>% 
  arrange(type, run, time) %>% 
  group_by(type, run) %>%
  mutate(ACF = SCmin+SCsub+SCcln, # Sum of all TB disease screened (Min+Sub+Cln)
         AllTB = SCmin+SCsub+SCcln+DXcln, # Sum of all TB disease diagnoses (ACF+BAU)
         AllTx = SCmin+SCsub+SCcln+DXcln+FPnds, # Sum of all diagnoses (ACF+BAU+FP)
         cumMor = cumsum(tMor), # Cumulative TB mortality
         cumFPnds = cumsum(FPnds), # Cumulative FP diagnoses
         cumFNdis = cumsum(FNdis), # Cumulative FN (diagnoses missed)
         cumNumSC = cumsum(NumSC), # Cumulative screening (ACF+FP)
         cumDXcln = cumsum(DXcln), # Cumulative BAU diagnoses 
         pFP = ifelse(FPnds == 0, 0, FPnds/(TPdis+FPnds))) %>% # FP treated (% over confirmed)
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(dMor = tMor - tMor[type == 'base'], # TB mortality difference to BAU
         dcumMor = cumMor - cumMor[type == 'base'], # Cumulative TB mortality difference to BAU
         pTBc = TBc/TBc[type == 'base'], # Proportional reduction TB prevalence to BAU
         pMor = rMor/rMor[type == 'base'], # Proportional reduction TB mortality to BAU
         dAllTB = AllTB - AllTB[type == 'base'], # All TB diagnoses difference to BAU
         dAllTx = AllTx - AllTx[type == 'base'], # All diagnoses (TB+FP) difference to BAU
         dBAU = cumDXcln - cumDXcln[type == 'base']) %>% # BAU TB diagnoses difference to BAU
  ungroup() %>% 
  group_by(type, run) %>% 
  mutate(cumAllTB = cumsum(dAllTB), cumAllTx = cumsum(dAllTx)) %>% 
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(NNS = 1/pTBc, # Number needed to screen
         NNT = 1/pMor) %>% # Number needed to treat to avert a TB death
  ungroup() %>% 
  group_by(type, run, time) %>% 
  mutate(pPRur = PRur/Pop, # Proportion rural
         pPUrb = PUrb/Pop, # Proportion urban
         pPHig = PHig/Pop, # Proportion high SES
         pPLow = PLow/Pop, # Proportion low SES
         pPRL = PRL/Pop, # Proportion rural - low SES
         pPRH = PRH/Pop, # Proportion rural - high SES
         pPUL = PUL/Pop, # Proportion urban - low SES
         pPUH = PUH/Pop) %>% # Proportion urban - high SES
  pivot_longer(cols = -c(time, type, run), names_to = "var", values_to = "values") %>% 
  group_by(time, type, var) %>% 
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE)) %>% 
  mutate(fill = ifelse(val < 0, "under", "over"))

# export(outs, here("outputs","runs.Rdata"))
# outs <- import(here("outputs","runs.Rdata"))

# 6. Plots ==========
prev_targets <- c(100, 50, 20)
dis_state <- factor(c("Clinical","Subclinical","Minimal"))
inf_dis <- c("Clinical", "Subclinical")
scenarios <- c("01: Xpert", "02: CXR->Xpert", "03: CXR", "Baseline")
treatments <- c("BAU: Clinical", "ACF: Clinical", "ACF: Minimal", "ACF: Subclinical")
types <- c("ACF", "BAU")
type_label <- labeller(type = c("acfa" = "01: Xpert", "acfb" = "02: CXR->Xpert", "acfc" = "03: CXR", "base" = "Baseline"))
dem_urbrur <- c("Rural", "Urban")
dem_ses <- c("High","Low")

# TB prevalence per scenarios 
#tiff(here("plots", "TBprev.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "TBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "TBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,250,25), limits = c(0,250), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_hline(yintercept = prev_targets, linetype = "dashed", color = "gray") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
#dev.off()

# Proportion reduction of TB prevalence
tiff(here("plots", "TBreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pTBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pTBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1.01,0.1), limits = c(0,1.01), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(min(acf_year),2050,5), limits = c(min(acf_year),2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(min(acf_year),2050)) + 
  labs(x = "Year", y = "Proportion reduction of TB prevalence") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_rect(aes(xmin = 2029 - 0.5, xmax = 2029 + 0.5, ymin = 0.210, ymax = 0.498), colour = "#1D2D5F", linetype = 'dashed', alpha = 0.2) + 
  geom_text(aes(x = 2028.5, y = 0.25, label = "ACT3 reduction", fontface = 'bold'), angle = 90, vjust = -0.5, hjust = 0, size = 3, color = "#1D2D5F") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Annual risk of infection 
tiff(here("plots", "ARI.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "ARI"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "ARI"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Annual risk of infection") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB mortality per scenarios 
tiff(here("plots", "TBmort.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "rMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "rMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB mortality rate (per 100K)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Total TB mortality 
tiff(here("plots", "totTBmort.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "tMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "tMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,15000,2500), expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Total TB mortality") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Proportion reduction of TB mortality
tiff(here("plots", "TBmorreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1.01,0.1), limits = c(0,1.01), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(min(acf_year),2050,5), limits = c(min(acf_year),2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(min(acf_year),2050)) + 
  labs(x = "Year", y = "Proportion reduction of TB mortality") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative TB deaths averted
tiff(here("plots", "cumTBdeathsavert.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dcumMor"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dcumMor"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative TB deaths averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB deaths averted
tiff(here("plots", "TBdeathsavert.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dMor"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dMor"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB deaths averted per year (compared to BAU)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB treatments
tiff(here("plots", "TBtreatment.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, scales = 'free_y', labeller = type_label) +
  geom_line(data = filter(outs, var %in% c("SCmin","SCsub","SCcln","DXcln")), aes(x = time, y = val, colour = var)) +
  scale_colour_manual(values = c("#000000","#F58B65","#4DC4CB","#FDC75D"), name = NULL, labels = treatments) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments (BAU vs ACF)
tiff(here("plots", "TBtreatBAUACF.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, scales = 'free_y', labeller = type_label) +
  geom_line(data = filter(outs, var %in% c("DXcln", "ACF")), aes(x = time, y = val, colour = var)) +
  scale_colour_manual(values = c("#CE2931","#000000"), name = NULL, labels = types) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# All TB treatments
tiff(here("plots", "TBtreatALLTB.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "AllTB"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "AllTB"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K"), limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050), ylim = c(0, 9e5)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# All treatments (TB+FP)
tiff(here("plots", "TBtreatALLTx.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "AllTx"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "AllTx"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K"), limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050), ylim = c(0, 9e5)) + 
  labs(x = "Year", y = "Number of treatments (including FP)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments averted
tiff(here("plots", "TBtreatmentaverted.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, labeller = type_label) +
  geom_line(data = filter(outs, var == "dAllTB" & type != 'base'), aes(x = time, y = val)) +
  geom_area(data = filter(outs, var == "dAllTB" & type != 'base'), aes(x = time, y = val, fill = fill), alpha = 0.5) +
  scale_fill_manual(values = c("over" = "#FF531A", "under" = "#1AC6FF"), name = NULL) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "none", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments averted (considering FP)
tiff(here("plots", "TBtreatmentaverted_withFP.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, labeller = type_label) +
  geom_line(data = filter(outs, var == "dAllTx" & type != 'base'), aes(x = time, y = val)) +
  geom_area(data = filter(outs, var == "dAllTx" & type != 'base'), aes(x = time, y = val, fill = fill), alpha = 0.5) +
  scale_fill_manual(values = c("over" = "#FF531A", "under" = "#1AC6FF"), name = NULL) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "none", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Cumulative TB treatments averted
tiff(here("plots", "cumTBtreatmentaverted_line.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dBAU"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dBAU"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Proportion prevalence (urban vs rural)
tiff(here("plots", "Prop_urbvrur.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "URt"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "URt"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB urban vs rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion prevalence (low vs high SES)
tiff(here("plots", "Prop_lovhi.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "HLt"), aes(x = time, y = (1-val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "HLt"), aes(x = time, ymin = (1-lo), ymax = (1-hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB low vs high") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Urban
tiff(here("plots", "Prop_UH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPUH"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPUH"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Rural
tiff(here("plots", "Prop_RH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPRH"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPRH"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Urban
tiff(here("plots", "Prop_UL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPUL"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPUL"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion Low SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Rural
tiff(here("plots", "Prop_RL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPRL"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPRL"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
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
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
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
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion disease state") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# False positives
tiff(here("plots", "Falsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "False positive diagnoses") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative false positives
tiff(here("plots", "CumFalsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumFPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumFPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative false positive diagnoses") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False negatives
tiff(here("plots", "Falseneg.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FNdis"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FNdis"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "False negative diagnoses (Missed diagnoses)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative false negatives
tiff(here("plots", "CumFalseneg.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumFNdis"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumFNdis"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative false negative diagnoses (Missed diagnoses)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False positive treated
tiff(here("plots", "FPtreated.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pFP"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pFP"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion FP (% of total screened)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Number needed to screen
tiff(here("plots", "NNS.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "NNS"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "NNS"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number needed to screen") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative number screened
tiff(here("plots", "NumScreen.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumNumSC"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumNumSC"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative number screened/treated") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()
