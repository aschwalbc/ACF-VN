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
library(progress) # Displays progress bar

# 1. Load data ==========
parms <- import(here("outputs","plaus_pts.Rdata"))
WPP <- import(here("data","pop","WPP.Rdata"))
WUP <- import(here("data","pop","WUP.Rdata"))
GDP <- import(here("data","pop","GDP.Rdata"))
base <- import(here("data","state_base.Rdata"))

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  select(-omega_ini, -iota_cln_ini, -phi_cln_ini, -rho_ini, -rho_fin, -sigma_ini, -sigma_fin) %>% 
  rename(omega = omega_fin, iota_cln = iota_cln_fin, phi_cln = phi_cln_fin) 

# 2.2 Intervention parameters
acf_year <- seq(2025,2027,1)

prop_target <- 1 # Proportion of population targeted for ACF
prop_reached <- 1 # Proportion of population participating in ACF

# 2.3 Xpert
xpert_fp_susinf <- 0.01 # (0.01-0.03)
xpert_fp_rec <- 0.01 # (0.01-0.03)
xpert_sens_min <- 0.01 # (0.01-0.03)
xpert_sens_sub <- 0.69 # (0.48-0.86)
xpert_sens_cln <- 0.85 # (0.79-0.90)
xpert_fp_tre <- 0.03 # (0.01-0.08)

# 2.4 Chest X-ray
cxr_fp_susinf <- 0.09 # (0.07-0.13)
cxr_fp_rec <- 0.16 # (0.14-0.27)
cxr_sens_min <- 0.85 # (0.80-0.90)
cxr_sens_sub <- 0.91 # (0.90-0.92)
cxr_sens_cln <- 0.91 # (0.90-0.92)
cxr_fp_tre <- 0.16 # (0.14-0.27)

# 2.3 ACF Scenario A - Xpert only, Xpert+ get treatment
acfa <- data.frame(parameter = character(), val = numeric()) %>% 
  add_row(parameter = "alpha_susinf", val = prop_target*prop_reached*xpert_fp_susinf) %>% 
  add_row(parameter = "alpha_rec", val = prop_target*prop_reached*xpert_fp_rec) %>% 
  add_row(parameter = "alpha_min", val = prop_target*prop_reached*xpert_sens_min) %>% 
  add_row(parameter = "alpha_sub", val = prop_target*prop_reached*xpert_sens_sub) %>%
  add_row(parameter = "alpha_cln", val = prop_target*prop_reached*xpert_sens_cln) %>% 
  add_row(parameter = "alpha_tre", val = prop_target*prop_reached*xpert_fp_tre)

# 2.4 ACF Scenario B - CXR only, CXR+ get Xpert, Xpert+ get treatment
acfb <- data.frame(parameter = character(), val = numeric()) %>% 
  add_row(parameter = "alpha_susinf", val = prop_target*prop_reached*xpert_fp_susinf*cxr_fp_susinf) %>% 
  add_row(parameter = "alpha_rec", val = prop_target*prop_reached*xpert_fp_rec*cxr_fp_rec) %>% 
  add_row(parameter = "alpha_min", val = prop_target*prop_reached*xpert_sens_min*cxr_sens_min) %>% 
  add_row(parameter = "alpha_sub", val = prop_target*prop_reached*xpert_sens_sub*cxr_sens_sub) %>%
  add_row(parameter = "alpha_cln", val = prop_target*prop_reached*xpert_sens_cln*cxr_sens_cln) %>% 
  add_row(parameter = "alpha_tre", val = prop_target*prop_reached*xpert_fp_tre*cxr_fp_tre)

# 2.5 ACF Scenario C - CXR only, CXR+ get treatment
acfc <- data.frame(parameter = character(), val = numeric()) %>% 
  add_row(parameter = "alpha_susinf", val = prop_target*prop_reached*cxr_fp_susinf) %>% 
  add_row(parameter = "alpha_rec", val = prop_target*prop_reached*cxr_fp_rec) %>% 
  add_row(parameter = "alpha_min", val = prop_target*prop_reached*cxr_sens_min) %>% 
  add_row(parameter = "alpha_sub", val = prop_target*prop_reached*cxr_sens_sub) %>%
  add_row(parameter = "alpha_cln", val = prop_target*prop_reached*cxr_sens_cln) %>% 
  add_row(parameter = "alpha_tre", val = prop_target*prop_reached*cxr_fp_tre)

# 2.6 Sociodemographic parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
chi <- approxfun(WPP$year, WPP$pop, method = 'linear', rule = 2)
rhoU <- approxfun(WUP$year, WUP$urbprop, method = 'linear', rule = 2)
rhoR <- approxfun(WUP$year, WUP$rurprop, method = 'linear', rule = 2)
sigL <- approxfun(GDP$year, GDP$lowses, method = 'linear', rule = 2)
sigH <- approxfun(GDP$year, (1-GDP$lowses), method = 'linear', rule = 2)

# 3. Models ==========
ode <- function(parms, base, interv = NULL, acf_times = NULL, end_time = 2050) {
  
  times <- seq(2020, end_time, by = 1)
  
  # Intervention parameters
  if(is.null(interv)){
    alpha_susinf <- 0
    alpha_rec <- 0
    alpha_min <- 0
    alpha_sub <- 0
    alpha_cln <- 0
    alpha_tre <- 0
  } else { 
    alpha_susinf <- interv[interv$parameter == 'alpha_susinf', 'val']
    alpha_rec <- interv[interv$parameter == 'alpha_rec', 'val']
    alpha_min <- interv[interv$parameter == 'alpha_min', 'val']
    alpha_sub <- interv[interv$parameter == 'alpha_sub', 'val']
    alpha_cln <- interv[interv$parameter == 'alpha_cln', 'val']
    alpha_tre <- interv[interv$parameter == 'alpha_tre', 'val']
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
  phi_min <- 0.0     # Treatment failure Minimal
  phi_sub <- 0.0     # Treatment failure Subclinical
  
  des <- function(times, state, parms) {
    with(as.list(c(times, state, parms)), {
      
      dN_RL  = nu(times)*rhoR(times)*sigL(times)*(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+
               N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH) - 
               ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - acf(floor(times))*alpha_susinf*N_RL - mu(times)*N_RL
      dSN_RL = acf(floor(times))*alpha_susinf*N_RL - delta*SN_RL - mu(times)*SN_RL
      dI_RL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf)+(P_RL*theta_recinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - acf(floor(times))*alpha_susinf*I_RL - mu(times)*I_RL
      dSI_RL = acf(floor(times))*alpha_susinf*I_RL - delta*SI_RL - mu(times)*SI_RL
      dO_RL  = gamma_infcle*I_RL + gamma_mincle*M_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf - acf(floor(times))*alpha_rec*O_RL + delta*(SN_RL+SI_RL+SO_RL) - mu(times)*O_RL
      dSO_RL = acf(floor(times))*alpha_rec*O_RL - delta*SO_RL - mu(times)*SO_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - acf(floor(times))*alpha_min*M_RL + phi_min*SM_RL + tau_min*P_RL - mu(times)*M_RL
      dSM_RL = acf(floor(times))*alpha_min*M_RL - phi_min*SM_RL - delta*SM_RL - mu(times)*SM_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - acf(floor(times))*alpha_sub*S_RL + phi_sub*SS_RL + tau_sub*P_RL - mu(times)*S_RL
      dSS_RL = acf(floor(times))*alpha_sub*S_RL - phi_sub*SS_RL - delta*SS_RL - mu(times)*SS_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - iota_cln*C_RL + phi_cln*RC_RL - acf(floor(times))*alpha_cln*C_RL + phi_cln*SC_RL - omega*C_RL - mu(times)*C_RL
      dRC_RL = iota_cln*C_RL - phi_cln*RC_RL - delta*RC_RL - mu(times)*RC_RL
      dSC_RL = acf(floor(times))*alpha_cln*C_RL - phi_cln*SC_RL - delta*SC_RL - mu(times)*SC_RL
      dP_RL  = delta*(SM_RL+SS_RL+RC_RL+SC_RL+SP_RL) - tau_min*P_RL - tau_sub*P_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - acf(floor(times))*alpha_tre*P_RL - mu(times)*P_RL
      dSP_RL = acf(floor(times))*alpha_tre*P_RL - delta*SP_RL - mu(times)*SP_RL

      dN_RH  = nu(times)*rhoR(times)*sigL(times)*(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+
               N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH) - 
               ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - acf(floor(times))*alpha_susinf*N_RH - mu(times)*N_RH
      dSN_RH = acf(floor(times))*alpha_susinf*N_RH - delta*SN_RH - mu(times)*SN_RH
      dI_RH  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(O_RH*theta_cleinf)+(P_RH*theta_recinf)) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - acf(floor(times))*alpha_susinf*I_RH - mu(times)*I_RH
      dSI_RH = acf(floor(times))*alpha_susinf*I_RH - delta*SI_RH - mu(times)*SI_RH
      dO_RH  = gamma_infcle*I_RH + gamma_mincle*M_RH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_cleinf - acf(floor(times))*alpha_rec*O_RH + delta*(SN_RH+SI_RH+SO_RH) - mu(times)*O_RH
      dSO_RH = acf(floor(times))*alpha_rec*O_RH - delta*SO_RH - mu(times)*SO_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH - acf(floor(times))*alpha_min*M_RH + phi_min*SM_RH + tau_min*P_RH - mu(times)*M_RH
      dSM_RH = acf(floor(times))*alpha_min*M_RH - phi_min*SM_RH - delta*SM_RH - mu(times)*SM_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH - acf(floor(times))*alpha_sub*S_RH + phi_sub*SS_RH + tau_sub*P_RH - mu(times)*S_RH
      dSS_RH = acf(floor(times))*alpha_sub*S_RH - phi_sub*SS_RH - delta*SS_RH - mu(times)*SS_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - iota_cln*C_RH + phi_cln*RC_RH - acf(floor(times))*alpha_cln*C_RH + phi_cln*SC_RH - omega*C_RH - mu(times)*C_RH
      dRC_RH = iota_cln*C_RH - phi_cln*RC_RH - delta*RC_RH - mu(times)*RC_RH
      dSC_RH = acf(floor(times))*alpha_cln*C_RH - phi_cln*SC_RH - delta*SC_RH - mu(times)*SC_RH
      dP_RH  = delta*(SM_RH+SS_RH+RC_RH+SC_RH+SP_RH) - tau_min*P_RH - tau_sub*P_RH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_recinf - acf(floor(times))*alpha_tre*P_RH - mu(times)*P_RH
      dSP_RH = acf(floor(times))*alpha_tre*P_RH - delta*SP_RH - mu(times)*SP_RH
      
      dN_UL  = nu(times)*rhoR(times)*sigL(times)*(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+
               N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH) - 
               ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - acf(floor(times))*alpha_susinf*N_UL - mu(times)*N_UL
      dSN_UL = acf(floor(times))*alpha_susinf*N_UL - delta*SN_UL - mu(times)*SN_UL
      dI_UL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(O_UL*theta_cleinf)+(P_UL*theta_recinf)) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - acf(floor(times))*alpha_susinf*I_UL - mu(times)*I_UL
      dSI_UL = acf(floor(times))*alpha_susinf*I_UL - delta*SI_UL - mu(times)*SI_UL
      dO_UL  = gamma_infcle*I_UL + gamma_mincle*M_UL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_cleinf - acf(floor(times))*alpha_rec*O_UL + delta*(SN_UL+SI_UL+SO_UL) - mu(times)*O_UL
      dSO_UL = acf(floor(times))*alpha_rec*O_UL - delta*SO_UL - mu(times)*SO_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL - acf(floor(times))*alpha_min*M_UL + phi_min*SM_UL + tau_min*P_UL - mu(times)*M_UL
      dSM_UL = acf(floor(times))*alpha_min*M_UL - phi_min*SM_UL - delta*SM_UL - mu(times)*SM_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL - acf(floor(times))*alpha_sub*S_UL + phi_sub*SS_UL + tau_sub*P_UL - mu(times)*S_UL
      dSS_UL = acf(floor(times))*alpha_sub*S_UL - phi_sub*SS_UL - delta*SS_UL - mu(times)*SS_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - iota_cln*C_UL + phi_cln*RC_UL - acf(floor(times))*alpha_cln*C_UL + phi_cln*SC_UL - omega*C_UL - mu(times)*C_UL
      dRC_UL = iota_cln*C_UL - phi_cln*RC_UL - delta*RC_UL - mu(times)*RC_UL
      dSC_UL = acf(floor(times))*alpha_cln*C_UL - phi_cln*SC_UL - delta*SC_UL - mu(times)*SC_UL
      dP_UL  = delta*(SM_UL+SS_UL+RC_UL+SC_UL+SP_UL) - tau_min*P_UL - tau_sub*P_UL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_recinf - acf(floor(times))*alpha_tre*P_UL - mu(times)*P_UL
      dSP_UL = acf(floor(times))*alpha_tre*P_UL - delta*SP_UL - mu(times)*SP_UL
      
      dN_UH  = nu(times)*rhoR(times)*sigL(times)*(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+
               N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH) - 
               ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - acf(floor(times))*alpha_susinf*N_UH - mu(times)*N_UH
      dSN_UH = acf(floor(times))*alpha_susinf*N_UH - delta*SN_UH - mu(times)*SN_UH
      dI_UH  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(O_UH*theta_cleinf)+(P_UH*theta_recinf)) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - acf(floor(times))*alpha_susinf*I_UH - mu(times)*I_UH
      dSI_UH = acf(floor(times))*alpha_susinf*I_UH - delta*SI_UH - mu(times)*SI_UH
      dO_UH  = gamma_infcle*I_UH + gamma_mincle*M_UH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_cleinf - acf(floor(times))*alpha_rec*O_UH + delta*(SN_UH+SI_UH+SO_UH) - mu(times)*O_UH
      dSO_UH = acf(floor(times))*alpha_rec*O_UH - delta*SO_UH - mu(times)*SO_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH - acf(floor(times))*alpha_min*M_UH + phi_min*SM_UH + tau_min*P_UH - mu(times)*M_UH
      dSM_UH = acf(floor(times))*alpha_min*M_UH - phi_min*SM_UH - delta*SM_UH - mu(times)*SM_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH - acf(floor(times))*alpha_sub*S_UH + phi_sub*SS_UH + tau_sub*P_UH - mu(times)*S_UH
      dSS_UH = acf(floor(times))*alpha_sub*S_UH - phi_sub*SS_UH - delta*SS_UH - mu(times)*SS_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - iota_cln*C_UH + phi_cln*RC_UH - acf(floor(times))*alpha_cln*C_UH + phi_cln*SC_UH - omega*C_UH - mu(times)*C_UH
      dRC_UH = iota_cln*C_UH - phi_cln*RC_UH - delta*RC_UH - mu(times)*RC_UH
      dSC_UH = acf(floor(times))*alpha_cln*C_UH - phi_cln*SC_UH - delta*SC_UH - mu(times)*SC_UH
      dP_UH  = delta*(SM_UH+SS_UH+RC_UH+SC_UH+SP_UH) - tau_min*P_UH - tau_sub*P_UH - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_recinf - acf(floor(times))*alpha_tre*P_UH - mu(times)*P_UH
      dSP_UH = acf(floor(times))*alpha_tre*P_UH - delta*SP_UH - mu(times)*SP_UH
      
      return(list(c(
        dN_RL, dN_RH, dN_UL, dN_UH, dSN_RL, dSN_RH, dSN_UL, dSN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dSI_RL, dSI_RH, dSI_UL, dSI_UH,
        dO_RL, dO_RH, dO_UL, dO_UH, dSO_RL, dSO_RH, dSO_UL, dSO_UH, dM_RL, dM_RH, dM_UL, dM_UH, dSM_RL, dSM_RH, dSM_UL, dSM_UH, 
        dS_RL, dS_RH, dS_UL, dS_UH, dSS_RL, dSS_RH, dSS_UL, dSS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH,
        dSC_RL, dSC_RH, dSC_UL, dSC_UH, dP_RL, dP_RH, dP_UL, dP_UH, dSP_RL, dSP_RH, dSP_UL, dSP_UH),
        Pop   = (N_RL+N_RH+N_UL+N_UH+SN_RL+SN_RH+SN_UL+SN_UH+I_RL+I_RH+I_UL+I_UH+SI_RL+SI_RH+SI_UL+SI_UH+
                   O_RL+O_RH+O_UL+O_UH+SO_RL+SO_RH+SO_UL+SO_UH+M_RL+M_RH+M_UL+M_UH+SM_RL+SM_RH+SM_UL+SM_UH+
                   S_RL+S_RH+S_UL+S_UH+SS_RL+SS_RH+SS_UL+SS_UH+C_RL+C_RH+C_UL+C_UH+RC_RL+RC_RH+RC_UL+RC_UH+
                   SC_RL+SC_RH+SC_UL+SC_UH+P_RL+P_RH+P_UL+P_UH+SP_RL+SP_RH+SP_UL+SP_UH), # Total population
        PRur  = (N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH), # Population rural
        PUrb  = (N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population urban
        PHig  = (N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population high SES
        PLow  = (N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL), # Population low SES
        PRL   = (N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL), # Population Rural - Low SES
        PRH   = (N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH), # Population Rural - High SES
        PUL   = (N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL), # Population Urban - Low SES
        PUH   = (N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH), # Population Urban - High SES
        Sub   = ((S_RL+S_RH+S_UL+S_UH)/chi(times)*1e5), # Subclinical TB (per 100k)
        SubLo = ((S_RL+S_UL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Subclinical TB in low SES (per 100k)
        SubHi = ((S_RH+S_UH)/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB in high SES (per 100k)
        SubUr = ((S_UH+S_UL)/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB in urban (per 100k)
        SubRu = ((S_RH+S_RL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Subclinical TB in urban (per 100k)
        SubRL = (S_RL/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Subclinical TB Rural - Low SES (per 100k)
        SubRH = (S_RH/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Subclinical TB Rural - High SES (per 100k)
        SubUL = (S_UL/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Subclinical TB Urban - Low SES (per 100k)
        SubUH = (S_UH/(N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Subclinical TB Urban - High SES (per 100k)
        Cln   = ((C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5), # Clinical TB (per 100k)
        ClnLo = ((C_RL+C_UL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Clinical TB in low SES (per 100k)
        ClnHi = ((C_RH+C_UH)/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB in high SES (per 100k)
        ClnUr = ((C_UH+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB in urban (per 100k)
        ClnRu = ((C_RH+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Clinical TB in urban (per 100k)
        ClnRL = (C_RL/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Clinical TB Rural - Low SES (per 100k)
        ClnRH = (C_RH/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Clinical TB Rural - High SES (per 100k)
        ClnUL = (C_UL/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Clinical TB Urban - Low SES (per 100k)
        ClnUH = (C_UH/(N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Clinical TB Urban - High SES (per 100k)
        TBc   = ((S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH)/chi(times)*1e5), # Infectious TB (per 100k)
        TBcLo = ((S_RL+S_UL+C_RL+C_UL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Infectious TB in low SES (per 100k)
        TBcHi = ((S_RH+S_UH+C_RH+C_UH)/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB in high SES (per 100k)
        TBcUr = ((S_UH+S_UL+C_UH+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL+N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB in urban (per 100k)
        TBcRu = ((S_RH+S_RL+C_RH+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL+N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Infectious TB in urban (per 100k)
        TBcRL = ((S_RL+C_RL)/(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+SM_RL+S_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL)*1e5), # Infectious TB Rural - Low SES (per 100k)
        TBcRH = ((S_RH+C_RH)/(N_RH+SN_RH+I_RH+SI_RH+O_RH+SO_RH+M_RH+SM_RH+S_RH+SS_RH+C_RH+RC_RH+SC_RH+P_RH+SP_RH)*1e5), # Infectious TB Rural - High SES (per 100k)
        TBcUL = ((S_UL+C_UL)/(N_UL+SN_UL+I_UL+SI_UL+O_UL+SO_UL+M_UL+SM_UL+S_UL+SS_UL+C_UL+RC_UL+SC_UL+P_UL+SP_UL)*1e5), # Infectious TB Urban - Low SES (per 100k)
        TBcUH = ((S_UH+C_UH)/(N_UH+SN_UH+I_UH+SI_UH+O_UH+SO_UH+M_UH+SM_UH+S_UH+SS_UH+C_UH+RC_UH+SC_UH+P_UH+SP_UH)*1e5), # Infectious TB Urban - High SES (per 100k)
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
        NumSC = ((acf(floor(times)))*((alpha_susinf*(N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH))+(alpha_rec*(O_RL+O_RH+O_UL+O_UH))+(alpha_min*(M_RL+M_RH+M_UL+M_UH))+(alpha_sub*(S_RL+S_RH+S_UL+S_UH))+(alpha_cln*(C_RL+C_RH+C_UL+C_UH))+(alpha_tre*(P_RL+P_RH+P_UL+P_UH)))),
        FPnds = (SN_RL+SI_RL+SO_RL+SP_RL+SN_RH+SI_RH+SO_RH+SP_RH+SN_UH+SI_UH+SO_UH+SP_UH+SN_UL+SI_UL+SO_UL+SP_UL),
        ARI   = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
    })
  }
  
  yini <- c(N_RL = base[, "N_RL"], N_RH = base[, "N_RH"], N_UL = base[, "N_UL"], N_UH = base[, "N_UH"],
            SN_RL = 0, SN_RH = 0, SN_UL = 0, SN_UH = 0,
            I_RL = base[, "I_RL"], I_RH = base[, "I_RH"], I_UL = base[, "I_UL"], I_UH = base[, "I_UH"],
            SI_RL = 0, SI_RH = 0, SI_UL = 0, SI_UH = 0,
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
outacfa_df <- do.call(rbind, outacfa)
outacfb_df <- do.call(rbind, outacfb)
outacfc_df <- do.call(rbind, outacfc)

# 5. Data curation ==========
outs <- rbind(outbase_df, outacfa_df, outacfb_df, outacfc_df) %>% 
  arrange(time) %>% 
  group_by(type) %>% 
  mutate(cumMor = cumsum(tMor), cumFPnds = cumsum(FPnds), cumNumSC = cumsum(NumSC)) %>% 
  group_by(time) %>% 
  mutate(dMor = tMor - tMor[type == 'base'], dcumMor = cumMor - cumMor[type == 'base'],
         pTBc = TBc/TBc[type == 'base']) %>% 
  mutate(pPRur = PRur/Pop, pPUrb = PUrb/Pop, pPHig = PHig/Pop, pPLow = PLow/Pop,
         pPRL = PRL/Pop, pPRH = PRH/Pop, pPUL = PUL/Pop, pPUH = PUH/Pop) %>% 
  pivot_longer(cols = -c(time, type, run), names_to = "var", values_to = "values") %>% 
  group_by(time, type, var) %>% 
  summarise(val = median(values), lo = quantile(values, 0.025), hi = quantile(values, 0.975))

# 6. Plots ==========
prev_targets <- c(50, 25, 10)
dis_state <- factor(c("Clinical","Subclinical","Minimal"))
inf_dis <- c("Clinical", "Subclinical")
scenarios <- c("01: Xpert", "02: CXR->Xpert", "03: CXR", "Baseline")
type_label <- labeller(type = c("acfa" = "01: Xpert", "acfb" = "02: CXR->Xpert", "acfc" = "03: CXR", "base" = "Baseline"))
dem_urbrur <- c("Rural", "Urban")
dem_ses <- c("High","Low")

# TB prevalence per scenarios 
tiff(here("plots", "TBprev.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "TBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "TBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
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

# Proportion reduction of TB prevalence
tiff(here("plots", "TBreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pTBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pTBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
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
  geom_ribbon(data = filter(outs, var == "URt"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
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
  geom_ribbon(data = filter(outs, var == "HLt"), aes(x = time, ymin = (1-lo), ymax = (1-hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
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
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
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
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
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
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050),expand = c(0,0)) +
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
  geom_ribbon(data = filter(outs, var == "dcumMor"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative TB deaths averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False positives
tiff(here("plots", "Falsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "False positive diagnoses") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False positives
tiff(here("plots", "Falsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
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
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-9, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative false positive diagnoses") +
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
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-9, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050),expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative number screened") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()