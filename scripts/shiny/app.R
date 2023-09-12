## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 00: Shiny Fit.R

# Load required libraries
library(shiny)
library(deSolve)
library(tidyverse)
library(data.table)
library(rsconnect)

# Baseline parameters
parms = c(
  beta = 8,            # Contact (per person/year) parameter
  kappa = 0.75,          # Relative infectiousness
  gamma_infcle = 1.83,  # REG: Infected -> Cleared
  lambda_infmin = 0.21, # PROG: Infected -> Minimal
  gamma_mincle = 0.16,  # REG: Minimal -> Cleared
  theta_cleinf = 0.88,  # REINF: Cleared -> Infected
  lambda_minsub = 0.25, # PROG: Minimal -> Subclinical
  lambda_infsub = 0.07, # PROG: Infected -> Subclinical
  gamma_submin = 1.58,  # REG: Subclinical -> Minimal
  lambda_subcln = 0.77, # PROG: Subclinical -> Clinical
  gamma_clnsub = 0.53,  # REG: Clinical -> Subclinical
  omega_ini = 0.30,     # Mortality (Initial)
  omega_fin = 0.23,     # Mortality (Final)
  iota_cln_ini = 0.44,  # Diagnosis Clinical (Initial)
  iota_cln_fin = 0.9,  # Diagnosis Clinical (Final)
  phi_cln_ini = 0.69,   # Treatment failure Clinical (Initial) 
  phi_cln_fin = 0.09,   # Treatment failure Clinical (Final)
  tau_min = 0.03,       # RELAP: Recovered -> Minimal 
  tau_sub = 0.04,       # RELAP: Recovered -> Subclinical
  rho_ini = 0.78,       # Proportion rural (Initial)
  rho_fin = 0.04,       # Proportion rural (Final)
  sigma_ini = 0.76,     # Proportion low SES (Initial)
  sigma_fin = 0.16)     # Proportion low SES (Final)

# Parameter ranges
ranges = list(
  beta = c(5,20),                 # Contact (per person/year) parameter (Horton et al. 2022 - PLOS GPH)
  kappa = c(0.5,1),               # Relative infectiousness (Emery et al. 2022)
  gamma_infcle = c(0.93,3.30),    # REG: Infected -> Cleared (Horton et al. 2023)
  lambda_infmin = c(0.04,0.23),   # PROG: Infected -> Minimal (Horton et al. 2023)
  gamma_mincle = c(0.14,0.23),    # REG: Minimal -> Cleared (Horton et al. 2023)
  theta_cleinf = c(0.7,1),        # REINF: Cleared -> Infected (Andrews et al 2012)
  lambda_minsub = c(0.21,0.28),   # PROG: Minimal -> Subclinical (Horton et al. 2023)
  lambda_infsub = c(0.01,0.10),   # PROG: Infected -> Subclinical (Horton et al. 2023)
  gamma_submin = c(1.24,2.03),    # REG: Subclinical -> Minimal (Horton et al. 2023)
  lambda_subcln = c(0.56,0.94),   # PROG: Subclinical -> Clinical (Horton et al. 2023)
  gamma_clnsub = c(0.46,0.72),    # REG: Clinical -> Subclinical (Horton et al. 2023)
  omega_ini = c(0.26,0.37),       # Mortality Clinical (Initial) (Richards et al. 2023 - Lancet GH)
  omega_fin = c(0,0.26),          # Mortality Clinical (Final)
  iota_cln_ini = c(0,1),          # Diagnosis Clinical (Initial)
  iota_cln_fin = c(0,1),          # Diagnosis Clinical (Final)
  phi_cln_ini = c(0.11,1),        # Treatment failure Clinical (Initial) 
  phi_cln_fin = c(0.07,0.11),     # Treatment failure Clinical (Final) (WHO treatment success rate)
  tau_min = c(0.02,0.05),         # RELAP: Recovered -> Minimal 
  tau_sub = c(0.02,0.05),         # RELAP: Recovered -> Subclinical
  rho_ini = c(0.7,0.85),          # Proportion rural (Initial)
  rho_fin = c(0,0.2),             # Proportion rural (Final)
  sigma_ini = c(0.7,0.85),        # Proportion low SES (Inital)
  sigma_fin = c(0,0.2))           # Proportion low SES (Final)

# Fitting targets
# unc <- 0.1 # 10% uncertainty
# 
# targets <- list(
#   TBc2008 = list(val = 199, sigma = unc*199),   # TB prevalence rate 1st survey in adults (Nguyen et al. Emerg Infect Dis 2021)
#   TBc2018 = list(val = 125, sigma = unc*125),   # TB prevalence rate 2nd survey in adults (Nguyen et al. Emerg Infect Dis 2021)
#   Mor2000 = list(val = 40.5, sigma = unc*40.5), # TB mortality 2000 in total population (WHO + UN WPP)
#   Mor2010 = list(val = 25.2, sigma = unc*25.2), # TB mortality 2010 in total population (WHO + UN WPP)
#   Dxs2010 = list(val = 59.7, sigma = unc*59.7), # Notifications 2010 in total population (WHO + UN WPP)
#   Dxs2020 = list(val = 56.2, sigma = unc*56.2), # Notifications 2020 in total population (WHO + UN WPP)
#   Spr2008 = list(val = 0.70, sigma = unc*0.70), # Proportion scTB 2008 in adults (Emery et al. medRxiv 2022)
#   Spr2018 = list(val = 0.66, sigma = unc*0.66), # Proportion scTB 2018 in adults (Emery et al. medRxiv 2022)
#   URs2008 = list(val = 0.24, sigma = unc*0.24), # Relative urban/rural in scTB 2008 (Surveys dataset)
#   URs2018 = list(val = 0.44, sigma = unc*0.44), # Relative urban/rural in scTB 2018 (Surveys dataset)
#   URc2008 = list(val = 0.28, sigma = unc*0.28), # Relative urban/rural in cTB 2008 (Surveys dataset)
#   URc2018 = list(val = 0.43, sigma = unc*0.43), # Relative urban/rural in cTB 2018 (Surveys dataset)
#   HLs2008 = list(val = 0.38, sigma = unc*0.38), # Relative high/low SES in scTB 2008 (Surveys dataset)
#   HLs2018 = list(val = 0.47, sigma = unc*0.47), # Relative high/low SES in scTB 2018 (Surveys dataset)
#   HLc2008 = list(val = 0.32, sigma = unc*0.32), # Relative high/low SES in cTB 2008 (Surveys dataset)
#   HLc2018 = list(val = 0.46, sigma = unc*0.46)) # Relative high/low SES in cTB 2018 (Surveys dataset)
# 
# targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
# targetsdb$var <- substr(rownames(targetsdb),1,3) # Create variable classification
# targetsdb$time <- as.numeric(substr(rownames(targetsdb),4,7)) # Create time variable
# targetsdb <- as.data.frame(cbind(targetsdb[grep("val", rownames(targetsdb)),c("time","var","V1")],V2 = targetsdb[grep("sigma", rownames(targetsdb)),c("V1")])) # Widen dataframe
# colnames(targetsdb)[c(3,4)] <- c("val","sigma") # Rename columns
# targetsdb$lo <- targetsdb$val-2*targetsdb$sigma # Lower bound
# targetsdb$hi <- targetsdb$val+2*targetsdb$sigma # Upper bound
# rownames(targetsdb) <- gsub(rownames(targetsdb),pattern = "\\.val", replacement = "")

targets <- list(
  TBc2008 = c(160,248),
  TBc2018 = c(98,159),
  Mor2000 = c(47.8,71.6),
  Mor2010 = c(26.8,40.2),
  Dxs2010 = c(63.5,95.3),
  Dxs2020 = c(58.6,87.8),
  Spr2008 = c(0.56,0.84),
  Spr2018 = c(0.53,0.79))

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,3) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),4,7)) # Create time variable
colnames(targetsdb) <- c("lo","hi","var","time") # Rename columns

# ODE model
ode <- function(parms, N = 100000, end_time = 2020) {
  
  # Static parameters  
  mu <- 1/70      # Age expectancy adult
  delta <- 1      # Treatment year
  
  # Inactive parameters (no ACF)
  iota_min <- 0       # Diagnosis Minimal
  iota_sub <- 0       # Diagnosis Subclinical
  phi_min <- 0        # Treatment failure Minimal
  phi_sub <- 0        # Treatment failure Subclinical
  theta_recinf <- 1   # REINF: Recovered -> Infected (No protection)
  
  # Function parameters
  forcer_omega <- matrix(c(1500, parms['omega_ini'], 1999, parms['omega_ini'], 2020, parms['omega_fin']), ncol = 2, byrow = TRUE)
  force_func_omega <- approxfun(x = forcer_omega[,1], y = forcer_omega[,2], method = "linear", rule = 2)
  
  forcer_iota <- matrix(c(1500, parms['iota_cln_ini'], 1999, parms['iota_cln_ini'], 2020, parms['iota_cln_fin']), ncol = 2, byrow = TRUE)
  force_func_iota <- approxfun(x = forcer_iota[,1], y = forcer_iota[,2], method = "linear", rule = 2)
  
  forcer_phi <- matrix(c(1500, parms['phi_cln_ini'], 1999, parms['phi_cln_ini'], 2020, parms['phi_cln_fin']), ncol = 2, byrow = TRUE)
  force_func_phi <- approxfun(x = forcer_phi[,1], y = forcer_phi[,2], method = "linear", rule = 2)
  
  forcer_rho <- matrix(c(1500, parms['rho_ini'], 1999, parms['rho_ini'], 2020, parms['rho_fin']), ncol = 2, byrow = TRUE)
  force_func_rhoR <- approxfun(x = forcer_rho[,1], y = forcer_rho[,2], method = "linear", rule = 2)
  force_func_rhoU <- approxfun(x = forcer_rho[,1], y = 1-forcer_rho[,2], method = "linear", rule = 2)
  
  forcer_sigma <- matrix(c(1500, parms['sigma_ini'], 1999, parms['sigma_ini'], 2020, parms['sigma_fin']), ncol = 2, byrow = TRUE)
  force_func_sigmaL <- approxfun(x = forcer_sigma[,1], y = forcer_sigma[,2], method = "linear", rule = 2)
  force_func_sigmaH <- approxfun(x = forcer_sigma[,1], y = 1-forcer_sigma[,2], method = "linear", rule = 2)
  
  des <- function(time, state, parms) {
    
    with(as.list(c(state, parms)), {
      
      dN_RL  = force_func_rhoR(time)*force_func_sigmaL(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - mu*N_RL
      dI_RL  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf)+(P_RL*theta_recinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - mu*I_RL
      dO_RL  = gamma_infcle*I_RL + gamma_mincle*M_RL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf - mu*O_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - iota_min*M_RL + phi_min*RM_RL + tau_min*P_RL - mu*M_RL
      dRM_RL = iota_min*M_RL - phi_min*RM_RL - delta*RM_RL - mu*RM_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - iota_sub*S_RL + phi_sub*RS_RL + tau_sub*P_RL - mu*S_RL
      dRS_RL = iota_sub*S_RL - phi_sub*RS_RL - delta*RS_RL - mu*RS_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - force_func_omega(time)*C_RL - mu*C_RL - force_func_iota(time)*C_RL + force_func_phi(time)*RC_RL
      dRC_RL = force_func_iota(time)*C_RL - force_func_phi(time)*RC_RL - delta*RC_RL - mu*RC_RL
      dP_RL  = delta*(RM_RL+RS_RL+RC_RL) - tau_min*P_RL - tau_sub*P_RL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu*P_RL
      
      dN_RH  = force_func_rhoR(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu*N_RH
      dI_RH  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(O_RH*theta_cleinf)+(P_RH*theta_recinf)) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu*I_RH
      dO_RH  = gamma_infcle*I_RH + gamma_mincle*M_RH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_cleinf - mu*O_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH - iota_min*M_RH + phi_min*RM_RH + tau_min*P_RH - mu*M_RH
      dRM_RH = iota_min*M_RH - phi_min*RM_RH - delta*RM_RH - mu*RM_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH - iota_sub*S_RH + phi_sub*RS_RH + tau_sub*P_RH - mu*S_RH
      dRS_RH = iota_sub*S_RH - phi_sub*RS_RH - delta*RS_RH - mu*RS_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - force_func_omega(time)*C_RH - mu*C_RH - force_func_iota(time)*C_RH + force_func_phi(time)*RC_RH
      dRC_RH = force_func_iota(time)*C_RH - force_func_phi(time)*RC_RH - delta*RC_RH - mu*RC_RH
      dP_RH  = delta*(RM_RH+RS_RH+RC_RH) - tau_min*P_RH - tau_sub*P_RH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_recinf - mu*P_RH
      
      dN_UL  = force_func_rhoU(time)*force_func_sigmaL(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu*N_UL
      dI_UL  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(O_UL*theta_cleinf)+(P_UL*theta_recinf)) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu*I_UL
      dO_UL  = gamma_infcle*I_UL + gamma_mincle*M_UL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_cleinf - mu*O_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL - iota_min*M_UL + phi_min*RM_UL + tau_min*P_UL - mu*M_UL
      dRM_UL = iota_min*M_UL - phi_min*RM_UL - delta*RM_UL - mu*RM_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL - iota_sub*S_UL + phi_sub*RS_UL + tau_sub*P_UL - mu*S_UL
      dRS_UL = iota_sub*S_UL - phi_sub*RS_UL - delta*RS_UL - mu*RS_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - force_func_omega(time)*C_UL - mu*C_UL - force_func_iota(time)*C_UL + force_func_phi(time)*RC_UL
      dRC_UL = force_func_iota(time)*C_UL - force_func_phi(time)*RC_UL - delta*RC_UL - mu*RC_UL
      dP_UL  = delta*(RM_UL+RS_UL+RC_UL) - tau_min*P_UL - tau_sub*P_UL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_recinf - mu*P_UL
      
      dN_UH  = force_func_rhoU(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu*N_UH
      dI_UH  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(O_UH*theta_cleinf)+(P_UH*theta_recinf)) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu*I_UH
      dO_UH  = gamma_infcle*I_UH + gamma_mincle*M_UH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_cleinf - mu*O_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH - iota_min*M_UH + phi_min*RM_UH + tau_min*P_UH - mu*M_UH
      dRM_UH = iota_min*M_UH - phi_min*RM_UH - delta*RM_UH - mu*RM_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH - iota_sub*S_UH + phi_sub*RS_UH + tau_sub*P_UH - mu*S_UH
      dRS_UH = iota_sub*S_UH - phi_sub*RS_UH - delta*RS_UH - mu*RS_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - force_func_omega(time)*C_UH - mu*C_UH - force_func_iota(time)*C_UH + force_func_phi(time)*RC_UH
      dRC_UH = force_func_iota(time)*C_UH - force_func_phi(time)*RC_UH - delta*RC_UH - mu*RC_UH
      dP_UH  = delta*(RM_UH+RS_UH+RC_UH) - tau_min*P_UH - tau_sub*P_UH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_recinf - mu*P_UH
      
      return(list(c( 
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH, dRM_RL, dRM_RH, dRM_UL, dRM_UH,
        dS_RL, dS_RH, dS_UL, dS_UH, dRS_RL, dRS_RH, dRS_UL, dRS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dP_RL, dP_RH, dP_UL, dP_UH),
        Pop = (N_RL+N_RH+N_UL+N_UH+I_RL+I_RH+I_UL+I_UH+O_RL+O_RH+O_UL+O_UH+M_RL+M_RH+M_UL+M_UH+dRM_RL+dRM_RH+dRM_UL+dRM_UH+S_RL+S_RH+S_UL+S_UH+RS_RL+RS_RH+RS_UL+RS_UH+C_RL+C_RH+C_UL+C_UH+RC_RL+RC_RH+RC_UL+RC_UH+P_RL+P_RH+P_UL+P_UH), # Total population
        Sub = (S_RL+S_RH+S_UL+S_UH), # Subclinical TB (per 100k)
        Cln = (C_RL+C_RH+C_UL+C_UH), # Clinical TB (per 100k)
        TBc = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # All TB (per 100k)
        Mor = (force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)), # Clinical TB mortality per time (per 100k)
        Dxs = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH)), # Notifications cTB per time in adults (per 100k)
        Spr = (S_RL+S_RH+S_UL+S_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB
        URs = (S_UL+S_UH)/(S_RL+S_RH), # Relative urban/rural in scTB
        URc = (C_UL+C_UH)/(C_RL+C_RH), # Relative urban/rural in cTB
        HLs = (S_RH+S_UH)/(S_RL+S_UL), # Relative high/low SES in scTB
        HLc = (C_RH+C_UH)/(C_RL+C_UL), # Relative high/low SES in cTB
        ARIsi = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+N_RH+N_UL+N_UH), # ARI: Susceptible -> Infected (%) 
        ARIoi = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(O_RL+O_RH+O_UL+O_UH)*theta_cleinf, # ARI: Cleared -> Infected (%) 
        ARIpi = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(P_RL+P_RH+P_UL+P_UH)*theta_recinf, # ARI: Recovered -> Infected (%) 
        ARI = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH))))) # ARI
    })
  }
  
  yini <- c(N_RL = 99000, N_RH = 0, N_UL = 0, N_UH = 0, I_RL = 0, I_RH = 0, I_UL = 0, I_UH = 0, O_RL = 0, O_RH = 0, O_UL = 0, O_UH = 0,
            M_RL = 0, M_RH = 0, M_UL = 0, M_UH = 0, RM_RL = 0, RM_RH = 0, RM_UL = 0, RM_UH = 0,
            S_RL = 0, S_RH = 0, S_UL = 0, S_UH = 0, RS_RL = 0, RS_RH = 0, RS_UL = 0, RS_UH = 0,
            C_RL = 1000, C_RH = 0, C_UL = 0, C_UH = 0, RC_RL = 0, RC_RH = 0, RC_UL = 0, RC_UH = 0,
            P_RL = 0, P_RH = 0, P_UL = 0, P_UH = 0)
  
  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

ui <- fluidPage(
  titlePanel("ACF-VN"),
  sidebarPanel(
    fluidRow(
      column(width = 4,
             lapply(names(parms)[1:8], function(param) {
               sliderInput(
                 param,
                 label = param,
                 min = ranges[[param]][1],
                 max = ranges[[param]][2],
                 value = parms[param]
               )
             })
      ),
      column(width = 4,
             lapply(names(parms)[9:16], function(param) {
               sliderInput(
                 param,
                 label = param,
                 min = ranges[[param]][1],
                 max = ranges[[param]][2],
                 value = parms[param]
               )
             })
      ),
      column(width = 4,
             lapply(names(parms)[17:23], function(param) {
               sliderInput(
                 param,
                 label = param,
                 min = ranges[[param]][1],
                 max = ranges[[param]][2],
                 value = parms[param]
               )
             })
      )
    ),
    br(),
  ),
  mainPanel(
    fluidRow(
      column(width = 3,
             plotOutput("plot1")
      ),
      column(width = 3,
             plotOutput("plot2")
      ),
      column(width = 3,
             plotOutput("plot3")
      ),
      column(width = 3,
             plotOutput("plot4")
      )
    ),
    fluidRow(
      column(width = 3,
             plotOutput("plot5")
      ),
      column(width = 3,
             plotOutput("plot6")
      ),
      column(width = 3,
             plotOutput("plot7")
      ),
      column(width = 3,
             plotOutput("plot8")
      )
    )
  )
)

server <- function(input, output) {
  # Reactive function to update plots based on slider inputs
  model_output <- reactive({
    parms <- c(
      beta = input$beta,
      kappa = input$kappa,
      gamma_infcle = input$gamma_infcle,  
      lambda_infmin = input$lambda_infmin, 
      gamma_mincle = input$gamma_mincle,  
      theta_cleinf = input$theta_cleinf,  
      lambda_minsub = input$lambda_minsub, 
      lambda_infsub = input$lambda_infsub, 
      gamma_submin = input$gamma_submin,  
      lambda_subcln = input$lambda_subcln, 
      gamma_clnsub = input$gamma_clnsub,  
      omega_ini = input$omega_ini,     
      omega_fin = input$omega_fin,     
      iota_cln_ini = input$iota_cln_ini,   
      iota_cln_fin = input$iota_cln_fin,   
      phi_cln_ini = input$phi_cln_ini,    
      phi_cln_fin = input$phi_cln_fin,   
      tau_min = input$tau_min,       
      tau_sub = input$tau_sub,       
      rho_ini = input$rho_ini,        
      rho_fin = input$rho_fin,        
      sigma_ini = input$sigma_ini,      
      sigma_fin = input$sigma_fin)
    
    ode_output <- ode(parms = parms, N = 100000, end_time = 2020)
    ode_output <- as.data.frame(ode_output) %>% 
      select(time, TBc, Spr, Mor, Dxs, URs, URc, HLs, HLc) %>% 
      pivot_longer(cols = c(TBc, Spr, Mor, Dxs, URs, URc, HLs, HLc), names_to = "var", values_to = "val")
    return(ode_output)
  })
  
  output$plot1 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "TBc" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "TBc"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "TBc"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 300, 50), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,300)) +
      labs(title = "",
           x = "Year",
           y = "TB prevalence rate") +
      theme_minimal()
  })
  
  output$plot2 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "Spr" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "Spr"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "Spr"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,1)) +
      labs(title = "",
           x = "Year",
           y = "Proportion scTB") +
      theme_minimal()
  })
  
  output$plot3 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "Mor" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "Mor"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "Mor"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 60, 5), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,60)) +
      labs(title = "",
           x = "Year",
           y = "TB mortality rate") +
      theme_minimal()
  })
  
  output$plot4 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "Dxs" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "Dxs"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "Dxs"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,100)) +
      labs(title = "",
           x = "Year",
           y = "Notification rate") +
      theme_minimal()
  })
  
  output$plot5 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "URs" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "URs"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "URs"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,2)) +
      labs(title = "",
           x = "Year",
           y = "Urban/Rural (Subclinical)") +
      theme_minimal()
  })
  
  output$plot6 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "URc" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "URc"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "URc"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,2)) +
      labs(title = "",
           x = "Year",
           y = "Urban/Rural (Clinical)") +
      theme_minimal()
  })
  
  output$plot7 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "HLs" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "HLs"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "HLs"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,2)) +
      labs(title = "",
           x = "Year",
           y = "High/Low (Subclinical)") +
      theme_minimal()
  })
  
  output$plot8 <- renderPlot({
    ode_output <- model_output()
    
    ggplot(subset(ode_output, var == "HLc" & time %in% seq(2000,2020)), aes(x=time)) +
      geom_line(aes(y = val), linewidth = 0.6, na.rm=TRUE, colour = "#CE2931") +
      #geom_point(data = subset(targetsdb, var == "HLc"), aes(y = val), size = 1, colour = "black") +
      geom_errorbar(data = subset(targetsdb, var == "HLc"), aes(ymin = lo, ymax = hi), linewidth = 0.6, width = 0.5) + 
      scale_x_continuous(breaks = seq(2000, 2020, 2), limits = c(1999.5, 2020.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0,2)) +
      labs(title = "",
           x = "Year",
           y = "High/Low (Clinical)") +
      theme_minimal()
  })

}

shinyApp(ui = ui, server = server)
