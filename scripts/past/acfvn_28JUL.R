## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 02: Fit.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(hmer) # Emulation and History Matching
library(lhs) # Latin Hypercube Samples
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(purrr) # Complete set of tools for functions and vectors
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame, allows use of j operator (:=)
library(beepr) # Sound cues when code has finished running

# 1. Parameters ==========
# 1.1 Baseline parameters (for testing)
parms = c(
  beta = 50,            # Contact (per person/year) parameter
  kappa = 0.5,          # Relative infectiousness
  gamma_infcle = 1.83,  # REG: Infected -> Cleared
  lambda_infmin = 0.10, # PROG: Infected -> Minimal
  gamma_mincle = 0.18,  # REG: Minimal -> Cleared
  theta_cleinf = 0.21,  # REINF: Cleared -> Infected
  lambda_minsub = 0.24, # PROG: Minimal -> Subclinical
  lambda_infsub = 0.04, # PROG: Infected -> Subclinical
  gamma_submin = 1.58,  # REG: Subclinical -> Minimal
  lambda_subcln = 0.72, # PROG: Subclinical -> Clinical
  gamma_clnsub = 0.57,  # REG: Clinical -> Subclinical
  omega_ini = 0.33,     # Mortality (Initial)
  omega_fin = 0.33,     # Mortality (Final)
  iota_cln_ini = 0.2,   # Diagnosis Clinical (Initial)
  iota_cln_fin = 0.4,   # Diagnosis Clinical (Final)
  phi_cln_ini = 0.5,    # Treatment failure Clinical (Initial) 
  phi_cln_fin = 0.08,   # Treatment failure Clinical (Final)
  tau_min = 0.03,       # RELAP: Recovered -> Minimal 
  tau_sub = 0.03,       # RELAP: Recovered -> Subclinical
  rho_ini = 0.7,        # Proportion rural (Initial)
  rho_fin = 0.5,        # Proportion rural (Final)
  sigma_ini = 0.8,      # Proportion low SES (Initial)
  sigma_fin = 0.2,      # Proportion low SES (Final)
  ratio_ca = 1/20)      # Ratio children-to-adult

# 1.2 Parameter ranges
ranges = list(
  beta = c(0,50),                 # Contact (per person/year) parameter (Horton et al. 2022 - PLOS GPH)
  kappa = c(0.5,1),               # Relative infectiousness (Emery et al. 2022)
  gamma_infcle = c(0.93,3.30),    # REG: Infected -> Cleared (Horton et al. 2023)
  lambda_infmin = c(0.04,0.23),   # PROG: Infected -> Minimal (Horton et al. 2023)
  gamma_mincle = c(0.14,0.23),    # REG: Minimal -> Cleared (Horton et al. 2023)
  theta_cleinf = c(0.14,0.30),    # REINF: Cleared -> Infected (Andrews et al 2012)
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
  rho_ini = c(0,1),               # Proportion rural (Initial)
  rho_fin = c(0,1),               # Proportion rural (Final)
  sigma_ini = c(0,1),             # Proportion low SES (Inital)
  sigma_fin = c(0,1),             # Proportion low SES (Final)
  ratio_ca = c(1/20,1/15))        # Child-to-adult ratio

# 2. Model ==========
ode <- function(parms, N = 100000, end_time = 2020) {
  
  # Static parameters  
  muc <- 1/70      # Age expectancy children
  mua <- 1/55      # Age expectancy adult
  alpha <- 1/15    # Age transition
  delta <- 1       # Treatment year
  
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
      
      dNc_RL = force_func_rhoR(time)*force_func_sigmaL(time)*(mua*(Na_RL+Ia_RL+Oa_RL+Ma_RL+RMa_RL+Sa_RL+RSa_RL+Ca_RL+RCa_RL+Pa_RL+Na_RH+Ia_RH+Oa_RH+Ma_RH+RMa_RH+Sa_RH+RSa_RH+Ca_RH+RCa_RH+Pa_RH+Na_UL+Ia_UL+Oa_UL+Ma_UL+RMa_UL+Sa_UL+RSa_UL+Ca_UL+RCa_UL+Pa_UL+Na_UH+Ia_UH+Oa_UH+Ma_UH+RMa_UH+Sa_UH+RSa_UH+Ca_UH+RCa_UH+Pa_UH) + 
        muc*(Nc_RL+Ic_RL+Oc_RL+Mc_RL+RMc_RL+Sc_RL+RSc_RL+Cc_RL+RCc_RL+Pc_RL+Nc_RH+Ic_RH+Oc_RH+Mc_RH+RMc_RH+Sc_RH+RSc_RH+Cc_RH+RCc_RH+Pc_RH+Nc_UL+Ic_UL+Oc_UL+Mc_UL+RMc_UL+Sc_UL+RSc_UL+Cc_UL+RCc_UL+Pc_UL+Nc_UH+Ic_UH+Oc_UH+Mc_UH+RMc_UH+Sc_UH+RSc_UH+Cc_UH+RCc_UH+Pc_UH) + force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)) - 
        beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Nc_RL/N) - muc*Nc_RL - alpha*Nc_RL
      dIc_RL  = beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Nc_RL/N)+((Oc_RL/N)*theta_cleinf)+((Pc_RL/N)*theta_recinf)) - gamma_infcle*Ic_RL*ratio_ca - lambda_infmin*Ic_RL*ratio_ca - lambda_infsub*Ic_RL*ratio_ca - muc*Ic_RL - alpha*Ic_RL
      dOc_RL  = gamma_infcle*Ic_RL*ratio_ca + gamma_mincle*Mc_RL*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oc_RL/N)*theta_cleinf) - muc*Oc_RL - alpha*Oc_RL  
      dMc_RL  = lambda_infmin*Ic_RL*ratio_ca + gamma_submin*Sc_RL*ratio_ca - gamma_mincle*Mc_RL*ratio_ca - lambda_minsub*Mc_RL*ratio_ca - iota_min*Mc_RL*ratio_ca + phi_min*RMc_RL*ratio_ca + tau_min*Pc_RL*ratio_ca - muc*Mc_RL - alpha*Mc_RL
      dRMc_RL = iota_min*Mc_RL*ratio_ca - phi_min*RMc_RL*ratio_ca - delta*RMc_RL - muc*RMc_RL - alpha*RMc_RL
      dSc_RL  = lambda_infsub*Ic_RL*ratio_ca + lambda_minsub*Mc_RL*ratio_ca + gamma_clnsub*Cc_RL*ratio_ca - gamma_submin*Sc_RL*ratio_ca - lambda_subcln*Sc_RL*ratio_ca - iota_sub*Sc_RL*ratio_ca + phi_sub*RSc_RL*ratio_ca + tau_sub*Pc_RL*ratio_ca - muc*Sc_RL - alpha*Sc_RL
      dRSc_RL = iota_sub*Sc_RL*ratio_ca - phi_sub*RSc_RL*ratio_ca - delta*RSc_RL - muc*RSc_RL - alpha*RSc_RL
      dCc_RL  = lambda_subcln*Sc_RL*ratio_ca - gamma_clnsub*Cc_RL*ratio_ca - force_func_omega(time)*Cc_RL - muc*Cc_RL - force_func_iota(time)*Cc_RL*ratio_ca + force_func_phi(time)*RCc_RL*ratio_ca - alpha*Cc_RL
      dRCc_RL = force_func_iota(time)*Cc_RL*ratio_ca - force_func_phi(time)*RCc_RL*ratio_ca - delta*RCc_RL - muc*RCc_RL - alpha*RCc_RL
      dPc_RL  = delta*(RMc_RL+RSc_RL+RCc_RL) - tau_min*Pc_RL*ratio_ca - tau_sub*Pc_RL*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pc_RL/N)*theta_recinf) - muc*Pc_RL - alpha*Pc_RL
      
      dNa_RL  = - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Na_RL/N) - mua*Na_RL + alpha*Nc_RL
      dIa_RL  = beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Na_RL/N)+((Oa_RL/N)*theta_cleinf)+((Pa_RL/N)*theta_recinf)) - gamma_infcle*Ia_RL - lambda_infmin*Ia_RL - lambda_infsub*Ia_RL - mua*Ia_RL + alpha*Ic_RL
      dOa_RL  = gamma_infcle*Ia_RL + gamma_mincle*Ma_RL - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oa_RL/N)*theta_cleinf) - mua*Oa_RL + alpha*Oc_RL  
      dMa_RL  = lambda_infmin*Ia_RL + gamma_submin*Sa_RL - gamma_mincle*Ma_RL - lambda_minsub*Ma_RL - iota_min*Ma_RL + phi_min*RMa_RL + tau_min*Pa_RL - mua*Ma_RL + alpha*Mc_RL
      dRMa_RL = iota_min*Ma_RL - phi_min*RMa_RL - delta*RMa_RL - mua*RMa_RL + alpha*RMc_RL
      dSa_RL  = lambda_infsub*Ia_RL + lambda_minsub*Ma_RL + gamma_clnsub*Ca_RL - gamma_submin*Sa_RL - lambda_subcln*Sa_RL - iota_sub*Sa_RL + phi_sub*RSa_RL + tau_sub*Pa_RL - mua*Sa_RL + alpha*Sc_RL
      dRSa_RL = iota_sub*Sa_RL - phi_sub*RSa_RL - delta*RSa_RL - mua*RSa_RL + alpha*RSc_RL
      dCa_RL  = lambda_subcln*Sa_RL - gamma_clnsub*Ca_RL - force_func_omega(time)*Ca_RL - mua*Ca_RL - force_func_iota(time)*Ca_RL + force_func_phi(time)*RCa_RL + alpha*Cc_RL
      dRCa_RL = force_func_iota(time)*Ca_RL - force_func_phi(time)*RCa_RL - delta*RCa_RL - mua*RCa_RL + alpha*RCc_RL
      dPa_RL  = delta*(RMa_RL+RSa_RL+RCa_RL) - tau_min*Pa_RL - tau_sub*Pa_RL - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pa_RL/N)*theta_recinf) - mua*Pa_RL + alpha*Pc_RL
      
      dNc_RH = force_func_rhoR(time)*force_func_sigmaH(time)*(mua*(Na_RL+Ia_RL+Oa_RL+Ma_RL+RMa_RL+Sa_RL+RSa_RL+Ca_RL+RCa_RL+Pa_RL+Na_RH+Ia_RH+Oa_RH+Ma_RH+RMa_RH+Sa_RH+RSa_RH+Ca_RH+RCa_RH+Pa_RH+Na_UL+Ia_UL+Oa_UL+Ma_UL+RMa_UL+Sa_UL+RSa_UL+Ca_UL+RCa_UL+Pa_UL+Na_UH+Ia_UH+Oa_UH+Ma_UH+RMa_UH+Sa_UH+RSa_UH+Ca_UH+RCa_UH+Pa_UH) + 
        muc*(Nc_RL+Ic_RL+Oc_RL+Mc_RL+RMc_RL+Sc_RL+RSc_RL+Cc_RL+RCc_RL+Pc_RL+Nc_RH+Ic_RH+Oc_RH+Mc_RH+RMc_RH+Sc_RH+RSc_RH+Cc_RH+RCc_RH+Pc_RH+Nc_UL+Ic_UL+Oc_UL+Mc_UL+RMc_UL+Sc_UL+RSc_UL+Cc_UL+RCc_UL+Pc_UL+Nc_UH+Ic_UH+Oc_UH+Mc_UH+RMc_UH+Sc_UH+RSc_UH+Cc_UH+RCc_UH+Pc_UH) + force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)) - 
        beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Nc_RH/N) - muc*Nc_RH - alpha*Nc_RH
      dIc_RH  = beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Nc_RH/N)+((Oc_RH/N)*theta_cleinf)+((Pc_RH/N)*theta_recinf)) - gamma_infcle*Ic_RH*ratio_ca - lambda_infmin*Ic_RH*ratio_ca - lambda_infsub*Ic_RH*ratio_ca - muc*Ic_RH - alpha*Ic_RH
      dOc_RH  = gamma_infcle*Ic_RH*ratio_ca + gamma_mincle*Mc_RH*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oc_RH/N)*theta_cleinf) - muc*Oc_RH - alpha*Oc_RH  
      dMc_RH  = lambda_infmin*Ic_RH*ratio_ca + gamma_submin*Sc_RH*ratio_ca - gamma_mincle*Mc_RH*ratio_ca - lambda_minsub*Mc_RH*ratio_ca - iota_min*Mc_RH*ratio_ca + phi_min*RMc_RH*ratio_ca + tau_min*Pc_RH*ratio_ca - muc*Mc_RH - alpha*Mc_RH
      dRMc_RH = iota_min*Mc_RH*ratio_ca - phi_min*RMc_RH*ratio_ca - delta*RMc_RH - muc*RMc_RH - alpha*RMc_RH
      dSc_RH  = lambda_infsub*Ic_RH*ratio_ca + lambda_minsub*Mc_RH*ratio_ca + gamma_clnsub*Cc_RH*ratio_ca - gamma_submin*Sc_RH*ratio_ca - lambda_subcln*Sc_RH*ratio_ca - iota_sub*Sc_RH*ratio_ca + phi_sub*RSc_RH*ratio_ca + tau_sub*Pc_RH*ratio_ca - muc*Sc_RH - alpha*Sc_RH
      dRSc_RH = iota_sub*Sc_RH*ratio_ca - phi_sub*RSc_RH*ratio_ca - delta*RSc_RH - muc*RSc_RH - alpha*RSc_RH
      dCc_RH  = lambda_subcln*Sc_RH*ratio_ca - gamma_clnsub*Cc_RH*ratio_ca - force_func_omega(time)*Cc_RH - muc*Cc_RH - force_func_iota(time)*Cc_RH*ratio_ca + force_func_phi(time)*RCc_RH*ratio_ca - alpha*Cc_RH
      dRCc_RH = force_func_iota(time)*Cc_RH*ratio_ca - force_func_phi(time)*RCc_RH*ratio_ca - delta*RCc_RH - muc*RCc_RH - alpha*RCc_RH
      dPc_RH  = delta*(RMc_RH+RSc_RH+RCc_RH) - tau_min*Pc_RH*ratio_ca - tau_sub*Pc_RH*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pc_RH/N)*theta_recinf) - muc*Pc_RH - alpha*Pc_RH
      
      dNa_RH  = - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Na_RH/N) - mua*Na_RH + alpha*Nc_RH
      dIa_RH  = beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Na_RH/N)+((Oa_RH/N)*theta_cleinf)+((Pa_RH/N)*theta_recinf)) - gamma_infcle*Ia_RH - lambda_infmin*Ia_RH - lambda_infsub*Ia_RH - mua*Ia_RH + alpha*Ic_RH
      dOa_RH  = gamma_infcle*Ia_RH + gamma_mincle*Ma_RH - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oa_RH/N)*theta_cleinf) - mua*Oa_RH + alpha*Oc_RH  
      dMa_RH  = lambda_infmin*Ia_RH + gamma_submin*Sa_RH - gamma_mincle*Ma_RH - lambda_minsub*Ma_RH - iota_min*Ma_RH + phi_min*RMa_RH + tau_min*Pa_RH - mua*Ma_RH + alpha*Mc_RH
      dRMa_RH = iota_min*Ma_RH - phi_min*RMa_RH - delta*RMa_RH - mua*RMa_RH + alpha*RMc_RH
      dSa_RH  = lambda_infsub*Ia_RH + lambda_minsub*Ma_RH + gamma_clnsub*Ca_RH - gamma_submin*Sa_RH - lambda_subcln*Sa_RH - iota_sub*Sa_RH + phi_sub*RSa_RH + tau_sub*Pa_RH - mua*Sa_RH + alpha*Sc_RH
      dRSa_RH = iota_sub*Sa_RH - phi_sub*RSa_RH - delta*RSa_RH - mua*RSa_RH + alpha*RSc_RH
      dCa_RH  = lambda_subcln*Sa_RH - gamma_clnsub*Ca_RH - force_func_omega(time)*Ca_RH - mua*Ca_RH - force_func_iota(time)*Ca_RH + force_func_phi(time)*RCa_RH + alpha*Cc_RH
      dRCa_RH = force_func_iota(time)*Ca_RH - force_func_phi(time)*RCa_RH - delta*RCa_RH - mua*RCa_RH + alpha*RCc_RH
      dPa_RH  = delta*(RMa_RH+RSa_RH+RCa_RH) - tau_min*Pa_RH - tau_sub*Pa_RH - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pa_RH/N)*theta_recinf) - mua*Pa_RH + alpha*Pc_RH
      
      dNc_UL = force_func_rhoU(time)*force_func_sigmaL(time)*(mua*(Na_RL+Ia_RL+Oa_RL+Ma_RL+RMa_RL+Sa_RL+RSa_RL+Ca_RL+RCa_RL+Pa_RL+Na_RH+Ia_RH+Oa_RH+Ma_RH+RMa_RH+Sa_RH+RSa_RH+Ca_RH+RCa_RH+Pa_RH+Na_UL+Ia_UL+Oa_UL+Ma_UL+RMa_UL+Sa_UL+RSa_UL+Ca_UL+RCa_UL+Pa_UL+Na_UH+Ia_UH+Oa_UH+Ma_UH+RMa_UH+Sa_UH+RSa_UH+Ca_UH+RCa_UH+Pa_UH) + 
        muc*(Nc_RL+Ic_RL+Oc_RL+Mc_RL+RMc_RL+Sc_RL+RSc_RL+Cc_RL+RCc_RL+Pc_RL+Nc_RH+Ic_RH+Oc_RH+Mc_RH+RMc_RH+Sc_RH+RSc_RH+Cc_RH+RCc_RH+Pc_RH+Nc_UL+Ic_UL+Oc_UL+Mc_UL+RMc_UL+Sc_UL+RSc_UL+Cc_UL+RCc_UL+Pc_UL+Nc_UH+Ic_UH+Oc_UH+Mc_UH+RMc_UH+Sc_UH+RSc_UH+Cc_UH+RCc_UH+Pc_UH) + force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)) - 
        beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Nc_UL/N) - muc*Nc_UL - alpha*Nc_UL
      dIc_UL  = beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Nc_UL/N)+((Oc_UL/N)*theta_cleinf)+((Pc_UL/N)*theta_recinf)) - gamma_infcle*Ic_UL*ratio_ca - lambda_infmin*Ic_UL*ratio_ca - lambda_infsub*Ic_UL*ratio_ca - muc*Ic_UL - alpha*Ic_UL
      dOc_UL  = gamma_infcle*Ic_UL*ratio_ca + gamma_mincle*Mc_UL*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oc_UL/N)*theta_cleinf) - muc*Oc_UL - alpha*Oc_UL  
      dMc_UL  = lambda_infmin*Ic_UL*ratio_ca + gamma_submin*Sc_UL*ratio_ca - gamma_mincle*Mc_UL*ratio_ca - lambda_minsub*Mc_UL*ratio_ca - iota_min*Mc_UL*ratio_ca + phi_min*RMc_UL*ratio_ca + tau_min*Pc_UL*ratio_ca - muc*Mc_UL - alpha*Mc_UL
      dRMc_UL = iota_min*Mc_UL*ratio_ca - phi_min*RMc_UL*ratio_ca - delta*RMc_UL - muc*RMc_UL - alpha*RMc_UL
      dSc_UL  = lambda_infsub*Ic_UL*ratio_ca + lambda_minsub*Mc_UL*ratio_ca + gamma_clnsub*Cc_UL*ratio_ca - gamma_submin*Sc_UL*ratio_ca - lambda_subcln*Sc_UL*ratio_ca - iota_sub*Sc_UL*ratio_ca + phi_sub*RSc_UL*ratio_ca + tau_sub*Pc_UL*ratio_ca - muc*Sc_UL - alpha*Sc_UL
      dRSc_UL = iota_sub*Sc_UL*ratio_ca - phi_sub*RSc_UL*ratio_ca - delta*RSc_UL - muc*RSc_UL - alpha*RSc_UL
      dCc_UL  = lambda_subcln*Sc_UL*ratio_ca - gamma_clnsub*Cc_UL*ratio_ca - force_func_omega(time)*Cc_UL - muc*Cc_UL - force_func_iota(time)*Cc_UL*ratio_ca + force_func_phi(time)*RCc_UL*ratio_ca - alpha*Cc_UL
      dRCc_UL = force_func_iota(time)*Cc_UL*ratio_ca - force_func_phi(time)*RCc_UL*ratio_ca - delta*RCc_UL - muc*RCc_UL - alpha*RCc_UL
      dPc_UL  = delta*(RMc_UL+RSc_UL+RCc_UL) - tau_min*Pc_UL*ratio_ca - tau_sub*Pc_UL*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pc_UL/N)*theta_recinf) - muc*Pc_UL - alpha*Pc_UL
      
      dNa_UL  = - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Na_UL/N) - mua*Na_UL + alpha*Nc_UL
      dIa_UL  = beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Na_UL/N)+((Oa_UL/N)*theta_cleinf)+((Pa_UL/N)*theta_recinf)) - gamma_infcle*Ia_UL - lambda_infmin*Ia_UL - lambda_infsub*Ia_UL - mua*Ia_UL + alpha*Ic_UL
      dOa_UL  = gamma_infcle*Ia_UL + gamma_mincle*Ma_UL - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oa_UL/N)*theta_cleinf) - mua*Oa_UL + alpha*Oc_UL  
      dMa_UL  = lambda_infmin*Ia_UL + gamma_submin*Sa_UL - gamma_mincle*Ma_UL - lambda_minsub*Ma_UL - iota_min*Ma_UL + phi_min*RMa_UL + tau_min*Pa_UL - mua*Ma_UL + alpha*Mc_UL
      dRMa_UL = iota_min*Ma_UL - phi_min*RMa_UL - delta*RMa_UL - mua*RMa_UL + alpha*RMc_UL
      dSa_UL  = lambda_infsub*Ia_UL + lambda_minsub*Ma_UL + gamma_clnsub*Ca_UL - gamma_submin*Sa_UL - lambda_subcln*Sa_UL - iota_sub*Sa_UL + phi_sub*RSa_UL + tau_sub*Pa_UL - mua*Sa_UL + alpha*Sc_UL
      dRSa_UL = iota_sub*Sa_UL - phi_sub*RSa_UL - delta*RSa_UL - mua*RSa_UL + alpha*RSc_UL
      dCa_UL  = lambda_subcln*Sa_UL - gamma_clnsub*Ca_UL - force_func_omega(time)*Ca_UL - mua*Ca_UL - force_func_iota(time)*Ca_UL + force_func_phi(time)*RCa_UL + alpha*Cc_UL
      dRCa_UL = force_func_iota(time)*Ca_UL - force_func_phi(time)*RCa_UL - delta*RCa_UL - mua*RCa_UL + alpha*RCc_UL
      dPa_UL  = delta*(RMa_UL+RSa_UL+RCa_UL) - tau_min*Pa_UL - tau_sub*Pa_UL - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pa_UL/N)*theta_recinf) - mua*Pa_UL + alpha*Pc_UL
      
      dNc_UH = force_func_rhoU(time)*force_func_sigmaH(time)*(mua*(Na_RL+Ia_RL+Oa_RL+Ma_RL+RMa_RL+Sa_RL+RSa_RL+Ca_RL+RCa_RL+Pa_RL+Na_RH+Ia_RH+Oa_RH+Ma_RH+RMa_RH+Sa_RH+RSa_RH+Ca_RH+RCa_RH+Pa_RH+Na_UL+Ia_UL+Oa_UL+Ma_UL+RMa_UL+Sa_UL+RSa_UL+Ca_UL+RCa_UL+Pa_UL+Na_UH+Ia_UH+Oa_UH+Ma_UH+RMa_UH+Sa_UH+RSa_UH+Ca_UH+RCa_UH+Pa_UH) + 
        muc*(Nc_RL+Ic_RL+Oc_RL+Mc_RL+RMc_RL+Sc_RL+RSc_RL+Cc_RL+RCc_RL+Pc_RL+Nc_RH+Ic_RH+Oc_RH+Mc_RH+RMc_RH+Sc_RH+RSc_RH+Cc_RH+RCc_RH+Pc_RH+Nc_UL+Ic_UL+Oc_UL+Mc_UL+RMc_UL+Sc_UL+RSc_UL+Cc_UL+RCc_UL+Pc_UL+Nc_UH+Ic_UH+Oc_UH+Mc_UH+RMc_UH+Sc_UH+RSc_UH+Cc_UH+RCc_UH+Pc_UH) + force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)) - 
        beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Nc_UH/N) - muc*Nc_UH - alpha*Nc_UH
      dIc_UH  = beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Nc_UH/N)+((Oc_UH/N)*theta_cleinf)+((Pc_UH/N)*theta_recinf)) - gamma_infcle*Ic_UH*ratio_ca - lambda_infmin*Ic_UH*ratio_ca - lambda_infsub*Ic_UH*ratio_ca - muc*Ic_UH - alpha*Ic_UH
      dOc_UH  = gamma_infcle*Ic_UH*ratio_ca + gamma_mincle*Mc_UH*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oc_UH/N)*theta_cleinf) - muc*Oc_UH - alpha*Oc_UH  
      dMc_UH  = lambda_infmin*Ic_UH*ratio_ca + gamma_submin*Sc_UH*ratio_ca - gamma_mincle*Mc_UH*ratio_ca - lambda_minsub*Mc_UH*ratio_ca - iota_min*Mc_UH*ratio_ca + phi_min*RMc_UH*ratio_ca + tau_min*Pc_UH*ratio_ca - muc*Mc_UH - alpha*Mc_UH
      dRMc_UH = iota_min*Mc_UH*ratio_ca - phi_min*RMc_UH*ratio_ca - delta*RMc_UH - muc*RMc_UH - alpha*RMc_UH
      dSc_UH  = lambda_infsub*Ic_UH*ratio_ca + lambda_minsub*Mc_UH*ratio_ca + gamma_clnsub*Cc_UH*ratio_ca - gamma_submin*Sc_UH*ratio_ca - lambda_subcln*Sc_UH*ratio_ca - iota_sub*Sc_UH*ratio_ca + phi_sub*RSc_UH*ratio_ca + tau_sub*Pc_UH*ratio_ca - muc*Sc_UH - alpha*Sc_UH
      dRSc_UH = iota_sub*Sc_UH*ratio_ca - phi_sub*RSc_UH*ratio_ca - delta*RSc_UH - muc*RSc_UH - alpha*RSc_UH
      dCc_UH  = lambda_subcln*Sc_UH*ratio_ca - gamma_clnsub*Cc_UH*ratio_ca - force_func_omega(time)*Cc_UH - muc*Cc_UH - force_func_iota(time)*Cc_UH*ratio_ca + force_func_phi(time)*RCc_UH*ratio_ca - alpha*Cc_UH
      dRCc_UH = force_func_iota(time)*Cc_UH*ratio_ca - force_func_phi(time)*RCc_UH*ratio_ca - delta*RCc_UH - muc*RCc_UH - alpha*RCc_UH
      dPc_UH  = delta*(RMc_UH+RSc_UH+RCc_UH) - tau_min*Pc_UH*ratio_ca - tau_sub*Pc_UH*ratio_ca - beta*ratio_ca*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pc_UH/N)*theta_recinf) - muc*Pc_UH - alpha*Pc_UH
      
      dNa_UH  = - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*(Na_UH/N) - mua*Na_UH + alpha*Nc_UH
      dIa_UH  = beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Na_UH/N)+((Oa_UH/N)*theta_cleinf)+((Pa_UH/N)*theta_recinf)) - gamma_infcle*Ia_UH - lambda_infmin*Ia_UH - lambda_infsub*Ia_UH - mua*Ia_UH + alpha*Ic_UH
      dOa_UH  = gamma_infcle*Ia_UH + gamma_mincle*Ma_UH - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Oa_UH/N)*theta_cleinf) - mua*Oa_UH + alpha*Oc_UH  
      dMa_UH  = lambda_infmin*Ia_UH + gamma_submin*Sa_UH - gamma_mincle*Ma_UH - lambda_minsub*Ma_UH - iota_min*Ma_UH + phi_min*RMa_UH + tau_min*Pa_UH - mua*Ma_UH + alpha*Mc_UH
      dRMa_UH = iota_min*Ma_UH - phi_min*RMa_UH - delta*RMa_UH - mua*RMa_UH + alpha*RMc_UH
      dSa_UH  = lambda_infsub*Ia_UH + lambda_minsub*Ma_UH + gamma_clnsub*Ca_UH - gamma_submin*Sa_UH - lambda_subcln*Sa_UH - iota_sub*Sa_UH + phi_sub*RSa_UH + tau_sub*Pa_UH - mua*Sa_UH + alpha*Sc_UH
      dRSa_UH = iota_sub*Sa_UH - phi_sub*RSa_UH - delta*RSa_UH - mua*RSa_UH + alpha*RSc_UH
      dCa_UH  = lambda_subcln*Sa_UH - gamma_clnsub*Ca_UH - force_func_omega(time)*Ca_UH - mua*Ca_UH - force_func_iota(time)*Ca_UH + force_func_phi(time)*RCa_UH + alpha*Cc_UH
      dRCa_UH = force_func_iota(time)*Ca_UH - force_func_phi(time)*RCa_UH - delta*RCa_UH - mua*RCa_UH + alpha*RCc_UH
      dPa_UH  = delta*(RMa_UH+RSa_UH+RCa_UH) - tau_min*Pa_UH - tau_sub*Pa_UH - beta*(kappa*(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)+(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH))*((Pa_UH/N)*theta_recinf) - mua*Pa_UH + alpha*Pc_UH
      
      return(list(c( 
        dNa_RL, dNa_RH, dNa_UL, dNa_UH, dNc_RL, dNc_RH, dNc_UL, dNc_UH,          # Susceptible
        dIa_RL, dIa_RH, dIa_UL, dIa_UH, dIc_RL, dIc_RH, dIc_UL, dIc_UH,          # Infected
        dOa_RL, dOa_RH, dOa_UL, dOa_UH, dOc_RL, dOc_RH, dOc_UL, dOc_UH,          # Cleared
        dMa_RL, dMa_RH, dMa_UL, dMa_UH, dMc_RL, dMc_RH, dMc_UL, dMc_UH,          # Minimal
        dRMa_RL, dRMa_RH, dRMa_UL, dRMa_UH, dRMc_RL, dRMc_RH, dRMc_UL, dRMc_UH,  # Treatment Minimal
        dSa_RL, dSa_RH, dSa_UL, dSa_UH, dSc_RL, dSc_RH, dSc_UL, dSc_UH,          # Subclinical
        dRSa_RL, dRSa_RH, dRSa_UL, dRSa_UH, dRSc_RL, dRSc_RH, dRSc_UL, dRSc_UH,  # Treatment Subclinical
        dCa_RL, dCa_RH, dCa_UL, dCa_UH, dCc_RL, dCc_RH, dCc_UL, dCc_UH,          # Clinical
        dRCa_RL, dRCa_RH, dRCa_UL, dRCa_UH, dRCc_RL, dRCc_RH, dRCc_UL, dRCc_UH,  # Treatment Clinical
        dPa_RL, dPa_RH, dPa_UL, dPa_UH, dPc_RL, dPc_RH, dPc_UL, dPc_UH),         # Post TB treatment
        Pa = (Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH), # Population adults
        Pc = (Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH), # Population children
        Sa = (Sa_RL+Sa_RH+Sa_UL+Sa_UH)/(Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH)*1e5, # All subclinical TB in adults (per 100k)
        Sc = (Sc_RL+Sc_RH+Sc_UL+Sc_UH)/(Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH)*1e5, # All subclinical TB in children (per 100k)
        St = (Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH), # All subclinical TB in all (per 100k)
        Ca = (Ca_RL+Ca_RH+Ca_UL+Ca_UH)/(Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH)*1e5, # All clinical TB in adults (per 100k)
        Cc = (Cc_RL+Cc_RH+Cc_UL+Cc_UH)/(Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH)*1e5, # All clinical TB in children (per 100k)
        Ct = (Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH), # All clinical TB in all (per 100k)
        Ta = (Sa_RL+Sa_RH+Sa_UL+Sa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH)/(Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH)*1e5, # All TB in adults (per 100k)
        Tc = (Sc_RL+Sc_RH+Sc_UL+Sc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)/(Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH)*1e5, # All TB in children (per 100k)
        Tt = (Sa_RL+Sa_RH+Sa_UL+Sa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH), # All TB in all (per 100k)
        Ua = (force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH))/(Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH)*1e5, # Clinical TB mortality per time in adults (per 100k)
        Uc = (force_func_omega(time)*(Cc_RL+Cc_RH+Cc_UL+Cc_UH))/(Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH)*1e5, # Clinical TB mortality per time in children (per 100k)
        Ut = (force_func_omega(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH)), # Clinical TB mortality per time in all (per 100k)
        Va = (force_func_iota(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH))/(Na_RL+Na_RH+Na_UL+Na_UH+Ia_RL+Ia_RH+Ia_UL+Ia_UH+Oa_RL+Oa_RH+Oa_UL+Oa_UH+Ma_RL+Ma_RH+Ma_UL+Ma_UH+dRMa_RL+dRMa_RH+dRMa_UL+dRMa_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH+RSa_RL+RSa_RH+RSa_UL+RSa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+RCa_RL+RCa_RH+RCa_UL+RCa_UH+Pa_RL+Pa_RH+Pa_UL+Pa_UH)*1e5, # Notifications cTB per time in adults (per 100k)
        Vc = (force_func_iota(time)*(Cc_RL+Cc_RH+Cc_UL+Cc_UH)*ratio_ca)/(Nc_RL+Nc_RH+Nc_UL+Nc_UH+Ic_RL+Ic_RH+Ic_UL+Ic_UH+Oc_RL+Oc_RH+Oc_UL+Oc_UH+Mc_RL+Mc_RH+Mc_UL+Mc_UH+dRMc_RL+dRMc_RH+dRMc_UL+dRMc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+RSc_RL+RSc_RH+RSc_UL+RSc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH+RCc_RL+RCc_RH+RCc_UL+RCc_UH+Pc_RL+Pc_RH+Pc_UL+Pc_UH)*1e5, # Notifications cTB per time in children (per 100k)
        Vt = (force_func_iota(time)*(Ca_RL+Ca_RH+Ca_UL+Ca_UH))+(force_func_iota(time)*Cc_RL+Cc_RH+Cc_UL+Cc_UH*ratio_ca), # Notifications cTB per time in total population
        Ra = (Sa_RL+Sa_RH+Sa_UL+Sa_UH)/(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH), # Proportion scTB in adults
        Rc = (Sc_RL+Sc_RH+Sc_UL+Sc_UH)/(Sc_RL+Sc_RH+Sc_UL+Sc_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH), # Proportion scTB in children
        Rt = (Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)/(Sa_RL+Sa_RH+Sa_UL+Sa_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH+Ca_RL+Ca_RH+Ca_UL+Ca_UH+Cc_RL+Cc_RH+Cc_UL+Cc_UH), # Proportion scTB in all
        Wa = (Sa_UL+Sa_UH)/(Sa_RL+Sa_RH), # Relative urban/rural in scTB (rS_UR) in adults
        Wc = (Sc_UL+Sc_UH)/(Sc_RL+Sc_RH), # Relative urban/rural in scTB (rS_UR) in children
        Wt = (Sa_UL+Sa_UH+Sc_UL+Sc_UH)/(Sa_RL+Sa_RH+Sc_RL+Sc_RH), # Relative urban/rural in scTB in all
        Ya = (Ca_UL+Ca_UH)/(Ca_RL+Ca_RH), # Relative urban/rural in cTB (rC_UR) in adults
        Yc = (Cc_UL+Cc_UH)/(Cc_RL+Cc_RH), # Relative urban/rural in cTB (rC_UR) in children
        Yt = (Ca_UL+Ca_UH+Cc_UL+Cc_UH)/(Ca_RL+Ca_RH+Cc_RL+Cc_RH), # Relative urban/rural in cTB in all
        Xa = (Sa_RH+Sa_UH)/(Sa_RL+Sa_UL), # Relative high/low SES in scTB (rS_HL) in adults
        Xc = (Sc_RH+Sc_UH)/(Sc_RL+Sc_UL), # Relative high/low SES in scTB (rS_HL) in children
        Xt = (Sa_RH+Sa_UH+Sc_RH+Sc_UH)/(Sa_RL+Sa_UL+Sc_RL+Sc_UL), # Relative high/low SES in scTB in all
        Za = (Ca_RH+Ca_UH)/(Ca_RL+Ca_UL), # Relative high/low SES in cTB (rC_HL) in adults
        Zc = (Cc_RH+Cc_UH)/(Cc_RL+Cc_UL), # Relative high/low SES in cTB (rC_HL) in children
        Zt = (Ca_RH+Ca_UH+Cc_RH+Cc_UH)/(Ca_RL+Ca_UL+Cc_RL+Cc_UL), # Relative high/low SES in cTB in all
        Qs = ((Sa_RL+Sa_RH+Sa_UL+Sa_UH)/(Sc_RL+Sc_RH+Sc_UL+Sc_UH)), # Adult-to-children ratio in scTB
        Qc = ((Ca_RL+Ca_RH+Ca_UL+Ca_UH)/(Cc_RL+Cc_RH+Cc_UL+Cc_UH)), # Adult-to-children ratio in cTB
        Qt = ((Ca_RL+Ca_RH+Ca_UL+Ca_UH+Sa_RL+Sa_RH+Sa_UL+Sa_UH)/(Cc_RL+Cc_RH+Cc_UL+Cc_UH+Sc_RL+Sc_RH+Sc_UL+Sc_UH)))) # Adult-to-children ratio in all TB
    })
  }
  
  yini <- c(Na_RL = 100000-1000, Na_RH = 0, Na_UL = 0, Na_UH = 0, Nc_RL = 0, Nc_RH = 0, Nc_UL = 0, Nc_UH = 0,
            Ia_RL = 0, Ia_RH = 0, Ia_UL = 0, Ia_UH = 0, Ic_RL = 0, Ic_RH = 0, Ic_UL = 0, Ic_UH = 0,
            Oa_RL = 0, Oa_RH = 0, Oa_UL = 0, Oa_UH = 0, Oc_RL = 0, Oc_RH = 0, Oc_UL = 0, Oc_UH = 0,
            Ma_RL = 0, Ma_RH = 0, Ma_UL = 0, Ma_UH = 0, Mc_RL = 0, Mc_RH = 0, Mc_UL = 0, Mc_UH = 0,
            RMa_RL = 0, RMa_RH = 0, RMa_UL = 0, RMa_UH = 0, RMc_RL = 0, RMc_RH = 0, RMc_UL = 0, RMc_UH = 0,
            Sa_RL = 0, Sa_RH = 0, Sa_UL = 0, Sa_UH = 0, Sc_RL = 0, Sc_RH = 0, Sc_UL = 0, Sc_UH = 0,
            RSa_RL = 0, RSa_RH = 0, RSa_UL = 0, RSa_UH = 0, RSc_RL = 0, RSc_RH = 0, RSc_UL = 0, RSc_UH = 0,
            Ca_RL = 1000, Ca_RH = 0, Ca_UL = 0, Ca_UH = 0, Cc_RL = 0, Cc_RH = 0, Cc_UL = 0, Cc_UH = 0,
            RCa_RL = 0, RCa_RH = 0, RCa_UL = 0, RCa_UH = 0, RCc_RL = 0, RCc_RH = 0, RCc_UL = 0, RCc_UH = 0,
            Pa_RL = 0, Pa_RH = 0, Pa_UL = 0, Pa_UH = 0, Pc_RL = 0, Pc_RH = 0, Pc_UL = 0, Pc_UH = 0)

  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

# Quick diagnostic for appropriate model run (Pop should always sum to 100k)
test <- as.data.frame(ode(parms)) %>% mutate(Pop = Pa + Pc)

# 3. HMER ==========
# 3.1 HMER-format results
hmer_res <- function(params, times, outputs) { 
  t_max = max(times) # Define max times
  all_res = ode(params, t_max) # Run model until max times 
  actual_res = all_res[all_res[,'time'] %in% times, c('time', outputs)] # Subset desired time and outputs
  shaped = reshape2::melt(actual_res[,outputs]) # Switch to long format
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = ""))) # Set for HMER-format
}

# 3.2 Fitting targets
unc <- 0.1 # 10% uncertainty

targets <- list(
  Ta2008 = list(val = 199, sigma = unc*199),   # TB prevalence rate 1st survey in adults (Nguyen et al. Emerg Infect Dis 2021)
  Ta2018 = list(val = 125, sigma = unc*125),   # TB prevalence rate 2nd survey in adults (Nguyen et al. Emerg Infect Dis 2021)
  Ut2000 = list(val = 40.5, sigma = unc*40.5), # TB mortality 2000 in total population (WHO + UN WPP)
  Ut2010 = list(val = 25.2, sigma = unc*25.2), # TB mortality 2010 in total population (WHO + UN WPP)
  Vt2010 = list(val = 59.7, sigma = unc*59.7), # Notifications 2010 in total population (WHO + UN WPP)
  Vt2020 = list(val = 56.2, sigma = unc*56.2), # Notifications 2020 in total population (WHO + UN WPP)
  Ra2008 = list(val = 0.70, sigma = unc*0.70), # Proportion scTB 2008 in adults (Emery et al. medRxiv 2022)
  Ra2018 = list(val = 0.66, sigma = unc*0.66), # Proportion scTB 2018 in adults (Emery et al. medRxiv 2022)
  Wa2008 = list(val = 0.24, sigma = unc*0.24), # Relative urban/rural in scTB 2008 (Surveys dataset)
  Wa2018 = list(val = 0.44, sigma = unc*0.44), # Relative urban/rural in scTB 2018 (Surveys dataset)
  Ya2008 = list(val = 0.28, sigma = unc*0.28), # Relative urban/rural in cTB 2008 (Surveys dataset)
  Ya2018 = list(val = 0.43, sigma = unc*0.43), # Relative urban/rural in cTB 2018 (Surveys dataset)
  Xa2008 = list(val = 0.38, sigma = unc*0.38), # Relative high/low SES in scTB 2008 (Surveys dataset)
  Xa2018 = list(val = 0.47, sigma = unc*0.47), # Relative high/low SES in scTB 2018 (Surveys dataset)
  Za2008 = list(val = 0.32, sigma = unc*0.32), # Relative high/low SES in cTB 2008 (Surveys dataset)
  Za2018 = list(val = 0.46, sigma = unc*0.46)) # Relative high/low SES in cTB 2018 (Surveys dataset)

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,2) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),3,6)) # Create time variable
targetsdb <- as.data.frame(cbind(targetsdb[grep("val", rownames(targetsdb)),c("time","var","V1")],V2 = targetsdb[grep("sigma", rownames(targetsdb)),c("V1")])) # Widen dataframe
colnames(targetsdb)[c(3,4)] <- c("val","sigma") # Rename columns
targetsdb$lo <- targetsdb$val-2*targetsdb$sigma # Lower bound
targetsdb$hi <- targetsdb$val+2*targetsdb$sigma # Upper bound
rownames(targetsdb) <- gsub(rownames(targetsdb),pattern = "\\.val", replacement = "")
export(targetsdb,here("data","targets.Rdata")) # Save data frame

# 3.3 Target check function
target_check <- function(wave = NULL, nwave = w, targets = targetsdb){
  # Prepare dataframes
  wave <- as.data.frame(wave)
  target <- rownames(targets)
  listdf <- list()
  
  # Perform target range check
  for (i in 1:length(target)){
    targ <- target[i]
    targval <- wave[, grep(colnames(wave), pattern = targ)]
    targlim <- targets[grep(rownames(targets), pattern = targ), c('lo','hi')]
    check <- targval >= targlim$lo & targval <= targlim$hi
    df <- data.frame("check" = check)
    colnames(df) <- targ
    listdf[[i]] <- df
  }
  targetchecks <- do.call("cbind",listdf)
  
  # Run check count
  runcounts <- rowSums(targetchecks)
  runcountsdf <- data.frame(freq = 0:(length(target)))
  
  for (i in 0:(length(target))) {
    runcountsdf$freq[i + 1] <- sum(runcounts == i)
  }
  
  rownames(runcountsdf) <- 0:(length(target))
  colnames(runcountsdf) <- paste0("wave",nwave)
  
  # Target check count
  targcounts <- colSums(targetchecks)
  targcountsdf <- data.frame(targcounts)
  colnames(targcountsdf) <- paste0("wave",nwave)
  
  # List output
  checkoutput <- list(runcountsdf, targcountsdf)
  names(checkoutput) <- c("runs","targets")
  return(checkoutput)
}

# 3.4 HMER initial run
ems <- list() # Empty list for emulators per wave
wave <- list() # Empty list for data used to train and validate per wave
checks <- list() # Empty list for target checks
wave_samp <- list() # Empty list for wave samples
wave_train <- list() # Empty list for wave training
wave_val <- list() # Empty list for wave validation
wave_res <- list() # Empty list for wave results
invalid <- list() # Empty list for invalid parameter sets
invalid_diag <- list() # Empty list for invalid parameter sets after diagnostics
non_imp_pts <- list() # Empty list for non-implausible points generated per wave

n_params <- length(ranges) # Number of parameters
n_points <- 10*n_params # HMER points = selection of parameters

ini_LHS_train <- lhs::maximinLHS(n_points, n_params)
ini_LHS_val <- lhs::maximinLHS(n_points, n_params)
ini_LHS <- rbind(ini_LHS_train, ini_LHS_val)
rm(ini_LHS_train, ini_LHS_val) # Clean objects

ini_pts <- setNames(data.frame(t(apply(ini_LHS, 1, function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + unlist(lapply(ranges, function(x) x[1]))))), names(ranges)) # Set random sets to create points
wave_res[[1]] <- data.frame(t(apply(ini_pts, 1, hmer_res, c(2000, 2008, 2010, 2018, 2020), c('Ta', 'Ut', 'Vt', 'Ra', 'Wa', 'Ya', 'Xa', 'Za'))))[,names(targets)] # Run ODE
checks[[1]] <- target_check(wave = wave_res[[1]], nwave = 1) # Perform target checks
wave[[1]] <- cbind(ini_pts, wave_res[[1]]) # Bind run points and results
rm(ini_pts, ini_LHS) # Clean objects

wave_samp[[1]] <- sample(1:nrow(wave[[1]]), n_points/2) # Sample half of runs
wave_train[[1]] <- wave[[1]][wave_samp[[1]],] # Designate half for emulator training
wave_val[[1]] <- wave[[1]][-wave_samp[[1]],] # Designate half for validation/diagnostics

ems[[1]] <- emulator_from_data(wave_train[[1]], names(targets), ranges) # Ranges for emulator

pdf(here("outputs", "val_w1.pdf"))
invalid[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()

for (j in 1:length(ems[[1]])) {
  misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  while(misclass > 0) {
    ems[[1]][[j]] <- ems[[1]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems[[1]])) {
  bad.model <- nrow(comparison_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  if (bad.model > floor(nrow(wave_val[[1]])/10)) {
    bad.ems <- c(bad.ems, j)
  }
}

ems[[1]] <- ems[[1]][!seq_along(ems[[1]]) %in% bad.ems]

pdf(here("outputs", "val_w1-diag.pdf"))
invalid_diag[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()

non_imp_pts[[1]] <- generate_new_design(ems[[1]], n_points, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[1]],here("outputs","pts_w1.Rdata")) # Save data frame
beepr::beep(2)

# 3.5 HMER loop runs
w <- 2 # Update wave run
wave_res[[w]] <- data.frame(t(apply(non_imp_pts[[w-1]], 1, hmer_res, c(2000, 2008, 2010, 2018, 2020), c('Ta', 'Ut', 'Vt', 'Ra', 'Wa', 'Ya', 'Xa', 'Za'))))[,names(targets)] # Run ODE
checks[[w]] <- target_check(wave = wave_res[[w]]) # Perform target checks
wave[[w]] <- cbind(non_imp_pts[[w-1]], wave_res[[w]])

wave_samp[[w]] <- sample(1:nrow(wave[[w]]), n_points/2) # Sample half of runs
wave_train[[w]] <- wave[[w]][wave_samp[[w]],] # Designate half for emulator training
wave_val[[w]] <- wave[[w]][-wave_samp[[w]],] # Designate half for validation/diagnostics

ems[[w]] <- emulator_from_data(wave_train[[w]], names(targets), ranges, check.ranges = TRUE) # Ranges for emulator

pdf(here("outputs", paste("val_w", w,".pdf", sep = "")))
invalid[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()

for (j in 1:length(ems[[w]])) {
  misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  while(misclass > 0) {
    ems[[w]][[j]] <- ems[[w]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems[[w]])) {
  bad.model <- nrow(comparison_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = TRUE))
  if (bad.model > floor(nrow(wave_val[[w]])/10)) {
    bad.ems <- c(bad.ems, j)
  }
}

ems[[w]] <- ems[[w]][!seq_along(ems[[w]]) %in% bad.ems]

pdf(here("outputs", paste("val_w", w,"-diag.pdf", sep = "")))
invalid_diag[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()

non_imp_pts[[w]] <- generate_new_design(c(rev(ems[1:w])), n_points, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[w]], here("outputs",paste("pts_w", w,".Rdata", sep = ""))) # Save data frame
beepr::beep(2)

runcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["runs"]]))
targetcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["targets"]]))

# 4. Results ==========
export(runcheck, here("outputs", "runcheck.Rdata")) # Save data frame
export(targetcheck, here("outputs", "targetcheck.Rdata")) # Save data frame

pts_fin <- import(here("outputs",paste("pts_w", w,".Rdata", sep = "")))

quants <- c(0.025,0.5,0.975) # Set quantiles

parameters <- apply(pts_fin, 2, quantile, probs = quants, na.rm = TRUE) # Set parameter quantiles
t_parameters <- data.table::transpose(as.data.frame(parameters)) # Transpose parameters
colnames(t_parameters) = rownames(parameters) # Set column names
rownames(t_parameters) = colnames(parameters) # Set row names
parameters <- t_parameters # Rename parameters
rm(t_parameters) # Clean objects

parameters$parameter <- c("beta","kappa","gamma_infcle","lambda_infmin","gamma_mincle","theta_cleinf",
                          "lambda_minsub","lambda_infsub","gamma_submin","lambda_subcln","gamma_clnsub",
                          "omega_ini","omega_fin", "iota_cln_ini", "iota_cln_fin", "phi_cln_ini", "phi_cln_fin",
                          "theta_recinf","tau_min","tau_sub","rho_ini","rho_fin","sigma_ini","sigma_fin","ratio_ca")
parameters[,c(1,2,3)] <- round(parameters[,c(1,2,3)],2) # Round to 2 decimal places
table <- data.frame(parameter = parameters$parameter, low  = parameters$`2.5%`, med = parameters$`50%`, hig = parameters$`97.5%`) # Output table
export(table,here("data","parameters.Rdata")) # Save data frame

results <- as.data.frame(apply(pts_fin, 1, ode))[-seq(1,521),] # Runs ODE for each set of points
results <- as.data.frame(t(apply(results, 1, quantile, probs = quants, na.rm = TRUE))) # Set parameter quantiles
comp <- colnames(ode(parms))[-1] # Compartment names (-time)
results$var <- c(t(replicate(521,comp))) # Compartment variable name
results$time <- rep(seq(1500,2020), length(comp)) # Time variable

export(results, here("outputs", paste("results_w", w,".Rdata", sep = ""))) # Save data frame
beepr::beep(2)
