## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 02: Fit.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(hmer) # Emulation and History Matching
library(lhs) # Latin Hypercube Samples
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(purrr) # Complete set of tools for functions and vectors
library(data.table) # Faster than data.frame
library(beepr) # Sound cues when code has finished running
library(progress) # Displays progress bar
library(tictoc) # Timing of R scripts

# 1. Parameters ==========
# 1.1 Baseline parameters (for testing)
parms = c(
  beta = 8,             # Contact (per person/year) parameter
  kappa = 0.75,         # Relative infectiousness
  gamma_infcle = 1.83,  # REG: Infected -> Cleared
  lambda_infmin = 0.21, # PROG: Infected -> Minimal
  gamma_mincle = 0.16,  # REG: Minimal -> Recovered
  theta_cleinf = 0.88,  # REINF: Cleared -> Infected
  lambda_minsub = 0.25, # PROG: Minimal -> Subclinical
  lambda_infsub = 0.07, # PROG: Infected -> Subclinical
  gamma_submin = 1.58,  # REG: Subclinical -> Minimal
  lambda_subcln = 0.77, # PROG: Subclinical -> Clinical
  gamma_clnsub = 0.53,  # REG: Clinical -> Subclinical
  omega_ini = 0.3,      # Mortality (Initial)
  omega_fin = 0.23,     # Mortality (Final)
  iota_cln_ini = 0.44,  # Diagnosis Clinical (Initial)
  iota_cln_fin = 0.9,   # Diagnosis Clinical (Final)
  phi_cln_ini = 0.69,   # Treatment failure Clinical (Initial) 
  phi_cln_fin = 0.09,   # Treatment failure Clinical (Final)
  tau_min = 0.005,      # RELAP: Treated -> Minimal 
  tau_sub = 0.005,      # RELAP: Treated -> Subclinical
  rho_ini = 0.78,       # Proportion rural (Initial)
  rho_fin = 0.04,       # Proportion rural (Final)
  sigma_ini = 0.76,     # Proportion low SES (Initial)
  sigma_fin = 0.16)     # Proportion low SES (Final)

# 1.2 Parameter ranges
ranges = list(
  beta = c(0,20),                 # Contact (per person/year) parameter (Horton et al. 2022 - PLOS GPH)
  kappa = c(0.5,1),               # Relative infectiousness (Emery et al. 2022)
  gamma_infcle = c(0.93,3.30),    # REG: Infected -> Cleared (Horton et al. 2023)
  lambda_infmin = c(0.04,0.23),   # PROG: Infected -> Minimal (Horton et al. 2023)
  gamma_mincle = c(0.14,0.23),    # REG: Minimal -> Recovered (Horton et al. 2023)
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
  tau_min = c(0.0,0.01),          # RELAP: Treated -> Minimal 
  tau_sub = c(0.0,0.01),          # RELAP: Treated -> Subclinical
  rho_ini = c(0,1),               # Proportion rural (Initial)
  rho_fin = c(0,1),               # Proportion rural (Final)
  sigma_ini = c(0,1),             # Proportion low SES (Inital)
  sigma_fin = c(0,1))             # Proportion low SES (Final)

# 2. Model ==========
ode <- function(parms, end_time = 2020) {
  
  # Static parameters  
  N <- 1e5 # Population size
  mu <- 1/70 # Age expectancy adult
  delta <- 1 # Treatment year
  theta_recinf <- 1 # REINF: Recovered -> Infected (No protection)
  theta_treinf <- 1 # REINF: Treated -> Infected (No protection)

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
      dI_RL  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(W_RL*theta_cleinf)+(O_RL*theta_recinf)+(P_RL*theta_treinf))) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - mu*I_RL
      dW_RL  = gamma_infcle*I_RL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_RL*theta_cleinf) - mu*W_RL
      dO_RL  = gamma_mincle*M_RL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_recinf) - mu*O_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL + tau_min*P_RL - mu*M_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL + tau_sub*P_RL - mu*S_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - force_func_omega(time)*C_RL - mu*C_RL - force_func_iota(time)*C_RL + force_func_phi(time)*RC_RL
      dRC_RL = force_func_iota(time)*C_RL - force_func_phi(time)*RC_RL - delta*RC_RL - mu*RC_RL
      dP_RL  = delta*RC_RL - tau_min*P_RL - tau_sub*P_RL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_treinf) - mu*P_RL
      
      dN_RH  = force_func_rhoR(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu*N_RH
      dI_RH  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(W_RH*theta_cleinf)+(O_RH*theta_recinf)+(P_RH*theta_treinf))) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu*I_RH
      dW_RH  = gamma_infcle*I_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_RH*theta_cleinf) - mu*W_RH
      dO_RH  = gamma_mincle*M_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_recinf) - mu*O_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH + tau_min*P_RH - mu*M_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH + tau_sub*P_RH - mu*S_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - force_func_omega(time)*C_RH - mu*C_RH - force_func_iota(time)*C_RH + force_func_phi(time)*RC_RH
      dRC_RH = force_func_iota(time)*C_RH - force_func_phi(time)*RC_RH - delta*RC_RH - mu*RC_RH
      dP_RH  = delta*RC_RH - tau_min*P_RH - tau_sub*P_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_treinf) - mu*P_RH
      
      dN_UL  = force_func_rhoU(time)*force_func_sigmaL(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu*N_UL
      dI_UL  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(W_UL*theta_cleinf)+(O_UL*theta_recinf)+(P_UL*theta_treinf))) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu*I_UL
      dW_UL  = gamma_infcle*I_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_UL*theta_cleinf) - mu*W_UL
      dO_UL  = gamma_mincle*M_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_recinf) - mu*O_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL + tau_min*P_UL - mu*M_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL + tau_sub*P_UL - mu*S_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - force_func_omega(time)*C_UL - mu*C_UL - force_func_iota(time)*C_UL + force_func_phi(time)*RC_UL
      dRC_UL = force_func_iota(time)*C_UL - force_func_phi(time)*RC_UL - delta*RC_UL - mu*RC_UL
      dP_UL  = delta*RC_UL - tau_min*P_UL - tau_sub*P_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_treinf) - mu*P_UL
      
      dN_UH  = force_func_rhoU(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu*N_UH
      dI_UH  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(W_UH*theta_cleinf)+(O_UH*theta_recinf)+(P_UH*theta_treinf))) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu*I_UH
      dW_UH  = gamma_infcle*I_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_UH*theta_cleinf) - mu*W_UH
      dO_UH  = gamma_mincle*M_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_recinf) - mu*O_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH + tau_min*P_UH - mu*M_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH + tau_sub*P_UH - mu*S_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - force_func_omega(time)*C_UH - mu*C_UH - force_func_iota(time)*C_UH + force_func_phi(time)*RC_UH
      dRC_UH = force_func_iota(time)*C_UH - force_func_phi(time)*RC_UH - delta*RC_UH - mu*RC_UH
      dP_UH  = delta*RC_UH - tau_min*P_UH - tau_sub*P_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_treinf) - mu*P_UH
      
      return(list(c( 
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dW_RL, dW_RH, dW_UL, dW_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH,
        dS_RL, dS_RH, dS_UL, dS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dP_RL, dP_RH, dP_UL, dP_UH),
        TBc   = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # All TB (per 100k)
        Mor   = (force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)), # Clinical TB mortality per time (per 100k)
        Dxs   = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH)), # Notifications cTB per time in adults (per 100k)
        Spr   = (S_RL+S_RH+S_UL+S_UH)/(S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH), # Proportion scTB
        URs   = (S_UL+S_UH)/(S_RL+S_RH), # Relative urban/rural in scTB
        URc   = (C_UL+C_UH)/(C_RL+C_RH), # Relative urban/rural in cTB
        HLs   = (S_RH+S_UH)/(S_RL+S_UL), # Relative high/low SES in scTB
        HLc   = (C_RH+C_UH)/(C_RL+C_UL))) # Relative high/low SES in cTB
    })
  }
  
  yini <- c(N_RL = 24750, N_RH = 24750, N_UL = 24750, N_UH = 24750, I_RL = 0, I_RH = 0, I_UL = 0, I_UH = 0, 
            W_RL = 0, W_RH = 0, W_UL = 0, W_UH = 0, O_RL = 0, O_RH = 0, O_UL = 0, O_UH = 0, 
            M_RL = 0, M_RH = 0, M_UL = 0, M_UH = 0, S_RL = 0, S_RH = 0, S_UL = 0, S_UH = 0, 
            C_RL = 250, C_RH = 250, C_UL = 250, C_UH = 250, RC_RL = 0, RC_RH = 0, RC_UL = 0, RC_UH = 0, 
            P_RL = 0, P_RH = 0, P_UL = 0, P_UH = 0)
  
  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

# Quick diagnostic
test <- as.data.frame(ode(parms))

# 3. HMER ==========
# 3.1 HMER-format results
hmer_res <- function(params, times, outputs) { 
  t_max = max(times) # Define max times
  all_res = ode(params, end_time = t_max) # Run model until max times 
  actual_res = all_res[all_res[,'time'] %in% times, c('time', outputs)] # Subset desired time and outputs
  shaped = reshape2::melt(actual_res[,outputs]) # Switch to long format
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = ""))) # Set for HMER-format
}

# 3.2 Fitting targets
targets <- list(
  TBc2008 = c(202,310),     # TB prevalence (Nguyen et al. Emerg Infect Dis 2021 - Revised)
  TBc2018 = c(177,290),     # TB prevalence (Nguyen et al. Emerg Infect Dis 2021 - Revised)
  Mor2000 = c(37.3,87.7),   # TB mortality (Over population aged >= 15yo)
  Mor2010 = c(22.8,45.7),   # TB mortality (Over population aged >= 15yo)
  Dxs2010 = c(63.5,95.3),   # TB notifications (Over population aged >= 15yo)
  Dxs2020 = c(58.5,87.8),   # TB notifications (Over population aged >= 15yo)
  URs2008 = c(0.192,0.288), # Urban/rural scTB
  URs2018 = c(0.352,0.528), # Urban/rural scTB
  URc2008 = c(0.224,0.336), # Urban/rural cTB
  URc2018 = c(0.344,0.516), # Urban/rural cTB
  HLs2008 = c(0.304,0.456), # High/low scTB
  HLs2018 = c(0.376,0.564), # High/low scTB
  HLc2008 = c(0.256,0.384), # High/low cTB
  HLc2018 = c(0.36,0.54),   # High/low cTB
  Spr2008 = c(0.56,0.84),   # Proportion scTB (Emery et al. medRxiv 2022)
  Spr2018 = c(0.53,0.79))   # Proportion scTB (Emery et al. medRxiv 2022)

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,3) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),4,7)) # Create time variable
colnames(targetsdb) <- c("lo","hi","var","time") # Rename columns

export(targetsdb,here("data","targets.Rdata")) # Save data frame

# 3.3 Target functions
# 3.3.1 Target check function
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

# 3.3.2 Target subset function
target_subset <- function(wave = NULL, nwave = w, sens = 0, targets = targetsdb, parms = parms){
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
  targetchecks <- as.data.frame(do.call("cbind",listdf))
  targetchecks <- as.data.frame(rowSums(targetchecks))
  colnames(targetchecks)[length(targetchecks)] <- "check"
  wave_check <- cbind(wave,targetchecks)
  wave_check <- wave_check[wave_check$check >= (length(target)-0),]
  wave_check <- wave_check[,1:length(wave_check)-1]
  return(wave_check)
}

# 3.4 HMER initial run
ems <- list() # Empty list for emulators per wave
wave <- list() # Empty list for data used to train and validate per wave
checks <- list() # Empty list for target checks
wave_train <- list() # Empty list for wave training
wave_val <- list() # Empty list for wave validation
wave_res <- list() # Empty list for wave results
wave_check <- list() # Empty list for wave checks
invalid_pts <- list() # Empty list for invalid parameter sets
invalid_diag <- list() # Empty list for invalid parameter sets after diagnostics
invalid_bad <- list() # Empty list for invalid parameter sets after diagnostics and removal of bad emulators
non_imp_pts <- list() # Empty list for non-implausible points generated per wave

ini_LHS_train <- lhs::maximinLHS(10*length(ranges), length(ranges))
ini_LHS_val <- lhs::maximinLHS(10*length(ranges), length(ranges))
ini_LHS <- rbind(ini_LHS_train, ini_LHS_val)
ini_pts <- setNames(data.frame(t(apply(ini_LHS, 1, function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + unlist(lapply(ranges, function(x) x[1]))))), names(ranges)) # Set random sets to create points
rm(ini_LHS_train, ini_LHS_val, ini_LHS) # Clean objects

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(ini_pts))
tmp <- list()
for (i in seq_len(nrow(ini_pts))) {
  res <- t(apply(ini_pts[i,], 1, hmer_res, c(2000, 2008, 2010, 2018, 2020), c('TBc', 'Mor', 'Dxs', 'Spr', 'URs', 'URc', 'HLs', 'HLc')))
  tmp[[i]] <- data.frame(res)[, names(targets)]
  pb$tick()  # Advance the progress bar
}
wave_res[[1]] <- do.call(rbind, tmp)
wave[[1]] <- cbind(ini_pts, wave_res[[1]]) # Bind run points and results
rm(pb, tmp, res, i) # Clean objects
print(pdf(here("outputs","simulator","w1_simulator.pdf")))
simulator_plot(wave_res, targets, normalize = TRUE, byhit = TRUE)
dev.off()
cat("Model runs completed\n")

checks[[1]] <- target_check(wave = wave_res[[1]], nwave = 1) # Perform target checks
runcheck <- checks[[1]][["runs"]]
targetcheck <- checks[[1]][["targets"]]
wave_check[[1]] <- target_subset(wave[[1]])
rm(ini_pts) # Clean objects
cat("Runs checked\n")

sample <- sample(1:nrow(wave[[1]]), round(length(wave[[1]][,1])/2))
wave_train[[1]] <- wave[[1]][sample,] 
wave_val[[1]] <- wave[[1]][-sample,]
rm(sample)
ems[[1]] <- emulator_from_data(wave_train[[1]], names(targets), ranges)
print(pdf(here("outputs", "activeparms", "w1_activeparms.pdf")))
plot_actives(ems[[1]])
dev.off()
cat("Emulators trained\n")

pdf(here("outputs", "diagnostics", "w1_diagnostics_pre.pdf"))
invalid_pts[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()
for (j in 1:length(ems[[1]])) {
  misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  while(misclass > 0) {
    ems[[1]][[j]] <- ems[[1]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  }
}
pdf(here("outputs", "diagnostics", "w1_diagnostics_post.pdf"))
invalid_diag[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()
bad.ems <- c()
for (j in 1:length(ems[[1]])) {
  bad.model <- nrow(comparison_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  if (bad.model > floor(nrow(wave_val[[1]])/10)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems[[1]] <- ems[[1]][!seq_along(ems[[1]]) %in% bad.ems]
cat("Available emulators:", length(ems[[1]]), "\n")
cat(names(ems[[1]]),"\n")
pdf(here("outputs", "diagnostics", "w1_diagnostics_post_badems.pdf"))
invalid_bad[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()
cat("Diagnostics performed\n")

non_imp_pts[[1]] <- generate_new_design(ems[[1]], (10*length(ranges))*2, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[1]],here("outputs", "pts", "w1_pts.Rdata")) # Save data frame
cat("New parameter points generated\n")
beepr::beep(2)

# 3.5 HMER loop runs
for (i in 13:40){
w <- i # Update wave run
tic()
cat("Running wave:", w, "\n")

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(non_imp_pts[[w-1]]))
tmp <- list()
for (i in seq_len(nrow(non_imp_pts[[w-1]]))) {
  res <- t(apply(non_imp_pts[[w-1]][i,], 1, hmer_res, c(2000, 2008, 2010, 2018, 2020), c('TBc', 'Mor', 'Dxs', 'Spr', 'URs', 'URc', 'HLs', 'HLc')))
  tmp[[i]] <- data.frame(res)[, names(targets)]
  pb$tick()  # Advance progress bar
}
wave_res[[w]] <- do.call(rbind, tmp)
wave[[w]] <- cbind(non_imp_pts[[w-1]], wave_res[[w]])
rm(pb, tmp, res, i) # Clean objects
pdf(here("outputs", "simulator", paste("w", w, "_simulator.pdf", sep = "")))
simulator_plot(wave_res[w], targets, normalize = TRUE, byhit = TRUE)
dev.off()
cat("Model runs completed\n")

checks[[w]] <- target_check(wave = wave_res[[w]]) # Perform target checks
runcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["runs"]]))
targetcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["targets"]]))
export(runcheck, here("outputs", "simulator", "runcheck.Rdata")) # Save data frame
export(targetcheck, here("outputs", "simulator", "targetcheck.Rdata")) # Save data frame
wave_check[[w]] <- target_subset(wave[[w]])
cat("Runs checked\n")

sample <- sample(1:nrow(wave[[w]]), round(length(wave[[w]][,1])/2))
wave_train[[w]] <- wave[[w]][sample,] 
wave_val[[w]] <- wave[[w]][-sample,]
rm(sample)
ems[[w]] <- emulator_from_data(wave_train[[w]], names(targets), ranges, check.ranges = TRUE)
pdf(here("outputs", "activeparms", paste("w", w, "_activeparms.pdf", sep = "")))
plot_actives(ems[[w]])
dev.off()
cat("Emulators trained\n")

pdf(here("outputs", "diagnostics", paste("w", w, "_diagnostics_pre.pdf", sep = "")))
invalid_pts[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()
for (j in 1:length(ems[[w]])) {
  misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  while(misclass > 0) {
    ems[[w]][[j]] <- ems[[w]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  }
}
pdf(here("outputs", "diagnostics", paste("w", w, "_diagnostics_post.pdf", sep = "")))
invalid_diag[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()
bad.ems <- c()
for (j in 1:length(ems[[w]])) {
  bad.model <- nrow(comparison_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  if (bad.model > floor(nrow(wave_val[[w]])/10)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems[[w]] <- ems[[w]][!seq_along(ems[[w]]) %in% bad.ems]
cat("Available emulators:", length(ems[[w]]), "\n")
cat(names(ems[[w]]),"\n")
pdf(here("outputs", "diagnostics", paste("w", w, "_diagnostics_post_badems.pdf", sep = "")))
invalid_bad[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()
export(ems[[w]], here("outputs", "ems", paste("w", w, "_ems.Rdata", sep = "")))
cat("Diagnostics performed\n")

non_imp_pts[[w]] <- generate_new_design(c(ems[1:w]), (10*length(ranges))*2, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[w]],here("outputs", "pts", paste("w", w, "_pts.Rdata", sep = ""))) # Save data frame
cat("New parameter points generated\n")
toc()
}
beepr::beep(2)

# 3.6 Reload data
# Non-implausible points
data_files <- list.files(here("outputs","pts"), pattern = "w[0-9]+_pts.Rdata", full.names = TRUE)
for (file in data_files) {
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  data <- import(file)
  non_imp_pts[[num]] <- data
}

# Waves
data_files <- list.files(here("outputs","waves"), pattern = "w[0-9]+_wave.Rdata", full.names = TRUE)
for (file in data_files) {
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  data <- import(file)
  wave[[num]] <- data
}

# 4. Results ==========
pts_fin <- non_imp_pts[[w]]

# Isolate best parameter sets
plaus_pts <- as.data.frame(do.call("rbind", wave_check))[1:length(parms)]
export(plaus_pts,here("outputs", "pts", "plaus_pts.Rdata")) # Save data frame

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
                          "tau_min","tau_sub","rho_ini","rho_fin","sigma_ini","sigma_fin")
parameters[,c(1,2,3)] <- round(parameters[,c(1,2,3)],2) # Round to 2 decimal places
table <- data.frame(parameter = parameters$parameter, low  = parameters$`2.5%`, med = parameters$`50%`, hig = parameters$`97.5%`) # Output table
export(table,here("data","parameters.Rdata")) # Save data frame

results <- as.data.frame(apply(pts_fin, 1, ode))[-seq(1,521),] # Runs ODE for each set of points
results <- as.data.frame(t(apply(results, 1, quantile, probs = quants, na.rm = TRUE))) # Set parameter quantiles
comp <- colnames(ode(parms))[-1] # Compartment names (-time)
results$var <- c(t(replicate(521,comp))) # Compartment variable name
results$time <- rep(seq(1500,2020), length(comp)) # Time variable

export(results, here("outputs", "results.Rdata")) # Save data frame
beepr::beep(2)
