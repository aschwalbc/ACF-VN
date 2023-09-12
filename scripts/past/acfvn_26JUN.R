## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 02: ACF-VN.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(hmer) # Emulation and History Matching
library(lhs) # Latin Hypercube Samples
library(deSolve) # Solvers for ordinary differential equations
library(ggplot2) # To build comparative plots
library(reshape2) # Reshaping data easily
library(purrr) # Complete set of tools for functions and vectors
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame, allows use of j operator (:=)
library(patchwork) # Plot composition
library(beepr) # Sound cues when code has finished running

# 1. Parameters ==========
# 1.1 Baseline parameters (for testing)
parms = c(
  beta = 100,           # Contact (per person/year) parameter
  alpha = 0.5,          # Relative infectiousness
  lambda_infmin = 0.08, # Prog: Infected -> Minimal
  lambda_infsub = 0.03, # Prog: Infected -> Subclinical
  lambda_minsub = 0.24, # Prog: Minimal -> Subclinical
  lambda_subcln = 0.72, # Prog: Subclinical -> Clinical
  gamma_infcle = 1.13,  # Reg: Infected -> Cleared
  gamma_mincle = 0.18,  # Reg: Minimal -> Cleared
  gamma_submin = 1.56,  # Reg: Subclinical -> Minimal
  gamma_clnsub = 0.58,  # Reg: Clinical -> Subclinical
  omega_ini = 0.33,     # Mortality (Initial)
  omega_fin = 0.33,     # Mortality (Final)
  iota_ini = 0.2,      # Diagnosis clinical (Initial)
  iota_fin = 0.4,       # Diagnosis clinical (Final)
  rho_ini = 0.7,        # Proportion rural (Initial)
  rho_fin = 0.5,        # Proportion rural (Final)
  sigma_ini = 0.8,      # Proportion low SES (Inital)
  sigma_fin = 0.2)      # Proportion low SES (Final)

# 1.2 Parameter ranges
ranges = list(
  beta = c(0,50),               # Contact (per person/year) parameter (Horton et al. 2022 - PLOS GPH)
  alpha = c(0.5,1),             # Relative infectiousness (Emery et al. 2022)
  lambda_infmin = c(0.03,0.15), # Prog: Infected -> Minimal (Horton et al. 2022)
  lambda_infsub = c(0.01,0.06), # Prog: Infected -> Subclinical (Horton et al. 2022)
  lambda_minsub = c(0.20,0.28), # Prog: Minimal -> Subclinical (Horton et al. 2022)
  lambda_subcln = c(0.57,0.94), # Prog: Subclinical -> Clinical (Horton et al. 2022)
  gamma_infcle = c(0.60,1.79),  # Reg: Infected -> Cleared (Horton et al. 2022)
  gamma_mincle = c(0.14,0.23),  # Reg: Minimal -> Cleared (Horton et al. 2022)
  gamma_submin = c(1.21,2.02),  # Reg: Subclinical -> Minimal (Horton et al. 2022)
  gamma_clnsub = c(0.46,0.73),  # Reg: Clinical -> Subclinical (Horton et al. 2022)
  omega_ini = c(0.26,0.37),     # Mortality (Initial) (Richards et al. 2023 - Lancet GH)
  omega_fin = c(0,1),           # Mortality (Final)
  iota_ini = c(0,1),            # Diagnosis clinical (Initial)
  iota_fin = c(0,1),            # Diagnosis clinical (Final)
  rho_ini = c(0,1),             # Proportion rural (Initial)
  rho_fin = c(0,1),             # Proportion rural (Final)
  sigma_ini = c(0,1),           # Proportion low SES (Inital)
  sigma_fin = c(0,1))           # Proportion low SES (Final)

# 2. Model ==========
ode <- function(parms, end_time = 2020) {
  
  P <- 100000 # Population (15/70* children)
  mu <- 1/70 # Age expectancy

  forcer_omega <- matrix(c(1500, parms['omega_ini'], 1999, parms['omega_ini'], 2020, parms['omega_fin']), ncol = 2, byrow = TRUE)
  force_func_omega <- approxfun(x = forcer_omega[,1], y = forcer_omega[,2], method = "linear", rule = 2)
  
  forcer_iota <- matrix(c(1500, parms['iota_ini'], 1999, parms['iota_ini'], 2020, parms['iota_fin']), ncol = 2, byrow = TRUE)
  force_func_iota <- approxfun(x = forcer_iota[,1], y = forcer_iota[,2], method = "linear", rule = 2)
  
  forcer_rho <- matrix(c(1500, parms['rho_ini'], 1999, parms['rho_ini'], 2020, parms['rho_fin']), ncol = 2, byrow = TRUE)
  force_func_rhoR <- approxfun(x = forcer_rho[,1], y = forcer_rho[,2], method = "linear", rule = 2)
  force_func_rhoU <- approxfun(x = forcer_rho[,1], y = 1-forcer_rho[,2], method = "linear", rule = 2)
  
  forcer_sigma <- matrix(c(1500, parms['sigma_ini'], 1999, parms['sigma_ini'], 2020, parms['sigma_fin']), ncol = 2, byrow = TRUE)
  force_func_sigmaL <- approxfun(x = forcer_sigma[,1], y = forcer_sigma[,2], method = "linear", rule = 2)
  force_func_sigmaH <- approxfun(x = forcer_sigma[,1], y = 1-forcer_sigma[,2], method = "linear", rule = 2)
  
  des <- function(time, state, parms) {
    
    with(as.list(c(state, parms)), {
      
      dN_RL = force_func_rhoR(time)*force_func_sigmaL(time)*(mu*P + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RL/P) - mu*N_RL
      dI_RL = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RL/P) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - mu*I_RL
      dO_RL = gamma_infcle*I_RL + gamma_mincle*M_RL - mu*O_RL    
      dM_RL = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - mu*M_RL
      dS_RL = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - mu*S_RL
      dC_RL = lambda_subcln*S_RL - gamma_clnsub*C_RL - force_func_omega(time)*C_RL - mu*C_RL - force_func_iota(time)*C_RL
      dD_RL = force_func_iota(time)*C_RL - mu*D_RL
      
      dN_RH = force_func_rhoR(time)*force_func_sigmaH(time)*(mu*P + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RH/P) - mu*N_RH
      dI_RH = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RH/P) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu*I_RH
      dO_RH = gamma_infcle*I_RH + gamma_mincle*M_RH - mu*O_RH    
      dM_RH = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH - mu*M_RH
      dS_RH = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH - mu*S_RH
      dC_RH = lambda_subcln*S_RH - gamma_clnsub*C_RH - force_func_omega(time)*C_RH - mu*C_RH - force_func_iota(time)*C_RH
      dD_RH = force_func_iota(time)*C_RH - mu*D_RH
      
      dN_UL = force_func_rhoU(time)*force_func_sigmaL(time)*(mu*P + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UL/P) - mu*N_UL
      dI_UL = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UL/P) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu*I_UL
      dO_UL = gamma_infcle*I_UL + gamma_mincle*M_UL - mu*O_UL    
      dM_UL = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL - mu*M_UL
      dS_UL = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL - mu*S_UL
      dC_UL = lambda_subcln*S_UL - gamma_clnsub*C_UL - force_func_omega(time)*C_UL - mu*C_UL - force_func_iota(time)*C_UL
      dD_UL = force_func_iota(time)*C_UL - mu*D_UL
      
      dN_UH = force_func_rhoU(time)*force_func_sigmaH(time)*(mu*P + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UH/P) - mu*N_UH
      dI_UH = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UH/P) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu*I_UH
      dO_UH = gamma_infcle*I_UH + gamma_mincle*M_UH - mu*O_UH    
      dM_UH = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH - mu*M_UH
      dS_UH = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH - mu*S_UH
      dC_UH = lambda_subcln*S_UH - gamma_clnsub*C_UH - force_func_omega(time)*C_UH - mu*C_UH - force_func_iota(time)*C_UH
      dD_UH = force_func_iota(time)*C_UH - mu*D_UH
      
      return(list(c(
        dN_RL, dN_RH, dN_UL, dN_UH,  # Susceptible
        dI_RL, dI_RH, dI_UL, dI_UH,  # Infected
        dO_RL, dO_RH, dO_UL, dO_UH,  # Cleared
        dM_RL, dM_RH, dM_UL, dM_UH,  # Minimal
        dS_RL, dS_RH, dS_UL, dS_UH,  # Subclinical
        dC_RL, dC_RH, dC_UL, dC_UH,  # Clinical
        dD_RL, dD_RH, dD_UL, dD_UH), # Diagnosed/treated
        W = (S_UL + S_UH)/(S_RL + S_RH), # Relative urban/rural in scTB (rS_UR)
        Y = (C_UL + C_UH)/(C_RL + C_RH), # Relative urban/rural in cTB (rC_UR)
        X = (S_RH + S_UH)/(S_RL + S_UL), # Relative high/low SES in scTB (rS_HL)
        Z = (C_RH + C_UH)/(C_RL + C_UL), # Relative high/low SES in cTB (rC_HL)
        S = S_RL + S_RH + S_UL + S_UH, # All subclinical TB
        C = C_RL + C_RH + C_UL + C_UH, # All clinical TB
        U = force_func_omega(time)*(C_RL + C_RH + C_UL + C_UH), # Per time
        V = force_func_iota(time)*(C_RL + C_RH + C_UL + C_UH))) # Per time
    })
  }
  
  yini <- c(N_RL = 100000-1000, N_RH = 0, N_UL = 0, N_UH = 0,
            I_RL = 0, I_RH = 0, I_UL = 0, I_UH = 0,
            O_RL = 0, O_RH = 0, O_UL = 0, O_UH = 0,
            M_RL = 0, M_RH = 0, M_UL = 0, M_UH = 0,
            S_RL = 0, S_RH = 0, S_UL = 0, S_UH = 0,
            C_RL = 1000, C_RH = 0, C_UL = 0, C_UH = 0,
            D_RL = 0, D_RH = 0, D_UL = 0, D_UH = 0)
  
  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

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
  S2007 = list(val = 144, sigma = unc*144), # Subclinical TB rate 1st survey
  S2017 = list(val = 128, sigma = unc*128), # Subclinical TB rate 2nd survey
  C2007 = list(val = 87, sigma = unc*87), # Clinical TB rate 1st survey
  C2017 = list(val = 73, sigma = unc*73), # Clinical TB rate 2nd survey
  U2000 = list(val = 40.5, sigma = unc*40.5), # Mortality 2000
  U2010 = list(val = 25.2, sigma = unc*25.2), # Mortality 2010
  U2020 = list(val = 10.3, sigma = unc*10.3), # Mortality 2020
  V2010 = list(val = 59.7, sigma = unc*59.7), # Diagnosed/treated 2010
  V2020 = list(val = 56.2, sigma = unc*56.2), # Diagnosed/treated 2020
  W2007 = list(val = 0.39, sigma = unc*0.39), # Relative urban/rural in scTB 2007 (rS_UR)
  W2017 = list(val = 0.84, sigma = unc*0.84), # Relative urban/rural in scTB 2017 (rS_UR)
  Y2007 = list(val = 0.34, sigma = unc*0.34), # Relative urban/rural in cTB 2007 (rC_UR)
  Y2017 = list(val = 0.67, sigma = unc*0.67), # Relative urban/rural in cTB 2017 (rC_UR)
  X2007 = list(val = 0.53, sigma = unc*0.53), # Relative high/low SES in scTB 2007 (rS_HL)
  X2017 = list(val = 0.83, sigma = unc*0.83), # Relative high/low SES in scTB 2017 (rS_HL)
  Z2007 = list(val = 0.50, sigma = unc*0.50), # Relative high/low SES in cTB 2007 (rC_HL)
  Z2017 = list(val = 0.91, sigma = unc*0.91)) # Relative high/low SES in cTB 2017 (rC_HL)

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,1) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),2,5)) # Create time variable
targetsdb <- as.data.frame(cbind(targetsdb[grep("val", rownames(targetsdb)),c("time","var","V1")],V2 = targetsdb[grep("sigma", rownames(targetsdb)),c("V1")])) # Widen dataframe
colnames(targetsdb)[c(3,4)] <- c("val","sigma") # Rename columns
targetsdb$lo <- targetsdb$val-2*targetsdb$sigma # Lower bound
targetsdb$hi <- targetsdb$val+2*targetsdb$sigma # Upper bound
export(targetsdb,here("data","targets.Rdata")) # Save data frame
rm(targetsdb)

# 3.3 HMER initial run
ems <- list() # Empty list for emulators per wave
wave <- list() # Empty list for data used to train and validate per wave
wave_samp <- list() # Empty list for wave samples
wave_train <- list() # Empty list for wave training
wave_val <- list() # Empty list for wave validation
wave_res <- list() # Empty list for wave results
invalid <- list() # Empty list for invalid parameter sets
invalid_diag <- list() # Empty list for invalid parameter sets after diagnostics
non_imp_pts <- list() # Empty list for non-implausible points generated per wave

n_params <- length(ranges) # Number of parameters
n_points <- 10*n_params # HMER points = selection of parameters (Quick = 10)
#n_points <- 20*n_params # HMER points = selection of parameters (Robust = 20)

ini_LHS_train <- lhs::maximinLHS(n_points, n_params)
ini_LHS_val <- lhs::maximinLHS(n_points, n_params)
ini_LHS <- rbind(ini_LHS_train, ini_LHS_val)
rm(ini_LHS_train, ini_LHS_val) # Clean objects

ini_pts <- setNames(data.frame(t(apply(ini_LHS, 1, function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + unlist(lapply(ranges, function(x) x[1]))))), names(ranges)) # Set random sets to create points
wave_res[[1]] <- data.frame(t(apply(ini_pts, 1, hmer_res, c(2000, 2007, 2010, 2017, 2020), c('S', 'C', 'U', 'V', 'W', 'Y', 'X', 'Z'))))[,names(targets)] # Run ODE
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
export(non_imp_pts[[1]],here("data","pts_w1.Rdata")) # Save data frame
beepr::beep(2)

# 3.4 HMER loop runs
w <- 10 # Update wave run
wave_res[[w]] <- data.frame(t(apply(non_imp_pts[[w-1]], 1, hmer_res, c(2000, 2007, 2010, 2017, 2020), c('S', 'C', 'U', 'V', 'W', 'Y', 'X', 'Z'))))[,names(targets)] # Run ODE
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
  bad.model <- nrow(comparison_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  if (bad.model > floor(nrow(wave_val[[w]])/10)) {
    bad.ems <- c(bad.ems, j)
  }
}

ems[[w]] <- ems[[w]][!seq_along(ems[[w]]) %in% bad.ems]

pdf(here("outputs", paste("val_w", w,"-diag.pdf", sep = "")))
invalid_diag[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()

non_imp_pts[[w]] <- generate_new_design(c(rev(ems[1:w])), n_points, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[w]], here("data",paste("pts_w", w,".Rdata", sep = ""))) # Save data frame
beepr::beep(2)

rm(list=ls()[!ls() %in% c("ode","parms", "w")]) # Clean objects (except "ode" and "parms")

# 4. Results ==========
w <- 10 # Update wave run
pts_fin <- import(here("data",paste("pts_w", w,".Rdata", sep = "")))

quants <- c(0.025,0.5,0.975) # Set quantiles

parameters <- apply(pts_fin, 2, quantile, probs = quants, na.rm = TRUE) # Set parameter quantiles
t_parameters <- data.table::transpose(as.data.frame(parameters)) # Transpose parameters
colnames(t_parameters) = rownames(parameters) # Set column names
rownames(t_parameters) = colnames(parameters) # Set row names
parameters <- t_parameters # Rename parameters
rm(t_parameters) # Clean objects

parameters$parameter <- c("beta","alpha","lambda_infmin","lambda_infsub","lambda_minsub","lambda_subcln",
                          "gamma_infcle","gamma_mincle","gamma_submin","gamma_clnsub","omega_ini","omega_fin",
                          "iota_ini","iota_fin", "rho_ini", "rho_fin", "sigma_ini", "sigma_fin")
parameters[,c(1,2,3)] <- round(parameters[,c(1,2,3)],2) # Round to 2 decimal places
table <- data.frame(parameter = parameters$parameter, low  = parameters$`2.5%`, med = parameters$`50%`, hig = parameters$`97.5%`) # Output table

results <- as.data.frame(apply(pts_fin, 1, ode))[-seq(1,521),] # Runs ODE for each set of points
results <- as.data.frame(t(apply(results, 1, quantile, probs = quants, na.rm = TRUE))) # Set parameter quantiles
comp <- colnames(ode(parms))[-1] # Compartment names (-time)
results$var <- c(t(replicate(521,comp))) # Compartment variable name
results$time <- rep(seq(1500,2020),36) # Time variable

export(results,here("data",paste("results_w", w,".Rdata", sep = ""))) # Save data frame
beepr::beep(2)
