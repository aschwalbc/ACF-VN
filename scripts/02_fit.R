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
  beta = 9,        # Contact (per person/year) parameter
  kappa = 0.75,     # Relative infectiousness
  infcle = 1.83,    # Infected -> Cleared
  infmin = 0.21,    # Infected -> Minimal
  minrec = 0.16,    # Minimal -> Recovered
  pi = 0.21,        # Protection from reinfection
  minsub = 0.25,    # Minimal -> Subclinical
  infsub = 0.07,    # Infected -> Subclinical
  submin = 1.58,    # Subclinical -> Minimal
  subcln = 0.77,    # Subclinical -> Clinical
  clnsub = 0.53,    # Clinical -> Subclinical
  mutb_ini = 0.3,   # TB mortality (Initial)
  mutb_fin = 0.23,  # TB mortality (Final)
  theta_ini = 0.44, # Diagnosis (Initial)
  theta_fin = 0.9,  # Diagnosis (Final)
  phi_ini = 0.69,   # Treatment failure (Initial) 
  phi_fin = 0.09,   # Treatment failure (Final)
  rho = 3.25)       # Risk of reinfection


# 1.2 Parameter ranges
ranges = list(
  beta = c(6,20),           # Contact (per person/year) parameter
  kappa = c(0.62,1),        # Relative infectiousness [Emery eLife 2023]
  infcle = c(0.93,3.30),    # Infected -> Cleared [Horton PNAS 2023]
  infmin = c(0.04,0.23),    # Infected -> Minimal [Horton PNAS 2023]
  minrec = c(0.14,0.23),    # Minimal -> Recovered [Horton PNAS 2023]
  pi = c(0.14,0.3),         # Protection from reinfection [Andrews et al 2012]
  minsub = c(0.21,0.28),    # Minimal -> Subclinical [Horton PNAS 2023]
  infsub = c(0.01,0.10),    # Infected -> Subclinical [Horton PNAS 2023]
  submin = c(1.24,2.03),    # Subclinical -> Minimal [Horton PNAS 2023]
  subcln = c(0.56,0.94),    # Subclinical -> Clinical [Horton PNAS 2023]
  clnsub = c(0.46,0.72),    # Clinical -> Subclinical [Horton PNAS 2023]
  mutb_ini = c(0.28,0.38),  # TB mortality (Initial) [Richards Lancet GH 2023]
  mutb_fin = c(0,0.28),     # TB mortality (Final)
  theta_ini = c(0,0.57),    # Diagnosis (Initial)
  theta_fin = c(0.57,0.77), # Diagnosis (Final) [WHO GTB]
  phi_ini = c(0.11,1),      # Treatment failure (Initial) 
  phi_fin = c(0.07,0.11),   # Treatment failure (Final) [WHO GTB]
  rho = c(2.14,4.27))       # Risk of reinfection [Verver AJRCCM 2005]

# 2. Model ==========
ode <- function(parms, end_time = 2020) {
  
  # Static parameters  
  N <- 1e5 # Population size
  mu <- 1/70 # Age expectancy adult
  delta <- 2 # Treatment duration (6 months)

  # Function parameters
  forcer_mutb <- matrix(c(1500, parms['mutb_ini'], 1999, parms['mutb_ini'], 2020, parms['mutb_fin']), ncol = 2, byrow = TRUE)
  force_func_mutb <- approxfun(x = forcer_mutb[,1], y = forcer_mutb[,2], method = "linear", rule = 2)
  
  forcer_theta <- matrix(c(1500, parms['theta_ini'], 1999, parms['theta_ini'], 2020, parms['theta_fin']), ncol = 2, byrow = TRUE)
  force_func_theta <- approxfun(x = forcer_theta[,1], y = forcer_theta[,2], method = "linear", rule = 2)
  
  forcer_phi <- matrix(c(1500, parms['phi_ini'], 1999, parms['phi_ini'], 2020, parms['phi_fin']), ncol = 2, byrow = TRUE)
  force_func_phi <- approxfun(x = forcer_phi[,1], y = forcer_phi[,2], method = "linear", rule = 2)
  
  des <- function(time, state, parms) {
    
    with(as.list(c(state, parms)), {
      
      dSUS = ((mu * N) + (force_func_mutb(time) * CLN)) - (((beta / N) * ((kappa * SUB) + CLN)) * SUS) - (mu * SUS)
      dINF = (((beta / N) * ((kappa * SUB) + CLN)) * (SUS + CLE + (pi * REC) + (rho * TRE))) - (infcle * INF) - (infmin * INF) - (infsub * INF) - (mu * INF)
      dCLE = (infcle * INF) - (((beta / N) * ((kappa * SUB) + CLN)) * CLE) - (mu * CLE)
      dREC = (minrec * MIN) - (((beta / N) * ((kappa * SUB) + CLN)) * (pi * REC)) - (mu * REC)
      dMIN = (infmin * INF) + (submin * SUB) - (minrec * MIN) - (minsub * MIN) - (mu * MIN)
      dSUB = (infsub * INF) + (minsub * MIN) + (clnsub * CLN) - (submin * SUB) - (subcln * SUB) - (mu * SUB)
      dCLN = (subcln * SUB) - (clnsub * CLN) - (force_func_theta(time) * CLN) + (force_func_phi(time) * TXT) - (force_func_mutb(time) * CLN) - (mu * CLN)
      dTXT = (force_func_theta(time) * CLN) - (force_func_phi(time) * TXT) - (delta * TXT) - (mu * TXT)
      dTRE = (delta * TXT) - (((beta / N) * ((kappa * SUB) + CLN)) * (rho * TRE)) - (mu * TRE)

      return(list(c( 
        dSUS, dINF, dCLE, dREC, dMIN, dSUB, dCLN, dTXT, dTRE),
        TBc   = (SUB + CLN), # TB prevalence (per 100k)
        Mor   = (force_func_mutb(time) * CLN), # TB mortality (per 100k)
        Dxs   = (force_func_theta(time) * CLN), # TB notifications (per 100k)
        Spr   = (SUB / (SUB + CLN)))) # Proportion subclinical TB (%)
    })
  }
  
  yini <- c(SUS = 1e5 - 1e3, INF = 0, CLE = 0, REC = 0,
            MIN = 0, SUB = 0, CLN = 1e3, TXT = 0, TRE = 0)
  
  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

# Quick diagnostic
tail(as.data.frame(ode(parms)))

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
  TBc2007 = c(202,310),     # TB prevalence [Nguyen Emerg Infect Dis 2021 - Revised]
  TBc2018 = c(177,290),     # TB prevalence [Nguyen Emerg Infect Dis 2021 - Revised]
  Mor2000 = c(37.3,87.7),   # TB mortality (Over population aged >= 15yo) [WHO GTB]
  Mor2010 = c(22.8,45.7),   # TB mortality (Over population aged >= 15yo) [WHO GTB]
  Dxs2010 = c(63.5,95.3),   # TB notifications (Over population aged >= 15yo) [WHO GTB]
  Dxs2020 = c(58.5,87.8),   # TB notifications (Over population aged >= 15yo) [WHO GTB]
  Spr2007 = c(0.56,0.84),   # Proportion scTB [Emery eLife 2023]
  Spr2018 = c(0.53,0.79))   # Proportion scTB [Emery eLife 2023]

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,3) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),4,7)) # Create time variable
colnames(targetsdb) <- c("lo","hi","var","time") # Rename columns

export(targetsdb,here("data","fit","targets.Rdata")) # Save data frame

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
target_subset <- function(wave = NULL, sens = 0, targets = targetsdb, parms = parms){
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
  wave_check <- wave_check[wave_check$check >= (length(target)-sens),]
  wave_check <- wave_check[,1:length(wave_check)-1]
  return(wave_check)
}

# 3.4 HMER initial run
wave_res <- list() # Empty list for wave results
wave <- list() # Empty list for data used to train and validate per wave
simulator <- list() # Empty list for simulator plots
checks <- list() # Empty list for target checks
wave_check <- list() # Empty list for wave checks
wave_train <- list() # Empty list for wave training
wave_val <- list() # Empty list for wave validation
ems <- list() # Empty list for emulators per wave
activeparms <- list() # Empty list for active parameters plots
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
  res <- t(apply(ini_pts[i,], 1, hmer_res, c(2000, 2007, 2010, 2018, 2020), c('TBc', 'Mor', 'Dxs', 'Spr')))
  tmp[[i]] <- data.frame(res)[, names(targets)]
  pb$tick()  # Advance the progress bar
}
wave_res[[1]] <- do.call(rbind, tmp)
export(wave_res[[1]], here("outputs", "wave_res", "w001_waveres.Rdata"))
wave[[1]] <- cbind(ini_pts, wave_res[[1]]) # Bind run points and results
export(wave[[1]], here("outputs", "waves", "w001_wave.Rdata"))
rm(pb, tmp, res, i) # Clean objects
simulator[[1]] <- simulator_plot(wave_res, targets, normalize = TRUE, byhit = TRUE)
pdf(here("outputs", "simulator", "w001_simulator.pdf"))
print(simulator[[1]])
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
activeparms[[1]] <- plot_actives(ems[[1]])
pdf(here("outputs", "activeparms", "w001_activeparms.pdf"))
print(activeparms[[1]])
dev.off()
cat("Emulators trained\n")

pdf(here("outputs", "diagnostics", "w001_diagnostics_pre.pdf"))
invalid_pts[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()
for (j in 1:length(ems[[1]])) {
  misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  while(misclass > 0) {
    ems[[1]][[j]] <- ems[[1]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[1]][[j]], targets, wave_val[[1]], plt = FALSE))
  }
}
pdf(here("outputs", "diagnostics", "w001_diagnostics_post.pdf"))
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
pdf(here("outputs", "diagnostics", "w001_diagnostics_post_badems.pdf"))
invalid_bad[[1]] <- validation_diagnostics(ems[[1]], validation = wave_val[[1]], targets = targets, plt = TRUE)
dev.off()
export(ems[[1]], here("outputs", "ems", "w001_ems.Rdata"))
cat("Diagnostics performed\n")

non_imp_pts[[1]] <- generate_new_design(ems[[1]], (10*length(ranges))*2, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[1]],here("outputs", "pts", "w001_pts.Rdata")) # Save data frame
cat("New parameter points generated\n")
beepr::beep(2)

# 3.5 HMER loop runs
# for (i in 180:180) {
w <- 180 # Update wave run
tic()
cat("Running wave:", w, "\n")

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(non_imp_pts[[w-1]]))
tmp <- list()
for (i in seq_len(nrow(non_imp_pts[[w-1]]))) {
  res <- t(apply(non_imp_pts[[w-1]][i,], 1, hmer_res, c(2000, 2007, 2010, 2018, 2020), c('TBc', 'Mor', 'Dxs', 'Spr')))
  tmp[[i]] <- data.frame(res)[, names(targets)]
  pb$tick()  # Advance progress bar
}

wave_res[[w]] <- do.call(rbind, tmp)
export(wave_res[[w]], here("outputs", "wave_res", sprintf("w%03d_waveres.Rdata", w)))
wave[[w]] <- cbind(non_imp_pts[[w-1]], wave_res[[w]])
export(wave[[w]], here("outputs", "waves", sprintf("w%03d_wave.Rdata", w)))
rm(pb, tmp, res, i) # Clean objects
simulator[[w]] <- simulator_plot(wave_res[w], targets, normalize = TRUE, byhit = TRUE)
pdf(here("outputs", "simulator", sprintf("w%03d_simulator.pdf", w)))
print(simulator[[w]])
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
activeparms[[w]] <- plot_actives(ems[[w]])
pdf(here("outputs", "activeparms", sprintf("w%03d_activeparms.pdf", w)))
print(activeparms[[w]])
dev.off()
cat("Emulators trained\n")

pdf(here("outputs", "diagnostics", sprintf("w%03d_diagnostics_pre.pdf", w)))
invalid_pts[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()
for (j in 1:length(ems[[w]])) {
  misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  while(misclass > 0) {
    ems[[w]][[j]] <- ems[[w]][[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems[[w]][[j]], targets, wave_val[[w]], plt = FALSE))
  }
}
pdf(here("outputs", "diagnostics", sprintf("w%03d_diagnostics_post.pdf", w)))
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
pdf(here("outputs", "diagnostics", sprintf("w%03d_diagnostics_post_badems.pdf", w)))
invalid_bad[[w]] <- validation_diagnostics(ems[[w]], validation = wave_val[[w]], targets = targets, plt = TRUE)
dev.off()
export(ems[[w]], here("outputs", "ems", sprintf("w%03d_ems.Rdata", w)))
cat("Diagnostics performed\n")

non_imp_pts[[w]] <- generate_new_design(c(ems[1:w]), (10*length(ranges))*2, targets, verbose = TRUE) # Generate new points
export(non_imp_pts[[w]], here("outputs", "pts", sprintf("w%03d_pts.Rdata", w))) # Save data frame
plaus_pts <- unique(as.data.frame(do.call("rbind", wave_check)))[1:length(parms)]
non_imp_pts[[w]] <- rbind(non_imp_pts[[w]], plaus_pts)
cat("New parameter points generated\n")
toc()
# }
beepr::beep(2)

# 3.6 Reload data
# Non-implausible points
data_files <- list.files(here("outputs","pts"), pattern = "w[0-9]+_pts.Rdata", full.names = TRUE)
for (file in data_files) {
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  data <- import(file)
  non_imp_pts[[num]] <- data
}

# Waves_res
data_files <- list.files(here("outputs","wave_res"), pattern = "w[0-9]+_waveres.Rdata", full.names = TRUE)
for (file in data_files) {
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  data <- import(file)
  wave_res[[num]] <- data
}

# Waves
data_files <- list.files(here("outputs","waves"), pattern = "w[0-9]+_wave.Rdata", full.names = TRUE)
for (file in data_files) {
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  data <- import(file)
  wave[[num]] <- data
}

# Emulators
data_files <- list.files(here("outputs","ems"), pattern = "w[0-9]+_ems.Rdata", full.names = TRUE)
for (file in data_files) {
  env <- new.env()
  load(file, env)
  objs <- mget(ls(env), env)
  num <- as.integer(gsub("[^0-9]", "", basename(file)))
  ems[[num]] <- objs
}

# Wave checks
for(w in 1:length(wave)) {
  checks[[w]] <- target_check(wave = wave_res[[w]])
  runcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["runs"]]))
  targetcheck <- do.call(cbind, lapply(1:w, function(i) checks[[i]][["targets"]]))
  wave_check[[w]] <- target_subset(wave[[w]])
}

rm(data_files, data, file, num, objs, env, w)

# 4. Results ==========
# Isolate best parameter sets
plaus_pts <- unique(as.data.frame(do.call("rbind", wave_check)))[1:length(parms)]
export(plaus_pts, here("outputs", "pts", "fitpts.Rdata")) # Save data frame

pts_fin <- plaus_pts

quants <- c(0.025,0.5,0.975) # Set quantiles
parameters <- apply(pts_fin, 2, quantile, probs = quants, na.rm = TRUE) # Set parameter quantiles
t_parameters <- data.table::transpose(as.data.frame(parameters)) # Transpose parameters
colnames(t_parameters) = rownames(parameters) # Set column names
rownames(t_parameters) = colnames(parameters) # Set row names
parameters <- t_parameters # Rename parameters
rm(t_parameters) # Clean objects

parameters$parameter <- c("beta", "kappa", "infcle", "infmin", "minrec", "pi", "minsub", "infsub", "submin", "subcln", "clnsub",
                          "mutb_ini", "mutb_fin", "theta_ini", "theta_fin", "phi_ini", "phi_fin", "rho")

parameters[,c(1,2,3)] <- round(parameters[,c(1,2,3)],2) # Round to 2 decimal places
table <- data.frame(parameter = parameters$parameter, low  = parameters$`2.5%`, med = parameters$`50%`, hig = parameters$`97.5%`) # Output table
export(table, here("data", "fit", "parms.Rdata")) # Save data frame

results <- as.data.frame(apply(pts_fin, 1, ode))[-seq(1,521),] # Runs ODE for each set of points
results <- as.data.frame(t(apply(results, 1, quantile, probs = quants, na.rm = TRUE))) # Set parameter quantiles
comp <- colnames(ode(parms))[-1] # Compartment names (-time)
results$var <- c(t(replicate(521,comp))) # Compartment variable name
results$time <- rep(seq(1500,2020), length(comp)) # Time variable

export(results, here("outputs", "results.Rdata")) # Save data frame
beepr::beep(2)
