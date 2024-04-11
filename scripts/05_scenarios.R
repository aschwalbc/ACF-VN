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
parms <- import(here("outputs","pts","fitpts.Rdata"))[1:1000,]
# parms <- parms %>% sample_frac(0.2)
WPP <- import(here("data","pop","WPP.Rdata"))
base <- import(here("data","fit","base.Rdata"))
DALYs <- import(here("docs","dalys","DALYs.xlsx"))
set.seed(251195) # Set seed for randomisation
options(scipen = 999) # Remove scientific notation

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  select(-contains('ini')) %>% 
  rename(mutb = mutb_fin, theta = theta_fin, phi = phi_fin) 

# 2.2 Intervention parameters
pop_target <- 1 # Proportion of population targeted for ACF
pop_reached <- 1 # Proportion of population participating in ACF

# 2.2.1 Xpert Ultra
xpert_fp_sic <- c(lo = 0.03, hi = 0.07)
xpert_fp_rec <- c(lo = 0.03, hi = 0.07)
xpert_sens_min <- c(lo = 0.03, hi = 0.07)
xpert_sens_sub <- c(lo = 0.68, hi = 0.86)
xpert_sens_cln <- c(lo = 0.86, hi = 0.95)
xpert_fp_tre <- c(lo = 0.03, hi = 0.30)

# 2.2.2 Chest X-ray
cxr_fp_sic <- c(lo = 0.07, hi = 0.13)
cxr_fp_rec <- c(lo = 0.07, hi = 0.13)
cxr_sens_min <- c(lo = 0.77, hi = 0.90)
cxr_sens_sub <- c(lo = 0.90, hi = 0.92)
cxr_sens_cln <- c(lo = 0.90, hi = 0.92)
cxr_fp_tre <- c(lo = 0.48, hi = 0.52)

# 2.2.3 Enhanced Xpert Ultra
enhxpert_fp_sic <- c(lo = 0.0016, hi = 0.0029)
enhxpert_fp_rec <- c(lo = 0.0016, hi = 0.0029)
enhxpert_sens_min <- c(lo = 0.0016, hi = 0.0029)
enhxpert_sens_sub <- c(lo = 0.68, hi = 0.86)
enhxpert_sens_cln <- c(lo = 0.86, hi = 0.95)
enhxpert_fp_tre <- c(lo = 0.01, hi = 0.08)

# 2.3 ACF scenarios
# 2.3.1 ACF Scenario A - Xpert only, Xpert+ get treatment
acfa <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = xpert_fp_sic['lo'], hi = xpert_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = xpert_fp_rec['lo'], hi = xpert_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = xpert_sens_min['lo'], hi = xpert_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = xpert_sens_sub['lo'], hi = xpert_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = xpert_sens_cln['lo'], hi = xpert_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = xpert_fp_tre['lo'], hi = xpert_fp_tre['hi']) 

# 2.3.2 ACF Scenario B - CXR only, CXR+ get Xpert, Xpert+ get treatment
acfb <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = xpert_fp_sic['lo']*cxr_fp_sic['lo'], hi = xpert_fp_sic['hi']*cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = xpert_fp_rec['lo']*cxr_fp_rec['lo'], hi = xpert_fp_rec['hi']*cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = xpert_sens_min['lo']*cxr_sens_min['lo'], hi = xpert_sens_min['hi']*cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = xpert_sens_sub['lo']*cxr_sens_sub['lo'], hi = xpert_sens_sub['hi']*cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = xpert_sens_cln['lo']*cxr_sens_cln['lo'], hi = xpert_sens_cln['hi']*cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = xpert_fp_tre['lo']*cxr_fp_tre['lo'], hi = xpert_fp_tre['hi']*cxr_fp_tre['hi'])

# 2.3.3 ACF Scenario C - CXR only, CXR+ get treatment
acfc <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = cxr_fp_sic['lo'], hi = cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = cxr_fp_rec['lo'], hi = cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = cxr_sens_min['lo'], hi = cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = cxr_sens_sub['lo'], hi = cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = cxr_sens_cln['lo'], hi = cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = cxr_fp_tre['lo'], hi = cxr_fp_tre['hi'])

# 2.3.4 ACF Scenario D - Enhanced Xpert only, Xpert+ get treatment *Same as ACF Scenario A*
acfd <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = enhxpert_fp_sic['lo'], hi = enhxpert_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = enhxpert_fp_rec['lo'], hi = enhxpert_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = enhxpert_sens_min['lo'], hi = enhxpert_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = enhxpert_sens_sub['lo'], hi = enhxpert_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = enhxpert_sens_cln['lo'], hi = enhxpert_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = enhxpert_fp_tre['lo'], hi = enhxpert_fp_tre['hi']) 

rm(list = ls(pattern = "^(xpert|cxr|enh)"))

# 2.7 Other parameters
# 2.7.1 Population parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
rm(WPP)

# 2.7.2 Cost data
screen_bau <- c(dstb = 104, drtb = 500) # Passive case-finding
screen_acf <- c(acfa = 8, acfb = 1.7, acfc = 1, acfd = 8) # ACF per algorithm
tb_rx <- c(dstb = 81, drtb = 973) # TB treatment

# 2.7.3 DALYs
daly_lo <- approxfun(DALYs$year, DALYs$lo, method = 'linear', rule = 2)
daly_hi <- approxfun(DALYs$year, DALYs$hi, method = 'linear', rule = 2)

# 2.7.4 MDR proportion
mdr_prop <- c(val = 0.049, lo = 0.028, hi = 0.069)  # DR-TB incidence over total incidence (2015-2021)

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
  
  # Proportion providing sputum in non-disease states
  if(is.null(acf_times)) {
    prop_sputum <- 1 
  } else {
    interv_name <- deparse(substitute(interv))
    
    if(interv_name == "acfa" | interv_name == "acfb" | interv_name == "acfd") {
      prop_sputum <- 0.6
    } else if(interv_name == "acfc") {
      prop_sputum <- 1
    }
  }

  # Active-case finding switch
  if(is.null(acf_times)) {
    acf <- function(times) 0
  } else { 
    values <- ifelse(floor(times) %in% acf_times, 1, 0)
    acf <- approxfun(times, values, rule = 2)
  }
  
  # Active-case finding costs
  if(is.null(interv)) {
    cm_screen_acf <- 0
  } else {
    interv_name <- deparse(substitute(interv))
    
    if(interv_name == "acfa") {
      cm_screen_acf <- screen_acf[["acfa"]]
    } else if(interv_name == "acfb") {
      cm_screen_acf <- screen_acf[["acfb"]]
    } else if(interv_name == "acfc") {
      cm_screen_acf <- screen_acf[["acfc"]]
    } else if(interv_name == "acfd") {
      cm_screen_acf <- screen_acf[["acfd"]]
    }
  }
  
  # MDR TB
  mdr <- runif(1, min = mdr_prop[["lo"]], max = mdr_prop[["hi"]])
  
  # DALY uncertainty
  daly <- function(times) {
    lo <- daly_lo(floor(times))
    hi <- daly_hi(floor(times))
    runif(1, min = lo, max = hi)
  }
  
  # Treatment duration 
  delta <- 2  # 6 months

  des <- function(times, state, parms) {
    with(as.list(c(times, state, parms)), {
      
      PopT  = (SUS + SUSx + INF + INFx + CLE + CLEx + REC + RECx + MIN + MINx + SUB + SUBx + CLN + CLNx + TXT + TRE + TREx)
      
      dSUS  = (nu(times) * (PopT)) - (((beta / PopT) * ((kappa * SUB) + CLN)) * SUS) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * SUS) + (delta * SUSx) - (mu(times) * SUS)
      dSUSx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * SUS) - (delta * SUSx) - (mu(times) * SUSx)
      dINF  = (((beta / PopT) * ((kappa * SUB) + CLN)) * (SUS + CLE + (pi * REC) + (rho * TRE))) - (infcle * INF) - (infmin * INF) - (infsub * INF) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * INF) - (mu(times) * INF)
      dINFx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * INF) - (delta * INFx) - (mu(times) * INFx)
      dCLE  = (infcle * INF) - (((beta / PopT) * ((kappa * SUB) + CLN)) * CLE) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * CLE) + (delta * (INFx + CLEx)) - (mu(times) * CLE)
      dCLEx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * CLE) - (delta * CLEx) - (mu(times) * CLEx)
      dREC  = (minrec * MIN) - (((beta / PopT) * ((kappa * SUB) + CLN)) * (pi * REC)) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_rec * REC) + (delta * RECx) - (mu(times) * REC)
      dRECx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_rec * REC) - (delta * RECx) - (mu(times) * RECx)
      dMIN  = (infmin * INF) + (submin * SUB) - (minrec * MIN) - (minsub * MIN) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_min * MIN) - (mu(times) * MIN)
      dMINx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_min * MIN) - (delta * MINx) - (mu(times) * MIN)
      dSUB  = (infsub * INF) + (minsub * MIN) + (clnsub * CLN) - (submin * SUB) - (subcln * SUB) - (acf(floor(times)) * pop_target * pop_reached * alpha_sub * SUB) - (mu(times) * SUB)
      dSUBx = (acf(floor(times)) * pop_target * pop_reached * alpha_sub * SUB) - (delta * SUBx) - (mu(times) * SUBx)
      dCLN  = (subcln * SUB) - (clnsub * CLN) - (theta * CLN) + (phi * (CLNx + TXT)) - (acf(floor(times)) * pop_target * pop_reached * alpha_cln * CLN) - (mutb * CLN) - (mu(times) * CLN)
      dCLNx = (acf(floor(times)) * pop_target * pop_reached * alpha_cln * CLN) - (phi * CLNx) - (delta * CLNx) - (mu(times) * CLNx)
      dTXT  = (theta * CLN) - (phi * TXT) - (delta * TXT) - (mu(times) * TXT)
      dTRE  = (delta * (MINx + SUBx + CLNx + TXT + TREx)) - (((beta / PopT) * ((kappa * SUB) + CLN)) * (rho * TRE)) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_tre * TRE) - (mu(times) * TRE)
      dTREx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_tre * TRE) - (delta * TREx) - (mu(times) * TRE)

      return(list(c(
        dSUS, dSUSx, dINF, dINFx, dCLE, dCLEx, dREC, dRECx, dMIN, dMINx, dSUB, dSUBx, dCLN, dCLNx, dTXT, dTRE, dTREx),
        rTBc    = ((SUB + CLN) / PopT * 1e5), # Infectious TB (per 100k)
        tTBc    = (SUB + CLN), # Total infectious TB
        rInc    = ((subcln * SUB)/ PopT * 1e5), # TB incidence (per 100k)
        tInc    = (subcln * SUB), # Total TB incidence
        rMor    = ((mutb * CLN) / PopT * 1e5), # Clinical TB mortality per time (per 100k)
        tMor    = (mutb * CLN), # Clinical TB mortality per time
        rDxs    = ((theta * CLN) / PopT * 1e5), # Notifications cTB per time in adults (per 100k)
        tDxs    = (theta * CLN), # Total notifications cTB per time in adults
        tInf    = (((beta / PopT) * ((kappa * SUB) + CLN)) * (SUS + CLE + (pi * REC) + (rho * TRE))), # Total infections
        tFP     = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * ((alpha_sic * (SUS + INF + CLE)) + (alpha_rec * REC) + (alpha_tre * (TRE)))), # FP diagnosed
        tTP     = ((acf(floor(times))) * ((pop_target * pop_reached * prop_sputum * (alpha_min * MIN)) + (pop_target * pop_reached * ((alpha_sub * SUB) + (alpha_cln * CLN))))), # TP diagnosed
        tTPM    = (acf(floor(times)) * (pop_target * pop_reached * prop_sputum * (alpha_min * MIN))), # TP minimal diagnosed
        tTPI    = (acf(floor(times)) * (pop_target * pop_reached * ((alpha_sub * SUB) + (alpha_cln * CLN)))), # TP subclinical + clinical diagnosed
        tScrn   = ((acf(floor(times))) * (((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + MIN + TRE)) + ((pop_target * pop_reached) * (SUB + CLN)))), # Total screened
        cBAUs   = (((1-mdr) * (theta * CLN)) * screen_bau[["dstb"]]), # Costs BAU DS-TB
        cBAUr   = (((mdr) * (theta * CLN)) * screen_bau[["drtb"]]), # Costs BAU DR-TB
        cACF    = (((acf(floor(times))) * (((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + MIN + TRE)) + ((pop_target * pop_reached) * (SUB + CLN)))) * cm_screen_acf), # Cost ACF
        cACFND  = (((acf(floor(times))) * ((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + TRE))) * cm_screen_acf), # Cost ACF non-disease state
        cACFM   = (((acf(floor(times))) * ((pop_target * pop_reached * prop_sputum * MIN))) * cm_screen_acf), # Cost ACF minimal state
        cACFI   = (((acf(floor(times))) * ((pop_target * pop_reached) * (SUB + CLN))) * cm_screen_acf), # Cost ACF infectious state
        cRxBAUs = (((1-mdr) * (theta * CLN)) * tb_rx[["dstb"]]), # Costs treatment BAU DS-TB
        cRxBAUr = (((mdr) * (theta * CLN)) * tb_rx[["drtb"]]), # Costs treatment BAU DR-TB
        cRxFPs  = (((1-mdr) * ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * ((alpha_sic * (SUS + INF + CLE)) + (alpha_rec * REC) + (alpha_tre * (TRE))))) * tb_rx[["dstb"]]), # Cost treatment false positives DS-TB
        cRxFPr  = (((mdr) * ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * ((alpha_sic * (SUS + INF + CLE)) + (alpha_rec * REC) + (alpha_tre * (TRE))))) * tb_rx[["drtb"]]), # Cost treatment false positives DR-TB
        cRxTPs  = (((1-mdr) * ((acf(floor(times)) * ((pop_target * pop_reached * prop_sputum * (alpha_min * MIN)) + ((pop_target * pop_reached) * ((alpha_sub * SUB) + (alpha_cln * CLN))))))) * tb_rx[["dstb"]]), # Cost treatment true positives DS-TB
        cRxTPr  = (((mdr) * ((acf(floor(times)) * ((pop_target * pop_reached * prop_sputum * (alpha_min * MIN)) + ((pop_target * pop_reached) * ((alpha_sub * SUB) + (alpha_cln * CLN))))))) * tb_rx[["drtb"]]), # Cost treatment true positives DR-TB
        cRxTPMs  = (((1-mdr) * (acf(floor(times)) * (pop_target * pop_reached * prop_sputum * (alpha_min * MIN)))) * tb_rx[["dstb"]]), # Cost treatment true positives minimal DS-TB
        cRxTPMr  = (((mdr) * (acf(floor(times)) * (pop_target * pop_reached * prop_sputum * (alpha_min * MIN)))) * tb_rx[["drtb"]]), # Cost treatment true positives minimal DR-TB
        cRxTPIs  = (((1-mdr) * (acf(floor(times)) * (pop_target * pop_reached * ((alpha_sub * SUB) + (alpha_cln * CLN))))) * tb_rx[["dstb"]]), # Cost treatment true infectious positives DS-TB
        cRxTPIr  = (((mdr) * (acf(floor(times)) * (pop_target * pop_reached * ((alpha_sub * SUB) + (alpha_cln * CLN))))) * tb_rx[["drtb"]]), # Cost treatment true infectious positives DR-TB
        DALYs   = (daly(times) * (subcln * SUB)), # DALY estimates
        ARI     = ((beta / PopT) * ((kappa * SUB) + CLN)))) # ARI
    })
  }
  
  yini <- c(SUS = base[, "SUS"], SUSx = 0, INF = base[, "INF"], INFx = 0, CLE = base[, "CLE"], CLEx = 0,
            REC = base[, "REC"], RECx = 0, MIN = base[, "MIN"], MINx = 0, SUB = base[, "SUB"], SUBx = 0,
            CLN = base[, "CLN"], CLNx = 0, TXT = base[, "TXT"], TRE = base[, "TRE"], TREx = 0)

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
  outbase[[i]] <- outbase[[i]] %>% mutate(type = 'base', run = i, round = "00")
  pb$tick()
}
r00_outbase_df <- do.call(rbind, outbase)
export(r00_outbase_df, here("outputs", "results", "r00_base.Rdata"))

# 4.2 Scenarios
acf_year <- list(
  "01" = 2025,
  "02" = seq(2025, 2026, 1),
  "03" = seq(2025, 2027, 1),
  "06" = seq(2025, 2030, 1),
  "07" = seq(2025, 2031, 1),
  "11" = seq(2025, 2035, 1),
  "12" = seq(2025, 2036, 1))

outacfa <- list()
outacfb <- list()
outacfc <- list()
outacfd <- list()

for (j in names(acf_year)) {
  rounds <- acf_year[[j]]
  print(names(acf_year[j]))
  
  if (j %in% c("03", "06", "11")) {
    outacfa[[j]] <- list()
  }
  if (j %in% c("03", "07", "12")) {
    outacfb[[j]] <- list()
  }
  if (j %in% c("01", "02", "03")) {
    outacfc[[j]] <- list()
  }
  if (j %in% c("03", "06", "11")) {
    outacfd[[j]] <- list()
  }
  
  pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))
  
  for (i in 1:nrow(parms)) {
    curr_parms <- as.data.frame(parms[i,])
    curr_base <- as.data.frame(base[i,-1])
    
    if (j %in% c("03", "06", "11")) {
      outacfa[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = acf_year[[j]]))
      outacfa[[j]][[i]] <- outacfa[[j]][[i]] %>% mutate(type = 'acfa', run = i, round = j)
    }
    
    if (j %in% c("03", "07", "12")) {
      outacfb[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfb, acf_times = acf_year[[j]]))
      outacfb[[j]][[i]] <- outacfb[[j]][[i]] %>% mutate(type = 'acfb', run = i, round = j)
    }
    
    if (j %in% c("01", "02", "03")) {
      outacfc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfc, acf_times = acf_year[[j]]))
      outacfc[[j]][[i]] <- outacfc[[j]][[i]] %>% mutate(type = 'acfc', run = i, round = j)
    }
    
    if (j %in% c("03", "06", "11")) {
      outacfd[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfd, acf_times = acf_year[[j]]))
      outacfd[[j]][[i]] <- outacfd[[j]][[i]] %>% mutate(type = 'acfd', run = i, round = j)
    }
    
    pb$tick()
  }
}

for (k in 1:3) { 
  acfa_name <- paste0("r", names(outacfa[k]), "_outacfa_df")
  assign(acfa_name, do.call(rbind, outacfa[[k]]))
  export(get(acfa_name), here("outputs", "results", paste0(acfa_name, ".Rdata")))
  
  acfb_name <- paste0("r", names(outacfb[k]), "_outacfb_df")
  assign(acfb_name, do.call(rbind, outacfb[[k]]))
  export(get(acfb_name), here("outputs", "results", paste0(acfb_name, ".Rdata")))
  
  acfc_name <- paste0("r", names(outacfc[k]), "_outacfc_df")
  assign(acfc_name, do.call(rbind, outacfc[[k]]))
  export(get(acfc_name), here("outputs", "results", paste0(acfc_name, ".Rdata")))
  
  acfd_name <- paste0("r", names(outacfd[k]), "_outacfd_df")
  assign(acfd_name, do.call(rbind, outacfd[[k]]))
  export(get(acfd_name), here("outputs", "results", paste0(acfd_name, ".Rdata")))
}

# 4.3 ACT3 comparison
act3 <- list()

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))

for (i in 1:nrow(parms)) {
  curr_parms <- as.data.frame(parms[i,])
  curr_base <- as.data.frame(base[i,-1])
  act3[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = seq(2025, 2027, 1))) 
  act3[[i]] <- act3[[i]] %>% mutate(type = 'acfa', run = i, round = '03')
  
  pb$tick()
}

act3 <- do.call(rbind, act3)
export(act3, here("outputs", "act3", "act3.Rdata"))

rm(list = ls())
