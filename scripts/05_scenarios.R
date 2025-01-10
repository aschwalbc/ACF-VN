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
# parms <- parms %>% sample_frac(0.1)
WPP <- import(here("data","pop","WPP.Rdata"))
base <- import(here("data","fit","base.Rdata"))
DALYs <- import(here("docs","dalys","DALYs.xlsx"))
set.seed(251195) # Set seed for randomisation

# 2. Strategies ==========
# 2.1 Baseline parameters
parms <- parms %>% 
  select(-contains('ini')) %>% 
  rename(mutb = mutb_fin, theta = theta_fin, phi = phi_fin) 

# 2.2 Intervention parameters
pop_target <- 1 # Proportion of population targeted for ACF
pop_reached <- 1 # Proportion of population participating in ACF

# 2.2.1 Xpert Ultra
ultra_fp_sic <- c(lo = 0.005, hi = 0.008)
ultra_fp_rec <- c(lo = 0.020, hi = 0.060)
ultra_sens_min <- c(lo = 0.026, hi = 0.070)
ultra_sens_sub <- c(lo = 0.676, hi = 0.856)
ultra_sens_cln <- c(lo = 0.862, hi = 0.947)
ultra_fp_tre <- c(lo = 0.020, hi = 0.060)

# 2.2.2 Chest X-ray
cxr_fp_sic <- c(lo = 0.069, hi = 0.134)
cxr_fp_rec <- c(lo = 0.481, hi = 0.524)
cxr_sens_min <- c(lo = 0.626, hi = 0.712)
cxr_sens_sub <- c(lo = 0.626, hi = 0.712)
cxr_sens_cln <- c(lo = 0.770, hi = 0.900)
cxr_fp_tre <- c(lo = 0.481, hi = 0.524)

# 2.2.3 Xpert MTB/RIF (SA)
xpert_fp_sic <- c(lo = 0.0016, hi = 0.0029)
xpert_fp_rec <- c(lo = 0.005, hi = 0.083)
xpert_sens_min <- c(lo = 0.007, hi = 0.030)
xpert_sens_sub <- c(lo = 0.484, hi = 0.717)
xpert_sens_cln <- c(lo = 0.786, hi = 0.899)
xpert_fp_tre <- c(lo = 0.005, hi = 0.083)

# 2.2.4 Chest X-ray (SA)
cxrsa_fp_sic <- c(lo = 0.069, hi = 0.134)
cxrsa_fp_rec <- c(lo = 0.481, hi = 0.524)
cxrsa_sens_min <- c(lo = 0.481, hi = 0.524)
cxrsa_sens_sub <- c(lo = 0.481, hi = 0.524)
cxrsa_sens_cln <- c(lo = 0.770, hi = 0.900)
cxrsa_fp_tre <- c(lo = 0.481, hi = 0.524)

# 2.3 ACF scenarios
# 2.3.1 ACF Scenario A - Ultra only, Ultra+ get treatment
acfa <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = ultra_fp_sic['lo'], hi = ultra_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = ultra_fp_rec['lo'], hi = ultra_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = ultra_sens_min['lo'], hi = ultra_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = ultra_sens_sub['lo'], hi = ultra_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = ultra_sens_cln['lo'], hi = ultra_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = ultra_fp_tre['lo'], hi = ultra_fp_tre['hi']) 

# 2.3.2 ACF Scenario B - CXR only, CXR+ get Ultra, Ultra+ get treatment
acfb <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = ultra_fp_sic['lo']*cxr_fp_sic['lo'], hi = ultra_fp_sic['hi']*cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = ultra_fp_rec['lo']*cxr_fp_rec['lo'], hi = ultra_fp_rec['hi']*cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = ultra_sens_min['lo']*cxr_sens_min['lo'], hi = ultra_sens_min['hi']*cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = ultra_sens_sub['lo']*cxr_sens_sub['lo'], hi = ultra_sens_sub['hi']*cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = ultra_sens_cln['lo']*cxr_sens_cln['lo'], hi = ultra_sens_cln['hi']*cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = ultra_fp_tre['lo']*cxr_fp_tre['lo'], hi = ultra_fp_tre['hi']*cxr_fp_tre['hi'])

# 2.3.3 ACF Scenario C - CXR only, CXR+ get treatment
acfc <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = cxr_fp_sic['lo'], hi = cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = cxr_fp_rec['lo'], hi = cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = cxr_sens_min['lo'], hi = cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = cxr_sens_sub['lo'], hi = cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = cxr_sens_cln['lo'], hi = cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = cxr_fp_tre['lo'], hi = cxr_fp_tre['hi'])

# 2.3.4 ACF SA #01 Scenario - MTB/RIF only, MTB/RIF+ get treatment
acfd <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = xpert_fp_sic['lo'], hi = xpert_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = xpert_fp_rec['lo'], hi = xpert_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = xpert_sens_min['lo'], hi = xpert_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = xpert_sens_sub['lo'], hi = xpert_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = xpert_sens_cln['lo'], hi = xpert_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = xpert_fp_tre['lo'], hi = xpert_fp_tre['hi'])

# 2.3.5 ACF SA #02 Scenario - Reduced NAAT unit cost
acfax <- acfa # Replicating for costing SA
acfbx <- acfb # Replicating for costing SA

# 2.3.6 ACF SA #03 Scenario - Clinical appraisal
acfac <- acfa # Replicating for clinical appraisal SA
acfbc <- acfb # Replicating for clinical appraisal SA
acfcc <- acfc # Replicating for clinical appraisal SA

clnap_parms <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "clnap_sic", lo = 0.047, hi = 0.074) %>% 
  add_row(parameter = "clnap_rec", lo = 0.089, hi = 0.162) %>% 
  add_row(parameter = "clnap_min", lo = 0.129, hi = 0.249) %>% 
  add_row(parameter = "clnap_sub", lo = 1.000, hi = 1.000) %>%
  add_row(parameter = "clnap_cln", lo = 1.000, hi = 1.000) %>% 
  add_row(parameter = "clnap_tre", lo = 0.545, hi = 0.581)

# 2.3.7 ACF SA #04 Scenario - CXR test positivity
acfbr <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = ultra_fp_sic['lo']*cxrsa_fp_sic['lo'], hi = ultra_fp_sic['hi']*cxrsa_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = ultra_fp_rec['lo']*cxrsa_fp_rec['lo'], hi = ultra_fp_rec['hi']*cxrsa_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = ultra_sens_min['lo']*cxrsa_sens_min['lo'], hi = ultra_sens_min['hi']*cxrsa_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = ultra_sens_sub['lo']*cxrsa_sens_sub['lo'], hi = ultra_sens_sub['hi']*cxrsa_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = ultra_sens_cln['lo']*cxrsa_sens_cln['lo'], hi = ultra_sens_cln['hi']*cxrsa_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = ultra_fp_tre['lo']*cxrsa_fp_tre['lo'], hi = ultra_fp_tre['hi']*cxrsa_fp_tre['hi'])

acfcr <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = cxrsa_fp_sic['lo'], hi = cxrsa_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = cxrsa_fp_rec['lo'], hi = cxrsa_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = cxrsa_sens_min['lo'], hi = cxrsa_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = cxrsa_sens_sub['lo'], hi = cxrsa_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = cxrsa_sens_cln['lo'], hi = cxrsa_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = cxrsa_fp_tre['lo'], hi = cxrsa_fp_tre['hi'])

rm(list = ls(pattern = "^(ultra|cxr|xpert)"))

# 2.7 Other parameters
# 2.7.1 Population parameters
mu <- approxfun(WPP$year, WPP$mortrate, method = 'linear', rule = 2)
nu <- approxfun(WPP$year, WPP$birthrate, method = 'linear', rule = 2)
rm(WPP)

# 2.7.2 Cost data
gamma_dist <- function(mean, sdpct = 0.2) {
  sd <- mean * sdpct
  varc <- sd^2
  shape <- mean^2 / varc
  scale <- varc / mean 
  return(c(shape = shape, scale = scale))
}

# 2.7.2.1 Passive case-finding
scr_bau_dstb <- gamma_dist(264)
scr_bau_drtb <- gamma_dist(1595)

# 2.7.2.2 Active case-finding
# 2.7.2.2.1 Main analysis
scr_acfa <- gamma_dist(8)
scr_acfb <- gamma_dist(1.7)
scr_acfc <- gamma_dist(1.2)

# 2.7.2.2.2 Sensitivity analyses
scr_acfd <- gamma_dist(8)

scr_acfax <- gamma_dist(3)
scr_acfbx <- gamma_dist(1.3)

scr_acfac <- gamma_dist(8)
scr_acfbc <- gamma_dist(1.7)
scr_acfcc <- gamma_dist(1.2)

scr_acfbr <- gamma_dist(1.7)
scr_acfcr <- gamma_dist(1.2)

# 2.7.2.3 TB treatment
tb_rx_dstb <- gamma_dist(81)
tb_rx_drtb <- gamma_dist(973)

# 2.7.3 DALYs
daly_lo <- approxfun(DALYs$year, DALYs$lo, method = 'linear', rule = 2)
daly_hi <- approxfun(DALYs$year, DALYs$hi, method = 'linear', rule = 2)

# 2.7.4 MDR proportion
mdr_prop <- c(val = 0.049, lo = 0.028, hi = 0.069)  # DR-TB incidence over total incidence (2015-2021)

# 3. Models ==========
ode <- function(parms, base, interv = NULL, acf_times = NULL, clnap = NULL, end_time = 2050) {
  
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
  
  # Clinical appraisal parameters
  if(is.null(clnap)){
    clnap_sic <- 1
    clnap_rec <- 1
    clnap_min <- 1
    clnap_sub <- 1
    clnap_cln <- 1
    clnap_tre <- 1
  } else { 
    clnap_sic <- runif(1, min = clnap[clnap$parameter == 'clnap_sic', 'lo'], max = clnap[clnap$parameter == 'clnap_sic', 'hi'])
    clnap_rec <- runif(1, min = clnap[clnap$parameter == 'clnap_rec', 'lo'], max = clnap[clnap$parameter == 'clnap_rec', 'hi'])
    clnap_min <- runif(1, min = clnap[clnap$parameter == 'clnap_min', 'lo'], max = clnap[clnap$parameter == 'clnap_min', 'hi'])
    clnap_sub <- runif(1, min = clnap[clnap$parameter == 'clnap_sub', 'lo'], max = clnap[clnap$parameter == 'clnap_sub', 'hi'])
    clnap_cln <- runif(1, min = clnap[clnap$parameter == 'clnap_cln', 'lo'], max = clnap[clnap$parameter == 'clnap_cln', 'hi'])
    clnap_tre <- runif(1, min = clnap[clnap$parameter == 'clnap_tre', 'lo'], max = clnap[clnap$parameter == 'clnap_tre', 'hi'])
  }
  
  # Proportion providing sputum in non-disease states
  if(is.null(acf_times)) {
    prop_sputum <- 1 
  } else {
    interv_name <- deparse(substitute(interv))
    
    if(interv_name %in% c('acfa', 'acfax', 'acfd', 'acfac')) {
      prop_sputum <- 0.6
    } else if(interv_name %in% c('acfb', 'acfc', 'acfbx', 'acfbc', 'acfcc', 'acfbr', 'acfcr')) {
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
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfa[['shape']], scale = scr_acfa[['scale']])
    } else if(interv_name == "acfb") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfb[['shape']], scale = scr_acfb[['scale']])
    } else if(interv_name == "acfc") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfc[['shape']], scale = scr_acfc[['scale']])
    } else if(interv_name == "acfd") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfd[['shape']], scale = scr_acfd[['scale']])
    } else if(interv_name == "acfax") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfax[['shape']], scale = scr_acfax[['scale']])
    } else if(interv_name == "acfbx") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfbx[['shape']], scale = scr_acfbx[['scale']])
    } else if(interv_name == "acfac") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfac[['shape']], scale = scr_acfac[['scale']])
    } else if(interv_name == "acfbc") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfbc[['shape']], scale = scr_acfbc[['scale']])
    } else if(interv_name == "acfcc") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfcc[['shape']], scale = scr_acfcc[['scale']])
    } else if(interv_name == "acfbr") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfbr[['shape']], scale = scr_acfbr[['scale']])
    } else if(interv_name == "acfcr") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfcr[['shape']], scale = scr_acfcr[['scale']])
    }
  }
  
  # BAU costs
  bau_dstb <- rgamma(n = 1, shape = scr_bau_dstb[['shape']], scale = scr_bau_dstb[['scale']])
  bau_drtb <- rgamma(n = 1, shape = scr_bau_drtb[['shape']], scale = scr_bau_drtb[['scale']])
  
  # Rx costs
  rx_dstb <- rgamma(n = 1, shape = tb_rx_dstb[['shape']], scale = tb_rx_dstb[['scale']])
  rx_drtb <- rgamma(n = 1, shape = tb_rx_drtb[['shape']], scale = tb_rx_drtb[['scale']])
  
  # Clinical appraisal costs
  if(is.null(clnap)){
    clnap_cost <- 0
  } else { 
    clnap_cost <- 2
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
      
      dSUS  = (nu(times) * (PopT)) - (((beta / PopT) * ((kappa * SUB) + CLN)) * SUS) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * SUS) + (delta * SUSx) - (mu(times) * SUS)
      dSUSx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * SUS) - (delta * SUSx) - (mu(times) * SUSx)

      dINF  = (((beta / PopT) * ((kappa * SUB) + CLN)) * (SUS + CLE + (pi * REC) + (rho * TRE))) - (infcle * INF) - (infmin * INF) - (infsub * INF) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * INF) - (mu(times) * INF)
      dINFx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * INF) - (delta * INFx) - (mu(times) * INFx)

      dCLE  = (infcle * INF) - (((beta / PopT) * ((kappa * SUB) + CLN)) * CLE) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * CLE) + (delta * (INFx + CLEx)) - (mu(times) * CLE)
      dCLEx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_sic * clnap_sic * CLE) - (delta * CLEx) - (mu(times) * CLEx)

      dREC  = (minrec * MIN) - (((beta / PopT) * ((kappa * SUB) + CLN)) * (pi * REC)) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_rec * clnap_rec * REC) + (delta * RECx) - (mu(times) * REC)
      dRECx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_rec * clnap_rec * REC) - (delta * RECx) - (mu(times) * RECx)

      dMIN  = (infmin * INF) + (submin * SUB) - (minrec * MIN) - (minsub * MIN) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_min * clnap_min * MIN) - (mu(times) * MIN)
      dMINx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_min * clnap_min * MIN) - (delta * MINx) - (mu(times) * MIN)

      dSUB  = (infsub * INF) + (minsub * MIN) + (clnsub * CLN) - (submin * SUB) - (subcln * SUB) - (acf(floor(times)) * pop_target * pop_reached * alpha_sub * clnap_sub * SUB) - (mu(times) * SUB)
      dSUBx = (acf(floor(times)) * pop_target * pop_reached * alpha_sub * clnap_sub * SUB) - (delta * SUBx) - (mu(times) * SUBx)

      dCLN  = (subcln * SUB) - (clnsub * CLN) - (theta * CLN) + (phi * (CLNx + TXT)) - (acf(floor(times)) * pop_target * pop_reached * alpha_cln * clnap_cln * CLN) - (mutb * CLN) - (mu(times) * CLN)
      dCLNx = (acf(floor(times)) * pop_target * pop_reached * alpha_cln * clnap_cln * CLN) - (phi * CLNx) - (delta * CLNx) - (mu(times) * CLNx)

      dTXT  = (theta * CLN) - (phi * TXT) - (delta * TXT) - (mu(times) * TXT)
      
      dTRE  = (delta * (MINx + SUBx + CLNx + TXT + TREx)) - (((beta / PopT) * ((kappa * SUB) + CLN)) * (rho * TRE)) - (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_tre * clnap_tre * TRE) - (mu(times) * TRE)
      dTREx = (acf(floor(times)) * pop_target * pop_reached * prop_sputum * alpha_tre * clnap_tre * TRE) - (delta * TREx) - (mu(times) * TRE)

      return(list(c(
        dSUS, dSUSx, dINF, dINFx, dCLE, dCLEx, dREC, dRECx, dMIN, dMINx, dSUB, dSUBx, dCLN, dCLNx, dTXT, dTRE, dTREx),
        rTBc      = ((SUB + CLN) / PopT * 1e5), # Infectious TB (per 100k)
        tTBc      = (SUB + CLN), # Total infectious TB
        rInc      = ((subcln * SUB)/ PopT * 1e5), # TB incidence (per 100k)
        tInc      = (subcln * SUB), # Total TB incidence
        rMor      = ((mutb * CLN) / PopT * 1e5), # Clinical TB mortality per time (per 100k)
        tMor      = (mutb * CLN), # Clinical TB mortality per time
        tScrn     = ((acf(floor(times)) * (((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + MIN + TRE)) + ((pop_target * pop_reached) * (SUB + CLN))))), # Total number screened
        tSUS      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * clnap_sic * SUS)), # Total susceptible diagnosed (SUS)
        tINF      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * clnap_sic * INF)), # Total infected diagnosed (INF)
        tCLE      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * clnap_sic * CLE)), # Total cleared diagnosed (CLE)
        tREC      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_rec * clnap_rec * REC)), # Total recovered diagnosed (REC)
        tMIN      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_min * clnap_min * MIN)), # Total minimal diagnosed (MIN)
        tSUB      = ((acf(floor(times)) * pop_target * pop_reached) * (alpha_sub * clnap_sub * SUB)), # Total subclinical diagnosed (SUB)
        tCLN      = ((acf(floor(times)) * pop_target * pop_reached) * (alpha_cln * clnap_cln * CLN)), # Total clinical diagnosed (CLN)
        tTRE      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_tre * clnap_tre * TRE)), # Total treated diagnosed (TRE)
        cBAUs     = (((1-mdr) * (theta * CLN)) * bau_dstb), # Costs BAU DS-TB
        cBAUr     = (((mdr) * (theta * CLN)) * bau_drtb), # Costs BAU DR-TB
        cACF      = ((acf(floor(times)) * (((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + MIN + TRE)) + ((pop_target * pop_reached) * (SUB + CLN)))) * cm_screen_acf), # Total cost ACF
        cACFSUS   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * SUS)) * cm_screen_acf), # Total cost ACF for susceptible (SUS)
        cACFINF   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * INF)) * cm_screen_acf), # Total cost ACF for infected (INF)
        cACFCLE   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * CLE)) * cm_screen_acf), # Total cost ACF for cleared (CLE)
        cACFREC   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * REC)) * cm_screen_acf), # Total cost ACF for recovered (REC)
        cACFMIN   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * MIN)) * cm_screen_acf), # Total cost ACF for minimal (MIN)
        cACFSUB   = ((acf(floor(times)) * (pop_target * pop_reached * SUB)) * cm_screen_acf), # Total cost ACF for subclinical (SUB)
        cACFCLN   = ((acf(floor(times)) * (pop_target * pop_reached * CLN)) * cm_screen_acf), # Total cost ACF for clinical (CLN)
        cACFTRE   = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * TRE)) * cm_screen_acf), # Total cost ACF for treated (TRE)
        cCLAP     = ((acf(floor(times)) * (((pop_target * pop_reached * prop_sputum) * (alpha_sic * (SUS + INF + CLE) + alpha_rec * REC + alpha_min * MIN + alpha_tre * TRE)) + ((pop_target * pop_reached) * (alpha_sub * SUB + alpha_cln * CLN)))) * clnap_cost), # Total cost clinical appraisal        
        cCLAPSUS  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_sic * SUS)) * clnap_cost), # Total cost clinical appraisal for susceptible (SUS)
        cCLAPINF  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_sic * INF)) * clnap_cost), # Total cost clinical appraisal for infected (INF)
        cCLAPCLE  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_sic * CLE)) * clnap_cost), # Total cost clinical appraisal for cleared (CLE)
        cCLAPREC  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_rec * REC)) * clnap_cost), # Total cost clinical appraisal for recovered (REC)
        cCLAPMIN  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_min * MIN)) * clnap_cost), # Total cost clinical appraisal for minimal (MIN)
        cCLAPSUB  = ((acf(floor(times)) * (pop_target * pop_reached * alpha_sub * SUB)) * clnap_cost), # Total cost clinical appraisal for subclinical (SUB)
        cCLAPCLN  = ((acf(floor(times)) * (pop_target * pop_reached * alpha_cln * CLN)) * clnap_cost), # Total cost clinical appraisal for clinical (CLN)
        cCLAPTRE  = ((acf(floor(times)) * (pop_target * pop_reached * prop_sputum * alpha_tre * TRE)) * clnap_cost), # Total cost clinical appraisal for treated (TRE)
        cRxBAUs   = (((1-mdr) * (theta * CLN)) * rx_dstb), # Costs treatment BAU DS-TB
        cRxBAUr   = (((mdr) * (theta * CLN)) * rx_drtb), # Costs treatment BAU DR-TB
        cRxSUS    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * clnap_sic * SUS)) * rx_dstb), # Cost DS-TB treatment for susceptible (SUS)
        cRxINF    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * clnap_sic * INF)) * rx_dstb), # Cost DS-TB treatment for infected (INF)
        cRxCLE    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * clnap_sic * CLE)) * rx_dstb), # Cost DS-TB treatment for cleared (CLE)
        cRxREC    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_rec * clnap_rec * REC)) * rx_dstb), # Cost DS-TB treatment for recovered (REC)
        cRxMIN    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_min * clnap_min * MIN)) * rx_dstb), # Cost DS-TB treatment for minimal (MIN)
        cRxSUBs   = ((1-mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_sub * clnap_sub * SUB)) * rx_dstb)), # Cost DS-TB treatment for subclinical (SUB)
        cRxSUBr   = ((mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_sub * clnap_sub * SUB)) * rx_drtb)), # Cost DR-TB treatment for subclinical (SUB)
        cRxCLNs   = ((1-mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_cln * clnap_cln * CLN)) * rx_dstb)), # Cost DS-TB treatment for clinical (CLN)
        cRxCLNr   = ((mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_cln * clnap_cln * CLN)) * rx_drtb)), # Cost DR-TB treatment for clinical (CLN)
        cRxTRE    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_tre * clnap_tre * TRE)) * rx_dstb), # Cost DS-TB treatment for minimal (TRE)
        ARI       = ((beta / PopT) * ((kappa * SUB) + CLN)), # Annual risk of infection
        DALYs     = (daly(times) * (subcln * SUB)))) # DALY estimates
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
  "02" = seq(2025, 2026, 1),
  "03" = seq(2025, 2027, 1),
  "04" = seq(2025, 2028, 1),
  "05" = seq(2025, 2029, 1),
  "06" = seq(2025, 2030, 1),
  "08" = seq(2025, 2032, 1),
  "09" = seq(2025, 2033, 1),
  "11" = seq(2025, 2035, 1),
  "14" = seq(2025, 2038, 1),
  "16" = seq(2025, 2040, 1))

# 4.2.1 Main analysis
outacfa <- list()
outacfb <- list()
outacfc <- list()

# 4.2.2 Sensitivity analyses
outacfd <- list()

outacfax <- list()
outacfbx <- list()

outacfac <- list()
outacfbc <- list()
outacfcc <- list()

outacfbr <- list()
outacfcr <- list()

for (j in names(acf_year)) {
  rounds <- acf_year[[j]]
  print(names(acf_year[j]))
  
  if (j %in% c("03", "06", "11")) {
    outacfa[[j]] <- list()
  }
  if (j %in% c("03", "08", "14")) {
    outacfb[[j]] <- list()
  }
  if (j %in% c("02", "03", "04")) {
    outacfc[[j]] <- list()
  }
  if (j %in% c("03", "08", "14")) {
    outacfd[[j]] <- list()
  }
  if (j %in% c("03", "06", "11")) {
    outacfax[[j]] <- list()
  }
  if (j %in% c("03", "08", "14")) {
    outacfbx[[j]] <- list()
  }
  if (j %in% c("03", "06", "11")) {
    outacfac[[j]] <- list()
  }
  if (j %in% c("03", "08", "14")) {
    outacfbc[[j]] <- list()
  }
  if (j %in% c("02", "03", "04")) {
    outacfcc[[j]] <- list()
  }
  if (j %in% c("04", "09", "16")) {
    outacfbr[[j]] <- list()
  }
  if (j %in% c("02", "03", "05")) {
    outacfcr[[j]] <- list()
  }
  
  pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))
  
  for (i in 1:nrow(parms)) {
    curr_parms <- as.data.frame(parms[i,])
    curr_base <- as.data.frame(base[i,-1])
    
    if (j %in% c("03", "06", "11")) { # Main analysis: NAAT-only
      outacfa[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = acf_year[[j]]))
      outacfa[[j]][[i]] <- outacfa[[j]][[i]] %>% mutate(type = 'acfa', run = i, round = j)
    }

    if (j %in% c("03", "08", "14")) { # Main analysis: CXR+NAAT
      outacfb[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfb, acf_times = acf_year[[j]]))
      outacfb[[j]][[i]] <- outacfb[[j]][[i]] %>% mutate(type = 'acfb', run = i, round = j)
    }

    if (j %in% c("02", "03", "04")) { # Main analysis: CXR
      outacfc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfc, acf_times = acf_year[[j]]))
      outacfc[[j]][[i]] <- outacfc[[j]][[i]] %>% mutate(type = 'acfc', run = i, round = j)
    }
    
    if (j %in% c("03", "08", "14")) { # Sensitivity analysis: Xpert MTB/RIF
      outacfd[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfd, acf_times = acf_year[[j]]))
      outacfd[[j]][[i]] <- outacfd[[j]][[i]] %>% mutate(type = 'acfd', run = i, round = j)
    }
    
    if (j %in% c("03", "06", "11")) { # Sensitivity analysis: Reduced NAAT unit cost
      outacfax[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfax, acf_times = acf_year[[j]]))
      outacfax[[j]][[i]] <- outacfax[[j]][[i]] %>% mutate(type = 'acfax', run = i, round = j)
    }
    
    if (j %in% c("03", "08", "14")) { # Sensitivity analysis: Reduced NAAT unit cost
      outacfbx[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfbx, acf_times = acf_year[[j]]))
      outacfbx[[j]][[i]] <- outacfbx[[j]][[i]] %>% mutate(type = 'acfbx', run = i, round = j)
    }
    
    if (j %in% c("03", "06", "11")) { # Sensitivity analysis: Clinical appraisal
      outacfac[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfac, acf_times = acf_year[[j]], clnap = clnap_parms))
      outacfac[[j]][[i]] <- outacfac[[j]][[i]] %>% mutate(type = 'acfac', run = i, round = j)
    }
    
    if (j %in% c("03", "08", "14")) { # Sensitivity analysis: Clinical appraisal
      outacfbc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfbc, acf_times = acf_year[[j]], clnap = clnap_parms))
      outacfbc[[j]][[i]] <- outacfbc[[j]][[i]] %>% mutate(type = 'acfbc', run = i, round = j)
    }
    
    if (j %in% c("02", "03", "04")) { # Sensitivity analysis: Clinical appraisal
      outacfcc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfcc, acf_times = acf_year[[j]], clnap = clnap_parms))
      outacfcc[[j]][[i]] <- outacfcc[[j]][[i]] %>% mutate(type = 'acfcc', run = i, round = j)
    }
    
    if (j %in% c("04", "09", "16")) { # Sensitivity analysis: CXR test positivity
      outacfbr[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfbr, acf_times = acf_year[[j]]))
      outacfbr[[j]][[i]] <- outacfbr[[j]][[i]] %>% mutate(type = 'acfbr', run = i, round = j)
    }
    
    if (j %in% c("02", "03", "05")) { # Sensitivity analysis: CXR test positivity
      outacfcr[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfcr, acf_times = acf_year[[j]]))
      outacfcr[[j]][[i]] <- outacfcr[[j]][[i]] %>% mutate(type = 'acfcr', run = i, round = j)
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
  
  acfax_name <- paste0("r", names(outacfax[k]), "_outacfax_df")
  assign(acfax_name, do.call(rbind, outacfax[[k]]))
  export(get(acfax_name), here("outputs", "results", paste0(acfax_name, ".Rdata")))
  
  acfbx_name <- paste0("r", names(outacfbx[k]), "_outacfbx_df")
  assign(acfbx_name, do.call(rbind, outacfbx[[k]]))
  export(get(acfbx_name), here("outputs", "results", paste0(acfbx_name, ".Rdata")))
  
  acfac_name <- paste0("r", names(outacfac[k]), "_outacfac_df")
  assign(acfac_name, do.call(rbind, outacfac[[k]]))
  export(get(acfac_name), here("outputs", "results", paste0(acfac_name, ".Rdata")))
  
  acfbc_name <- paste0("r", names(outacfbc[k]), "_outacfbc_df")
  assign(acfbc_name, do.call(rbind, outacfbc[[k]]))
  export(get(acfbc_name), here("outputs", "results", paste0(acfbc_name, ".Rdata")))
  
  acfcc_name <- paste0("r", names(outacfcc[k]), "_outacfcc_df")
  assign(acfcc_name, do.call(rbind, outacfcc[[k]]))
  export(get(acfcc_name), here("outputs", "results", paste0(acfcc_name, ".Rdata")))
  
  acfbr_name <- paste0("r", names(outacfbr[k]), "_outacfbr_df")
  assign(acfbr_name, do.call(rbind, outacfbr[[k]]))
  export(get(acfbr_name), here("outputs", "results", paste0(acfbr_name, ".Rdata")))
  
  acfcr_name <- paste0("r", names(outacfcr[k]), "_outacfcr_df")
  assign(acfcr_name, do.call(rbind, outacfcr[[k]]))
  export(get(acfcr_name), here("outputs", "results", paste0(acfcr_name, ".Rdata")))
}

# 4.3 ACT3 comparison
act3_ultra <- list()
act3_mtbrif <- list()

pb <- progress_bar$new(format = "[:bar] :percent :eta", total = nrow(parms))

for (i in 1:nrow(parms)) {
  curr_parms <- as.data.frame(parms[i,])
  curr_base <- as.data.frame(base[i,-1])
  
  act3_ultra[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfa, acf_times = seq(2025, 2027, 1))) 
  act3_ultra[[i]] <- act3_ultra[[i]] %>% mutate(type = 'acfa', run = i, round = '03')
  
  act3_mtbrif[[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfd, acf_times = seq(2025, 2027, 1))) 
  act3_mtbrif[[i]] <- act3_mtbrif[[i]] %>% mutate(type = 'acfd', run = i, round = '03')
  
  pb$tick()
}

act3_ultra <- do.call(rbind, act3_ultra)
act3_mtbrif <- do.call(rbind, act3_mtbrif)

export(act3_ultra, here("outputs", "act3", "act3_ultra.Rdata"))
export(act3_mtbrif, here("outputs", "act3", "act3_mtbrif.Rdata"))

rm(list = ls())
