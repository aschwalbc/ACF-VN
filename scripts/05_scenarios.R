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
ultra_fp_rec <- c(lo = 0.026, hi = 0.070)
ultra_sens_min <- c(lo = 0.026, hi = 0.070)
ultra_sens_sub <- c(lo = 0.676, hi = 0.856)
ultra_sens_cln <- c(lo = 0.862, hi = 0.947)
ultra_fp_tre <- c(lo = 0.020, hi = 0.060)

# 2.2.2 Chest X-ray
cxr_fp_sic <- c(lo = 0.069, hi = 0.134)
cxr_fp_rec <- c(lo = 0.069, hi = 0.134)
cxr_sens_min <- c(lo = 0.770, hi = 0.900)
cxr_sens_sub <- c(lo = 0.900, hi = 0.920)
cxr_sens_cln <- c(lo = 0.900, hi = 0.920)
cxr_fp_tre <- c(lo = 0.481, hi = 0.524)

# 2.2.3 Xpert MTB/RIF (SA)
xpert_fp_sic <- c(lo = 0.0016, hi = 0.0029)
xpert_fp_rec <- c(lo = 0.007, hi = 0.030)
xpert_sens_min <- c(lo = 0.007, hi = 0.030)
xpert_sens_sub <- c(lo = 0.484, hi = 0.717)
xpert_sens_cln <- c(lo = 0.786, hi = 0.899)
xpert_fp_tre <- c(lo = 0.005, hi = 0.083)

# 2.3 ACF scenarios
# 2.3.1 ACF Scenario A - Ultra only, Ultra+ get treatment
acfa <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = ultra_fp_sic['lo'], hi = ultra_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = ultra_fp_rec['lo'], hi = ultra_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = ultra_sens_min['lo'], hi = ultra_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = ultra_sens_sub['lo'], hi = ultra_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = ultra_sens_cln['lo'], hi = ultra_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = ultra_fp_tre['lo'], hi = ultra_fp_tre['hi']) 
acfax <- acfa # Replicating for costing sensitivity analysis

# 2.3.2 ACF Scenario B - CXR only, CXR+ get Ultra, Ultra+ get treatment
acfb <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = ultra_fp_sic['lo']*cxr_fp_sic['lo'], hi = ultra_fp_sic['hi']*cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = ultra_fp_rec['lo']*cxr_fp_rec['lo'], hi = ultra_fp_rec['hi']*cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = ultra_sens_min['lo']*cxr_sens_min['lo'], hi = ultra_sens_min['hi']*cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = ultra_sens_sub['lo']*cxr_sens_sub['lo'], hi = ultra_sens_sub['hi']*cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = ultra_sens_cln['lo']*cxr_sens_cln['lo'], hi = ultra_sens_cln['hi']*cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = ultra_fp_tre['lo']*cxr_fp_tre['lo'], hi = ultra_fp_tre['hi']*cxr_fp_tre['hi'])
acfbx <- acfb # Replicating for costing sensitivity analysis

# 2.3.3 ACF Scenario C - CXR only, CXR+ get treatment
acfc <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = cxr_fp_sic['lo'], hi = cxr_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = cxr_fp_rec['lo'], hi = cxr_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = cxr_sens_min['lo'], hi = cxr_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = cxr_sens_sub['lo'], hi = cxr_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = cxr_sens_cln['lo'], hi = cxr_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = cxr_fp_tre['lo'], hi = cxr_fp_tre['hi'])

# 2.3.4 ACF SA Scenario - MTB/RIF only, MTB/RIF+ get treatment
acfd <- data.frame(parameter = character(), lo = numeric(), hi = numeric()) %>% 
  add_row(parameter = "alpha_sic", lo = xpert_fp_sic['lo'], hi = xpert_fp_sic['hi']) %>% 
  add_row(parameter = "alpha_rec", lo = xpert_fp_rec['lo'], hi = xpert_fp_rec['hi']) %>% 
  add_row(parameter = "alpha_min", lo = xpert_sens_min['lo'], hi = xpert_sens_min['hi']) %>% 
  add_row(parameter = "alpha_sub", lo = xpert_sens_sub['lo'], hi = xpert_sens_sub['hi']) %>%
  add_row(parameter = "alpha_cln", lo = xpert_sens_cln['lo'], hi = xpert_sens_cln['hi']) %>% 
  add_row(parameter = "alpha_tre", lo = xpert_fp_tre['lo'], hi = xpert_fp_tre['hi'])

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

# 2.7.2.2 Active case-finding per algorithm
scr_acfa <- gamma_dist(8)
scr_acfax <- c(usd1 = 3.0, usd2 = 3.8)
scr_acfb <- gamma_dist(1.7)
scr_acfbx <- c(usd1 = 1.3, usd2 = 1.4)
scr_acfc <- gamma_dist(1.2)
scr_acfd <- gamma_dist(8)

# 2.7.2.3 TB treatment
tb_rx_dstb <- gamma_dist(81)
tb_rx_drtb <- gamma_dist(973)

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
    
    if(interv_name %in% c('acfa', 'acfax', 'acfb', 'acfbx', 'acfd')) {
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
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfa[['shape']], scale = scr_acfa[['scale']])
    } else if(interv_name == "acfax") {
      cm_screen_acf <- runif(n = 1, min = scr_acfax[['usd1']], max = scr_acfax[['usd2']])
    } else if(interv_name == "acfb") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfb[['shape']], scale = scr_acfb[['scale']])
    } else if(interv_name == "acfbx") {
      cm_screen_acf <- runif(n = 1, min = scr_acfbx[['usd1']], max = scr_acfbx[['usd2']])
    } else if(interv_name == "acfc") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfc[['shape']], scale = scr_acfc[['scale']])
    } else if(interv_name == "acfd") {
      cm_screen_acf <- rgamma(n = 1, shape = scr_acfd[['shape']], scale = scr_acfd[['scale']])
    }
  }
  
  # BAU costs
  bau_dstb <- rgamma(n = 1, shape = scr_bau_dstb[['shape']], scale = scr_bau_dstb[['scale']])
  bau_drtb <- rgamma(n = 1, shape = scr_bau_drtb[['shape']], scale = scr_bau_drtb[['scale']])
  
  # Rx costs
  rx_dstb <- rgamma(n = 1, shape = tb_rx_dstb[['shape']], scale = tb_rx_dstb[['scale']])
  rx_drtb <- rgamma(n = 1, shape = tb_rx_drtb[['shape']], scale = tb_rx_drtb[['scale']])

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
        rTBc      = ((SUB + CLN) / PopT * 1e5), # Infectious TB (per 100k)
        tTBc      = (SUB + CLN), # Total infectious TB
        rInc      = ((subcln * SUB)/ PopT * 1e5), # TB incidence (per 100k)
        tInc      = (subcln * SUB), # Total TB incidence
        rMor      = ((mutb * CLN) / PopT * 1e5), # Clinical TB mortality per time (per 100k)
        tMor      = (mutb * CLN), # Clinical TB mortality per time
        tScrn     = ((acf(floor(times)) * (((pop_target * pop_reached * prop_sputum) * (SUS + INF + CLE + REC + MIN + TRE)) + ((pop_target * pop_reached) * (SUB + CLN))))), # Total number screened
        tSUS      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * SUS)), # Total susceptible diagnosed (SUS)
        tINF      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * INF)), # Total infected diagnosed (INF)
        tCLE      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_sic * CLE)), # Total cleared diagnosed (CLE)
        tREC      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_rec * REC)), # Total recovered diagnosed (REC)
        tMIN      = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum) * (alpha_min * MIN)), # Total minimal diagnosed (MIN)
        tSUB      = ((acf(floor(times)) * pop_target * pop_reached) * (alpha_sub * SUB)), # Total subclinical diagnosed (SUB)
        tCLN      = ((acf(floor(times)) * pop_target * pop_reached) * (alpha_cln * CLN)), # Total clinical diagnosed (CLN)
        tTRE      = ((acf(floor(times)) * pop_target * pop_reached) * (alpha_tre * TRE)), # Total treated diagnosed (TRE)
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
        cRxBAUs   = (((1-mdr) * (theta * CLN)) * rx_dstb), # Costs treatment BAU DS-TB
        cRxBAUr   = (((mdr) * (theta * CLN)) * rx_drtb), # Costs treatment BAU DR-TB
        cRxSUS    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * SUS)) * rx_dstb), # Cost DS-TB treatment for susceptible (SUS)
        cRxINF    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * INF)) * rx_dstb), # Cost DS-TB treatment for infected (INF)
        cRxCLE    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_sic * CLE)) * rx_dstb), # Cost DS-TB treatment for cleared (CLE)
        cRxREC    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_rec * REC)) * rx_dstb), # Cost DS-TB treatment for recovered (REC)
        cRxMIN    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_min * MIN)) * rx_dstb), # Cost DS-TB treatment for minimal (MIN)
        cRxSUBs   = ((1-mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_sub * SUB)) * rx_dstb)), # Cost DS-TB treatment for subclinical (SUB)
        cRxSUBr   = ((mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_sub * SUB)) * rx_drtb)), # Cost DR-TB treatment for subclinical (SUB)
        cRxCLNs   = ((1-mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_cln * CLN)) * rx_dstb)), # Cost DS-TB treatment for clinical (CLN)
        cRxCLNr   = ((mdr) * ((acf(floor(times)) * pop_target * pop_reached * (alpha_cln * CLN)) * rx_drtb)), # Cost DR-TB treatment for clinical (CLN)
        cRxTRE    = ((acf(floor(times)) * pop_target * pop_reached * prop_sputum * (alpha_tre * TRE)) * rx_dstb), # Cost DS-TB treatment for minimal (TRE)
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
  "01" = 2025,
  "02" = seq(2025, 2026, 1),
  "03" = seq(2025, 2027, 1),
  "06" = seq(2025, 2030, 1),
  "07" = seq(2025, 2031, 1),
  "08" = seq(2025, 2032, 1), 
  "11" = seq(2025, 2035, 1),
  "12" = seq(2025, 2036, 1),
  "14" = seq(2025, 2038, 1))

outacfa <- list()
outacfax <- list()
outacfb <- list()
outacfbx <- list()
outacfc <- list()
outacfd <- list()

for (j in names(acf_year)) {
  rounds <- acf_year[[j]]
  print(names(acf_year[j]))
  
  if (j %in% c("03", "06", "11")) {
    outacfa[[j]] <- list()
  }
  if (j %in% c("03", "06", "11")) {
    outacfax[[j]] <- list()
  }
  if (j %in% c("03", "07", "12")) {
    outacfb[[j]] <- list()
  }
  if (j %in% c("03", "07", "12")) {
    outacfbx[[j]] <- list()
  }
  if (j %in% c("01", "02", "03")) {
    outacfc[[j]] <- list()
  }
  if (j %in% c("03", "08", "14")) {
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

    if (j %in% c("03", "06", "11")) {
      outacfax[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfax, acf_times = acf_year[[j]]))
      outacfax[[j]][[i]] <- outacfax[[j]][[i]] %>% mutate(type = 'acfax', run = i, round = j)
    }

    if (j %in% c("03", "07", "12")) {
      outacfb[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfb, acf_times = acf_year[[j]]))
      outacfb[[j]][[i]] <- outacfb[[j]][[i]] %>% mutate(type = 'acfb', run = i, round = j)
    }

    if (j %in% c("03", "07", "12")) {
      outacfbx[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfbx, acf_times = acf_year[[j]]))
      outacfbx[[j]][[i]] <- outacfbx[[j]][[i]] %>% mutate(type = 'acfbx', run = i, round = j)
    }

    if (j %in% c("01", "02", "03")) {
      outacfc[[j]][[i]] <- as.data.frame(ode(parms = curr_parms, base = curr_base, interv = acfc, acf_times = acf_year[[j]]))
      outacfc[[j]][[i]] <- outacfc[[j]][[i]] %>% mutate(type = 'acfc', run = i, round = j)
    }
    
    if (j %in% c("03", "08", "14")) {
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

  acfax_name <- paste0("r", names(outacfax[k]), "_outacfax_df")
  assign(acfax_name, do.call(rbind, outacfax[[k]]))
  export(get(acfax_name), here("outputs", "results", paste0(acfax_name, ".Rdata")))

  acfb_name <- paste0("r", names(outacfb[k]), "_outacfb_df")
  assign(acfb_name, do.call(rbind, outacfb[[k]]))
  export(get(acfb_name), here("outputs", "results", paste0(acfb_name, ".Rdata")))

  acfbx_name <- paste0("r", names(outacfbx[k]), "_outacfbx_df")
  assign(acfbx_name, do.call(rbind, outacfbx[[k]]))
  export(get(acfbx_name), here("outputs", "results", paste0(acfbx_name, ".Rdata")))

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
