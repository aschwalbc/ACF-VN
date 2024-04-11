## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 04: Baseline.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(janitor) # Cleaning dirty data
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame
library(progress) # Displays progress bar

# 1. Load data ==========
parms <- import(here("outputs", "pts","fitpts.Rdata"))[1:1000,] # Fitted parameters
WPP <- import(here("data","pop","WPP.Rdata")) # Population data

# 3. Model ==========
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

basemod <- data.frame()
pb <- progress_bar$new(format = "[:bar] :percent", total = nrow(parms))

for(i in 1:nrow(parms)) {
  output <- as.data.frame(ode((unlist(parms[i,]))))
  output <- output %>% mutate(run = i)
  basemod <- bind_rows(basemod, output)
  pb$tick()
}
rm(parms, ode, output, i, pb)

# 4. Results ==========
# Baseline population
base <- basemod %>% 
  filter(time == 2020) %>% # Filter data for 2020
  within(rm(time, TBc, Mor, Dxs, Spr)) %>% # Remove time variable
  mutate(across(-run, ~ . / 100000)) %>% 
  mutate(POP = as.numeric(WPP[WPP$year == 2020, 'pop'])) %>%  # Add adult population in 2020
  mutate(across(-c(run, POP), ~ . * POP)) %>% # Calculate population per compartment
  select(-POP) %>% # Remove variable
  relocate(run)

export(base, here("data", "fit", "base.Rdata")) # Save data frame
rm(basemod, base, WPP)
