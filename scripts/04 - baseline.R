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

# 1. Load data ==========
par <- import(here("data","parameters_manual.Rdata")) # Manual fit parameters
# par <- import(here("data","parameters.Rdata")) # HMER fit parameters
WPP <- import(here("data","pop","WPP_Pop_1950-2100.csv")) # Population size 1950-2100
WPPb <- import(here("data","pop","WPP_Births_1950-2100.csv")) # Births 1950-2100
WPPda <- import(here("data","pop","WPP_Deaths_1950-2021.csv")) # Deaths 1950-2021
WPPdb <- import(here("data","pop","WPP_Deaths_2022-2100.csv")) # Deaths 2022-2100
WUP <- import(here("data","pop","WUP_Urban_1950-2050.csv")) # Proportion Urban 1950-2050
WDI <- import(here("data","pop","WB_WDI.csv")) # GDP World Development Indicators
WEO <- import(here("data","pop","WEO.csv")) # World Economic Outputs

# 2. Model prep ==========
# 2.1 Fitted parameters (Posterior)
val <- c('med','low',"hig")

parms <- c(
  beta = par[par$parameter == 'beta', val[1]],                   # Contact (per person/year) parameter
  kappa = par[par$parameter == 'kappa', val[1]],                 # Relative infectiousness
  gamma_infcle = par[par$parameter == 'gamma_infcle', val[1]],   # REG: Infected -> Cleared
  lambda_infmin = par[par$parameter == 'lambda_infmin', val[1]], # PROG: Infected -> Minimal
  gamma_mincle = par[par$parameter == 'gamma_mincle', val[1]],   # REG: Minimal -> Cleared
  theta_cleinf = par[par$parameter == 'theta_cleinf', val[1]],   # REINF: Cleared -> Infected
  lambda_minsub = par[par$parameter == 'lambda_minsub', val[1]], # PROG: Minimal -> Subclinical
  lambda_infsub = par[par$parameter == 'lambda_infsub', val[1]], # PROG: Infected -> Subclinical
  gamma_submin = par[par$parameter == 'gamma_submin', val[1]],   # REG: Subclinical -> Minimal
  lambda_subcln = par[par$parameter == 'lambda_subcln', val[1]], # PROG: Subclinical -> Clinical
  gamma_clnsub = par[par$parameter == 'gamma_clnsub', val[1]],   # REG: Clinical -> Subclinical
  omega_ini = par[par$parameter == 'omega_ini', val[1]],         # Mortality (Initial)
  omega_fin = par[par$parameter == 'omega_fin', val[1]],         # Mortality (Final)
  iota_cln_ini = par[par$parameter == 'iota_cln_ini', val[1]],   # Diagnosis Clinical (Initial)
  iota_cln_fin = par[par$parameter == 'iota_cln_fin', val[1]],   # Diagnosis Clinical (Final)
  phi_cln_ini = par[par$parameter == 'phi_cln_ini', val[1]],     # Treatment failure Clinical (Initial) 
  phi_cln_fin = par[par$parameter == 'phi_cln_fin', val[1]],     # Treatment failure Clinical (Final)
  tau_min = par[par$parameter == 'tau_min', val[1]],             # RELAP: Recovered -> Minimal 
  tau_sub = par[par$parameter == 'tau_sub', val[1]],             # RELAP: Recovered -> Subclinical
  rho_ini = par[par$parameter == 'rho_ini', val[1]],             # Proportion rural (Initial)
  rho_fin = par[par$parameter == 'rho_fin', val[1]],             # Proportion rural (Final)
  sigma_ini = par[par$parameter == 'sigma_ini', val[1]],         # Proportion low SES (Initial)
  sigma_fin = par[par$parameter == 'sigma_fin', val[1]])         # Proportion low SES (Final)
rm(par, val)

# 2.2 Population and demographics
WPP <- WPP %>% # Population data
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrp, pop = PopTotal) %>%
  select(iso3, year, agegp, pop) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020 & year <= 2050) %>% 
  filter(!agegp %in% c('0-4','5-9','10-14')) %>%
  group_by(iso3, year) %>% 
  summarise(pop = sum(pop)*1e3)

WPPb <- WPPb %>% # Birth data
  rename(iso3 = ISO3_code, year = Time, births = Births) %>%
  select(iso3, year, births) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020 & year <= 2050) %>%
  mutate(births = births*1e3)

WPPda <- WPPda %>% # Mortality data (1950-2021)
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrpStart, mort = DeathTotal) %>%
  select(iso3, year, agegp, mort) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020) %>% 
  filter(agegp >= 15) %>%
  group_by(iso3, year) %>%
  summarise(mort = sum(mort)*1e3)

WPPdb <- WPPdb %>% # Mortality data (2022-2100)
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrpStart, mort = DeathTotal) %>%
  select(iso3, year, agegp, mort) %>%
  filter(iso3 == "VNM") %>%
  filter(year <= 2050) %>% 
  filter(agegp >= 15) %>%
  group_by(iso3, year) %>%
  summarise(mort = sum(mort)*1e3)

WPPd <- rbind(WPPda,WPPdb)
rm(WPPda,WPPdb)

WPP <- WPP %>%
  left_join(WPPb, by=c("iso3","year")) %>%
  inner_join(WPPd, by=c("iso3","year")) %>%
  mutate(mortrate = mort/pop, birthrate = births/pop) 
rm(WPPb,WPPd)

export(WPP,here("data","pop","WPP.Rdata")) # Save data frame

WUP <- clean_names(WUP) %>% # Urbanisation data
  filter(index == 115) %>% 
  select(starts_with("x")) %>%
  rename_all(~gsub("^x", "", .)) %>% 
  mutate_all(~ . / 100) %>% 
  pivot_longer(cols = everything(), names_to = "year", values_to = "urbprop") %>% 
  filter(year >= 2020) %>% 
  mutate(rurprop = 1-urbprop, iso3 = "VNM") %>% 
  select(iso3, year, urbprop, rurprop)

export(WUP,here("data","pop","WUP.Rdata")) # Save data frame

WEO <- WEO %>% # World Economic Output (GDP)
  setNames(WEO[1,]) %>%
  slice(2:n()) %>% 
  clean_names() %>% 
  rename_all(~gsub("^x", "", .)) %>% 
  rename(iso3 = iso, var = subject_descriptor) %>% 
  select(iso3, var, matches("^\\d")) %>% 
  mutate_at(vars(-iso3, -var), ~as.numeric(gsub(",", "", gsub("\\.", "", ., fixed = TRUE)))) %>% 
  filter(var == "Gross domestic product per capita, constant prices") %>% 
  select(-var) %>% 
  pivot_longer(cols = -iso3, names_to = "year", values_to = "gdp") %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(lowses = NA) %>% 
  mutate(lowses = case_when(year == 2008 ~ 0.5946, year == 2018 ~ 0.6034, TRUE ~ lowses))
  
gdpses <- WEO %>% 
  filter(!is.na(lowses)) %>% 
  lm(lowses ~ gdp, data = .)

WEO <- WEO %>% 
  mutate(lowses = predict(gdpses, newdata = .))

sesmodel <- lm(lowses ~ year, data = WEO)

WEO_exp <- data.frame(year = 2029:2050) %>% 
  mutate(lowses = predict(sesmodel, newdata = .))

GDP <- rbind(select(WEO, c(year, lowses)), WEO_exp) %>% 
  filter(year >= 2020)

# 3. Model ==========
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

basemod <- as.data.frame(ode(parms)) # Run baseline model
rm(parms, ode)

# 4. Results ==========
# Baseline population
basepop <- basemod %>% 
  filter(time == 2020) %>% # Filter data for 2020
  select(contains('_')) %>% # Remove time variable
  mutate(across(everything(), ~ . / 100000)) %>% 
  mutate(Pop = as.numeric(WPP[WPP$year == 2020, 'pop'])) %>%  # Add adult population in 2020
  mutate(across(contains('_'), ~ . * Pop)) %>% # Calculate population per compartment
  select(-Pop) # Remove variable

export(basepop, here("data","basepop.Rdata")) # Save data frame
rm(basemod, basepop, WPP)
