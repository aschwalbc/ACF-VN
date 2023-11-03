## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 00: Extra.R

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
parms <- import(here("outputs", "pts","fitpts.Rdata")) %>% rename(gamma_minrec = gamma_mincle) # Fitted parameters
WPP <- import(here("data","pop","WPP_Pop_1950-2100.csv")) # Population size 1950-2100
WPPb <- import(here("data","pop","WPP_Births_1950-2100.csv")) # Births 1950-2100
WPPd <- import(here("data","pop","WPP_Deaths_1950-2021.csv")) # Deaths 1950-2021

# 2. Data curation ==========
WPP <- WPP %>% # Population data
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrp, pop = PopTotal) %>%
  select(iso3, year, agegp, pop) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2000 & year <= 2020) %>% 
  filter(!agegp %in% c('0-4','5-9','10-14')) %>%
  group_by(iso3, year) %>% 
  summarise(pop = sum(pop)*1e3)

WPPb <- WPPb %>% # Birth data
  rename(iso3 = ISO3_code, year = Time, births = Births) %>%
  select(iso3, year, births) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2000 & year <= 2020) %>% 
  mutate(births = births*1e3)

WPPd <- WPPd %>% # Mortality data (1950-2021)
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrpStart, mort = DeathTotal) %>%
  select(iso3, year, agegp, mort) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2000) %>% 
  filter(agegp >= 15) %>%
  group_by(iso3, year) %>%
  summarise(mort = sum(mort)*1e3)

WPP <- WPP %>%
  left_join(WPPb, by=c("iso3","year")) %>%
  inner_join(WPPd, by=c("iso3","year")) %>%
  mutate(mortrate = mort/pop, birthrate = births/pop) 
rm(WPPb,WPPd)

# 3. Model ==========
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
      dO_RL  = gamma_minrec*M_RL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_recinf) - mu*O_RL
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_minrec*M_RL - lambda_minsub*M_RL + tau_min*P_RL - mu*M_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL + tau_sub*P_RL - mu*S_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - force_func_omega(time)*C_RL - mu*C_RL - force_func_iota(time)*C_RL + force_func_phi(time)*RC_RL
      dRC_RL = force_func_iota(time)*C_RL - force_func_phi(time)*RC_RL - delta*RC_RL - mu*RC_RL
      dP_RL  = delta*RC_RL - tau_min*P_RL - tau_sub*P_RL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_treinf) - mu*P_RL
      
      dN_RH  = force_func_rhoR(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu*N_RH
      dI_RH  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(W_RH*theta_cleinf)+(O_RH*theta_recinf)+(P_RH*theta_treinf))) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu*I_RH
      dW_RH  = gamma_infcle*I_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_RH*theta_cleinf) - mu*W_RH
      dO_RH  = gamma_minrec*M_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_recinf) - mu*O_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_minrec*M_RH - lambda_minsub*M_RH + tau_min*P_RH - mu*M_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH + tau_sub*P_RH - mu*S_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - force_func_omega(time)*C_RH - mu*C_RH - force_func_iota(time)*C_RH + force_func_phi(time)*RC_RH
      dRC_RH = force_func_iota(time)*C_RH - force_func_phi(time)*RC_RH - delta*RC_RH - mu*RC_RH
      dP_RH  = delta*RC_RH - tau_min*P_RH - tau_sub*P_RH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_treinf) - mu*P_RH
      
      dN_UL  = force_func_rhoU(time)*force_func_sigmaL(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu*N_UL
      dI_UL  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(W_UL*theta_cleinf)+(O_UL*theta_recinf)+(P_UL*theta_treinf))) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu*I_UL
      dW_UL  = gamma_infcle*I_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_UL*theta_cleinf) - mu*W_UL
      dO_UL  = gamma_minrec*M_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_recinf) - mu*O_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_minrec*M_UL - lambda_minsub*M_UL + tau_min*P_UL - mu*M_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL + tau_sub*P_UL - mu*S_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - force_func_omega(time)*C_UL - mu*C_UL - force_func_iota(time)*C_UL + force_func_phi(time)*RC_UL
      dRC_UL = force_func_iota(time)*C_UL - force_func_phi(time)*RC_UL - delta*RC_UL - mu*RC_UL
      dP_UL  = delta*RC_UL - tau_min*P_UL - tau_sub*P_UL - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_treinf) - mu*P_UL
      
      dN_UH  = force_func_rhoU(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu*N_UH
      dI_UH  = (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(W_UH*theta_cleinf)+(O_UH*theta_recinf)+(P_UH*theta_treinf))) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu*I_UH
      dW_UH  = gamma_infcle*I_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*W_UH*theta_cleinf) - mu*W_UH
      dO_UH  = gamma_minrec*M_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_recinf) - mu*O_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_minrec*M_UH - lambda_minsub*M_UH + tau_min*P_UH - mu*M_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH + tau_sub*P_UH - mu*S_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - force_func_omega(time)*C_UH - mu*C_UH - force_func_iota(time)*C_UH + force_func_phi(time)*RC_UH
      dRC_UH = force_func_iota(time)*C_UH - force_func_phi(time)*RC_UH - delta*RC_UH - mu*RC_UH
      dP_UH  = delta*RC_UH - tau_min*P_UH - tau_sub*P_UH - (((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_treinf) - mu*P_UH
  
      return(list(c( 
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dW_RL, dW_RH, dW_UL, dW_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH,
        dS_RL, dS_RH, dS_UL, dS_UH, dC_RL, dC_RH, dC_UL, dC_UH, dRC_RL, dRC_RH, dRC_UL, dRC_UH, dP_RL, dP_RH, dP_UL, dP_UH),
        TBc   = (S_RL+S_RH+S_UL+S_UH+C_RL+C_RH+C_UL+C_UH),
        Dxs   = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH)),
        Dxout  = (force_func_iota(time)),
        PDRcln = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH))/(C_RL+C_RH+C_UL+C_UH),
        PDRinf = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH))/(C_RL+C_RH+C_UL+C_UH+S_RL+S_RH+S_UL+S_UH),
        PDRall = (force_func_iota(time)*(C_RL+C_RH+C_UL+C_UH))/(C_RL+C_RH+C_UL+C_UH+S_RL+S_RH+S_UL+S_UH+M_RL+M_RH+M_UL+M_UH)))
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
resi <- basemod %>% 
  select(time, run, TBc, Dxs, Dxout, PDRcln, PDRinf, PDRall) %>% 
  pivot_longer(cols = -c(time, run), names_to = "var", values_to = "values") %>% 
  group_by(time, var) %>% 
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE)) %>% 
  filter(time %in% c(2010, 2019))

# TB prev reduction
prevrate2008 <- c(val = 250, lo = 202, hi = 310)
prevrate2018 <- c(val = 227, lo = 177, hi = 290)

n <- 10000

prevrate2008_samples <- rnorm(n, mean = prevrate2008["val"], sd = (prevrate2008["hi"] - prevrate2008["lo"])/3.92)
prevrate2018_samples <- rnorm(n, mean = prevrate2018["val"], sd = (prevrate2018["hi"] - prevrate2018["lo"])/3.92)

df <- data.frame(run = 1:n, prevrate2008 = prevrate2008_samples, prevrate2018 = prevrate2018_samples)

red <- df %>% 
  mutate(prevred = prevrate2018/prevrate2008) %>% 
  summarise(val = median(prevred), lo = quantile(prevred, 0.025), hi = quantile(prevred, 0.975))
