## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 00: Manual Fit.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(deSolve) # Solvers for ordinary differential equations
library(reshape2) # Reshaping data easily
library(purrr) # Complete set of tools for functions and vectors
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame, allows use of j operator (:=)
library(patchwork) # Plot composition
library(gridExtra) # Miscellaneous functions for grid graphics
library(cowplot) # Publication quality figures

# 1. Parameters ==========
# 1.1 Baseline parameters (for testing)
parms = c(
  beta = 9,             # Contact (per person/year) parameter
  kappa = 0.75,         # Relative infectiousness
  gamma_infcle = 1.83,  # REG: Infected -> Cleared
  lambda_infmin = 0.21, # PROG: Infected -> Minimal
  gamma_mincle = 0.16,  # REG: Minimal -> Cleared
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
  tau_min = 0.005,      # RELAP: Recovered -> Minimal 
  tau_sub = 0.005,      # RELAP: Recovered -> Subclinical
  rho_ini = 0.78,       # Proportion rural (Initial)
  rho_fin = 0.04,       # Proportion rural (Final)
  sigma_ini = 0.76,     # Proportion low SES (Initial)
  sigma_fin = 0.16)     # Proportion low SES (Final)

# 1.2 Parameter ranges
ranges = list(
  beta = c(5,10),                 # Contact (per person/year) parameter (Horton et al. 2022 - PLOS GPH)
  kappa = c(0.5,1),               # Relative infectiousness (Emery et al. 2022)
  gamma_infcle = c(0.93,3.30),    # REG: Infected -> Cleared (Horton et al. 2023)
  lambda_infmin = c(0.04,0.23),   # PROG: Infected -> Minimal (Horton et al. 2023)
  gamma_mincle = c(0.14,0.23),    # REG: Minimal -> Cleared (Horton et al. 2023)
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
  tau_min = c(0.0,0.01),          # RELAP: Recovered -> Minimal 
  tau_sub = c(0.0,0.01),          # RELAP: Recovered -> Subclinical
  rho_ini = c(0.7,0.85),          # Proportion rural (Initial)
  rho_fin = c(0,0.2),             # Proportion rural (Final)
  sigma_ini = c(0.7,0.85),        # Proportion low SES (Inital)
  sigma_fin = c(0,0.2))           # Proportion low SES (Final)

# 1.3 Fitting targets
targets <- list(
  TBc2008 = c(159.2,238.8), # TB prevalence (Nguyen et al. Emerg Infect Dis 2021)
  # TBc2018 = c(98,159),      # TB prevalence (Nguyen et al. Emerg Infect Dis 2021)
  TBc2018 = c(115,231.4),   # TB decline in men (-18%; 95%CI: -42.2 to +16.3) from first survey estimate (199)
  Mor2000 = c(37.3,87.7),   # TB mortality (Over population aged >= 15yo)
  # Mor2000 = c(25.4,59.8),   # TB mortality (Over all population)
  Mor2010 = c(22.8,45.7),   # TB mortality (Over population aged >= 15yo)
  # Mor2010 = c(17.3,34.5),   # TB mortality (Over all population)
  Dxs2010 = c(63.5,95.3),   # TB notifications (Over population aged >= 15yo)
  # Dxs2010 = c(48,72),       # TB notifications (Over all population)
  Dxs2020 = c(58.5,87.8),   # TB notifications (Over population aged >= 15yo)
  # Dxs2020 = c(45.2,67.8),   # TB notifications (Over all population)
  # URs2008 = c(0.192,0.288), # Urban/rural scTB
  # URs2018 = c(0.352,0.528), # Urban/rural scTB
  # URc2008 = c(0.224,0.336), # Urban/rural cTB
  # URc2018 = c(0.344,0.516), # Urban/rural cTB
  # HLs2008 = c(0.304,0.456), # High/low scTB
  # HLs2018 = c(0.376,0.564), # High/low scTB
  # HLc2008 = c(0.256,0.384), # High/low cTB
  # HLc2018 = c(0.36,0.54),   # High/low cTB
  Spr2008 = c(0.56,0.84),   # Proportion scTB (Emery et al. medRxiv 2022)
  Spr2018 = c(0.53,0.79))   # Proportion scTB (Emery et al. medRxiv 2022)

targetsdb <- as.data.frame(t(as.data.frame(targets))) # Create dataframe with fitting targets
targetsdb$var <- substr(rownames(targetsdb),1,3) # Create variable classification
targetsdb$time <- as.numeric(substr(rownames(targetsdb),4,7)) # Create time variable
colnames(targetsdb) <- c("lo","hi","var","time") # Rename columns

# 2. Model ==========
ode <- function(parms, end_time = 2020) {
  
  # Static parameters  
  N <- 1e5 # Population size
  mu <- 1/70 # Age expectancy adult
  delta <- 1 # Treatment year
  theta_recinf <- 1 # REINF: Recovered -> Infected (No protection)
  
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
      dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL + tau_min*P_RL - mu*M_RL
      dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL + tau_sub*P_RL - mu*S_RL
      dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - force_func_omega(time)*C_RL - mu*C_RL - force_func_iota(time)*C_RL + force_func_phi(time)*RC_RL
      dRC_RL = force_func_iota(time)*C_RL - force_func_phi(time)*RC_RL - delta*RC_RL - mu*RC_RL
      dP_RL  = delta*RC_RL - tau_min*P_RL - tau_sub*P_RL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu*P_RL
      
      dN_RH  = force_func_rhoR(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RH - mu*N_RH
      dI_RH  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RH+(O_RH*theta_cleinf)+(P_RH*theta_recinf)) - gamma_infcle*I_RH - lambda_infmin*I_RH - lambda_infsub*I_RH - mu*I_RH
      dO_RH  = gamma_infcle*I_RH + gamma_mincle*M_RH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RH*theta_cleinf - mu*O_RH
      dM_RH  = lambda_infmin*I_RH + gamma_submin*S_RH - gamma_mincle*M_RH - lambda_minsub*M_RH + tau_min*P_RH - mu*M_RH
      dS_RH  = lambda_infsub*I_RH + lambda_minsub*M_RH + gamma_clnsub*C_RH - gamma_submin*S_RH - lambda_subcln*S_RH + tau_sub*P_RH - mu*S_RH
      dC_RH  = lambda_subcln*S_RH - gamma_clnsub*C_RH - force_func_omega(time)*C_RH - mu*C_RH - force_func_iota(time)*C_RH + force_func_phi(time)*RC_RH
      dRC_RH = force_func_iota(time)*C_RH - force_func_phi(time)*RC_RH - delta*RC_RH - mu*RC_RH
      dP_RH  = delta*RC_RH - tau_min*P_RH - tau_sub*P_RH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RH*theta_recinf - mu*P_RH
      
      dN_UL  = force_func_rhoU(time)*force_func_sigmaL(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UL - mu*N_UL
      dI_UL  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UL+(O_UL*theta_cleinf)+(P_UL*theta_recinf)) - gamma_infcle*I_UL - lambda_infmin*I_UL - lambda_infsub*I_UL - mu*I_UL
      dO_UL  = gamma_infcle*I_UL + gamma_mincle*M_UL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UL*theta_cleinf - mu*O_UL
      dM_UL  = lambda_infmin*I_UL + gamma_submin*S_UL - gamma_mincle*M_UL - lambda_minsub*M_UL + tau_min*P_UL - mu*M_UL
      dS_UL  = lambda_infsub*I_UL + lambda_minsub*M_UL + gamma_clnsub*C_UL - gamma_submin*S_UL - lambda_subcln*S_UL + tau_sub*P_UL - mu*S_UL
      dC_UL  = lambda_subcln*S_UL - gamma_clnsub*C_UL - force_func_omega(time)*C_UL - mu*C_UL - force_func_iota(time)*C_UL + force_func_phi(time)*RC_UL
      dRC_UL = force_func_iota(time)*C_UL - force_func_phi(time)*RC_UL - delta*RC_UL - mu*RC_UL
      dP_UL  = delta*RC_UL - tau_min*P_UL - tau_sub*P_UL - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UL*theta_recinf - mu*P_UL
      
      dN_UH  = force_func_rhoU(time)*force_func_sigmaH(time)*(mu*N + force_func_omega(time)*(C_RL+C_RH+C_UL+C_UH)) - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_UH - mu*N_UH
      dI_UH  = ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_UH+(O_UH*theta_cleinf)+(P_UH*theta_recinf)) - gamma_infcle*I_UH - lambda_infmin*I_UH - lambda_infsub*I_UH - mu*I_UH
      dO_UH  = gamma_infcle*I_UH + gamma_mincle*M_UH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_UH*theta_cleinf - mu*O_UH
      dM_UH  = lambda_infmin*I_UH + gamma_submin*S_UH - gamma_mincle*M_UH - lambda_minsub*M_UH + tau_min*P_UH - mu*M_UH
      dS_UH  = lambda_infsub*I_UH + lambda_minsub*M_UH + gamma_clnsub*C_UH - gamma_submin*S_UH - lambda_subcln*S_UH + tau_sub*P_UH - mu*S_UH
      dC_UH  = lambda_subcln*S_UH - gamma_clnsub*C_UH - force_func_omega(time)*C_UH - mu*C_UH - force_func_iota(time)*C_UH + force_func_phi(time)*RC_UH
      dRC_UH = force_func_iota(time)*C_UH - force_func_phi(time)*RC_UH - delta*RC_UH - mu*RC_UH
      dP_UH  = delta*RC_UH - tau_min*P_UH - tau_sub*P_UH - ((beta/N)*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_UH*theta_recinf - mu*P_UH
      
      return(list(c( 
        dN_RL, dN_RH, dN_UL, dN_UH, dI_RL, dI_RH, dI_UL, dI_UH, dO_RL, dO_RH, dO_UL, dO_UH, dM_RL, dM_RH, dM_UL, dM_UH,
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
            O_RL = 0, O_RH = 0, O_UL = 0, O_UH = 0, M_RL = 0, M_RH = 0, M_UL = 0, M_UH = 0,
            S_RL = 0, S_RH = 0, S_UL = 0, S_UH = 0, C_RL = 250, C_RH = 250, C_UL = 250, C_UH = 250, 
            RC_RL = 0, RC_RH = 0, RC_UL = 0, RC_UH = 0, P_RL = 0, P_RH = 0, P_UL = 0, P_UH = 0)
  
  times <- seq(1500, end_time, by = 1)
  out <- deSolve::ode(yini, times, des, parms)
  return(out)
}

# Run ode
oderes <- as.data.frame(ode(parms))

# 3. Plots ==========
# 3.1 Target plots
results <- oderes %>% 
  select(time, TBc, Spr, Mor, Dxs, URs, URc, HLs, HLc) %>% 
  pivot_longer(cols = c(TBc, Spr, Mor, Dxs, URs, URc, HLs, HLc), names_to = "var", values_to = "val")

title_size = 13
axis_title_size = 13
x_axis_text_size = 8
y_axis_text_size = 10
title_align = 0.5
title_face = "bold"
line_thickness = 0.6
colour_line = "#CE2931"
point_size = 1
alpha = 0.2
patches_height = 1.6*16
height = 16
aspect_ratio = 1
error_bar_width = 0.5
date_line_thickness = 0.4
date_line_colour = "gray35"
prior_colour = "deepskyblue2"

breaks = seq(2000, 2020, 2)
limits = c(1999.5, 2020.5)
from_to = seq(2000,2020)

a = ggplot(subset(results, var == "TBc" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "TBc"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "TBc"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) + 
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, 50), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,300)) +
  labs(title = "",
       x = "Year",
       y = "TB prevalence rate") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

b = ggplot(subset(results, var == "Spr" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "Spr"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "Spr"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "",
       x = "Year",
       y = "Proportion scTB") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

c = ggplot(subset(results, var == "Mor" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "Mor"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "Mor"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 60, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,60)) +
  labs(title = "",
       x = "Year",
       y = "TB mortality rate") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

d = ggplot(subset(results, var == "Dxs" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "Dxs"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "Dxs"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,100)) +
  labs(title = "",
       x = "Year",
       y = "Notification rate") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

e = ggplot(subset(results, var == "URs" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "URs"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "URs"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,2)) +
  labs(title = "",
       x = "Year",
       y = "Urban/Rural (Subclinical)") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

f = ggplot(subset(results, var == "URc" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "URc"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "URc"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,2)) +
  labs(title = "",
       x = "Year",
       y = "Urban/Rural (Clinical)") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

g = ggplot(subset(results, var == "HLs" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "HLs"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "HLs"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,2)) +
  labs(title = "",
       x = "Year",
       y = "High/Low (Subclinical)") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

h = ggplot(subset(results, var == "HLc" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = val), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  #geom_point(data = subset(targetsdb, var == "HLc"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(targetsdb, var == "HLc"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,2)) +
  labs(title = "",
       x = "Year",
       y = "High/Low (Clinical)") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

plotx = (a|b|c|d)/(e|f|g|h)
plotx
aspect_ratio = 2

# 3.2 State plots
states <- oderes %>% 
  filter(time %in% c(2008,2018)) %>% 
  pivot_longer(cols = -time, names_to = "var", values_to = "val") %>% 
  mutate(var = gsub("_[A-Z]+", "", var)) %>% 
  filter(var %in% c("N", "I", "O", "M", "S", "C", "RC", "P")) %>% 
  mutate(var = case_when(var == "N" ~ "SUS", var == "I" ~ "INF", var == "O" ~ "CLE",
                         var == "M" ~ "MIN", var == "S" ~ "SUB", var == "C" ~ "CLN",
                         var == "RC" ~ "RxCLN", var == "P" ~ "REC")) %>% 
  group_by(var, time) %>% 
  summarise(val = sum(val)) %>% 
  mutate(pval = (val/100000)) %>% 
  mutate(var = factor(var, levels = c("SUS", "INF", "CLE", "MIN", "SUB", "CLN", "RxCLN", "REC")))

x1 = ggplot(filter(states, time == 2008 & var %in% c("SUS", "INF", "CLE"))) +
  geom_col(aes(x = var, y = pval)) +
  labs(title = '2008', x = 'State', y = 'Percentage') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme(legend.position = "bottom")

x2 = ggplot(filter(states, time == 2018 & var %in% c("SUS", "INF", "CLE"))) +
  geom_col(aes(x = var, y = pval)) +
  labs(title = '2018', x = 'State', y = 'Percentage') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme(legend.position = "bottom")

y1 = ggplot(filter(states, time == 2008 & var %in% c("MIN", "SUB", "CLN", "RxCLN", "REC"))) +
  geom_col(aes(x = var, y = pval)) +
  labs(title = '2008', x = 'State', y = 'Percentage') +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed", size = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme(legend.position = "bottom")

y2 = ggplot(filter(states, time == 2018 & var %in% c("MIN", "SUB", "CLN", "RxCLN", "REC"))) +
  geom_col(aes(x = var, y = pval)) +
  labs(title = '2018', x = 'State', y = 'Percentage') +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed", size = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme(legend.position = "bottom")

ploty = (x1 | x2) / (y1 | y2)

# 3.3 ARI plots
ari <- oderes %>% 
  filter(time >= 2000) %>% 
  select(time, starts_with("ARI")) %>% 
  mutate(ARI = ARIsi + ARIoi + ARIpi) %>% 
  mutate(ARI = ARI/100000, ARIsi = ARIsi/100000, ARIoi = ARIoi/100000, ARIpi = ARIpi/100000) %>% 
  pivot_longer(cols = -time, names_to = "var", values_to = "val") %>% 
  mutate(var = case_when(var == "ARI" ~ "Total ARI", var == "ARIoi" ~ "CLE -> INF", 
                         var == "ARIpi" ~ "REC -> INF", var == "ARIsi" ~ "SUS -> INF")) %>% 
  mutate(var = factor(var, levels = c("Total ARI", "SUS -> INF", "CLE -> INF", "REC -> INF")))

plotz = ggplot(ari) + 
  geom_line(aes(x = time, y = val, colour = var)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = 'Time (year)', y = 'ARI')

arit <- oderes %>% 
  filter(time >= 2000) %>% 
  select(time, ARI) %>% 
  pivot_longer(cols = -time, names_to = "var", values_to = "val")

plotz = ggplot(arit) + 
  geom_line(aes(x = time, y = val)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0,max(arit$val))) +
  labs(x = 'Time (year)', y = 'ARI')

plotyz <- plot_grid(ploty, plotz, ncol = 2, rel_widths = c(3, 2))

plotxyz <- plot_grid(plotx, plotyz, nrow = 2, rel_heights = c(1,1))

allplot <- plot_grid(plotxyz, table, ncol = 2, rel_widths = c(3, 1))

save_plot <- function(plot, filename, path) {
  timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  filename <- paste0(filename, "_", timestamp, ".tiff")
  tiff(file.path(path, filename), width = 18, height = 10, units = 'in', res = 100)
  print(plot)
  dev.off()
}

save_plot(allplot, "plot", here("scripts", "manual fits"))
