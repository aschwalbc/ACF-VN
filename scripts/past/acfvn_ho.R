

library(hmer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
library(data.table)
library(patchwork)

ode_results = function(parms, end_time = 2020) {
  
  P = 100000
  mu = 1/70

  forcer_i = matrix(c(1500, parms['i1'], 1999, parms['i1'], 2020, parms['i2']), ncol = 2, byrow = TRUE)
  force_func_i = approxfun(x = forcer_i[,1], y = forcer_i[,2], method = "linear", rule = 2)
  
  forcer = matrix(c(1500, parms['j'], 1999, parms['j'], 2020, parms['k']), ncol = 2, byrow = TRUE)
  force_func = approxfun(x = forcer[,1], y = forcer[,2], method = "linear", rule = 2)
  
  forcer_p_R = matrix(c(1500, parms['l'], 1999, parms['l'], 2020, parms['m']), ncol = 2, byrow = TRUE)
  force_func_p_R = approxfun(x = forcer_p_R[,1], y = forcer_p_R[,2], method = "linear", rule = 2)
  force_func_p_U = approxfun(x = forcer_p_R[,1], y = 1-forcer_p_R[,2], method = "linear", rule = 2)
  
  forcer_p_L = matrix(c(1500, parms['n'], 1999, parms['n'], 2020, parms['o']), ncol = 2, byrow = TRUE)
  force_func_p_L = approxfun(x = forcer_p_L[,1], y = forcer_p_L[,2], method = "linear", rule = 2)
  force_func_p_H = approxfun(x = forcer_p_L[,1], y = 1-forcer_p_L[,2], method = "linear", rule = 2)
  
  des = function(time, state, parms) {
    
    with(as.list(c(state, parms)), {
      
      dN_RL = force_func_p_R(time)*force_func_p_L(time)*(mu*P + force_func_i(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RL/P) - mu*N_RL
      dI_RL = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RL/P) - a*I_RL - b*I_RL - c*I_RL - mu*I_RL
      dO_RL = a*I_RL + d*M_RL - mu*O_RL    
      dM_RL = b*I_RL + f*S_RL - d*M_RL - e*M_RL - mu*M_RL
      dS_RL = c*I_RL + e*M_RL + h*C_RL - f*S_RL - g*S_RL - mu*S_RL
      dC_RL = g*S_RL - h*C_RL - force_func_i(time)*C_RL - mu*C_RL - force_func(time)*C_RL
      dD_RL = force_func(time)*C_RL - mu*D_RL
      
      dN_RH = force_func_p_R(time)*force_func_p_H(time)*(mu*P + force_func_i(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RH/P) - mu*N_RH
      dI_RH = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_RH/P) - a*I_RH - b*I_RH - c*I_RH - mu*I_RH
      dO_RH = a*I_RH + d*M_RH - mu*O_RH    
      dM_RH = b*I_RH + f*S_RH - d*M_RH - e*M_RH - mu*M_RH
      dS_RH = c*I_RH + e*M_RH + h*C_RH - f*S_RH - g*S_RH - mu*S_RH
      dC_RH = g*S_RH - h*C_RH - force_func_i(time)*C_RH - mu*C_RH - force_func(time)*C_RH
      dD_RH = force_func(time)*C_RH - mu*D_RH
      
      dN_UL = force_func_p_U(time)*force_func_p_L(time)*(mu*P + force_func_i(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UL/P) - mu*N_UL
      dI_UL = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UL/P) - a*I_UL - b*I_UL - c*I_UL - mu*I_UL
      dO_UL = a*I_UL + d*M_UL - mu*O_UL    
      dM_UL = b*I_UL + f*S_UL - d*M_UL - e*M_UL - mu*M_UL
      dS_UL = c*I_UL + e*M_UL + h*C_UL - f*S_UL - g*S_UL - mu*S_UL
      dC_UL = g*S_UL - h*C_UL - force_func_i(time)*C_UL - mu*C_UL - force_func(time)*C_UL
      dD_UL = force_func(time)*C_UL - mu*D_UL
      
      dN_UH = force_func_p_U(time)*force_func_p_H(time)*(mu*P + force_func_i(time)*(C_RL+C_RH+C_UL+C_UH)) - beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UH/P) - mu*N_UH
      dI_UH = beta*(alpha*(S_RL+S_RH+S_UL+S_UH)+(C_RL+C_RH+C_UL+C_UH))*(N_UH/P) - a*I_UH - b*I_UH - c*I_UH - mu*I_UH
      dO_UH = a*I_UH + d*M_UH - mu*O_UH    
      dM_UH = b*I_UH + f*S_UH - d*M_UH - e*M_UH - mu*M_UH
      dS_UH = c*I_UH + e*M_UH + h*C_UH - f*S_UH - g*S_UH - mu*S_UH
      dC_UH = g*S_UH - h*C_UH - force_func_i(time)*C_UH - mu*C_UH - force_func(time)*C_UH
      dD_UH = force_func(time)*C_UH - mu*D_UH
      
      return(list(c(dN_RL, dN_RH, dN_UL, dN_UH,
                    dI_RL, dI_RH, dI_UL, dI_UH,
                    dO_RL, dO_RH, dO_UL, dO_UH,
                    dM_RL, dM_RH, dM_UL, dM_UH,
                    dS_RL, dS_RH, dS_UL, dS_UH,
                    dC_RL, dC_RH, dC_UL, dC_UH,
                    dD_RL, dD_RH, dD_UL, dD_UH),
                    W = (S_UL + S_UH)/(S_RL + S_RH), #rS_UR
                    X = (C_UL + C_UH)/(C_RL + C_RH), #rC_UR
                    Y = (S_RH + S_UH)/(S_RL + S_UL), #rS_HL
                    Z = (C_RH + C_UH)/(C_RL + C_UL), #rC_HL
                    S = S_RL + S_RH + S_UL + S_UH,
                    C = C_RL + C_RH + C_UL + C_UH,
                    H = force_func_i(time)*(C_RL + C_RH + C_UL + C_UH), 
                    T = force_func(time)*(C_RL + C_RH + C_UL + C_UH)))
    })
  }
  
  yini = c(N_RL = 100000-1000, N_RH = 0, N_UL = 0, N_UH = 0,
           I_RL = 0, I_RH = 0, I_UL = 0, I_UH = 0,
           O_RL = 0, O_RH = 0, O_UL = 0, O_UH = 0,
           M_RL = 0, M_RH = 0, M_UL = 0, M_UH = 0,
           S_RL = 0, S_RH = 0, S_UL = 0, S_UH = 0,
           C_RL = 1000, C_RH = 0, C_UL = 0, C_UH = 0,
           D_RL = 0, D_RH = 0, D_UL = 0, D_UH = 0)
  
  times = seq(1500, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms)
  return(out)
}


get_results = function(params, times, outputs) {
  t_max = max(times)
  all_res = ode_results(params, t_max)
  actual_res = all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped = reshape2::melt(actual_res[,outputs])
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = "")))
}


ranges = list(
  beta = c(0,50),     #contact parameter
  alpha = c(0.5,1.0), #relative infectiousness 
  a = c(0.64,1.84),   #infout
  b = c(0.04,0.17),   #infmin
  c = c(0.01,0.06),   #infsub
  d = c(0.14,0.23),   #minout
  e = c(0.21,0.28),   #minsub
  f = c(1.22,2.01),   #sub_min
  g = c(0.56,0.95),   #sub_clin
  h = c(0.46,0.73),   #clin_sub
  i1 = c(0.28,0.38),  #mort
  i2 = c(0,1),        #mort
  j = c(0,1),         #clin_diag
  k = c(0,1),         #clin_diag
  l = c(0,1),         #p_R initial
  m = c(0,1),         #p_R final
  n = c(0,1),         #P_l initial
  o = c(0,1))         #P_l final 

unc = 0.1

targets = list(
  S2007 = list(val = 144, sigma = unc*144),
  S2017 = list(val = 128, sigma = unc*128),
  C2007 = list(val = 87, sigma = unc*87),
  C2017 = list(val = 73, sigma = unc*73),
  H2000 = list(val = 32.9, sigma = unc*32.9),
  H2010 = list(val = 22.6, sigma = unc*22.6),
  H2020 = list(val = 10.3, sigma = unc*10.3),
  T2010 = list(val = 53.6, sigma = unc*53.6),
  T2020 = list(val = 55.8, sigma = unc*55.8),
  W2007 = list(val = 0.39, sigma = unc*0.39), #rS_UR2007
  W2017 = list(val = 0.84, sigma = unc*0.84), #rS_UR2017
  X2007 = list(val = 0.34, sigma = unc*0.34), #rC_UR2007
  X2017 = list(val = 0.67, sigma = unc*0.67), #rC_UR2017
  Y2007 = list(val = 0.53, sigma = unc*0.53), #rS_HL2007
  Y2017 = list(val = 0.83, sigma = unc*0.83), #rS_HL2017
  Z2007 = list(val = 0.50, sigma = unc*0.50), #rC_HL2007
  Z2017 = list(val = 0.91, sigma = unc*0.91)  #rC_HL2017
  )



data = as.data.frame(t(as.data.frame(targets)))
data$var = substr(rownames(data),1,1)
data$time = as.numeric(substr(rownames(data),2,5))
data = as.data.frame(cbind(data[grep("val", rownames(data)),c("time","var","V1")],V2 = data[grep("sigma", rownames(data)),c("V1")]))
colnames(data)[c(3,4)] = c("val","sigma")
data$lo = data$val-2*data$sigma
data$hi = data$val+2*data$sigma



number_params = length(ranges)
number_points = 10*number_params
#number_points = 20*number_params



initial_LHS = lhs::randomLHS(number_points, number_params)
initial_points = setNames(data.frame(t(apply(initial_LHS, 1, function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

initial_results <- data.frame(t(apply(initial_points, 1, get_results, c(2000, 2007, 2010, 2017, 2020), c('S', 'C', 'H', 'T', 'W', 'X', 'Y', 'Z'))))[,names(targets)]

wave0 = cbind(initial_points, initial_results)

t_sample <- sample(1:nrow(wave0), number_points)
training <- wave0[t_sample,]
validation <- wave0[-t_sample,]

ems_wave1 <- emulator_from_data(training, names(targets), ranges)



new_points <- generate_new_runs(ems_wave1, number_points, targets, verbose = TRUE)

min_val <- list()
max_val <- list()
new_ranges <- list()
for (i in 1:length(ranges)) {
  par <- names(ranges)[[i]]
  min_val[[par]] <- max(min(new_points[,par])-0.05*diff(range(new_points[,par])), 
                        ranges[[par]][1])
  max_val[[par]] <- min(max(new_points[,par])+0.05*diff(range(new_points[,par])),
                        ranges[[par]][2])
  new_ranges[[par]] <- c(min_val[[par]], max_val[[par]])
}


new_initial_results <- data.frame(t(apply(new_points, 1, get_results, c(2000, 2007, 2010, 2017, 2020), c('S', 'C', 'H', 'T', 'W', 'X', 'Y', 'Z'))))[,names(targets)]



wave1 <- cbind(new_points, new_initial_results)

new_t_sample <- sample(1:nrow(wave1), number_points)
new_training <- wave1[new_t_sample,]
new_validation <- wave1[-new_t_sample,]

ems_wave2 <- emulator_from_data(new_training, names(targets), new_ranges)



new_new_points <- generate_new_runs(c(ems_wave2, ems_wave1), number_points, targets, verbose=TRUE)

min_val <- list()
max_val <- list()
new_new_ranges <- list()
for (i in 1:length(ranges)) {
  par <- names(ranges)[[i]]
  min_val[[par]] <- max(min(new_new_points[,par])-0.05*diff(range(new_new_points[,par])),
                        ranges[[par]][1])
  max_val[[par]] <- min(max(new_new_points[,par])+0.05*diff(range(new_new_points[,par])),
                        ranges[[par]][2])
  new_new_ranges[[par]] <- c(min_val[[par]], max_val[[par]])
}


new_new_initial_results <- data.frame(t(apply(new_new_points, 1, get_results, c(2000, 2007, 2010, 2017, 2020), c('S', 'C', 'H', 'T', 'W', 'X', 'Y', 'Z'))))[,names(targets)]

wave2 <- cbind(new_new_points, new_new_initial_results)

new_new_t_sample <- sample(1:nrow(wave1), number_points)
new_new_training <- wave1[new_new_t_sample,]
new_new_validation <- wave1[-new_new_t_sample,]

ems_wave3 <- emulator_from_data(new_new_training, names(targets), new_new_ranges)



new_new_new_points <- generate_new_runs(c(ems_wave3, ems_wave2, ems_wave1), number_points, targets, verbose=TRUE)



#final_points = new_points
#final_points = new_new_points
final_points = new_new_new_points

quants <- c(0.025,0.5,0.975)

parameters = apply(final_points, 2, quantile, probs = quants, na.rm = TRUE)
parameters_transposed = transpose(as.data.frame(parameters))
colnames(parameters_transposed) = rownames(parameters)
rownames(parameters_transposed) = colnames(parameters)
parameters = parameters_transposed


results = as.data.frame(apply(final_points, 1, ode_results))[-seq(1,521),]
results = as.data.frame(t(apply(results, 1, quantile, probs = quants, na.rm = TRUE)))
x = colnames(ode_results(parms))[-1]
results$var = c(t(replicate(521,x)))
results$time = rep(seq(1500,2020),36)



parameters$parameter = c("beta","alpha","inf_out","inf_min","inf_sub","min_out","min_sub","sub_min","sub_clin","clin_sub","mort1","mort2","clin_diag1","clin_diag2", "p_R1", "p_R2", "p_L1", "p_L2")
parameters[,c(1,2,3)] = round(parameters[,c(1,2,3)],2)
table = data.frame(parameter = parameters$parameter, low  = parameters$`2.5%`, med = parameters$`50%`, hig = parameters$`97.5%`)



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



breaks = seq(1999, 2020, 1)
limits = c(1999, 2020.5)
from_to = seq(1999,2020)



a = ggplot(subset(results, var == "S" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "S"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "S"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) + 
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 250, 25), limits = c(0, 250), expand = c(0, 0)) +
  labs(title = "",
       x = "Year",
       y = "Subclinical") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(), #element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(), #element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

b = ggplot(subset(results, var == "C" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "C"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "C"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 250, 25), limits = c(0, 250), expand = c(0, 0)) +
  labs(title = "",
       x = "Year",
       y = "Clinical") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(), #element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(), #element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

c = ggplot(subset(results, var == "H" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "H"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "H"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 60, 5), limits = c(0, 60), expand = c(0, 0)) +
  labs(title = "",
       x = "Year",
       y = "Mortality") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(), #element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(), #element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

d = ggplot(subset(results, var == "T" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "T"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "T"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "",
       x = "Year",
       y = "Notifications") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_blank(), #element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_blank(), #element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")




e = ggplot(subset(results, var == "W" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "W"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "W"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2), expand = c(0, 0)) +
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

f = ggplot(subset(results, var == "X" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "X"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "X"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2), expand = c(0, 0)) +
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

g = ggplot(subset(results, var == "Y" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "Y"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "Y"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2), expand = c(0, 0)) +
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

h = ggplot(subset(results, var == "Z" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = subset(data, var == "Z"), aes(y = val), size = point_size, colour = "black") +
  geom_errorbar(data = subset(data, var == "Z"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = breaks, limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.25), limits = c(0, 2), expand = c(0, 0)) +
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

plot = (a|b|c|d)/(e|f|g|h)

aspect_ratio = 2

ggsave("plot.png" , width = aspect_ratio*height, height = height, units = "cm", dpi = 300)


parms = c(
  beta = 100,
  alpha = 0.5,
  a = 1.13,
  b = 0.08,
  c = 0.03,
  d = 0.18,
  e = 0.24,
  f = 1.56,
  g = 0.72,
  h = 0.58,
  i1 = 0.33,
  i2 = 0.66,
  #i2 = 0.33,
  j = 0.20,
  k = 0.4,
  l = 0.7,
  m = 0.5,
  n = 0.8,
  o = 0.2)





