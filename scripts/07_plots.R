## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 07: Plots.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(ggplot2) # To build plots
library(RColorBrewer) # Colour palletes
library(extrafont) # Add specific fonts

# 1. Load data ==========
outs <- import(here("outputs", "outs", "outs.Rdata"))
icers <- import(here("outputs", "outs", "icers.Rdata"))

# 2. Plots ==========
prev_targets <- c(100, 50, 20)
dis_state <- factor(c("Symptomatic", "Asymptomatic", "Unconfirmed"))
inf_dis <- c("Symptomatic", "Asymptomatic")
scenarios <- c("NAAT", "CXR+NAAT", "CXR", "BAU")
type_labels <- c(acfa = "NAAT", acfb = "CXR+NAAT", acfc = "CXR", base = 'BAU',
                 acfax = "NAAT (1USD)", acfbx = "CXR+NAAT (1USD)")
goal_labels <- c('100' = '100 per 100k', '50' = '50 per 100k', '20' = '20 per 100k')
facet_labeller <- labeller(type = type_labels, goal = goal_labels)

# TB prevalence rate
acftbprev <- filter(outs, type %in% c('acfa', 'acfb', 'acfc') & var == 'rTBc') %>% mutate(db = 'acf')
bautbprev <- filter(outs, type == 'base' & var == 'rTBc') %>% mutate(db = 'bau') %>% within(rm(type, goal))

png(here("plots", "01_TBprev.png"), width = 10, height = 5, units = 'in', res = 1000) 
ggplot() +
  geom_hline(data = acftbprev, aes(yintercept = as.numeric(as.character(goal))*1.1, group = interaction(type, goal)), color = "darkgrey", linetype = "dashed") +
  geom_hline(data = acftbprev, aes(yintercept = as.numeric(as.character(goal))*0.9, group = interaction(type, goal)), color = "darkgrey", linetype = "dashed") +
  facet_grid(cols = vars(goal), labeller = facet_labeller) +
  geom_line(data = acftbprev, aes(x = time, y = val, colour = type)) +
  geom_line(data = bautbprev, aes(x = time, y = val, colour = "base")) +
  geom_ribbon(data = acftbprev, aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  geom_ribbon(data = bautbprev, aes(x = time, ymin = lo, ymax = hi, fill = "base"), alpha = 0.2) +
  scale_color_manual(values = c("acfa" = "#CE2931", "acfb" = "#2984CE", "acfc" = "#FFBC47", "base" = "#828282"),
                     labels = type_labels, name = 'Algorithms') +
  scale_fill_manual(values = c("acfa" = "#CE2931", "acfb" = "#2984CE", "acfc" = "#FFBC47", "base" = "#828282"), 
                    labels = type_labels, name = 'Algorithms') +
  scale_y_continuous(breaks = seq(0, 300, 100), limits = c(0, 300), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(2025, 2050, 5), limits = c(2023, 2050), expand = c(0, 0)) +
  labs(x = "Year", y = "TB prevalence (per 100K)") +
  theme_bw() +
  theme(legend.position = "bottom", plot.margin = margin(10, 15, 5, 10, "pt"),
        panel.spacing = unit(0.8, "lines"), text = element_text(family = "Open Sans"))
dev.off()

# Number treated
png(here("plots", "02_Ntreated.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var %in% c('cumTPinf', 'cumMIN', 'cumFP') & goal == '50' & time == 2050 & res == 'main'),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c('cumTPinf', 'cumMIN', 'cumFP'))), stat = "identity", position = 'stack') +
  # geom_bar(data = filter(outs, var %in% c('cumTPinf', 'cumMIN', 'cumFP') & goal == '50' & time == 2050 & type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumTPinf", "cumMIN", "cumFP"))), stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumFP" = "#D3D3D3", "cumMIN" = "#900C3F", "cumTPinf" = "#FF5733"),
                    labels = c("cumFP" = "Non-disease", "cumMIN" = "Non-infectious TB", "cumTPinf" = "Infectious TB")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = 'M'), breaks = seq(0,20e6,5e6)) +
  labs(y = "Number treated", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Number treated (disease-only)
png(here("plots", "02_Ntreated_dis.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var %in% c('cumTPinf', 'cumMIN') & goal == '50' & time == 2050 & res == 'main'),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c('cumTPinf', 'cumMIN'))), stat = "identity", position = 'stack') +
  # geom_bar(data = filter(outs, var %in% c('cumTPinf', 'cumMIN') & goal == '50' & time == 2050 & type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumTPinf", "cumMIN"))), stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumMIN" = "#900C3F", "cumTPinf" = "#FF5733"),
                    labels = c("cumMIN" = "Non-infectious TB", "cumTPinf" = "Infectious TB")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = 'M'), breaks = seq(0,1e6,0.5e6)) +
  labs(y = "Number treated", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Total cost of all
png(here("plots", "07_Totalcost.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, res == 'main', var == 'cumcAll' & !goal == 'none' & time == 2050), aes(y = abs(val), x = goal, fill = goal), stat = "identity") +
  scale_fill_manual(values = c("#9F86C0","#5E548E","#231942")) +
  geom_errorbar(data = filter(outs, res == 'main', var == 'cumcAll' & goal != 'none' & time == 2050), aes(ymin = abs(lo), ymax = abs(hi), x = goal), colour = 'darkgrey', width = 0.2, stat = "identity") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B'), breaks = seq(0,7e9,1e9)) +
  coord_cartesian(ylim = c(0, 6.5e9)) +
  labs(y = "Total costs (USD)", fill = 'TB prevalence threshold (per 100k)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Total cost of all
png(here("plots", "07_Totalcost_50t.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  # geom_bar(data = filter(outs, res == 'main', var %in% c('cumcRxACF', 'cumcACF') & goal == '50' & time == 2050),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcRxACF", "cumcACF"))), stat = "identity", position = 'stack') +
  geom_bar(data = filter(outs, var %in% c('cumcRxACF', 'cumcACF') & goal == '50' & time == 2050, type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcRxACF", "cumcACF"))), stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumcRxACF" = "#9F86C0", "cumcACF" = "#5E548E"),
                    labels = c("cumcRxACF" = "Treatment", "cumcACF" = "Screening")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B')) +
  labs(y = "Total costs of intervention (USD)", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Total cost of ACF screening stratified
png(here("plots", "03_ACFscreeningcost.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  # geom_bar(data = filter(outs, res == 'main', var %in% c('cumcACFFP','cumcACFMIN','cumcACFTPinf') & goal == '50' & time == 2050),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcACFTPinf", "cumcACFMIN", "cumcACFFP"))),
  #          stat = "identity", position = 'stack') +
  geom_bar(data = filter(outs, var %in% c('cumcACFFP','cumcACFMIN','cumcACFTPinf') & goal == '50' & time == 2050, type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcACFTPinf", "cumcACFMIN", "cumcACFFP"))),
           stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumcACFTPinf" = "#FF5733", "cumcACFMIN" = "#900C3F", "cumcACFFP" = "#D3D3D3"),
                    labels = c("cumcACFTPinf" = "Infectious TB", "cumcACFMIN" = "Non-infectious TB", "cumcACFFP" = "Non-disease")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B'), breaks = seq(0,6e9,1e9)) +
  coord_cartesian(ylim = c(0, 3e9)) +
  labs(y = "Total costs of ACF screening (USD)", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        text = element_text(family = "Open Sans"))
dev.off()

# Total cost of ACF screening stratified (disease-only)
png(here("plots", "03_ACFscreeningcost_dis.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  # geom_bar(data = filter(outs, res == 'main', var %in% c('cumcACFMIN','cumcACFTPinf') & goal == '50' & time == 2050),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c('cumcACFTPinf', 'cumcACFMIN'))), stat = "identity", position = 'stack') +
  geom_bar(data = filter(outs, var %in% c('cumcACFMIN','cumcACFTPinf') & goal == '50' & time == 2050, type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c('cumcACFTPinf', 'cumcACFMIN'))), stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumcACFTPinf" = "#FF5733", "cumcACFMIN" = "#900C3F"),
                    labels = c("cumcACFTPinf" = "Infectious TB", "cumcACFMIN" = "Non-infectious TB")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = 'M'), breaks = seq(0,25e6,5e6)) +
  coord_cartesian(ylim = c(0, 25e6)) +
  labs(y = "Total costs of ACF screening (USD)", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()
 
# Total cost of ACF treatment stratified
png(here("plots", "04_ACFtreatmentcost.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  # geom_bar(data = filter(outs, res == 'main', var %in% c('cumcRxACFFP','cumcRxMIN','cumcRxACFTPinf') & goal == '50' & time == 2050),
  #          aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcRxACFTPinf", "cumcRxMIN", "cumcRxACFFP"))),
  #          stat = "identity", position = 'stack') +
  geom_bar(data = filter(outs, var %in% c('cumcRxACFFP','cumcRxMIN','cumcRxACFTPinf') & goal == '50' & time == 2050, type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
           aes(y = abs(val), x = goal, fill = factor(var, levels = c("cumcRxACFTPinf", "cumcRxMIN", "cumcRxACFFP"))),
           stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumcRxACFFP" = "#D3D3D3", "cumcRxMIN" = "#900C3F", "cumcRxACFTPinf" = "#FF5733"),
                    labels = c("cumcRxACFFP" = "Non-disease", "cumcRxMIN" = "Non-infectious TB", "cumcRxACFTPinf" = "Infectious TB")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B'), breaks = seq(0,6e9,0.5e9)) +
  coord_cartesian(ylim = c(0, 1.5e9)) +
  labs(y = "Total costs of ACF treatment (USD)", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Total cost of ACF treatment stratified (disease-only)
png(here("plots", "04_ACFtreatmentcost_dis.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  # geom_bar(data = filter(outs, res == 'main', var %in% c('cumcRxMIN', 'cumcRxACFTPinf') & goal == '50' & time == 2050),
  #          aes(y = abs(val), x = goal, fill = var), stat = "identity", position = 'stack') +
  geom_bar(data = filter(outs, var %in% c('cumcRxMIN', 'cumcRxACFTPinf') & goal == '50' & time == 2050, type %in% c('acfa', 'acfb', 'acfc', 'acfd')),
           aes(y = abs(val), x = goal, fill = var), stat = "identity", position = 'stack') +
  scale_fill_manual(values = c("cumcRxACFTPinf" = "#FF5733", "cumcRxMIN" = "#900C3F"),
                    labels = c("cumcRxACFTPinf" = "Infectious TB", "cumcRxMIN" = "Non-infectious TB")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = 'M'), breaks = seq(0,150e6,50e6)) +
  coord_cartesian(ylim = c(0, 110e6)) +
  labs(y = "Total costs of ACF treatment (USD)", fill = 'State') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(family = "Open Sans"))
dev.off()

# Cost-effectiveness plane
cua <- data.frame(scenario = c('Baseline', 'MTB/RIF', 'Ultra', 'CXR+Ultra', 'CXR'),
                    dalys = c(0, 4281265, 4356480, 4414209, 5593902),
                    costs = c(0, 2952030989, 2271612723, 355788865, 1159573412))

cua_sa <- data.frame(scenario = 'CXR+Ultra (1USD)', dalys = 4396945, costs = 234876119)

png(here("plots", "10_WHOCEplane.png"), width = 10, height = 5, units = 'in', res = 1000)
ggplot(data = filter(icers, type %in% c('acfa', 'acfb', 'acfc', 'acfd'))) +
  # geom_point(aes(x = abs(dfcumDALYs), y = dfcumcAll, colour = type, fill = type), shape = 21, alpha = 0.6) +
  # scale_colour_discrete(labels = c('acfa' = 'Ultra', 'acfb' = 'CXR+Ultra', 'acfc' = 'CXR', 'acfd' = 'MTB/RIF'), name = NULL) +
  geom_segment(aes(x = -6e6, xend = 6e6, y = 1602 * -6e6, yend = 1602 * 6e6), colour = '#FFBF00', linetype = 'dashed', linewidth = 1) +
  geom_segment(aes(x = 0, xend = 4414209, y = 0, yend = 355788865), colour = '#1D2D5F') +
  geom_segment(aes(x = 0, xend = 4396945, y = 0, yend = 234876119), colour = '#CE2931', linetype = 'dashed') +
  geom_segment(aes(x = 4414209, xend = 5593902, y = 355788865, yend = 1159573412), colour = '#1D2D5F') +
  geom_point(cua, mapping = aes(x = dalys, y = costs), colour = '#1D2D5F') +
  geom_point(cua_sa, mapping = aes(x = dalys, y = costs), colour = '#CE2931') +
  geom_label(cua, mapping = aes(x = dalys, y = costs, label = scenario), 
             fill = '#1D2D5F', colour = '#FFFFFF', fontface = 'bold', nudge_y = 2e8) +
  geom_label(cua_sa, mapping = aes(x = dalys, y = costs, label = scenario), 
             fill = '#CE2931', colour = '#FFFFFF', fontface = 'bold', nudge_y = -2e8) +
  scale_x_continuous(breaks = seq(0, 6e6, 1e6), labels = scales::label_number(scale = 1e-6, suffix = 'M')) +
  scale_y_continuous(breaks = seq(0, 4.5e9, 0.5e9), labels = scales::label_number(scale = 1e-9, suffix = 'B')) +
  coord_cartesian(xlim = c(0, 6e6), ylim = c(0, 4.5e9)) +
  labs(x = 'DALYs averted', y = 'Incremental costs (USD)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        text = element_text(family = "Open Sans"))
dev.off()

png(here("plots", "15_Costsaving_yr.png"), width = 10, height = 5, units = 'in', res = 1000)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#828282") +
  geom_line(data = filter(outs, var %in% c('dfcBAU', 'dfcRxBAU', 'cACF', 'cRxACF'), goal == '50', 
                          type %in% c('acfb', 'acfbx')), aes(x = time, y = val, colour = var)) +
  geom_ribbon(data = filter(outs, var %in% c('dfcBAU', 'dfcRxBAU', 'cACF', 'cRxACF'), goal == '50', 
                            type %in% c('acfb', 'acfbx')), aes(x = time, ymin = lo, ymax = hi, fill = var), alpha = 0.1) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL,
                     labels = c("ACF screening", "ACF treatment", "BAU screening", "BAU treatment")) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL,
                    labels = c("ACF screening", "ACF treatment", "BAU screening", "BAU treatment")) +
  scale_y_continuous(breaks = seq(-450e6, 850e6, 20e6), labels = scales::label_number(scale = 1e-6, suffix = 'M')) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2024, 2050), ylim = c(-15e6,200e6)) +
  labs(x = "Year", y = "Incremental costs (USD)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", panel.spacing = unit(0.8, "lines"),
        plot.margin = margin(10,15,5,10,"pt"), text = element_text(family = "Open Sans"))
dev.off()
