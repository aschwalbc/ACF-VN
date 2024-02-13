## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 07: Plots.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(ggplot2) # To build plots

# 1. Load data ==========
outs <- import(here("outputs", "outs", "outs.Rdata"))

# 2. Plots ==========
prev_targets <- c(100, 50, 20)
dis_state <- factor(c("Clinical","Subclinical","Minimal"))
inf_dis <- c("Clinical", "Subclinical")
scenarios <- c("Xpert", "CXR+Xpert", "CXR", "BAU")
treatments <- c("BAU: Clinical", "ACF: Clinical", "ACF: Minimal", "ACF: Subclinical")
types <- c("ACF", "BAU")
type_labels <- c(acfa = "Xpert", acfb = "CXR+Xpert", acfc = "CXR")

# TB prevalence rate per scenarios 
png(here("plots", "00_TBprev01_100.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_hline(yintercept = prev_targets[1], linetype = "dashed", color = "gray") +
  geom_line(data = filter(outs, var == 'rTBc' & goal %in% c('100', 'none')), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == 'rTBc' & goal %in% c('100', 'none')), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,250,25), limits = c(0,250), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)", title = 'TB prevalence threshold: 100/100,000') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

png(here("plots", "00_TBprev02_50.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_hline(yintercept = prev_targets[2], linetype = "dashed", color = "gray") +
  geom_line(data = filter(outs, var == 'rTBc' & goal %in% c('50', 'none')), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == 'rTBc' & goal %in% c('50', 'none')), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,250,25), limits = c(0,250), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)", title = 'TB prevalence threshold: 50/100,000') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

png(here("plots", "00_TBprev03_20.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_hline(yintercept = prev_targets[3], linetype = "dashed", color = "gray") +
  geom_line(data = filter(outs, var == 'rTBc' & goal %in% c('20', 'none')), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == 'rTBc' & goal %in% c('20', 'none')), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,250,25), limits = c(0,250), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)", title = 'TB prevalence threshold: 20/100,000') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Incident TB averted
png(here("plots", "01_IncTBaverted.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var == 'dfcumInc' & !goal == 'none' & time == 2050), aes(y = abs(val), x = goal, fill = goal), stat = "identity") +
  scale_fill_manual(values = c("#9F86C0","#5E548E","#231942")) +
  geom_errorbar(data = filter(outs, var == 'dfcumInc' & goal != 'none' & time == 2050), aes(ymin = abs(lo), ymax = abs(hi), x = goal), colour = 'darkgrey', width = 0.2, stat = "identity") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = 'M')) + 
  labs(y = "Incident TB averted by 2050", fill = 'TB prevalence threshold (per 100k)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

# TB mortality averted
png(here("plots", "02_MortTBaverted.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var == 'dfcumMor' & !goal == 'none' & time == 2050), aes(y = abs(val), x = goal, fill = goal), stat = "identity") +
  scale_fill_manual(values = c("#9F86C0","#5E548E","#231942")) +
  geom_errorbar(data = filter(outs, var == 'dfcumMor' & goal != 'none' & time == 2050), aes(ymin = abs(lo), ymax = abs(hi), x = goal), colour = 'darkgrey', width = 0.2, stat = "identity") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = 'K')) + 
  labs(y = "TB mortality averted by 2050", fill = 'TB prevalence threshold (per 100k)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

# Total cost of ACF screening
png(here("plots", "03_ACFscreening_cost.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var == 'cumcACF' & !goal == 'none' & time == 2050), aes(y = abs(val), x = goal, fill = goal), stat = "identity") +
  scale_fill_manual(values = c("#9F86C0","#5E548E","#231942")) +
  geom_errorbar(data = filter(outs, var == 'cumcACF' & goal != 'none' & time == 2050), aes(ymin = abs(lo), ymax = abs(hi), x = goal), colour = 'darkgrey', width = 0.2, stat = "identity") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B')) + 
  labs(y = "Total costs of ACF screening (USD)", fill = 'TB prevalence threshold (per 100k)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

# Total cost of ACF treatment
png(here("plots", "04_ACFtreatment_cost.png"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_grid(~type, labeller = labeller(type = type_labels)) +
  geom_bar(data = filter(outs, var == 'cumcRxACF' & !goal == 'none' & time == 2050), aes(y = abs(val), x = goal, fill = goal), stat = "identity") +
  scale_fill_manual(values = c("#9F86C0","#5E548E","#231942")) +
  geom_errorbar(data = filter(outs, var == 'cumcRxACF' & goal != 'none' & time == 2050), aes(ymin = abs(lo), ymax = abs(hi), x = goal), colour = 'darkgrey', width = 0.2, stat = "identity") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = 'B')) + 
  labs(y = "Total costs of ACF treatment (USD)", fill = 'TB prevalence threshold (per 100k)') +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(10,15,5,10,"pt"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()
