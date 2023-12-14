## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 06: Plots.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(ggplot2) # To build plots

# 1. Load data ==========


# 6. Plots ==========
prev_targets <- c(100, 50, 20)
dis_state <- factor(c("Clinical","Subclinical","Minimal"))
inf_dis <- c("Clinical", "Subclinical")
scenarios <- c("01: Xpert", "02: CXR->Xpert", "03: CXR", "Baseline")
treatments <- c("BAU: Clinical", "ACF: Clinical", "ACF: Minimal", "ACF: Subclinical")
types <- c("ACF", "BAU")
type_label <- labeller(type = c("acfa" = "01: Xpert", "acfb" = "02: CXR->Xpert", "acfc" = "03: CXR", "base" = "Baseline"))
dem_urbrur <- c("Rural", "Urban")
dem_ses <- c("High","Low")

# TB prevalence per scenarios 
#tiff(here("plots", "TBprev.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "TBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "TBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,250,25), limits = c(0,250), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB prevalence rate (per 100K)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_hline(yintercept = prev_targets, linetype = "dashed", color = "gray") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
#dev.off()

# Proportion reduction of TB prevalence
tiff(here("plots", "TBreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pTBc"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pTBc"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1.01,0.1), limits = c(0,1.01), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(min(acf_year),2050,5), limits = c(min(acf_year),2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(min(acf_year),2050)) + 
  labs(x = "Year", y = "Proportion reduction of TB prevalence") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  geom_rect(aes(xmin = 2029 - 0.5, xmax = 2029 + 0.5, ymin = 0.210, ymax = 0.498), colour = "#1D2D5F", linetype = 'dashed', alpha = 0.2) + 
  geom_text(aes(x = 2028.5, y = 0.25, label = "ACT3 reduction", fontface = 'bold'), angle = 90, vjust = -0.5, hjust = 0, size = 3, color = "#1D2D5F") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Annual risk of infection 
tiff(here("plots", "ARI.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "ARI"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "ARI"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Annual risk of infection") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB mortality per scenarios 
tiff(here("plots", "TBmort.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "rMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "rMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB mortality rate (per 100K)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Total TB mortality 
tiff(here("plots", "totTBmort.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "tMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "tMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,15000,2500), expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Total TB mortality") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Proportion reduction of TB mortality
tiff(here("plots", "TBmorreduct.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "pMor"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pMor"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1.01,0.1), limits = c(0,1.01), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(min(acf_year),2050,5), limits = c(min(acf_year),2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(min(acf_year),2050)) + 
  labs(x = "Year", y = "Proportion reduction of TB mortality") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative TB deaths averted
tiff(here("plots", "cumTBdeathsavert.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dcumMor"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dcumMor"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative TB deaths averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB deaths averted
tiff(here("plots", "TBdeathsavert.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dMor"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dMor"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "TB deaths averted per year (compared to BAU)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# TB treatments
tiff(here("plots", "TBtreatment.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, scales = 'free_y', labeller = type_label) +
  geom_line(data = filter(outs, var %in% c("SCmin","SCsub","SCcln","DXcln")), aes(x = time, y = val, colour = var)) +
  scale_colour_manual(values = c("#000000","#F58B65","#4DC4CB","#FDC75D"), name = NULL, labels = treatments) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments (BAU vs ACF)
tiff(here("plots", "TBtreatBAUACF.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, scales = 'free_y', labeller = type_label) +
  geom_line(data = filter(outs, var %in% c("DXcln", "ACF")), aes(x = time, y = val, colour = var)) +
  scale_colour_manual(values = c("#CE2931","#000000"), name = NULL, labels = types) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# All TB treatments
tiff(here("plots", "TBtreatALLTB.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "AllTB"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "AllTB"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K"), limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050), ylim = c(0, 9e5)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# All treatments (TB+FP)
tiff(here("plots", "TBtreatALLTx.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "AllTx"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "AllTx"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K"), limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050), ylim = c(0, 9e5)) + 
  labs(x = "Year", y = "Number of treatments (including FP)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments averted
tiff(here("plots", "TBtreatmentaverted.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, labeller = type_label) +
  geom_line(data = filter(outs, var == "dAllTB" & type != 'base'), aes(x = time, y = val)) +
  geom_area(data = filter(outs, var == "dAllTB" & type != 'base'), aes(x = time, y = val, fill = fill), alpha = 0.5) +
  scale_fill_manual(values = c("over" = "#FF531A", "under" = "#1AC6FF"), name = NULL) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "none", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# TB treatments averted (considering FP)
tiff(here("plots", "TBtreatmentaverted_withFP.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  facet_wrap(~type, labeller = type_label) +
  geom_line(data = filter(outs, var == "dAllTx" & type != 'base'), aes(x = time, y = val)) +
  geom_area(data = filter(outs, var == "dAllTx" & type != 'base'), aes(x = time, y = val, fill = fill), alpha = 0.5) +
  scale_fill_manual(values = c("over" = "#FF531A", "under" = "#1AC6FF"), name = NULL) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments") +
  theme_bw() +
  theme(legend.position = "none", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Cumulative TB treatments averted
tiff(here("plots", "cumTBtreatmentaverted_line.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "dBAU"), aes(x = time, y = abs(val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "dBAU"), aes(x = time, ymin = abs(lo), ymax = abs(hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number of TB treatments averted") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Proportion prevalence (urban vs rural)
tiff(here("plots", "Prop_urbvrur.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "URt"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "URt"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB urban vs rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion prevalence (low vs high SES)
tiff(here("plots", "Prop_lovhi.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "HLt"), aes(x = time, y = (1-val), colour = type)) +
  geom_ribbon(data = filter(outs, var == "HLt"), aes(x = time, ymin = (1-lo), ymax = (1-hi), fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion TB low vs high") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Urban
tiff(here("plots", "Prop_UH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPUH"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPUH"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: High SES - Rural
tiff(here("plots", "Prop_RH.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPRH"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPRH"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion High SES - Rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Urban
tiff(here("plots", "Prop_UL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPUL"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPUL"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion Low SES - Urban") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion population: Low SES - Rural
tiff(here("plots", "Prop_RL.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pPRL"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pPRL"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion Low SES - Rural") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Proportion subclinical
tiff(here("plots", "Prop_scTB.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  facet_wrap(~type, labeller = type_label) +
  geom_area(data = filter(outs, var %in% c("Cpr","Spr")), aes(x = time, y = val, fill = var), position = "fill") +
  scale_fill_manual(values = c("#F58B65","#FDC75D"), name = NULL, labels = inf_dis) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion infectious TB") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# Disease states
tiff(here("plots", "Prop_disstates.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  facet_wrap(~type, labeller = type_label) +
  geom_area(data = filter(outs, var %in% c("ClnIf","SubIf","MinIf")), aes(x = time, y = val, fill = reorder(var, -val, decreasing = TRUE)), position = "fill") +
  scale_fill_manual(values = c("#F58B65","#FDC75D","#4DC4CB"), name = "State:", labels = dis_state) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion disease state") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"), panel.spacing.x = unit(10, "mm"))
dev.off()

# False positives
tiff(here("plots", "Falsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "False positive diagnoses") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative false positives
tiff(here("plots", "CumFalsepos.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumFPnds"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumFPnds"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative false positive diagnoses") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False negatives
tiff(here("plots", "Falseneg.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "FNdis"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "FNdis"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-3, suffix = "K")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "False negative diagnoses (Missed diagnoses)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative false negatives
tiff(here("plots", "CumFalseneg.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumFNdis"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumFNdis"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative false negative diagnoses (Missed diagnoses)") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# False positive treated
tiff(here("plots", "FPtreated.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot(data = outs) +
  geom_line(data = filter(outs, var == "pFP"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "pFP"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(2020,2050,10), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Proportion FP (% of total screened)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
dev.off()

# Number needed to screen
tiff(here("plots", "NNS.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "NNS"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "NNS"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Number needed to screen") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()

# Cumulative number screened
tiff(here("plots", "NumScreen.tiff"), width = 6, height = 5, units = 'in', res = 150)
ggplot() +
  geom_line(data = filter(outs, var == "cumNumSC"), aes(x = time, y = val, colour = type)) +
  geom_ribbon(data = filter(outs, var == "cumNumSC"), aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_fill_manual(values = c("#CE2931","#2984CE","#FFBC47","#1D2D5F"), name = NULL, labels = scenarios) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(breaks = seq(2020,2050,5), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2020,2050)) + 
  labs(x = "Year", y = "Cumulative number screened/treated") +
  geom_rect(data = data.frame(acf_year), aes(xmin = acf_year, xmax = acf_year + 1, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(10,15,5,10,"pt"))
dev.off()
