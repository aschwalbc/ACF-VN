## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 03: Plots.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(ggplot2) # To build comparative plots
library(tidyverse) # To use tidyverse
library(reshape2) # Reshaping data easily
library(data.table) # Faster than data.frame, allows use of j operator (:=)
library(patchwork) # Plot composition
library(extrafont) # Add specific fonts

# 1. Load data ==========
results <- as.data.table(import(here("outputs","results.Rdata")))
targets <- as.data.table(import(here("data","fit","targets.Rdata")))

# 2. Plot ==========
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

# 3. Fit plots ==========
a = ggplot(subset(results, var == "TBc" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_errorbar(data = subset(targets, var == "TBc"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) + 
  scale_x_continuous(breaks = seq(2000,2020,2), limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 350, 50), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,350)) +
  labs(title = "",
       x = "Year",
       y = "TB prevalence rate per 100k") +
  theme(text = element_text(family = "Open sans"),
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

b = ggplot(subset(results, var == "Spr" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_errorbar(data = subset(targets, var == "Spr"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = seq(2000,2020,2), limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "",
       x = "Year",
       y = "Proportion asymptomatic TB") +
  theme(text = element_text(family = "Open sans"),
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

c = ggplot(subset(results, var == "Mor" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_errorbar(data = subset(targets, var == "Mor"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = seq(2000,2020,2), limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 90, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,90)) +
  labs(title = "",
       x = "Year",
       y = "TB mortality rate per 100k") +
  theme(text = element_text(family = "Open sans"),
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

d = ggplot(subset(results, var == "Dxs" & time %in% from_to), aes(x=time)) +
  geom_line(aes(y = `50%`), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_errorbar(data = subset(targets, var == "Dxs"), aes(ymin = lo, ymax = hi), size = line_thickness, width = error_bar_width) +
  scale_x_continuous(breaks = seq(2000,2020,2), limits = limits, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,100)) +
  labs(title = "",
       x = "Year",
       y = "Notification rate per 100k") +
  theme(text = element_text(family = "Open sans"),
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

plot = (a|b|c|d)
aspect_ratio = 2

png(here("plots", "00_fit.png"),  width = 14, height = 5, units = 'in', res = 1000)
print(plot)
dev.off()
