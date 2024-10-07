# ACT3 TB Prevalence

library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
set.seed(123)

act3 <- import(here("outputs","act3","ACT3.csv"))
act3$group <- factor(act3$group, levels = c('int', 'ctrl'))

png(here("outputs","act3","ACT3.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot(act3, aes(x = year, y = rate, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, 
                position = position_dodge(width = 0.9)) +
  labs(x = 'Year', y = 'TB prevalence rate (per 100k)',
       title = 'ACT3', fill = 'Group') +
  scale_fill_manual(values = c("int" = "salmon", "ctrl" = "turquoise"),
                    labels = c("int" = "Intervention", "ctrl" = "Control")) +
  theme_bw() +
  theme(legend.position = 'bottom')
dev.off()

# TB prevalence reduction (Int Y1 vs Int Y4)
filter(act3, group == 'int' & year %in% c(1, 4))

rate1 <- 389.2; ci_lower1 <- 333.8 ; ci_upper1 <- 453.5
sd1 <- (ci_upper1 - ci_lower1) / (2 * 1.96)

rate2 <- 125.7; ci_lower2 <- 95.1; ci_upper2 <- 165.8
sd2 <- (ci_upper2 - ci_lower2) / (2 * 1.96)

# Simplified approach
1 - (rate2 / rate1)
1 - (ci_lower2 / ci_upper1)
1 - (ci_upper2 / ci_lower1)

# Sampling approach
smp1 <- rnorm(100000, mean = rate1, sd = sd1)
smp2 <- rnorm(100000, mean = rate2, sd = sd2)

prpdf <- 1 - (smp2 / smp1)

summary(prpdf)

inty1inty4 <- data.frame(type = 'IntY1-IntY4',time = 2028, val = abs(summary(prpdf)[[3]]),
                         lo = abs(summary(prpdf)[[5]]), hi = abs(summary(prpdf)[[2]]))

# TB prevalence decline (Int Y1 vs Ctrl Y4)
filter(act3, (group == 'int' & year == 1) | (group == 'ctrl' & year == 4))

rate1 <- 389.2; ci_lower1 <- 333.8 ; ci_upper1 <- 453.5
sd1 <- (ci_upper1 - ci_lower1) / (2 * 1.96)

rate2 <- 225.5; ci_lower2 <- 183.3; ci_upper2 <- 277.2
sd2 <- (ci_upper2 - ci_lower2) / (2 * 1.96)

# Simplified approach
1 - (rate2 / rate1)
1 - (ci_lower2 / ci_upper1)
1 - (ci_upper2 / ci_lower1)

# Sampling approach
smp1 <- rnorm(100000, mean = rate1, sd = sd1)
smp2 <- rnorm(100000, mean = rate2, sd = sd2)

prpdf <- 1 - (smp2 / smp1)

summary(prpdf)

inty1ctrly4 <- data.frame(type = 'IntY1-CtrlY4',time = 2028, val = abs(summary(prpdf)[[3]]),
                         lo = abs(summary(prpdf)[[5]]), hi = abs(summary(prpdf)[[2]]))

# Comparison graph
act <- import(here("outputs","act3","act3.Rdata"))
base <- import(here("outputs","results","r00_base.Rdata"))
comp <- rbind(inty1inty4, inty1ctrly4)

df <- act %>% 
  rbind(base) %>% 
  select(time, type, round, run, tTBc) %>% 
  group_by(run, time) %>%
  mutate(prTBc = tTBc / tTBc[type == 'base'],
         redTBc = 1 - (tTBc / tTBc[type == 'base'])) %>% 
  within(rm(tTBc)) %>% 
  pivot_longer(cols = -c(time, type, run, round), names_to = "var", values_to = "values") %>%
  group_by(time, type, round, var) %>%
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

filter(df, time == 2028 & type == 'acfa')

png(here("plots","S00_act3_int1vint4.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  geom_line(filter(df, var == 'prTBc'), mapping = aes(x = time, y = val, colour = type)) + 
  geom_ribbon(filter(df, var == 'prTBc'), mapping = aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931", "#1D2D5F"), labels = c('Xpert-only', 'BAU')) +
  scale_fill_manual(values = c("#CE2931", "#1D2D5F"), labels = c('Xpert-only', 'BAU')) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 100), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,1), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2025, 2035), ylim = c(0, 1)) + 
  labs(x = "Year", y = "TB prevalence proportional reduction") +
  geom_rect(aes(xmin = 2025, xmax = 2028, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) + 
  geom_point(aes(x = 2028, y = 1 - 0.6754), shape = 17) +
  geom_errorbar(aes(x = 2028, ymin = 1 - 0.7073, ymax = 1 - 0.6397), width = 0.5) + 
  theme_bw() +
  theme(legend.position = 'none', plot.margin = margin(10,15,5,10,"pt"), panel.grid.minor.x = element_blank())
dev.off()

png(here("plots","S00_act3_int1vctrl4.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  geom_line(filter(df, var == 'prTBc'), mapping = aes(x = time, y = val, colour = type)) + 
  geom_ribbon(filter(df, var == 'prTBc'), mapping = aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931", "#1D2D5F"), labels = c('Xpert-only', 'BAU')) +
  scale_fill_manual(values = c("#CE2931", "#1D2D5F"), labels = c('Xpert-only', 'BAU')) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 100), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,1), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2025, 2035), ylim = c(0, 1)) + 
  labs(x = "Year", y = "TB prevalence proportional reduction") +
  geom_rect(aes(xmin = 2025, xmax = 2028, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) + 
  geom_point(aes(x = 2028, y = 1 - 0.4166), shape = 17) +
  geom_errorbar(aes(x = 2028, ymin = 1 - 0.4672, ymax = 1 - 0.3615), width = 0.5) + 
  theme_bw() +
  theme(legend.position = 'none', plot.margin = margin(10,15,5,10,"pt"), panel.grid.minor.x = element_blank())
dev.off()

png(here("plots","S00_act3.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  geom_line(filter(df, var == 'prTBc'), mapping = aes(x = time, y = val, colour = type)) + 
  geom_ribbon(filter(df, var == 'prTBc'), mapping = aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931", "#1D2D5F","#636363","#182C25"), labels = c('Xpert-only', 'BAU','IntY1-CtrlY4','IntY1-IntY4')) +
  scale_fill_manual(values = c("#CE2931", "#1D2D5F"), labels = c('Xpert-only', 'BAU')) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 100), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,1), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2025, 2035), ylim = c(0, 1)) + 
  labs(x = "Year", y = "TB prevalence proportional reduction") +
  geom_rect(aes(xmin = 2025, xmax = 2028, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) + 
  geom_point(data = comp, aes(x = time, y = 1 - val, colour = type), shape = 17) +
  geom_errorbar(data = comp, aes(x = time, ymin = 1 - hi, ymax = 1 - lo, colour = type), width = 0.5) + 
  theme_bw() +
  theme(legend.position = 'none', plot.margin = margin(10,15,5,10,"pt"), panel.grid.minor.x = element_blank())
dev.off()          
