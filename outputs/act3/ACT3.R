# ACT3 TB Prevalence

library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(extrafont) # Add specific fonts
set.seed(123)

act3 <- import(here("outputs","act3","ACT3.csv"))
act3$group <- factor(act3$group, levels = c('int', 'ctrl'))

png(here("outputs","act3","ACT3.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot(act3, aes(x = year, y = rate, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, 
                position = position_dodge(width = 0.9)) +
  labs(x = 'Year', y = 'TB prevalence rate (per 100k)', fill = 'Group') +
  scale_fill_manual(values = c("int" = "salmon", "ctrl" = "turquoise"),
                    labels = c("int" = "Intervention", "ctrl" = "Control")) +
  theme_bw() +
  theme(legend.position = 'bottom', text = element_text(family = "Open Sans"))
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

inty1inty4 <- data.frame(type = 'Intervention Y1 v Intervention Y4',time = 2028, val = abs(summary(prpdf)[[3]]),
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

inty1ctrly4 <- data.frame(type = 'Intervention Y1 v Control Y4',time = 2028, val = abs(summary(prpdf)[[3]]),
                         lo = abs(summary(prpdf)[[5]]), hi = abs(summary(prpdf)[[2]]))

# Comparison graph
act <- import(here("outputs","act3","act3_ultra.Rdata"))
# act <- import(here("outputs","act3","act3_mtbrif.Rdata"))
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
            hi = quantile(values, 0.975, na.rm = TRUE)) %>% 
  filter(type != 'base')

filter(df, time == 2028 & type == 'acfa')

png(here("plots","00_act3.png"), width = 6, height = 5, units = 'in', res = 1000)
ggplot() +
  geom_line(filter(df, var == 'prTBc'), mapping = aes(x = time, y = val, colour = type)) + 
  geom_ribbon(filter(df, var == 'prTBc'), mapping = aes(x = time, ymin = lo, ymax = hi, fill = type), alpha = 0.2) +
  scale_color_manual(values = c("#CE2931","#828282","#828282"), labels = c('Xpert-only', 'BAU','IntY1-CtrlY4','IntY1-IntY4')) +
  scale_fill_manual(values = c("#CE2931"), labels = c('Xpert-only', 'BAU')) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 100), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2020,2050,1), limits = c(2020,2050), expand = c(0,0)) +
  coord_cartesian(xlim = c(2025, 2035), ylim = c(0, 1)) + 
  labs(x = "Year", y = "Proportional reduction in TB prevalence") +
  geom_rect(aes(xmin = 2025, xmax = 2028, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2) + 
  geom_point(data = comp, aes(x = time, y = 1 - val, colour = type, shape = type), size = 2) +
  scale_shape_manual(values = c(15, 17)) +
  geom_errorbar(data = comp, aes(x = time, ymin = 1 - hi, ymax = 1 - lo, colour = type), width = 0.5) + 
  theme_bw() +
  theme(legend.position = 'bottom', plot.margin = margin(10,15,5,10,"pt"), panel.grid.minor.x = element_blank(),
        text = element_text(family = "Open Sans"), axis.title.x = element_text(margin = margin(t = 10))) + 
  guides(shape = guide_legend(title = "Outcome"), colour = 'none', fill = 'none')
dev.off()          
