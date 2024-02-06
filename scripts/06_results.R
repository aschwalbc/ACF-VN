## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 06: Results.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame
options(scipen = 999)

# 1. Load data ==========
folder <- here("outputs", "results")
files <- list.files(folder, full.names = TRUE)

for (i in files) {
  file <- tools::file_path_sans_ext(basename(i))
  dt <- setDT(import(i))
  assign(file, dt)
  print(file)
}

rm(list = setdiff(ls(), ls(pattern = "outacf|base")))

# 2. Data curation ==========
df_names <- ls(pattern = "^r")

outs <- rbindlist(lapply(df_names, function(df_name) {
  df <- get(df_name)
  
  result <- df %>% 
    select(type, run, round, time, everything()) %>% 
    select(-contains("_"))
  
  return(result)
}))

rm(list = setdiff(ls(), "outs"))

outs <- outs %>% 
  mutate(cBAU = cBAUs + cBAUr, # Cost of BAU (DSTB + DRTB)
         cRxBAU = cRxBAUs + cRxBAUr, # Cost of treatment BAU (DSTB + DRTB)
         cRxTP = cRxTPs + cRxTPr,  # Cost of treatment ACF TP (DSTB + DRTB)
         cRxFP = cRxFPs + cRxFPr,  # Cost of treatment ACF FP (DSTB + DRTB)
         cRxACF = cRxTPs + cRxTPr + cRxFPs + cRxFPr, # Cost of treatment ACF
         cCF = cBAUs + cBAUr + cACF, # Cost of case finding
         cRx = cRxBAUs + cRxBAUr + cRxTPs + cRxTPr + cRxFPs + cRxFPr) %>% # Cost of treatment
  mutate(goal = case_when(type == 'acfa' & round == '02' ~ '100', type == 'acfa' & round == '05' ~ '50', type == 'acfa' & round == '10' ~ '20',
                          type == 'acfb' & round == '02' ~ '100', type == 'acfb' & round == '05' ~ '50', type == 'acfb' & round == '12' ~ '20',
                          type == 'acfc' & round == '01' ~ '100', type == 'acfc' & round == '02' ~ '50', type == 'acfc' & round == '03' ~ '20',
                          type == 'base' ~ 'none'))
  
outs_m <- outs %>%
  arrange(type, run, time) %>%
  group_by(run, time) %>%
  mutate(prTBc = tTBc / tTBc[type == 'base'],  # Proportional TB prevalence to BAU
         redTBc = 1 - (tTBc / tTBc[type == 'base']), # Proportional TB prevalence reduction 
         dfTBc = tTBc - tTBc[type == 'base'],  # TB prevalence difference to BAU
         prInc = tInc / tInc[type == 'base'],  # Proportional TB incidence to BAU
         dfInc = tInc - tInc[type == 'base'],  # TB incidence difference to BAU
         prMor = rMor / rMor[type == 'base'],  # Proportional TB mortality to BAU
         dfMor = tMor - tMor[type == 'base']) %>% # TB mortality difference to BAU
  ungroup() %>%
  group_by(type, run) %>% 
  mutate(prFP = ifelse(tFPos == 0, NA, tFPos / (tTPos + tFPos))) %>% # Proportion FP
  ungroup() %>% 
  pivot_longer(cols = -c(time, type, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE),
            lo = quantile(values, 0.025, na.rm = TRUE),
            hi = quantile(values, 0.975, na.rm = TRUE))
  
outs_yr <- outs %>% 
  arrange(type, run, time) %>% 
  filter(time == floor(time)) %>% 
  group_by(type, run, round) %>%
  mutate(cumTBc = cumsum(tTBc), # Cumulative TB prevalence
         cumInc = cumsum(tInc), # Cumulative TB incidence
         cumMor = cumsum(tMor), # Cumulative TB mortality
         cumTPos = cumsum(tTPos), # Cumulative TP diagnoses
         cumFPos = cumsum(tFPos), # Cumulative FP diagnoses
         cumScrn = cumsum(tScrn), # Cumulative population screened
         cumcACF = cumsum(cACF), # Cumulative cost of ACF
         cumcRxTP = cumsum(cRxTP), # Cumulative TP treatment
         cumcRxFP = cumsum(cRxFP), # Cumulative FP treatment
         cumcRxACF = cumsum(cRxACF), # Cumulative ACF treatment
         cumcCF = cumsum(cCF), # Cumulative cost of case finding
         cumcRx = cumsum(cRx), # Cumulative cost of all treatments
         cumDALYs = cumsum(DALYs)) %>% #Cumulative DALYs
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(dfcumInc = cumInc - cumInc[type == 'base'], # TB incidence averted
         dfcumMor = cumMor - cumMor[type == 'base'], # TB mortality averted
         dfcumcCF = cumcCF - cumcCF[type == 'base'], # CF cost averted
         dfcumcRx = cumcRx - cumcRx[type == 'base'], # Rx cost averted
         dfcumDALYs = cumDALYs - cumDALYs[type == 'base']) %>% # DALYs averted
  ungroup() %>% 
  group_by(type, run) %>% 
  mutate(cumprFP = ifelse(cumFPos == 0, NA, cumFPos / (cumTPos + cumFPos))) %>% 
  ungroup() %>%
  select(time, type, run, round, goal, contains('cum')) %>%
  pivot_longer(cols = -c(time, type, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

outs <- outs_m %>% 
  rbind(outs_yr) %>% 
  arrange(factor(type, levels = c("base", "acfa", "acfb", "acfc")), round, time) 
rm(outs_m, outs_yr)

export(outs, here("outputs", "outs", "outs.Rdata"))

# 3. Extract results ==========
outs <- import(here("outputs", "outs", "outs.Rdata"))

outini <- outs %>% 
  filter((type == 'acfa' & round == '02' & time == 2027) |
           (type == 'acfa' & round == '05' & time == 2030) |
           (type == 'acfa' & round == '10' & time == 2035) |
           (type == 'acfb' & round == '02' & time == 2027) |
           (type == 'acfb' & round == '05' & time == 2030) |
           (type == 'acfb' & round == '12' & time == 2037) |
           (type == 'acfc' & round == '01' & time == 2026) |
           (type == 'acfc' & round == '02' & time == 2027) |
           (type == 'acfc' & round == '03' & time == 2028))

filter(outini, var == 'rTBc') # TB prevalence rate
filter(outini, var == 'redTBc') # TB prevalence rate reduction
filter(outini, var == 'rInc') # TB incidence rate
filter(outini, var == 'rMor') # TB mortality rate

outfin <- outs %>% 
  filter(time == 2050 & 
           ((type == 'acfa' & round == '02') | (type == 'acfa' & round == '05') | (type == 'acfa' & round == '10') |
            (type == 'acfb' & round == '02') | (type == 'acfb' & round == '05') | (type == 'acfb' & round == '12') | 
            (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03')))

filter(outfin, var == 'dfcumInc') # TB incidence averted 
filter(outfin, var == 'dfcumMor') # TB mortality averted
filter(outfin, var == 'cumcACF') # Total costs of ACF screening
filter(outfin, var == 'cumcCF') # Total costs of case finding
filter(outfin, var == 'cumTPos') # Total true positives in ACF
filter(outfin, var == 'cumFPos') # Total true positives in ACF
filter(outfin, var == 'cumprFP') # Proportion FP treated
filter(outfin, var == 'cumcRxACF') # Total costs of ACF treatment
filter(outfin, var == 'cumcRx') # Total costs of treatment
filter(outfin, var == 'dfcumDALYs') # DALYs averted

icer <- outs %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
          (type == 'acfa' & round == '02') | (type == 'acfa' & round == '05') | (type == 'acfa' & round == '10') |
          (type == 'acfb' & round == '02') | (type == 'acfb' & round == '05') | (type == 'acfb' & round == '12') | 
          (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03'))) %>%
  filter(var == 'dfcumDALYs' | var == 'cumcCF' | var == 'cumcRx') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  mutate(val_COST = val_cumcCF + val_cumcRx, lo_COST = lo_cumcCF + lo_cumcRx, hi_COST = hi_cumcCF + hi_cumcRx) %>% 
  rename(val_DALY = val_dfcumDALYs, lo_DALY = lo_dfcumDALYs, hi_DALY = hi_dfcumDALYs) %>% 
  select(type, round, val_COST, lo_COST, hi_COST, val_DALY, lo_DALY, hi_DALY) %>% 
  mutate(cost_val = (val_COST - val_COST[type == 'base']), 
         cost_lo = (lo_COST - lo_COST[type == 'base']),
         cost_hi = (hi_COST - hi_COST[type == 'base']),
         daly_val = (abs(val_DALY) - abs(val_DALY[type == 'base'])),
         daly_lo = (abs(lo_DALY) - lo_DALY[type == 'base']),
         daly_hi = (abs(hi_DALY) - hi_DALY[type == 'base'])) %>% 
  select(type, round, starts_with('cost'), starts_with('daly')) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = cost_val / daly_val, icer_lo = cost_lo / daly_lo, icer_hi = cost_hi / daly_hi) %>% 
  select(type, round, starts_with('icer')) %>% 
  mutate(var = 'icer') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(type, round, var, val, lo, hi)
