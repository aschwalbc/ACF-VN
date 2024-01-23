## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 06: Results.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame

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

outs_m <- outs %>%
  arrange(type, run, time) %>%
  mutate(cPCF = cPCFs + cPCFr,  # Cost of PCF (DSTB + DRTB)
         cRxPCF = cRxPCFs + cRxPCFr,  # Cost of treatment PCF (DSTB + DRTB)
         cRxTP = cRxTPs + cRxTPr,  # Cost of treatment ACF TP (DSTB + DRTB)
         cRxFP = cRxFPs + cRxFPr,  # Cost of treatment ACF FP (DSTB + DRTB)
         cCF = cPCF + cACF, # Cost of case finding
         cRx = cRxPCF + cRxTP + cRxFP) %>% # Cost of treatment
  group_by(run, time) %>%
  mutate(prTBc = tTBc / tTBc[type == 'base'],  # Proportional TB prevalence to BAU
         dfTBc = abs(tTBc - tTBc[type == 'base']),  # TB prevalence difference to BAU
         prInc = tInc / tInc[type == 'base'],  # Proportional TB incidence to BAU
         dfInc = abs(tInc - tInc[type == 'base']),  # TB incidence difference to BAU
         prMor = rMor / rMor[type == 'base'],  # Proportional TB mortality to BAU
         dfMor = abs(tMor - tMor[type == 'base']), # TB mortality difference to BAU
         redTBc = 1 - (tTBc / tTBc[type == 'base'])) %>% # Proportional TB prevalence reduction 
  ungroup() %>%
  group_by(type, run) %>% 
  mutate(prFP = ifelse(tFPos == 0, NA, tFPos / (tTPos + tFPos))) %>% # Proportion FP
  ungroup() %>% 
  pivot_longer(cols = -c(time, type, run, round), names_to = "var", values_to = "values") %>%
  group_by(time, type, round, var) %>%
  summarise(val = median(values, na.rm = TRUE),
            lo = quantile(values, 0.025, na.rm = TRUE),
            hi = quantile(values, 0.975, na.rm = TRUE))
  
outs_yr <- outs %>% 
  arrange(type, run, time) %>% 
  filter(time == floor(time)) %>% 
  group_by(type, run, round) %>%
  mutate(cumTBc = cumsum(tTBc), # Cumulative TB prevalence
         cumMor = cumsum(tMor), # Cumulative TB mortality
         cumInc = cumsum(tInc), # Cumulative TB incidence
         cumFPos = cumsum(tFPos), # Cumulative FP diagnoses
         cumTPos = cumsum(tTPos), # Cumulative TP diagnoses
         cumRx = cumsum(tFPos) + cumsum(tTPos)) %>% # Cumulative treated
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(dfcumInc = abs(cumInc - cumInc[type == 'base']), # TB incidence averted
         dfcumMor = abs(cumMor - cumMor[type == 'base'])) %>% # TB mortality averted
  ungroup() %>% 
  group_by(type, run) %>% 
  mutate(cumprFP = ifelse(cumFPos == 0, NA, cumFPos / cumRx)) %>% 
  ungroup() %>%
  select(time, type, run, round, contains('cum')) %>%
  pivot_longer(cols = -c(time, type, run, round), names_to = "var", values_to = "values") %>%
  group_by(time, type, round, var) %>%
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
filter(outini, var == 'rMor') # TB mortality rate

outfin <- outs %>% 
  filter(time == 2050 & 
           ((type == 'acfa' & round == '02') | (type == 'acfa' & round == '05') | (type == 'acfa' & round == '10') |
            (type == 'acfb' & round == '02') | (type == 'acfb' & round == '05') | (type == 'acfb' & round == '12') | 
            (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03')))

filter(outfin, var == 'dfcumInc') # TB incidence averted 
filter(outfin, var == 'dfcumMor') # TB mortality averted
filter(outfin, var == 'cumFPos') # Cumulative false positives
filter(outfin, var == 'cumTPos') # Cumulative true positives
filter(outfin, var == 'cumRx') # Cumulative treated
filter(outfin, var == 'cumprFP') # Proportion FP treated


