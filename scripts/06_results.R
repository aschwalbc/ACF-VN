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

options(scipen = 999)

outs_m <- outs %>% 
  arrange(type, run, time) %>% 
  mutate(cPCF = cPCFs + cPCFr, # Cost of PCF (DSTB + DRTB)
         cRxPCF = cRxPCFs + cRxPCFr, # Cost of treatment PCF (DSTB + DRTB)
         cRxTP = cRxTPs + cRxTPr, # Cost of treatment ACF TP (DSTB + DRTB)
         cRxFP = cRxFPs + cRxFPr) %>% # Cost of treatment ACF FP (DSTB + DRTB)
  mutate(cCF = cPCF + cACF,
         cRx = cRxPCF + cRxTP + cRxFP) %>% 
  group_by(run, time) %>% 
  mutate(prTBc = tTBc/tTBc[type == 'base'], # Proportional TB prevalence to BAU
         dfTBc = abs(tTBc - tTBc[type == 'base']), # TB prevalence difference to BAU
         prMor = rMor/rMor[type == 'base'], # Proportional TB mortality to BAU
         dfMor = abs(tMor - tMor[type == 'base'])) %>%  # TB mortality difference to BAU
  ungroup() %>% 
  mutate(redTBc = 1-prTBc) %>% # Proportional TB prevalence reduction
  group_by(type, run) %>% 
  mutate(prFP = ifelse(tFPos == 0, 0, tFPos/(tTPos+tFPos))) %>% 
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
         cumFPos = cumsum(tFPos), # Cumulative FP diagnoses
         cumTPos = cumsum(tTPos), # Cumulative TP diagnoses
         cumRx = cumsum(tFPos) + cumsum(tTPos)) %>% # Cumulative treated
  ungroup() %>% 
  group_by(type, run) %>% 
  mutate(cumprFP = ifelse(cumFPos == 0, 0, cumFPos/cumRx)) %>% 
  select(time, type, run, round, cumTBc, cumMor, cumFPos, cumTPos, cumRx, cumprFP) %>% 
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
