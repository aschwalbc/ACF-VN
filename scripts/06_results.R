## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 06: Results.R

# Packages ==========
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(data.table) # Faster than data.frame

# 1. Load data ==========
runs <- list()

folder <- here("outputs", "results")
files <- list.files(folder, full.names = TRUE)

for (i in files) {
  file <- tools::file_path_sans_ext(basename(i))
  dt <- setDT(import(i))
  prefix <- substring(file, 1, 3)
  if (!prefix %in% names(runs)) {
    runs[[prefix]] <- list()
  }
  runs[[prefix]] <- append(runs[[prefix]], list(dt))
  print(file)
}

rm(list = setdiff(ls(), "runs"))

for (i in 1:length(runs)) {
  name <- paste0(names(runs)[i], "outs")
  if (i == 1) {
    dt <- rbindlist(list(runs[[i]][[1]]))
  } else {
    dt <- rbindlist(list(runs[[i]][[1]], runs[[i]][[2]], runs[[i]][[3]]))
  }
  assign(name, dt)
  print(name)
}

rm(dt, runs, i, name)

# 2. Data curation ==========
df_names <- ls(pattern = "outs")

outs <- rbindlist(lapply(df_names, function(df_name) {
  df <- get(df_name)
  df %>% select(type, run, time, round, Pop, Sub, Cln, TBc, rMor, tMor, FPnds, TPdis)
}))

rm(list = setdiff(ls(), "outs"))

outs_m <- outs %>% 
  arrange(type, run, time) %>% 
  group_by(run, time) %>% 
  mutate(pTBc = TBc/TBc[type == 'base'], # Proportional reduction TB prevalence to BAU
         dMor = tMor - tMor[type == 'base'], # TB mortality difference to BAU
         pMor = rMor/rMor[type == 'base']) %>%  # Proportional reduction TB mortality to BAU
  ungroup() %>% 
  mutate(redpTBc = 1-pTBc) %>% # Proportional reduction TB mortality to BAU
  group_by(type, run) %>% 
  mutate(pFP = ifelse(FPnds == 0, 0, FPnds/(TPdis+FPnds))) %>% 
  pivot_longer(cols = -c(time, type, run, round), names_to = "var", values_to = "values") %>% 
  group_by(time, type, var, round) %>% 
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

outs_yr <- outs %>% 
  arrange(type, run, time) %>% 
  filter(time == floor(time)) %>% 
  group_by(type, run, round) %>%
  mutate(cumMor = cumsum(tMor), # Cumulative TB mortality
         cumFPnds = cumsum(FPnds), # Cumulative FP diagnoses
         cumTPdis = cumsum(TPdis)) %>% # Cumulative TP diagnoses
  ungroup() %>% 
  group_by(type, run) %>% 
  mutate(cumpFP = ifelse(cumFPnds == 0, 0, cumFPnds/(cumTPdis+cumFPnds))) %>% 
  pivot_longer(cols = -c(time, type, run, round), names_to = "var", values_to = "values") %>% 
  group_by(time, type, var, round) %>% 
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

export(outs, here("outputs", "outs", "outs.Rdata"))
export(outs_m, here("outputs", "outs", "outs_m.Rdata"))
export(outs_yr, here("outputs", "outs", "outs_yr.Rdata"))
