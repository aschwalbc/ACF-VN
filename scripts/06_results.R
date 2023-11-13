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

# 2. Data curation ==========
outs <- list()

for (i in 2:length(runs)) {
  outs[[i]] <- rbind(runs[[1]], runs[[i]][[1]], runs[[i]][[1]], runs[[i]][[1]]) %>% 
    arrange(type, run, time) %>% 
    group_by(type, run) %>%
    mutate(ACF = SCmin+SCsub+SCcln, # Sum of all TB disease screened (Min+Sub+Cln)
           AllTB = SCmin+SCsub+SCcln+DXcln, # Sum of all TB disease diagnoses (ACF+BAU)
           AllTx = SCmin+SCsub+SCcln+DXcln+FPnds, # Sum of all diagnoses (ACF+BAU+FP)
           cumMor = cumsum(tMor), # Cumulative TB mortality
           cumFPnds = cumsum(FPnds), # Cumulative FP diagnoses
           cumFNdis = cumsum(FNdis), # Cumulative FN (diagnoses missed)
           cumNumSC = cumsum(NumSC), # Cumulative screening (ACF+FP)
           cumDXcln = cumsum(DXcln), # Cumulative BAU diagnoses 
           pFP = ifelse(FPnds == 0, 0, FPnds/(TPdis+FPnds))) %>% # FP treated (% over confirmed)
    ungroup() %>% 
    group_by(run, time) %>% 
    mutate(dMor = tMor - tMor[type == 'base'], # TB mortality difference to BAU
           dcumMor = cumMor - cumMor[type == 'base'], # Cumulative TB mortality difference to BAU
           pTBc = TBc/TBc[type == 'base'], # Proportional reduction TB prevalence to BAU
           pMor = rMor/rMor[type == 'base'], # Proportional reduction TB mortality to BAU
           dAllTB = AllTB - AllTB[type == 'base'], # All TB diagnoses difference to BAU
           dAllTx = AllTx - AllTx[type == 'base'], # All diagnoses (TB+FP) difference to BAU
           dBAU = cumDXcln - cumDXcln[type == 'base']) %>% # BAU TB diagnoses difference to BAU
    ungroup() %>% 
    group_by(type, run) %>% 
    mutate(cumAllTB = cumsum(dAllTB), cumAllTx = cumsum(dAllTx)) %>% 
    ungroup() %>% 
    group_by(run, time) %>% 
    mutate(NNS = 1/pTBc, # Number needed to screen
           NNT = 1/pMor) %>% # Number needed to treat to avert a TB death
    ungroup() %>% 
    group_by(type, run, time) %>% 
    mutate(pPRur = PRur/Pop, # Proportion rural
           pPUrb = PUrb/Pop, # Proportion urban
           pPHig = PHig/Pop, # Proportion high SES
           pPLow = PLow/Pop, # Proportion low SES
           pPRL = PRL/Pop, # Proportion rural - low SES
           pPRH = PRH/Pop, # Proportion rural - high SES
           pPUL = PUL/Pop, # Proportion urban - low SES
           pPUH = PUH/Pop) %>% # Proportion urban - high SES
    pivot_longer(cols = -c(time, type, run), names_to = "var", values_to = "values") %>% 
    group_by(time, type, var) %>% 
    summarise(val = median(values, na.rm = TRUE), 
              lo = quantile(values, 0.025, na.rm = TRUE), 
              hi = quantile(values, 0.975, na.rm = TRUE)) %>% 
    mutate(fill = ifelse(val < 0, "under", "over"))
}

r1outs <- rbind(runs[[1]][[1]], runs[[2]][[1]])

for (i in 2:length(runs)) {
  outs[[i]] <- rbindlist(list(runs[[i]][[1]], runs[[i]][[2]], runs[[i]][[3]]))
}

outs <- rbindlist(c(list(runs[[1]][[1]]), outs))

# export(outs, here("outputs","runs.Rdata"))
# outs <- import(here("outputs","runs.Rdata"))
