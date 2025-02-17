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
    within(rm(SUS, SUSx, INF, INFx, CLE, CLEx, REC, RECx, MIN, MINx, SUB, SUBx, CLN, CLNx, TXT, TRE, TREx))
  
  return(result)
}))

rm(list = setdiff(ls(), "outs"))

outs <- outs %>% 
  mutate(tFP = tSUS + tINF + tCLE + tREC + tTRE, # Total FP
         tTP = tMIN + tSUB + tCLN, # Total TP
         tTPinf = tSUB + tCLN, # Total infectious TP
         tACF = tSUS + tINF + tCLE + tREC + tTRE + tMIN + tSUB + tCLN, # Total ACF
         cBAU = cBAUs + cBAUr, # Cost of BAU (DSTB + DRTB)
         cACFFP = cACFSUS + cACFINF + cACFCLE + cACFREC + cACFTRE + cCLAPSUS + cCLAPINF + cCLAPCLE + cCLAPREC + cCLAPTRE, # Cost of ACF FP
         cACFTP = cACFMIN + cACFSUB + cACFCLN + cCLAPMIN + cCLAPSUB + cCLAPCLN, # Cost of ACF TP
         cACFTPinf = cACFSUB + cACFCLN + cCLAPSUB + cCLAPCLN, # Cost of ACF infectious TP
         cRxBAU = cRxBAUs + cRxBAUr, # Cost of treatment BAU (DSTB + DRTB)
         cRxACFFP = cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE, # Cost of treatment ACF FP 
         cRxACFTP = cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of treatment ACF TP (DSTB + DRTB)
         cRxACFTPinf = cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of treatment ACF infectious TP (DSTB + DRTB)
         cCF = cBAUs + cBAUr + cACF, # Cost of case finding (PCF + ACF)
         cRxACF = cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE + cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of treatment ACF
         cRx = cRxBAUs + cRxBAUr + cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE + cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of treatment (PCF + ACF)
         cIntv = cACF + cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE + cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of intervention (-BAU)
         cAll = cBAUs + cBAUr + cACF + cCLAP + cRxBAUs + cRxBAUr + cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE + cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr) %>% # All costs
  mutate(goal = case_when(type == 'acfa' & round == '03' ~ '100', type == 'acfa' & round == '06' ~ '50', type == 'acfa' & round == '11' ~ '20',
                          type == 'acfb' & round == '03' ~ '100', type == 'acfb' & round == '08' ~ '50', type == 'acfb' & round == '14' ~ '20',
                          type == 'acfc' & round == '02' ~ '100', type == 'acfc' & round == '03' ~ '50', type == 'acfc' & round == '04' ~ '20',
                          type == 'acfd' & round == '03' ~ '100', type == 'acfd' & round == '08' ~ '50', type == 'acfd' & round == '14' ~ '20',
                          type == 'acfax' & round == '03' ~ '100', type == 'acfax' & round == '06' ~ '50', type == 'acfax' & round == '11' ~ '20',
                          type == 'acfbx' & round == '03' ~ '100', type == 'acfbx' & round == '08' ~ '50', type == 'acfbx' & round == '14' ~ '20',
                          type == 'acfac' & round == '03' ~ '100', type == 'acfac' & round == '06' ~ '50', type == 'acfac' & round == '11' ~ '20',
                          type == 'acfbc' & round == '03' ~ '100', type == 'acfbc' & round == '08' ~ '50', type == 'acfbc' & round == '14' ~ '20',
                          type == 'acfcc' & round == '02' ~ '100', type == 'acfcc' & round == '03' ~ '50', type == 'acfcc' & round == '04' ~ '20',
                          type == 'acfbr' & round == '04' ~ '100', type == 'acfbr' & round == '09' ~ '50', type == 'acfbr' & round == '16' ~ '20',
                          type == 'acfcr' & round == '02' ~ '100', type == 'acfcr' & round == '03' ~ '50', type == 'acfcr' & round == '05' ~ '20',
                          type == 'base' ~ 'none')) %>% 
  mutate(goal = factor(goal, levels = c('100', '50', '20', 'none'))) %>% 
  mutate(res = case_when(type == 'acfa' ~ 'main', type == 'acfb' ~ 'main', type == 'acfc' ~ 'main', 
                         type == 'acfd' ~ 'sa: mtb/rif', type == 'acfax' ~ 'sa: 1usd', type == 'acfbx' ~ 'sa: 1usd',
                         type == 'acfac' ~ 'sa: clnap', type == 'acfbc' ~ 'sa: clnap', type == 'acfcc' ~ 'sa: clnap', 
                         type == 'acfbr' ~ 'sa: cxrpos', type == 'acfcr' ~ 'sa: cxrpos', type == 'base' ~ 'main'))
  
outs_m <- outs %>%
  arrange(type, run, time) %>%
  group_by(run, time) %>%
  mutate(prTBc = tTBc / tTBc[type == 'base'], # Proportional TB prevalence to BAU
         redTBc = 1 - (tTBc / tTBc[type == 'base']), # Proportional TB prevalence reduction 
         dfTBc = tTBc - tTBc[type == 'base'], # TB prevalence difference to BAU
         prInc = tInc / tInc[type == 'base'], # Proportional TB incidence to BAU
         dfInc = tInc - tInc[type == 'base'], # TB incidence difference to BAU
         prMor = rMor / rMor[type == 'base'], # Proportional TB mortality to BAU
         dfMor = tMor - tMor[type == 'base'], # TB mortality difference to BAU
         dfcBAU = cBAU - cBAU[type == 'base'], # BAU diagnosis costs difference to BAU
         dfcRxBAU = cRxBAU - cRxBAU[type == 'base']) %>% # BAU treatment costs difference to BAU
  ungroup() %>%
  pivot_longer(cols = -c(time, type, res, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, res, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE),
            lo = quantile(values, 0.025, na.rm = TRUE),
            hi = quantile(values, 0.975, na.rm = TRUE))

outs_yr <- outs %>% 
  arrange(type, res, run, time) %>% 
  filter(time == floor(time)) %>% 
  filter(time >= 2025) %>% 
  group_by(type, res, run, round) %>%
  mutate(cumTBc = cumsum(tTBc), # Cumulative TB prevalence
         cumInc = cumsum(tInc), # Cumulative TB incidence
         cumMor = cumsum(tMor), # Cumulative TB mortality
         cumScrn = cumsum(tScrn), # Cumulative screened
         cumSUS = cumsum(tSUS), # Cumulative diagnosed SUS
         cumINF = cumsum(tINF), # Cumulative diagnosed INF
         cumCLE = cumsum(tCLE), # Cumulative diagnosed CLE
         cumREC = cumsum(tREC), # Cumulative diagnosed REC
         cumMIN = cumsum(tMIN), # Cumulative diagnosed MIN
         cumSUB = cumsum(tSUB), # Cumulative diagnosed SUB
         cumCLN = cumsum(tCLN), # Cumulative diagnosed CLN
         cumTRE = cumsum(tTRE), # Cumulative diagnosed TRE
         cumFP = cumsum(tFP), # Cumulative diagnosed FP
         cumTP = cumsum(tTP), # Cumulative diagnosed TP
         cumTPinf = cumsum(tTPinf), # Cumulative diagnosed infectious TP
         cumACF = cumsum(tACF), # Cumulative diagnosed ACF
         cumcBAU = cumsum(cBAU), # Cumulative cost of BAU
         cumcACFSUS = cumsum(cACFSUS + cCLAPSUS), # Cumulative cost of ACF SUS
         cumcACFINF = cumsum(cACFINF + cCLAPINF), # Cumulative cost of ACF INF
         cumcACFCLE = cumsum(cACFCLE + cCLAPCLE), # Cumulative cost of ACF CLE
         cumcACFREC = cumsum(cACFREC + cCLAPREC), # Cumulative cost of ACF REC
         cumcACFMIN = cumsum(cACFMIN + cCLAPMIN), # Cumulative cost of ACF MIN
         cumcACFSUB = cumsum(cACFSUB + cCLAPSUB), # Cumulative cost of ACF SUB
         cumcACFCLN = cumsum(cACFCLN + cCLAPCLN), # Cumulative cost of ACF CLN
         cumcACFTRE = cumsum(cACFTRE + cCLAPTRE), # Cumulative cost of ACF TRE
         cumcACFFP = cumsum(cACFFP + cCLAPSUS + cCLAPINF + cCLAPCLE + cCLAPREC + cCLAPTRE), # Cumulative cost of ACF FP
         cumcACFTP = cumsum(cACFTP + cCLAPMIN + cCLAPSUB + cCLAPCLN), # Cumulative cost of ACF TP
         cumcACFTPinf = cumsum(cACFTPinf + cCLAPSUB + cCLAPCLN), # Cumulative cost of ACF infectious TP
         cumcACF = cumsum(cACF + cCLAP), # Cumulative cost of ACF
         cumcRxBAU = cumsum(cRxBAU), # Cumulative cost of BAU treatment
         cumcRxSUS = cumsum(cRxSUS), # Cumulative cost of treatment SUS
         cumcRxINF = cumsum(cRxINF), # Cumulative cost of treatment INF
         cumcRxCLE = cumsum(cRxCLE), # Cumulative cost of treatment CLE
         cumcRxREC = cumsum(cRxREC), # Cumulative cost of treatment REC
         cumcRxMIN = cumsum(cRxMIN), # Cumulative cost of treatment MIN
         cumcRxSUB = cumsum(cRxSUBs + cRxSUBr), # Cumulative cost of treatment SUB (DSTB + DRTB)
         cumcRxCLN = cumsum(cRxCLNs + cRxCLNr), # Cumulative cost of treatment CLN (DSTB + DRTB)
         cumcRxTRE = cumsum(cRxTRE), # Cumulative cost of treatment TRE
         cumcRxACFFP = cumsum(cRxACFFP), # Cumulative cost of treatment FP
         cumcRxACFTP = cumsum(cRxACFTP), # Cumulative cost of treatment TP
         cumcRxACFTPinf = cumsum(cRxACFTPinf), # Cumulative cost of treatment infectious TP
         cumcRxACF = cumsum(cRxACF), # Cumulative cost of treatment ACF
         cumcIntv = cumsum(cIntv), # Cumulative intervention costs
         cumcPCF = cumsum(cBAU + cRxBAU), # Cumulative BAU costs (diagnosis + treatment)
         cumcCF = cumsum(cCF), # Cumulative cost of case finding
         cumcRx = cumsum(cRx), # Cumulative cost of all treatments
         cumDALYs = cumsum(DALYs), # Cumulative DALYs
         cumcAll = cumsum(cAll)) %>% # Cumulative costs of all 
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(dfcumInc = cumInc - cumInc[type == 'base'], # TB incidence averted
         dfcumMor = cumMor - cumMor[type == 'base'], # TB mortality averted
         dfcumcBAU = cumcBAU - cumcBAU[type == 'base'], # BAU cost averted
         dfcumcCF = cumcCF - cumcCF[type == 'base'], # CF cost averted
         dfcumcRxBAU = cumcRxBAU - cumcRxBAU[type == 'base'], # BAU Rx cost averted
         dfcumcRx = cumcRx - cumcRx[type == 'base'], # Rx cost averted
         dfcumcPCF = cumcPCF - cumcPCF[type == 'base'], # BAU (diagnosis + treatment) averted
         dfcumcAll = cumcAll - cumcAll[type == 'base'], # All costs averted
         dfcumDALYs = cumDALYs - cumDALYs[type == 'base'], # DALYs averted
         prpcACF = cumcACF / (cumcACF + cumcRxACF), # Proportion costs of ACF screening over total costs
         prpcBAU = (cumcBAU + cumcRxBAU) / cumcAll) %>% # Proportion costs due to BAU
  ungroup() %>% 
  mutate(cumprFP = ifelse(cumFP == 0, NA, cumFP / (cumACF)), # Proportion FP
         cumrtFP = ifelse(cumFP == 0, NA, cumFP / cumTP)) %>% # Ratio TP:FP)
  select(time, type, res, run, round, goal, contains('cum'), contains('prp'), contains('cp')) %>%
  pivot_longer(cols = -c(time, type, res, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, res, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

icers <- outs_yr %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
                           (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
                           (type == 'acfb' & round == '03') | (type == 'acfb' & round == '08') | (type == 'acfb' & round == '14') | 
                           (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') | (type == 'acfc' & round == '04') |
                           (type == 'acfd' & round == '03') | (type == 'acfd' & round == '08') | (type == 'acfd' & round == '14') |
                           (type == 'acfax' & round == '03') | (type == 'acfax' & round == '06') | (type == 'acfax' & round == '11') |
                           (type == 'acfbx' & round == '03') | (type == 'acfbx' & round == '08') | (type == 'acfbx' & round == '14') | 
                           (type == 'acfac' & round == '03') | (type == 'acfac' & round == '06') | (type == 'acfac' & round == '11') |
                           (type == 'acfbc' & round == '03') | (type == 'acfbc' & round == '08') | (type == 'acfbc' & round == '14') | 
                           (type == 'acfcc' & round == '02') | (type == 'acfcc' & round == '03') | (type == 'acfcc' & round == '04') |
                           (type == 'acfbr' & round == '04') | (type == 'acfbr' & round == '09') | (type == 'acfbr' & round == '16') | 
                           (type == 'acfcr' & round == '02') | (type == 'acfcr' & round == '03') | (type == 'acfcr' & round == '05'))) %>%
  filter(var == 'dfcumDALYs' | var == 'dfcumcAll') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  rename(val_COST = val_dfcumcAll, lo_COST = lo_dfcumcAll, hi_COST = hi_dfcumcAll,
         val_DALY = val_dfcumDALYs, lo_DALY = lo_dfcumDALYs, hi_DALY = hi_dfcumDALYs) %>% 
  select(time, type, res, round, goal, val_COST, lo_COST, hi_COST, val_DALY, lo_DALY, hi_DALY) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = abs(val_COST / val_DALY), icer_lo = abs(lo_COST / lo_DALY), icer_hi = abs(hi_COST / hi_DALY)) %>% 
  select(time, type, res, round, goal, starts_with('icer')) %>% 
  mutate(var = 'icer') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(time, type, res, round, goal, var, val, lo, hi)

daly <- outs %>%
  arrange(type, res, run, time) %>%
  filter(time == floor(time)) %>%
  filter(time >= 2025, type != 'base') %>%
  group_by(type, res, run, round) %>%
  mutate(cumDALYs = cumsum(DALYs)) %>%
  ungroup() %>%
  filter(time == 2050) %>%
  select(type, goal, run, cumDALYs) %>%
  pivot_wider(names_from = type, values_from = cumDALYs) %>%
  pivot_longer(cols = -c(goal, run), names_to = "type_1", values_to = "value_1") %>%
  expand_grid(type_2 = unique(outs$type)) %>%
  filter(type_2 != 'base', type_1 != type_2) %>% 
  left_join(outs %>%
              filter(time == floor(time)) %>%
              filter(time >= 2025, type != 'base') %>%
              group_by(type, res, run, round) %>%
              mutate(cumDALYs = cumsum(DALYs)) %>%
              ungroup() %>%
              filter(time == 2050, type != 'base') %>%
              select(type, goal, run, cumDALYs) %>%
              rename(value_2 = cumDALYs, type_2 = type), by = c("goal", "run", "type_2")) %>%
  mutate(diff = value_1 - value_2,
         comparison = paste(type_1, "v", type_2, sep = "")) %>%
  select(goal, run, comparison, diff) %>%
  group_by(goal, comparison) %>%
  summarise(val = median(diff, na.rm = TRUE), 
            lo = quantile(diff, 0.025, na.rm = TRUE), 
            hi = quantile(diff, 0.975, na.rm = TRUE),
            .groups = "drop")

cost <- outs %>%
  arrange(type, res, run, time) %>%
  filter(time == floor(time)) %>%
  filter(time >= 2025, type != 'base') %>%
  group_by(type, res, run, round) %>%
  mutate(cumcAll = cumsum(cAll)) %>%
  ungroup() %>%
  filter(time == 2050) %>%
  select(type, goal, run, cumcAll) %>%
  pivot_wider(names_from = type, values_from = cumcAll) %>%
  pivot_longer(cols = -c(goal, run), names_to = "type_1", values_to = "value_1") %>%
  expand_grid(type_2 = unique(outs$type)) %>%
  filter(type_2 != 'base', type_1 != type_2) %>% 
  left_join(outs %>%
              filter(time == floor(time)) %>%
              filter(time >= 2025, type != 'base') %>%
              group_by(type, res, run, round) %>%
              mutate(cumcAll = cumsum(cAll)) %>%
              ungroup() %>%
              filter(time == 2050, type != 'base') %>%
              select(type, goal, run, cumcAll) %>%
              rename(value_2 = cumcAll, type_2 = type), by = c("goal", "run", "type_2")) %>%
  mutate(diff = value_1 - value_2,
         comparison = paste(type_1, "v", type_2, sep = "")) %>%
  select(goal, run, comparison, diff) %>%
  group_by(goal, comparison) %>%
  summarise(val = median(diff, na.rm = TRUE), 
            lo = quantile(diff, 0.025, na.rm = TRUE), 
            hi = quantile(diff, 0.975, na.rm = TRUE),
            .groups = "drop")

outs <- rbind(outs_m, outs_yr, icers) %>% 
  arrange(factor(type, levels = c("base", "acfa", "acfb", "acfc", "acfd", "acfax", "acfbx", "acfac", "acfbc", "acfcc", "acfbr", "acfcr")), round, time) 
rm(outs_m, outs_yr, icers)

export(outs, here("outputs", "outs", "outs.Rdata"))
export(daly, here("outputs", "outs", "daly.Rdata"))
export(cost, here("outputs", "outs", "cost.Rdata"))

rm(outs, daly, cost)

# 3. Extract results ==========
outs <- import(here("outputs", "outs", "outs.Rdata"))
daly <- import(here("outputs", "outs", "daly.Rdata"))
cost <- import(here("outputs", "outs", "cost.Rdata"))

outini <- outs %>% 
  filter((type == 'acfa' & round == '03' & time == 2028) |
           (type == 'acfa' & round == '06' & time == 2031) |
           (type == 'acfa' & round == '11' & time == 2036) |
           (type == 'acfb' & round == '03' & time == 2028) |
           (type == 'acfb' & round == '08' & time == 2033) |
           (type == 'acfb' & round == '14' & time == 2039) |
           (type == 'acfc' & round == '02' & time == 2027) |
           (type == 'acfc' & round == '03' & time == 2028) |
           (type == 'acfc' & round == '04' & time == 2029) |
           (type == 'acfd' & round == '03' & time == 2028) |
           (type == 'acfd' & round == '08' & time == 2033) |
           (type == 'acfd' & round == '14' & time == 2039) |
           (type == 'acfax' & round == '03' & time == 2028) |
           (type == 'acfax' & round == '06' & time == 2031) |
           (type == 'acfax' & round == '11' & time == 2036) |
           (type == 'acfbx' & round == '03' & time == 2028) |
           (type == 'acfbx' & round == '08' & time == 2033) |
           (type == 'acfbx' & round == '14' & time == 2039) |
           (type == 'acfac' & round == '03' & time == 2028) |
           (type == 'acfac' & round == '06' & time == 2031) |
           (type == 'acfac' & round == '11' & time == 2036) |
           (type == 'acfbc' & round == '03' & time == 2028) |
           (type == 'acfbc' & round == '08' & time == 2033) |
           (type == 'acfbc' & round == '14' & time == 2039) |
           (type == 'acfcc' & round == '02' & time == 2027) |
           (type == 'acfcc' & round == '03' & time == 2028) |
           (type == 'acfcc' & round == '04' & time == 2029) |
           (type == 'acfbr' & round == '03' & time == 2028) |
           (type == 'acfbr' & round == '08' & time == 2033) |
           (type == 'acfbr' & round == '14' & time == 2039) |
           (type == 'acfcr' & round == '02' & time == 2027) |
           (type == 'acfcr' & round == '03' & time == 2028) |
           (type == 'acfcr' & round == '04' & time == 2029))

outfin <- outs %>% 
  filter(goal %in% c('50', 'none')) %>% # CHANGE THRESHOLD HERE
  filter(time == 2050) %>% 
  filter((type == 'base') |
           (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
           (type == 'acfb' & round == '03') | (type == 'acfb' & round == '08') | (type == 'acfb' & round == '14') | 
           (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') | (type == 'acfc' & round == '04') |
           (type == 'acfd' & round == '03') | (type == 'acfd' & round == '08') | (type == 'acfd' & round == '14') |
           (type == 'acfax' & round == '03') | (type == 'acfax' & round == '06') | (type == 'acfax' & round == '11') |
           (type == 'acfbx' & round == '03') | (type == 'acfbx' & round == '08') | (type == 'acfbx' & round == '14') | 
           (type == 'acfac' & round == '03') | (type == 'acfac' & round == '06') | (type == 'acfac' & round == '11') |
           (type == 'acfbc' & round == '03') | (type == 'acfbc' & round == '08') | (type == 'acfbc' & round == '14') | 
           (type == 'acfcc' & round == '02') | (type == 'acfcc' & round == '03') | (type == 'acfcc' & round == '04') |
           (type == 'acfbr' & round == '04') | (type == 'acfbr' & round == '09') | (type == 'acfbr' & round == '16') | 
           (type == 'acfcr' & round == '02') | (type == 'acfcr' & round == '03') | (type == 'acfcr' & round == '05')) %>% 
  mutate(val = case_when(var == 'cumcIntv' ~ as.numeric(val) / as.numeric(round), TRUE ~ as.numeric(val)),
         lo = case_when(var == 'cumcIntv' ~ as.numeric(lo) / as.numeric(round), TRUE ~ as.numeric(lo)),
         hi = case_when(var == 'cumcIntv' ~ as.numeric(hi) / as.numeric(round), TRUE ~ as.numeric(hi))) %>%
  mutate(val = case_when(var == 'dfcumcPCF' ~ as.numeric(val) / 25, TRUE ~ as.numeric(val)),
         lo = case_when(var == 'dfcumcPCF' ~ as.numeric(lo) / 25, TRUE ~ as.numeric(lo)),
         hi = case_when(var == 'dfcumcPCF' ~ as.numeric(hi) / 25, TRUE ~ as.numeric(hi))) %>%
  mutate(val = format(val, big.mark = ","), lo = format(lo, big.mark = ","), hi = format(hi, big.mark = ","))

# Results
# Main table
filter(outfin, var == 'cumInc') # Cumulative TB incidence 2050
filter(outfin, var == 'cumMor') # Cumulative TB deaths 2050
filter(outfin, var == 'cumDALYs') # Cumulative DALYs 2050

filter(outfin, var == 'cumTP') # Cumulative TP
filter(outfin, var == 'cumFP') # Cumulative FP

filter(outfin, var == 'cumcCF') # Cumulative case-finding costs 2050
filter(outfin, var == 'cumcRx') # Cumulative treatment costs 2050
filter(outfin, var == 'cumcAll') # Cumulative costs 2050

filter(outfin, var == 'cumcIntv') # Cumulative ACF costs 2050
filter(outfin, var == 'dfcumcPCF') # Cumulative costs 2050
filter(outfin, var == 'icer') # ICER

# Main text
filter(outfin, var == 'cumMIN') # Cumulative minimal
filter(outfin, var == 'cumTPinf') # Cumulative infectious TP
filter(outfin, var == 'cumrtFP') # Ratio TP:FP

filter(outfin, var == 'prpcACF') # Proportion costs ACF
filter(outfin, var == 'prpcBAU') # Proportion costs BAU

# CUA tables
filter(outfin, var == 'dfcumDALYs')
filter(outfin, var == 'dfcumcAll')

rm(list = ls())