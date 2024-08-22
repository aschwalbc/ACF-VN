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
         cACFFP = cACFSUS + cACFINF + cACFCLE + cACFREC + cACFTRE, # Cost of ACF FP
         cACFTP = cACFMIN + cACFSUB + cACFCLN, # Cost of ACF TP
         cACFTPinf = cACFSUB + cACFCLN, # Cost of ACF infectious TP
         cRxBAU = cRxBAUs + cRxBAUr, # Cost of treatment BAU (DSTB + DRTB)
         cRxACFFP = cRxSUS + cRxINF + cRxCLE + cRxREC + cRxTRE, # Cost of treatment ACF FP 
         cRxACFTP = cRxMIN + cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr, # Cost of treatment ACF TP (DSTB + DRTB)
         cRxACFTPinf = cRxSUBs + cRxSUBr + cRxCLNs + cRxCLNr) %>%  # Cost of treatment ACF infectious TP (DSTB + DRTB)
  mutate(cCF = cBAU + cACF, # Cost of case finding (PCF + ACF)
         cRxACF = cRxACFFP + cRxACFTP, # Cost of treatment ACF
         cRx = cRxBAU + cRxACFFP + cRxACFTP, # Cost of treatment (PCF + ACF)
         cAll = cBAU + cACF + cRxBAU + cRxACFFP + cRxACFTP) %>% # All costs
  mutate(goal = case_when(type == 'acfa' & round == '03' ~ '100', type == 'acfa' & round == '06' ~ '50', type == 'acfa' & round == '11' ~ '20',
                          type == 'acfax' & round == '03' ~ '100', type == 'acfax' & round == '06' ~ '50', type == 'acfax' & round == '11' ~ '20',
                          type == 'acfb' & round == '03' ~ '100', type == 'acfb' & round == '07' ~ '50', type == 'acfb' & round == '12' ~ '20',
                          type == 'acfbx' & round == '03' ~ '100', type == 'acfbx' & round == '07' ~ '50', type == 'acfbx' & round == '12' ~ '20',
                          type == 'acfc' & round == '01' ~ '100', type == 'acfc' & round == '02' ~ '50', type == 'acfc' & round == '03' ~ '20',
                          type == 'base' ~ 'none')) %>% 
  mutate(goal = factor(goal, levels = c('100', '50', '20', 'none'))) %>% 
  mutate(res = case_when(type == 'acfa' ~ 'main', type == 'acfax' ~ 'sens',
                         type == 'acfb' ~ 'main', type == 'acfbx' ~ 'sens',
                         type == 'acfc' ~ 'main', type == 'base' ~ 'main'))
  
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
  pivot_longer(cols = -c(time, type, res, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, res, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE),
            lo = quantile(values, 0.025, na.rm = TRUE),
            hi = quantile(values, 0.975, na.rm = TRUE))

outs_yr <- outs %>% 
  arrange(type, run, time) %>% 
  filter(time == floor(time)) %>% 
  filter(time >= 2025) %>% 
  group_by(type, run, round) %>%
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
         cumcACFSUS = cumsum(cACFSUS), # Cumulative cost of ACF SUS
         cumcACFINF = cumsum(cACFINF), # Cumulative cost of ACF INF
         cumcACFCLE = cumsum(cACFCLE), # Cumulative cost of ACF CLE
         cumcACFREC = cumsum(cACFREC), # Cumulative cost of ACF REC
         cumcACFMIN = cumsum(cACFMIN), # Cumulative cost of ACF MIN
         cumcACFSUB = cumsum(cACFSUB), # Cumulative cost of ACF SUB
         cumcACFCLN = cumsum(cACFCLN), # Cumulative cost of ACF CLN
         cumcACFTRE = cumsum(cACFTRE), # Cumulative cost of ACF TRE
         cumcACFFP = cumsum(cACFFP), # Cumulative cost of ACF FP
         cumcACFTP = cumsum(cACFTP), # Cumulative cost of ACF TP
         cumcACFTPinf = cumsum(cACFTPinf), # Cumulative cost of ACF infectious TP
         cumcACF = cumsum(cACF), # Cumulative cost of ACF
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
         cumcCF = cumsum(cCF), # Cumulative cost of case finding
         cumcRx = cumsum(cRx), # Cumulative cost of all treatments
         cumDALYs = cumsum(DALYs), # Cumulative DALYs
         cumcAll = cumsum(cAll)) %>% # Cumulative costs of all 
  ungroup() %>% 
  group_by(run, time) %>% 
  mutate(dfcumInc = cumInc - cumInc[type == 'base'], # TB incidence averted
         dfcumMor = cumMor - cumMor[type == 'base'], # TB mortality averted
         dfcumcCF = cumcCF - cumcCF[type == 'base'], # CF cost averted
         dfcumcRx = cumcRx - cumcRx[type == 'base'], # Rx cost averted
         dfcumcAll = cumcAll - cumcAll[type == 'base'], # All costs averted
         dfcumDALYs = cumDALYs - cumDALYs[type == 'base'], # DALYs averted
         prpcACF = cumcACF / (cumcACF + cumcRxACF)) %>% # Proportion costs of ACF screening over total costs
  ungroup() %>% 
  mutate(cumprFP = ifelse(cumFP == 0, NA, cumFP / (cumACF)), # Proportion FP
         cumrtFP = ifelse(cumFP == 0, NA, cumFP / cumTP), # Ratio TP:FP
         cumrtFPScr = ifelse(cumScrn == 0, NA, cumScrn / cumFP), # Ratio FP:Scr
         cumrtTPScr = ifelse(cumScrn == 0, NA, cumScrn / cumTP), # Ratio TP:Scr
         cumrtScr = ifelse(cumScrn == 0, NA, cumScrn / (cumTP + cumFP)), # Ratio Rx:Scr
         cpIncACF = ifelse((cumcACF + cumcRxACF) == 0, NA, (cumcACF + cumcRxACF) / abs(dfcumInc)), # Cost per incident TB averted (ACF costs)
         cpIncCF = ifelse((cumcCF + cumcRx) == 0, NA, (cumcCF + cumcRx) / abs(dfcumInc)), # Cost per incident TB averted (All costs)
         cpMorACF = ifelse((cumcACF + cumcRxACF) == 0, NA, (cumcACF + cumcRxACF) / abs(dfcumMor)), # Cost per TB death averted (ACF costs)
         cpMorCF = ifelse((cumcCF + cumcRx) == 0, NA, (cumcCF + cumcRx) / abs(dfcumMor)), # Cost per TB death averted (All costs)
         cpDALYsACF = ifelse((cumcACF + cumcRxACF) == 0, NA, (cumcACF + cumcRxACF) / abs(dfcumDALYs)), # Cost per DALYs averted (ACF costs) 
         cpDALYsCF = ifelse((cumcCF + cumcRx) == 0, NA, (cumcCF + cumcRx) / abs(dfcumDALYs))) %>% # Cost per DALYs averted (All costs) 
  select(time, type, res, run, round, goal, contains('cum'), contains('prp'), contains('cp')) %>%
  pivot_longer(cols = -c(time, type, res, run, round, goal), names_to = "var", values_to = "values") %>%
  group_by(time, type, res, round, goal, var) %>%
  summarise(val = median(values, na.rm = TRUE), 
            lo = quantile(values, 0.025, na.rm = TRUE), 
            hi = quantile(values, 0.975, na.rm = TRUE))

outs <- rbind(outs_m, outs_yr) %>% 
  arrange(factor(type, levels = c("base", "acfa", "acfb", "acfc", "acfax", "acfbx")), round, time) 
rm(outs_m, outs_yr)

export(outs, here("outputs", "outs", "outs.Rdata"))
rm(outs)

# 3. Extract results ==========
outs <- import(here("outputs", "outs", "outs.Rdata"))

outini <- outs %>% 
  filter((type == 'acfa' & round == '03' & time == 2028) |
           (type == 'acfa' & round == '06' & time == 2031) |
           (type == 'acfa' & round == '11' & time == 2036) |
           (type == 'acfb' & round == '03' & time == 2028) |
           (type == 'acfb' & round == '07' & time == 2032) |
           (type == 'acfb' & round == '12' & time == 2037) |
           (type == 'acfc' & round == '01' & time == 2026) |
           (type == 'acfc' & round == '02' & time == 2027) |
           (type == 'acfc' & round == '03' & time == 2028) |
           (type == 'acfd' & round == '03' & time == 2028) |
           (type == 'acfd' & round == '06' & time == 2031) |
           (type == 'acfd' & round == '11' & time == 2036))

filter(outini, var == 'rTBc') # TB prevalence rate
filter(outini, var == 'redTBc') # TB prevalence rate reduction
filter(outini, var == 'rInc') # TB incidence rate
filter(outini, var == 'rMor') # TB mortality rate

outfin <- outs %>% 
  filter(time == 2050) %>% 
  filter((type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
           (type == 'acfb' & round == '03') | (type == 'acfb' & round == '07') | (type == 'acfb' & round == '12') | 
           (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') |
           (type == 'acfd' & round == '03') | (type == 'acfd' & round == '06') | (type == 'acfd' & round == '11'))

filter(outfin, var == 'dfcumInc') # TB incidence averted 
filter(outfin, var == 'dfcumMor') # TB mortality averted
filter(outfin, var == 'dfcumInf') # Mtb infections averted
filter(outfin, var == 'cumDALYs') # Total DALYs
filter(outfin, var == 'cumcAll') # All costs
filter(outfin, var == 'cumcACF') # Total costs of ACF screening
filter(outfin, var == 'cumcCF') # Total costs of case finding
filter(outfin, var == 'cumTP') # Total true positives in ACF
filter(outfin, var == 'cumFP') # Total true positives in ACF
filter(outfin, var == 'cumprFP') # Proportion FP treated
filter(outfin, var == 'cumcRxACF') # Total costs of ACF treatment
filter(outfin, var == 'cumcRx') # Total costs of treatment
filter(outfin, var == 'prpcACF') # Proportion of costs for ACF over all costs
filter(outfin, var == 'prpcCF') # Proportion of costs for CF over all costs
filter(outfin, var == 'dfcumDALYs') # DALYs averted

icerDALY <- outs %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
          (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
          (type == 'acfb' & round == '03') | (type == 'acfb' & round == '07') | (type == 'acfb' & round == '12') | 
          (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') |
          (type == 'acfd' & round == '03') | (type == 'acfd' & round == '06') | (type == 'acfd' & round == '11'))) %>%
  filter(var == 'dfcumDALYs' | var == 'cumcCF' | var == 'cumcRx') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  mutate(val_COST = val_cumcCF + val_cumcRx, lo_COST = lo_cumcCF + lo_cumcRx, hi_COST = hi_cumcCF + hi_cumcRx) %>% 
  rename(val_DALY = val_dfcumDALYs, lo_DALY = lo_dfcumDALYs, hi_DALY = hi_dfcumDALYs) %>% 
  select(time, type, round, goal, val_COST, lo_COST, hi_COST, val_DALY, lo_DALY, hi_DALY) %>% 
  mutate(cost_val = (val_COST - val_COST[type == 'base']), 
         cost_lo = (lo_COST - lo_COST[type == 'base']),
         cost_hi = (hi_COST - hi_COST[type == 'base']),
         daly_val = (abs(val_DALY) - abs(val_DALY[type == 'base'])),
         daly_lo = (abs(lo_DALY) - lo_DALY[type == 'base']),
         daly_hi = (abs(hi_DALY) - hi_DALY[type == 'base'])) %>% 
  select(time, type, round, goal, starts_with('cost'), starts_with('daly')) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = cost_val / daly_val, icer_lo = cost_lo / daly_lo, icer_hi = cost_hi / daly_hi) %>% 
  select(time, type, round, goal, starts_with('icer')) %>% 
  mutate(var = 'icerDALY') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(time, type, round, goal, var, val, lo, hi)

icerInc <- outs %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
          (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
          (type == 'acfb' & round == '03') | (type == 'acfb' & round == '07') | (type == 'acfb' & round == '12') | 
          (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') |
          (type == 'acfd' & round == '03') | (type == 'acfd' & round == '06') | (type == 'acfd' & round == '11'))) %>%
  filter(var == 'dfcumInc' | var == 'cumcCF' | var == 'cumcRx') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  mutate(val_COST = val_cumcCF + val_cumcRx, lo_COST = lo_cumcCF + lo_cumcRx, hi_COST = hi_cumcCF + hi_cumcRx) %>% 
  rename(val_Inc = val_dfcumInc, lo_Inc = lo_dfcumInc, hi_Inc = hi_dfcumInc) %>% 
  select(time, type, round, goal, val_COST, lo_COST, hi_COST, val_Inc, lo_Inc, hi_Inc) %>% 
  mutate(cost_val = (val_COST - val_COST[type == 'base']), 
         cost_lo = (lo_COST - lo_COST[type == 'base']),
         cost_hi = (hi_COST - hi_COST[type == 'base']),
         inc_val = (abs(val_Inc) - abs(val_Inc[type == 'base'])),
         inc_lo = (abs(lo_Inc) - lo_Inc[type == 'base']),
         inc_hi = (abs(hi_Inc) - hi_Inc[type == 'base'])) %>% 
  select(time, type, round, goal, starts_with('cost'), starts_with('inc')) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = cost_val / inc_val, icer_lo = cost_lo / inc_lo, icer_hi = cost_hi / inc_hi) %>% 
  select(time, type, round, goal, starts_with('icer')) %>% 
  mutate(var = 'icerInc') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(time, type, round, goal, var, val, lo, hi)

icerMor <- outs %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
                           (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
                           (type == 'acfb' & round == '03') | (type == 'acfb' & round == '07') | (type == 'acfb' & round == '12') | 
                           (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') |
                           (type == 'acfd' & round == '03') | (type == 'acfd' & round == '06') | (type == 'acfd' & round == '11'))) %>%
  filter(var == 'dfcumMor' | var == 'cumcCF' | var == 'cumcRx') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  mutate(val_COST = val_cumcCF + val_cumcRx, lo_COST = lo_cumcCF + lo_cumcRx, hi_COST = hi_cumcCF + hi_cumcRx) %>% 
  rename(val_Mor = val_dfcumMor, lo_Mor = lo_dfcumMor, hi_Mor = hi_dfcumMor) %>% 
  select(time, type, round, goal, val_COST, lo_COST, hi_COST, val_Mor, lo_Mor, hi_Mor) %>% 
  mutate(cost_val = (val_COST - val_COST[type == 'base']), 
         cost_lo = (lo_COST - lo_COST[type == 'base']),
         cost_hi = (hi_COST - hi_COST[type == 'base']),
         mor_val = (abs(val_Mor) - abs(val_Mor[type == 'base'])),
         mor_lo = (abs(lo_Mor) - lo_Mor[type == 'base']),
         mor_hi = (abs(hi_Mor) - hi_Mor[type == 'base'])) %>% 
  select(time, type, round, goal, starts_with('cost'), starts_with('mor')) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = cost_val / mor_val, icer_lo = cost_lo / mor_lo, icer_hi = cost_hi / mor_hi) %>% 
  select(time, type, round, goal, starts_with('icer')) %>% 
  mutate(var = 'icerMor') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(time, type, round, goal, var, val, lo, hi)

icerInf <- outs %>%
  filter(time == 2050 & ((type == 'base' & round == '00') |
                           (type == 'acfa' & round == '03') | (type == 'acfa' & round == '06') | (type == 'acfa' & round == '11') |
                           (type == 'acfb' & round == '03') | (type == 'acfb' & round == '07') | (type == 'acfb' & round == '12') | 
                           (type == 'acfc' & round == '01') | (type == 'acfc' & round == '02') | (type == 'acfc' & round == '03') |
                           (type == 'acfd' & round == '03') | (type == 'acfd' & round == '06') | (type == 'acfd' & round == '11'))) %>%
  filter(var == 'dfcumInf' | var == 'cumcCF' | var == 'cumcRx') %>% 
  pivot_wider(names_from = var, values_from = c(val, lo, hi)) %>% 
  mutate(val_COST = val_cumcCF + val_cumcRx, lo_COST = lo_cumcCF + lo_cumcRx, hi_COST = hi_cumcCF + hi_cumcRx) %>% 
  rename(val_Inf = val_dfcumInf, lo_Inf = lo_dfcumInf, hi_Inf = hi_dfcumInf) %>% 
  select(time, type, round, goal, val_COST, lo_COST, hi_COST, val_Inf, lo_Inf, hi_Inf) %>% 
  mutate(cost_val = (val_COST - val_COST[type == 'base']), 
         cost_lo = (lo_COST - lo_COST[type == 'base']),
         cost_hi = (hi_COST - hi_COST[type == 'base']),
         inf_val = (abs(val_Inf) - abs(val_Inf[type == 'base'])),
         inf_lo = (abs(lo_Inf) - lo_Inf[type == 'base']),
         inf_hi = (abs(hi_Inf) - hi_Inf[type == 'base'])) %>% 
  select(time, type, round, goal, starts_with('cost'), starts_with('inf')) %>% 
  filter(!type == 'base') %>% 
  mutate(icer_val = cost_val / inf_val, icer_lo = cost_lo / inf_lo, icer_hi = cost_hi / inf_hi) %>% 
  select(time, type, round, goal, starts_with('icer')) %>% 
  mutate(var = 'icerInf') %>% 
  rename(val = icer_val, lo = icer_lo, hi = icer_hi) %>% 
  select(time, type, round, goal, var, val, lo, hi)

outs <- outs %>% 
  rbind(icerDALY, icerInc, icerMor, icerInf)
export(outs, here("outputs", "outs", "outs.Rdata"))
rm(list = ls())
