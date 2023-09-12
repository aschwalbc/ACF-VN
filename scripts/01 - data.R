## Analysis code for ACF-VN
## Distributed under CC BY 4.0
## RScript 01: Data.R

# Packages ==========
library(data.table) # Faster than data.frame, allows use of j operator (:=)
library(rio) # Facilitates importing and exporting
library(here) # Building file paths
library(tidyverse) # To use tidyverse
library(naniar) # Visualise missingness

# 1. Load data ==========
data <- as.data.table(import(here("data","data.csv")))

# 2. Data curation ==========
VN <- data %>% 
  select(survey_id, NewID_survey, sex_final, age_final, agegroup, stratum_id, quartile_Nicola_All, tbcase, starts_with("cough"),
         symptoms_any_NF, fever, weight_loss, night_sweats, sputum_c, blood_sputum_c) %>% 
  rename(id = survey_id, survey_n = NewID_survey, sex = sex_final, age = age_final, agegp = agegroup, stratum = stratum_id, ses = quartile_Nicola_All,
         cough_surv = cough_clean, cough2w_surv = coughgt2w, cough_nf = cough_cleannf, cough2w_nf = coughgt2wnf, sympt_nf = symptoms_any_NF) %>% 
  mutate(survey_n = case_when(survey_n == 0 ~ 'surv1', survey_n == 1 ~ 'surv2'), # Surveys: surv1 = 2006-2007 and surv2 = 2017-2018
         stratum = case_when(stratum %in% c('Remote','Rural') ~ 'rural', stratum == "Urban" ~ 'urban'),
         ses = case_when(ses %in% c(1,2) ~ 'low', ses %in% c(3,4) ~ 'high'),
         tbcase = case_when(tbcase == 0 ~ 'no', tbcase == 1 ~ 'yes'),
         cough_nf = case_when(cough_nf == 0 ~ 'no', cough_nf == 1 ~ 'yes', cough_nf %in% c(2,9) ~ NA), # Cough: Khong (0) = no, Co (1) = yes, Khongco TT (9) = NA, 2 (2) = NA
         cough2w_nf = case_when(cough2w_nf == 1 ~ 'more than 2w', cough2w_nf == 2 ~ 'less than 2w'),
         cough_surv = case_when(cough_surv == 'Khong' ~ 'no', cough_surv == 'Co' ~ 'yes'),
         cough2w_surv = case_when(cough2w_surv == 'Tren 2 tuan' ~ 'more than 2w', cough2w_surv == 'Duoi 2 tuan' ~ 'less than 2w'),
         sympt_nf = case_when(sympt_nf == 0 ~ 'no', sympt_nf == 1 ~ 'yes'),
         sympt_surv = case_when(cough == 1 |fever == 1 | weight_loss == 1 | night_sweats == 1 | sputum_c == 1 | blood_sputum_c == 1 ~ 'yes', TRUE ~ 'no')) %>% 
  #mutate(cough2w_surv = ifelse(cough_surv == 'yes' & is.na(cough2w_surv), 'less than 2w', cough2w_surv)) %>% # fix missingness of cough duration
  select(id, survey_n, sex, age, agegp, stratum, ses, tbcase, cough_surv, cough_nf, cough2w_surv, cough2w_nf, sympt_surv, sympt_nf)

# 3. Data exploration ==========
# Missingness of cough duration per stratum and SES status
table(VN$cough2w_surv, VN$stratum, useNA = 'always')

gg_miss_var(VN, show_pct = TRUE, facet = stratum) + ylim(0,100)
gg_miss_var(VN, show_pct = TRUE, facet = ses) + ylim(0,100)

# 4. Results ==========
N = 100000 # Per 100,000 inhabitants

# Survey 1: 2007-2008
VNs1 <- (filter(VN, survey_n == 'surv1'))
popsurv1 <- nrow(VNs1)

nrow(filter(VNs1, tbcase == 'yes')) # Total TB cases
nrow(filter(VNs1, tbcase == 'yes'))/popsurv1*N # # Total TB cases (per 100k)

# Cough duration (>2w vs <2w)
nrow(filter(VNs1, tbcase == 'yes' & cough2w_surv == 'more than 2w')) #cTB
nrow(filter(VNs1, tbcase == 'yes' & cough2w_surv == 'more than 2w'))/popsurv1*N # cTB (per 100k)

nrow(filter(VNs1, tbcase == 'yes' & cough2w_surv == 'less than 2w')) # scTB
nrow(filter(VNs1, tbcase == 'yes' & cough2w_surv == 'less than 2w'))/popsurv1*N # scTB (per 100k)

# Presence of cough ***
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes')) #cTB
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes'))/popsurv1*N # cTB (per 100k)

nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no')) # scTB
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no'))/popsurv1*N # scTB (per 100k)

# Presence of symptoms
nrow(filter(VNs1, tbcase == 'yes' & sympt_surv == 'yes')) #cTB
nrow(filter(VNs1, tbcase == 'yes' & sympt_surv == 'yes'))/popsurv1*N # cTB (per 100k)

nrow(filter(VNs1, tbcase == 'yes' & sympt_surv == 'no')) # scTB
nrow(filter(VNs1, tbcase == 'yes' & sympt_surv == 'no'))/popsurv1*N # scTB (per 100k)

nrow(filter(VNs1, tbcase == 'yes' & sympt_nf == 'yes')) #cTB
nrow(filter(VNs1, tbcase == 'yes' & sympt_nf == 'yes'))/popsurv1*N # cTB (per 100k)

nrow(filter(VNs1, tbcase == 'yes' & sympt_nf == 'no')) # scTB
nrow(filter(VNs1, tbcase == 'yes' & sympt_nf == 'no'))/popsurv1*N # scTB (per 100k)

# Relative stratum (urban/rural)
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no' & stratum == 'urban'))/nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no' & (stratum == 'urban' | stratum == 'rural')))
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes' & stratum == 'urban'))/nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes' & (stratum == 'urban' | stratum == 'rural')))

# Relative stratum (high/low)
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no' & ses == 'high'))/nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'no' & (ses == 'high' | ses == 'low')))
nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes' & ses == 'high'))/nrow(filter(VNs1, tbcase == 'yes' & cough_surv == 'yes' & (ses == 'high' | ses == 'low')))

# Survey 2: 2017-2018
VNs2 <- (filter(VN, survey_n == 'surv2'))
popsurv2 <- nrow(VNs2)

nrow(filter(VNs2, tbcase == 'yes')) # Total TB cases
nrow(filter(VNs2, tbcase == 'yes'))/popsurv1*N # # Total TB cases (per 100k)

# Cough duration (>2w vs <2w)
nrow(filter(VNs2, tbcase == 'yes' & cough2w_surv == 'more than 2w')) #cTB
nrow(filter(VNs2, tbcase == 'yes' & cough2w_surv == 'more than 2w'))/popsurv2*N # cTB (per 100k)

nrow(filter(VNs2, tbcase == 'yes' & cough2w_surv == 'less than 2w')) # scTB
nrow(filter(VNs2, tbcase == 'yes' & cough2w_surv == 'less than 2w'))/popsurv2*N # scTB (per 100k)

# Presence of cough ***
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes')) #cTB
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes'))/popsurv2*N # cTB (per 100k)

nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no')) # scTB
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no'))/popsurv2*N # scTB (per 100k)

# Presence of symptoms
nrow(filter(VNs2, tbcase == 'yes' & sympt_surv == 'yes')) #cTB
nrow(filter(VNs2, tbcase == 'yes' & sympt_surv == 'yes'))/popsurv2*N # cTB (per 100k)

nrow(filter(VNs2, tbcase == 'yes' & sympt_surv == 'no')) # scTB
nrow(filter(VNs2, tbcase == 'yes' & sympt_surv == 'no'))/popsurv2*N # scTB (per 100k)

nrow(filter(VNs2, tbcase == 'yes' & sympt_nf == 'yes')) #cTB
nrow(filter(VNs2, tbcase == 'yes' & sympt_nf == 'yes'))/popsurv2*N # cTB (per 100k)

nrow(filter(VNs2, tbcase == 'yes' & sympt_nf == 'no')) # scTB
nrow(filter(VNs2, tbcase == 'yes' & sympt_nf == 'no'))/popsurv2*N # scTB (per 100k)

# Relative stratum (urban/rural)
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no' & stratum == 'urban'))/nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no' & (stratum == 'urban' | stratum == 'rural')))
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes' & stratum == 'urban'))/nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes' & (stratum == 'urban' | stratum == 'rural')))

# Relative stratum (high/low)
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no' & ses == 'high'))/nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'no' & (ses == 'high' | ses == 'low')))
nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes' & ses == 'high'))/nrow(filter(VNs2, tbcase == 'yes' & cough_surv == 'yes' & (ses == 'high' | ses == 'low')))

