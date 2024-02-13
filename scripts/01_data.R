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
data <- as.data.table(import(here("data","tbs","data.csv")))
WPP <- import(here("data","pop","WPP_Pop_1950-2100.csv")) # Population size 1950-2100
WPPb <- import(here("data","pop","WPP_Births_1950-2100.csv")) # Births 1950-2100
WPPda <- import(here("data","pop","WPP_Deaths_1950-2021.csv")) # Deaths 1950-2021
WPPdb <- import(here("data","pop","WPP_Deaths_2022-2100.csv")) # Deaths 2022-2100
WUP <- import(here("data","pop","WUP_Urban_1950-2050.csv")) # Proportion Urban 1950-2050
WDI <- import(here("data","pop","WB_WDI.csv")) # GDP World Development Indicators
WEO <- import(here("data","pop","WEO.csv")) # World Economic Outputs

# 2. TB prevalence survey data ==========
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

# Missingness of cough duration per stratum and SES status
table(VN$cough2w_surv, VN$stratum, useNA = 'always')

gg_miss_var(VN, show_pct = TRUE, facet = stratum) + ylim(0,100)
gg_miss_var(VN, show_pct = TRUE, facet = ses) + ylim(0,100)

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

# 3. Population and demographics ==========
WPP <- WPP %>% # Population data
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrp, pop = PopTotal) %>%
  select(iso3, year, agegp, pop) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020 & year <= 2050) %>% 
  filter(!agegp %in% c('0-4','5-9','10-14')) %>%
  group_by(iso3, year) %>% 
  summarise(pop = sum(pop)*1e3)

WPPb <- WPPb %>% # Birth data
  rename(iso3 = ISO3_code, year = Time, births = Births) %>%
  select(iso3, year, births) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020 & year <= 2050) %>%
  mutate(births = births*1e3)

WPPda <- WPPda %>% # Mortality data (1950-2021)
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrpStart, mort = DeathTotal) %>%
  select(iso3, year, agegp, mort) %>%
  filter(iso3 == "VNM") %>%
  filter(year >= 2020) %>% 
  filter(agegp >= 15) %>%
  group_by(iso3, year) %>%
  summarise(mort = sum(mort)*1e3)

WPPdb <- WPPdb %>% # Mortality data (2022-2100)
  rename(iso3 = ISO3_code, year = Time, agegp = AgeGrpStart, mort = DeathTotal) %>%
  select(iso3, year, agegp, mort) %>%
  filter(iso3 == "VNM") %>%
  filter(year <= 2050) %>% 
  filter(agegp >= 15) %>%
  group_by(iso3, year) %>%
  summarise(mort = sum(mort)*1e3)

WPPd <- rbind(WPPda,WPPdb)
rm(WPPda,WPPdb)

WPP <- WPP %>%
  left_join(WPPb, by=c("iso3","year")) %>%
  inner_join(WPPd, by=c("iso3","year")) %>%
  mutate(mortrate = mort/pop, birthrate = births/pop) 
rm(WPPb,WPPd)

export(WPP,here("data","ses","WPP.Rdata")) # Save data frame

WUP <- clean_names(WUP) %>% # Urbanisation data
  filter(index == 115) %>% 
  select(starts_with("x")) %>%
  rename_all(~gsub("^x", "", .)) %>% 
  mutate_all(~ . / 100) %>% 
  pivot_longer(cols = everything(), names_to = "year", values_to = "urbprop") %>% 
  filter(year >= 2020) %>% 
  mutate(rurprop = 1-urbprop, iso3 = "VNM") %>% 
  select(iso3, year, urbprop, rurprop)

export(WUP,here("data","ses","WUP.Rdata")) # Save data frame

WEO <- WEO %>% # World Economic Output (GDP)
  setNames(WEO[1,]) %>%
  slice(2:n()) %>% 
  clean_names() %>% 
  rename_all(~gsub("^x", "", .)) %>% 
  rename(iso3 = iso, var = subject_descriptor) %>% 
  select(iso3, var, matches("^\\d")) %>% 
  mutate_at(vars(-iso3, -var), ~as.numeric(gsub(",", "", gsub("\\.", "", ., fixed = TRUE)))) %>% 
  filter(var == "Gross domestic product per capita, constant prices") %>% 
  select(-var) %>% 
  pivot_longer(cols = -iso3, names_to = "year", values_to = "gdp") %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(lowses = NA) %>% 
  mutate(lowses = case_when(year == 2008 ~ 0.5946, year == 2018 ~ 0.6034, TRUE ~ lowses))

gdpses <- WEO %>% 
  filter(!is.na(lowses)) %>% 
  lm(lowses ~ gdp, data = .)

WEO <- WEO %>% 
  mutate(lowses = predict(gdpses, newdata = .))

sesmodel <- lm(lowses ~ year, data = WEO)

WEO_exp <- data.frame(year = 2029:2050) %>% 
  mutate(lowses = predict(sesmodel, newdata = .))

GDP <- rbind(select(WEO, c(year, lowses)), WEO_exp) %>% 
  filter(year >= 2020)

export(GDP, here("data","ses","GDP.Rdata")) # Save data frame

