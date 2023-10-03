dN_RL  = nu(times)*(N_RL+SN_RL+I_RL+SI_RL+O_RL+SO_RL+M_RL+RM_RL+SM_RL+S_RL+RS_RL+SS_RL+C_RL+RC_RL+SC_RL+P_RL+SP_RL) - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*N_RL - acf(floor(times))*alpha_nds*N_RL - mu(times)*N_RL
dSN_RL = acf(floor(times))*alpha_nds*N_RL - delta*SN_RL - mu(times)*SN_RL
dI_RL  = ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*(N_RL+(O_RL*theta_cleinf)+(P_RL*theta_recinf)) - gamma_infcle*I_RL - lambda_infmin*I_RL - lambda_infsub*I_RL - acf(floor(times))*alpha_nds*I_RL - mu(times)*I_RL
dSI_RL = acf(floor(times))*alpha_nds*I_RL - delta*SI_RL - mu(times)*SI_RL
dO_RL  = gamma_infcle*I_RL + gamma_mincle*M_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*O_RL*theta_cleinf - acf(floor(times))*alpha_nds*O_RL + delta(SN_RL+SI_RL+SO_RL) - mu(times)*O_RL
dSO_RL = acf(floor(times))*alpha_nds*O_RL - delta*SO_RL - mu(times)*SO_RL
dM_RL  = lambda_infmin*I_RL + gamma_submin*S_RL - gamma_mincle*M_RL - lambda_minsub*M_RL - iota_min*M_RL + phi_min*RM_RL - acf(floor(times))*alpha_min*M_RL + phi_min*SM_RL + tau_min*P_RL - mu(times)*M_RL
dRM_RL = iota_min*M_RL - phi_min*RM_RL - delta*RM_RL - mu(times)*RM_RL
dSM_RL = acf(floor(times))*alpha_min*M_RL - phi_min*SM_RL - delta*SM_RL - mu(times)*SM_RL
dS_RL  = lambda_infsub*I_RL + lambda_minsub*M_RL + gamma_clnsub*C_RL - gamma_submin*S_RL - lambda_subcln*S_RL - iota_sub*S_RL + phi_sub*RS_RL - acf(floor(times))*alpha_sub*S_RL + phi_sub*SS_RL + tau_sub*P_RL - mu(times)*S_RL
dRS_RL = iota_sub*S_RL - phi_sub*RS_RL - delta*RS_RL - mu(times)*RS_RL
dSS_RL = acf(floor(times))*alpha_sub*S_RL - phi_sub*SS_RL - delta*SS_RL - mu(times)*SS_RL
dC_RL  = lambda_subcln*S_RL - gamma_clnsub*C_RL - iota_cln*C_RL + phi_cln*RC_RL - acf(floor(times))*alpha_cln*C_RL + phi_cln*SC_RL - omega*C_RL - mu(times)*C_RL
dRC_RL = iota_cln*C_RL - phi_cln*RC_RL - delta*RC_RL - mu(times)*RC_RL
dSC_RL = acf(floor(times))*alpha_cln*C_RL - phi_cln*SC_RL - delta*SC_RL - mu(times)*SC_RL
dP_RL  = delta*(RM_RL+SM_RL+RS_RL+SS_RL+RC_RL+SC_RL+SP_RL) - tau_min*P_RL - tau_sub*P_RL - ((beta/chi(times))*((kappa*(S_RL+S_RH+S_UL+S_UH))+(C_RL+C_RH+C_UL+C_UH)))*P_RL*theta_recinf - mu(times)*P_RL
dSP_RL = acf(floor(times))*alpha_nds*P_RL - delta*SP_RL - mu(times)*SP_RL

WDI <- clean_names(WDI) %>% # World Development Indicators (GDP)
  slice(1:(n()-5)) %>% 
  select(-country_name, - series_code) %>% 
  rename(iso3 = country_code, series = series_name) %>% 
  rename_with(~gsub("x(\\d+)_yr(\\d+)", "\\2", .x), starts_with("x")) %>% 
  mutate_if(is.numeric, as.character) %>% 
  mutate_all(~ifelse(str_detect(., "\\.{2}"), "", .)) %>% 
  mutate_at(vars(-iso3, -series), as.numeric) %>% 
  filter(series == "GDP growth (annual %)") %>% 
  pivot_longer(cols = c(-iso3, -series), names_to = "year", values_to = "gdp") %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year >= 1985)

gdpmodel <- lm(gdp ~ year, data = WDI)
pred_yr <- data.frame(year = 2023:2050)
pred_gdp <- predict(gdpmodel, newdata = pred_yr)
pred_data <- data.frame(year = 2023:2050, gdp = pred_gdp)
GDP <- rbind(select(WDI, c(year, gdp)), pred_data)
rownames(GDP) <- NULL