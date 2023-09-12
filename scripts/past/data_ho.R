

data = as.data.frame(read.csv("~/Desktop/data.csv"))

data[is.na(data$coughgt2wnf),"coughgt2wnf"] = 0


N = 100000

dim(subset(data, survey == "2nd survey" & tbcase == 1))[1]
dim(subset(data, survey == "2nd survey"))[1]

dim(subset(data, survey == "2nd survey" & tbcase == 1 & cough == 0))[1]/dim(subset(data, survey == "2nd survey"))[1]

N*dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1]/dim(subset(data, survey == "1st survey"))[1]
N*dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1]/dim(subset(data, survey == "2nd survey"))[1]

N*dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1))[1]/dim(subset(data, survey == "1st survey"))[1]
N*dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1))[1]/dim(subset(data, survey == "2nd survey"))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & stratum_id == "Urban"))[1]/dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & (stratum_id == "Rural" | stratum_id == "Remote")))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & stratum_id == "Urban"))[1]/dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & (stratum_id == "Rural" | stratum_id == "Remote")))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1 & stratum_id == "Urban"))[1]/dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1 & (stratum_id == "Rural" | stratum_id == "Remote")))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1 & stratum_id == "Urban"))[1]/dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1 & (stratum_id == "Rural" | stratum_id == "Remote")))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & quartile_Nicola_All %in% c(3,4)))[1]/dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & quartile_Nicola_All %in% c(1,2)))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & quartile_Nicola_All %in% c(3,4)))[1]/dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2) & quartile_Nicola_All %in% c(1,2)))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1 & quartile_Nicola_All %in% c(3,4)))[1]/dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1 & quartile_Nicola_All %in% c(1,2)))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1 & quartile_Nicola_All %in% c(3,4)))[1]/dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1 & quartile_Nicola_All %in% c(1,2)))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1]
dim(subset(data, survey == "1st survey" & tbcase == 1 & cough == 0))[1]
dim(subset(data, survey == "1st survey"))[1]

dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & cough == 0))[1]
dim(subset(data, survey == "2nd survey"))[1]

dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1))[1]
dim(subset(data, survey == "1st survey" & tbcase == 1 & cough == 1))[1]
dim(subset(data, survey == "1st survey"))[1]

dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1))[1]
dim(subset(data, survey == "2nd survey" & tbcase == 1 & cough == 1))[1]
dim(subset(data, survey == "2nd survey"))[1]

# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough_cleannf == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 0))[1]
# 
# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough_c == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough_2w == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & sputum_c == 0))[1]
# 
# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough == 0))[1]
# 
# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough_clean == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2w == 0))[1]
# dim(subset(data, survey == "1st survey" & tbcase == 1 & cough_CAT == 0))[1]
# 
# subset(data, survey == "1st survey" & tbcase == 1)$coughgt2w
# subset(data, survey == "1st survey" & tbcase == 1)$coughgt2wnf
# 
# clin1 = dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf == 1))[1] #clinical
# sub1 = dim(subset(data, survey == "1st survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1] #subclinical
# 
# clin2 = dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf == 1))[1] #clinical
# sub2 = dim(subset(data, survey == "2nd survey" & tbcase == 1 & coughgt2wnf %in% c(0,2)))[1] #subclinical



















