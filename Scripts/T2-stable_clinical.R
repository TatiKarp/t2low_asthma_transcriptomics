# This script will check clinical characteristics for T2-low and T2-high stable groups 
library(dplyr)
library(readr)
library(edgeR) 
library(tableone)
setwd("/Users/tatiana/Work/RP2/ATLANTIS")
# master.Table <- read_csv("./th_high_th_low/master_table_th2.csv")

#### define T2-high and T2-low phenotype ####
big_master_table <-  read_csv("./atlantis_patient_data.csv") %>%
  dplyr::select (c(PT, ASTHEA, VISIT, SYS_COR, BIO, LABEOSV, FENRES)) %>%
  filter(ASTHEA == "A") %>%
  mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES')) %>%
  filter(include == 'YES') %>%
  mutate(group_th = case_when(
    ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
    ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
    ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
    ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low',
    TRUE ~ 'undeterm')) %>%
  mutate(group_th_num = case_when(
    ((LABEOSV > 0.3) & (FENRES > 25)) ~ 1,
    ((LABEOSV < 0.15) & (FENRES < 25)) ~ -1,
    ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 1,
    ((LABEOSV < 0.15) & (is.na(FENRES))) ~ -1,
    TRUE ~ 0))

big_master_table_wide <- big_master_table %>%
  mutate(PT = as.character(PT)) %>%
  dplyr::select(c(PT, VISIT, group_th)) %>%
  tidyr::pivot_wider (names_from = VISIT, values_from= group_th)

num_big_master_table_wide <- big_master_table %>%
  mutate(PT = as.character(PT)) %>%
  dplyr::select(c(PT, VISIT, group_th_num)) %>%
  tidyr::pivot_wider (names_from = VISIT, values_from = group_th_num) %>%
  mutate(at_least_2 = `VISIT 1` + `VISIT 2` + `VISIT 3`)

big_master_table_wide <- big_master_table_wide %>%
  mutate(stable_low = if_else(`VISIT 1` == "low" & `VISIT 2` == "low" & `VISIT 3` == "low", "YES","NO"),
         stable_high = if_else(`VISIT 1` == "high" & `VISIT 2` == "high" & `VISIT 3` == "high", "YES","NO"))

master.Table <- read.csv("./th_high_th_low/master_table_th2.csv") %>%
  mutate(PT = as.character(PT)) %>%
  left_join(big_master_table_wide %>%
              dplyr::select(c(PT, stable_low, stable_high)), by = c("PT" = "PT")) %>%
  mutate(stable_low = if_else(group_th == "healthy", "healthy", stable_low),
         stable_high = if_else(group_th == "healthy", "healthy", stable_high)) %>%
  left_join(num_big_master_table_wide %>%
              dplyr::select(c(PT, at_least_2)), by = c("PT" = "PT")) %>%
  mutate(at_least_2 = as.character(at_least_2)) %>%
  mutate(at_least_2 = if_else(group_th == "healthy", "healthy", at_least_2)) 

#write.csv(master.Table, './th_high_th_low/master_table_stable_th2.csv', row.names = F)
############### prepare clinical data ###########################
big_master_table <-  read.csv("./atlantis_patient_data.csv", header =TRUE, na.strings=c("","NA"))
additional_data <- read.csv("./Umi_dedup/Dif_expr/Asthma_groups/clinicaldata_bhr_feno_mld.csv", header =TRUE, na.strings=c("","NA"))%>%
  dplyr::select(-c(X,FENRES, MLD_ratio))
clinical_table <- big_master_table[,c('PT','PACKNO','BMI','PHADRES', 'DUR_DIS','AGE_DIAG','B_FEV1F','GINA', 'NUM_EX',"MORE1EX",
                                      'LAMA','BIO','SYS_COR','ICS_MEAN_DOSE','ICS_DDOSE_EQ',
                                      'acq6_score','LABEOSV','LABNEUV','LABMACV','BRONCHP', 'LYMPHOP', 'EOSP', 'MACROP', 'NEUTROP',
                                      'B_TLCPNVF', 'B_RVTLCPNVF', 'B_FEV1PNVG','B_FEV1FPNVG', 'FENRES', 'PCD',
                                      'B_R520PNVR', 'B_SCONDPNVF','B_SACINPNVF', 'B_F50PNVR', 'LA', 'WA','TA',
                                      'Pi10', 'WA_TA100', 'VI_856', 'VI_950', 'lung_ratio', 'MLD_ratio',
                                      'ICS', 'ICS_LABA', 'ICS_DDOSE_EQ', 'ICS_LABA_DDOSE_EQ')]
clinical_table <- clinical_table[!duplicated(clinical_table$PT), ] %>%
  mutate (GINA = as.factor(GINA),
          NUM_EX = as.factor(NUM_EX))

master.Table <- master.Table %>%
  filter(QC_check == 'YES') %>%
  mutate(gender = as.factor(gender),
         smoking.status = as.factor(smoking.status),
         group_th = as.factor(group_th)) %>%
  left_join(clinical_table %>%
              mutate(PT = as.character(PT)), by = c('PT'='PT')) %>%
  left_join(additional_data %>%
              mutate(PT = as.character(PT)), by = c('PT'='PT')) %>%
  mutate(any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes")) %>%
  mutate(ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                                coalesce(ICS_DDOSE_EQ,0) + coalesce(ICS_LABA_DDOSE_EQ,0)))

## ICS_dose_sum = ICS_DDOSE_EQ + ICS_LABA_DDOSE_EQ
## any_ICS is one of the ICS or ICS_LABA == Yes

non_normally <- c('DUR_DIS', 'AGE_DIAG','PACKNO','acq6_score', 'LABEOSV', 'BRONCHP', 'LYMPHOP',
                  'MACROP', 'NEUTROP', 'B_R520','VI_856', 'VI_950','FENRES', 'EOSP','B_R520PNVR', 
                  'B_SCONDPNVF','B_SACINPNVF', "NUM_EX", 'ICS_dose_sum')

var_for_table <- c('gender','age','smoking.status','PACKNO','BMI','PHADRES', 'DUR_DIS','AGE_DIAG','B_FEV1F','GINA', 'NUM_EX','LAMA','BIO',
                   'SYS_COR', 'any_ICS', 'ICS_dose_sum',
                   'acq6_score','LABEOSV','LABNEUV','LABMACV','BRONCHP', 'LYMPHOP', 'EOSP', 'MACROP', 'NEUTROP',
                   'B_TLCPNVF', 'B_FEV1PNVG','B_FEV1FPNVG', 'FENRES', 'PCD','bhr',
                   'B_RVTLCPNVF',  'B_R520PNVR', 'B_SCONDPNVF','B_SACINPNVF', 'B_F50PNVR',
                   'VI_856', 'VI_950', 'lung_ratio','MLD_ratio', "MORE1EX")

## check NAs 
high_not_na <- master.Table %>% 
  group_by(stable_high) %>%
  summarise(across(everything(), ~sum(!is.na(.)), .names = "not_NA_{.col}")) %>%
  t()

low_not_na <- master.Table %>% 
  group_by(stable_low) %>%
  summarise(across(everything(), ~sum(!is.na(.)), .names = "not_NA_{.col}")) %>%
  t()


## compare high_stable vs healthy ####
master.Table.high <- master.Table%>%
  filter(stable_high %in% c('YES', 'healthy'))%>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "stable_high" , data = master.Table.high)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_high_stable_healthy_comparisons.csv")


## compare low_stable vs healthy ####
master.Table.low <- master.Table%>%
  filter(stable_low %in% c('YES', 'healthy'))%>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "stable_low" , data = master.Table.low)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_low_stable_healthy_comparisons.csv")

## compare low stable vs high stable 
master.Table.high.low <- master.Table%>%
  filter(stable_high == "YES" | stable_low == "YES") %>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "stable_low" , data = master.Table.high.low)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_low_stable_high_stable_comparisons.csv")
