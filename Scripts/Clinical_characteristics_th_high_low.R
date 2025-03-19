### this script will create clinical table for th2 high / th2-low asthma
library(ggplot2)
library(dplyr)
library(tableone)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")

############### prepare the data ###########################
## upload all data
master.Table <- read.csv("./th_high_th_low/master_table_th2.csv")

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

## add seasons 
master.Table.seasons <- read.csv("./Season/Season_new_date/ATLANTIS_master_table_seasons_new.csv")

# combine all
master.Table <- master.Table %>%
  filter(QC_check == 'YES') %>%
  mutate(gender = as.factor(gender),
         smoking.status = as.factor(smoking.status),
         group_th = as.factor(group_th)) %>%
  left_join(clinical_table, by = c('PT'='PT')) %>%
  left_join(additional_data, by = c('PT'='PT')) %>%
  mutate(any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes")) %>%
  mutate(ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                                coalesce(ICS_DDOSE_EQ,0) + coalesce(ICS_LABA_DDOSE_EQ,0))) %>%
  left_join(master.Table.seasons %>%
              dplyr::select(c(GenomeScan_ID, season.Maaike, season_2_Maaike)), by = c("GenomeScan_ID" = "GenomeScan_ID"))


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
                   'VI_856', 'VI_950', 'lung_ratio','MLD_ratio', "MORE1EX", "season.Maaike", "season_2_Maaike")

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th" , data = master.Table)

Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)

write.csv(Final_statistics, file = "./th_high_th_low/th2_high_low_baseline_statistics_all_comparisons.csv")
          
#### no undetermined samples ####
master.Table.no.undet <- master.Table%>%
  filter(as.character(group_th) != 'undeterm')%>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th" , data = master.Table.no.undet)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_high_low_baseline_statistics_no_undeterm_comparisons.csv")

## compare high vs healthy ####
master.Table.high <- master.Table%>%
  filter(!(group_th %in% c('low', 'undeterm')))%>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th" , data = master.Table.high)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_high_healthy_comparisons.csv")

## compare low vs healthy ####
master.Table.low <- master.Table %>%
  filter(!(group_th %in% c('high', 'undeterm'))) %>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th" , data = master.Table.low)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/th2_low_healthy_comparisons.csv")

## compare undetermined vs healthy 
master.Table.undetermined <- master.Table %>%
  filter(group_th %in% c("undeterm", "healthy")) %>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th", data = master.Table.undetermined)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/Intermediate_healthy_comparisons.csv")

## compare high vs low 
master.Table.high.low <- master.Table %>%
  filter(group_th %in% c("high", "low")) %>%
  droplevels()

tabtotal <- CreateTableOne(vars = var_for_table, strata = "group_th", data = master.Table.high.low)
Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)
write.csv(Final_statistics, file = "./th_high_th_low/high_low_comparisons.csv")

## chi square exacerbations

master.Table.high.low <- master.Table %>%
  filter(group_th %in% c("high", "low"))
table(master.Table.high.low$group_th, master.Table.high.low$PHADRES)
chisq.test(master.Table.high.low$group_th, master.Table.high.low$PHADRES, correct=FALSE)


# chi square check : 
master.Table.filtered <- master.Table %>%
  filter(group_th %in% c("high", "low")) %>%
  droplevels()
table(master.Table.filtered$group_th, master.Table.filtered$PHADRES)
chisq.test(master.Table.filtered$group_th, master.Table.filtered$PHADRES, correct=TRUE)


