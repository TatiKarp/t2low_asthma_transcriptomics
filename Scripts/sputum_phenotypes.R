# Perform de analysis between healthy controls and sputum neutrophilic and paucigranulocytic phenotypes

library(edgeR)
library(dplyr)

setwd("~/Work/RP2/ATLANTIS")

# 1. Load meta and clinical data
big_master_table <-  read.csv("./atlantis_patient_data.csv", header =TRUE, na.strings=c("","NA")) %>%
  dplyr::select(c('PT', 'VISIT', 'LABEOSV','LABNEUV','LABMACV','BRONCHP',
                  'LYMPHOP', 'EOSP', 'MACROP', 'NEUTROP','FENRES')) %>%
  filter(VISIT == "VISIT 1")

master.Table <- read.csv("./th_high_th_low/master_table_th2.csv") %>%
  left_join(big_master_table, by = c("PT" = "PT"))

filtered_master_table <- master.Table %>%
  filter(!is.na(EOSP)) %>%
  filter(!is.na(NEUTROP)) %>%
  filter(!(group_th == "healthy")) %>%
  filter(group_th == "low") %>%
  mutate(sputum_group = case_when(
    (EOSP < 2 & NEUTROP >= 50) ~ "sp_neu",
    (EOSP < 2 & NEUTROP < 50) ~ "sp_pau",
    .default = "none"
  ))

# table(filtered_master_table$sputum_group)

# sp_neu sp_pau 
# 21     18 

# 2. Create :
# -sputum paugranulocytic group: <2% eosinophils + <50% neutrophils
# -sputum neutrophilic group: <2% eosinophils + ≥50% neutrophils

master.Table.sp <- master.Table %>%
  left_join(filtered_master_table %>%
              dplyr::select(c(GenomeScan_ID, sputum_group)), by = c("GenomeScan_ID" = "GenomeScan_ID")) %>%
  mutate(sputum_group = case_when(
    group_th == "healthy" ~ "healthy",
    group_th == "low" ~ sputum_group)) %>%
  filter(!is.na(sputum_group))

table(master.Table.sp$group_th, master.Table.sp$sputum_group)

# 3. Dif expression only subset  
# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table.sp$GenomeScan_ID)) %>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + sputum_group + age + gender + smoking.status, data = master.Table.sp)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
# FALSE  TRUE 
# 34698 16873 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")

DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgse 4.0 and above- to reproduce the results of edgeR 3.0

contrasts <- limma::makeContrasts(
  sp_neu = sputum_groupsp_neu - sputum_grouphealthy,
  sp_pau = sputum_groupsp_pau - sputum_grouphealthy,
  levels = design
)

qlf_sp_neu  <- edgeR::glmQLFTest(fit, contrast = contrasts[,"sp_neu"])
qlf_sp_pau <- edgeR::glmQLFTest(fit, contrast = contrasts[,"sp_pau"])

summary(decideTests(qlf_sp_neu)) 
summary(decideTests(qlf_sp_pau)) 


## 4. Dif expr all samples:

master.Table.all <- master.Table %>%
  left_join(filtered_master_table %>%
              dplyr::select(c(GenomeScan_ID, sputum_group)), by = c("GenomeScan_ID" = "GenomeScan_ID")) %>%
  mutate(group_th = case_when(
    sputum_group == "sp_pau" ~ "sp_pau",
    sputum_group == "sp_neu" ~ "sp_neu",
    .default = group_th)) 

table(master.Table.all$group_th)

# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table.all$GenomeScan_ID)) %>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + group_th + age + gender + smoking.status, data = master.Table.all)

# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
table(keep)
# FALSE  TRUE 
# 33580 17991 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)

# Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgse 4.0 and above- to reproduce the results of edgeR 3.0

contrasts <- limma::makeContrasts(
  sp_neu = group_thsp_neu - group_thhealthy,
  sp_pau = group_thsp_pau - group_thhealthy,
  levels = design
)

qlf_sp_neu  <- edgeR::glmQLFTest(fit, contrast = contrasts[,"sp_neu"])
qlf_sp_pau <- edgeR::glmQLFTest(fit, contrast = contrasts[,"sp_pau"])

summary(decideTests(qlf_sp_neu)) 
summary(decideTests(qlf_sp_pau)) 
