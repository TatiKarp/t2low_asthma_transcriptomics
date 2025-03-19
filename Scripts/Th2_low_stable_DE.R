# this script will select stable Th2-low group and perform DGE analysis comparing that group vs healthy

library(dplyr)
library(readr)
library(ggplot2)
library(biomaRt)
library(edgeR) 
setwd("/Users/tatiana/Work/RP2/ATLANTIS")
# master.Table <- read_csv("./th_high_th_low/master_table_th2.csv")

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
# column at_least_2 == -3 if every 3 visits is Th2 low,
# column at_least_2 == -2 if at least 2 visits is Th2 low,


## plot only for samples with RNA-seq data

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

write.csv(master.Table, "./th_high_th_low/master_table_th2_stable.csv", row.names = FALSE)

## 22 samples with "stable" Th2-low

master.Table.stable <- master.Table %>%
  filter(!is.na(stable_low))
# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table.stable$GenomeScan_ID))%>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + stable_low + age + gender + smoking.status, data = master.Table.stable)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# Fit the model
DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL,design)

contrasts <- limma::makeContrasts(
  lowstable_healthy = stable_lowYES - stable_lowhealthy,
  levels = design
)

qlf  <- edgeR::glmQLFTest(fit, contrast = contrasts[,"lowstable_healthy"])

## add gene names
ensembl38 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 108)

dataprobes <- getBM(filters = "ensembl_gene_id", attributes = c(
  "ensembl_gene_id",
  "external_gene_name"),
  values = as.character(rownames(expression.data)),
  mart = ensembl38
)

de.results <- edgeR::topTags(
  qlf,
  n=nrow(DGEL))$table %>% 
  tibble::rownames_to_column("Gene") %>%
  left_join(dataprobes, by = c("Gene" = "ensembl_gene_id"))

###### Define Th2-low stable as at least 2 visits as Th2-low #######

master.Table.stable.2 <- master.Table %>%
  filter(!is.na(at_least_2)) %>%
  mutate(groups_at_least_2 = case_when(
    ((at_least_2 == "-2") | (at_least_2 == "-3")) ~ "stable_low",
    at_least_2 == "healthy" ~ "healthy", 
  TRUE	~ "rest"))
# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table.stable.2$GenomeScan_ID))%>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + groups_at_least_2 + age + gender + smoking.status, data = master.Table.stable.2)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# Fit the model
DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design)

contrasts <- limma::makeContrasts(
  lowstable_healthy = groups_at_least_2stable_low - groups_at_least_2healthy,
  levels = design
)

qlf  <- edgeR::glmQLFTest(fit, contrast = contrasts[,"lowstable_healthy"])
summary(decideTests(qlf)) 

# add gene names
ensembl38 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 108)

dataprobes <- getBM(filters = "ensembl_gene_id", attributes = c(
  "ensembl_gene_id",
  "external_gene_name"),
  values = as.character(rownames(expression.data)),
  mart = ensembl38
)

de.results <- edgeR::topTags(
  qlf,
  n=nrow(DGEL))$table %>% 
  tibble::rownames_to_column("Gene") %>%
  left_join(dataprobes, by = c("Gene" = "ensembl_gene_id"))

#------------------
# plot volcano plot
#------------------

library(ggplot2)
library(ggrepel)

source ("/Users/tatiana/Work/RP2/ATLANTIS_project/source_scripts/DE_scripts.R")
# th2 low
png("./th_high_th_low/plots/TH2_low_stable_genes_volcano.png",width=900, height=900,res = 150)
volcano_low <- volcano.plot(de.results, de.results$logFC, de.results$FDR, "Th2 low stable vs healthy")

print(volcano_low)
dev.off()
