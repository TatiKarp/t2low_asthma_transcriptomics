library(dplyr)
library(edgeR)
library(biomaRt)

setwd('/Users/tatiana/Work/RP2/ATLANTIS')
master.Table <- read.csv('./Umi_dedup/Dif_expr/ATLANTIS_master_table_QC.csv',header = TRUE)%>%
  filter (QC_check == 'YES')
# add additional data 
big_master_table <-  read.csv('./atlantis_patient_data.csv', header =TRUE)%>%
  dplyr::select(c(PT, SYS_COR, BIO, LABEOSV, FENRES))

big_master_table <- big_master_table[!duplicated(big_master_table$PT), ]

master.Table.TH <- master.Table%>%
  filter(asthma.status == 'A')%>%
  left_join(big_master_table, by = c('PT'='PT'))%>%
  mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES'))%>%
  mutate(group_th = case_when(
    ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
    ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
    ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
    ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low'))%>%
  mutate(group_th = if_else(is.na(group_th), 'undeterm', group_th))%>%
  mutate(group_th = if_else(include == 'NO', 'cor_bio', group_th))

table(master.Table.TH$group_th)
# cor_bio     high      low undeterm 
# 22       63       82      194 

master.Table <- master.Table%>%
  left_join(master.Table.TH%>%
              dplyr::select(c(GenomeScan_ID, group_th)))%>%
  mutate(group_th = if_else(is.na(group_th), 'healthy', group_th))%>%
  filter(!(group_th  == 'cor_bio'))%>%
  mutate(group_th = as.factor(group_th))

write.table(master.Table, './th_high_th_low/master_table_th2.csv', sep = ',', row.names = F, quote = F)
# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + group_th + age + gender + smoking.status, data = master.Table)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# Fit the model
DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgse 4.0 and above- to reproduce the results of edgeR 3.0

contrasts <- limma::makeContrasts(
  high_healthy = group_thhigh - group_thhealthy,
  low_healthy = group_thlow - group_thhealthy,
  high_low = group_thhigh - group_thlow,
  levels = design
)

qlf_high   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_healthy"])
qlf_low <- edgeR::glmQLFTest(fit, contrast = contrasts[,"low_healthy"])
qlf_high_low <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_low"])

summary(decideTests(qlf_high)) 
summary(decideTests(qlf_low)) 

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

# create the result table
de.results.high <- edgeR::topTags(
  qlf_high,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = dataprobes,
    by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
              PValue = signif(PValue, digits = 3),
              FDR = signif(FDR, digits = 3))
# write.table(de.results.high,'./th_high_th_low/DE_genes_THhigh_healthy.csv', sep = ',', row.names = F, quote = F)
# de.results.high.sign <- de.results.high %>%
#   filter(FDR < 0.05) %>%
#   arrange(-logFC) 
# write.table(de.results.high.sign,'./th_high_th_low/DE_genes_THhigh_healthy_sign.csv', sep = ',', row.names = F, quote = F)

de.results.low <- edgeR::topTags(
  qlf_low,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = dataprobes,
    by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
          PValue = signif(PValue, digits = 3),
          FDR = signif(FDR, digits = 3))
# write.table(de.results.low,'./th_high_th_low/DE_genes_THlow_healthy.csv', sep = ',', row.names = F, quote = F)

de.results.high.low <- edgeR::topTags(
  qlf_high_low,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = dataprobes,
    by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
          PValue = signif(PValue, digits = 3),
          FDR = signif(FDR, digits = 3))
# write.table(de.results.high.low,'./th_high_th_low/DE_genes_THhigh_THlow.csv', sep = ',', row.names = F, quote = F)

