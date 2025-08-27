library(dplyr) 

## 1. Prepare the data 
setwd("~/Work/RP2/ATLANTIS")


# load T2-high genes 
de.results.high <- read.csv("./th_high_th_low/DE_genes_THhigh_healthy.csv") %>%
  filter(FDR < 0.05)

# load expression data (logCPM)
cpm_expr <- read.csv("./th_high_th_low/logcpm_norm_expr.csv") %>%
  tibble::column_to_rownames("X")

# load metadata
clinical_table <-  read.csv("./atlantis_patient_data.csv", header =TRUE, na.strings=c("","NA")) %>%
  dplyr::select(c('PT', 'VISIT', 'LABEOSV', 'FENRES')) %>%
  filter(VISIT == "VISIT 1")

master.Table <- read.csv("./th_high_th_low/master_table_th2.csv") %>%
  left_join(clinical_table, by = c("PT" = "PT")) 
master.Table.high <- master.Table %>%
  filter(group_th == "high")

## 2. Run correlation 
# 2.1. prepare data for correlation - only T2high
expr_sign_genes <- cpm_expr[de.results.high$Gene, master.Table.high$GenomeScan_ID] %>%
  t() %>%
  as.data.frame()

blood_eos_cor <- psych::corr.test(master.Table.high$LABEOSV, 
                                  expr_sign_genes, 
                                  method = "pearson",
                                  adjust = "fdr")
sum(blood_eos_cor$p.adj < 0.05)
sum(blood_eos_cor$p < 0.05)
blood_eos_cor_df <- data.frame(cor_p_bec = signif(t(blood_eos_cor$p), 2),
                               cor_p_adj_bec = signif(t(blood_eos_cor$p.adj), 2))
                                  

feno_cor <- psych::corr.test(master.Table.high$FENRES, 
                             expr_sign_genes, 
                             method = "pearson",
                             adjust = "fdr")
sum(feno_cor$p.adj < 0.05)
sum(feno_cor$p < 0.05)

feno_cor_df <- data.frame(cor_p_feno = signif(t(feno_cor$p), 2),
                          cor_p_adj_feno = signif(t(feno_cor$p.adj), 2))

# Combine results: 
genes_vs_bec_feno <- de.results.high %>%
  left_join(blood_eos_cor_df %>%
              tibble::rownames_to_column("Gene"), by = c("Gene" = "Gene")) %>%
  left_join(feno_cor_df %>%
              tibble::rownames_to_column("Gene"), by = c("Gene" = "Gene")) 
  

# 2.2. prepare data for correlation - all samples
expr_sign_genes_all <- cpm_expr[de.results.high$Gene, master.Table$GenomeScan_ID] %>%
  t() %>%
  as.data.frame()

blood_eos_cor_all <- psych::corr.test(master.Table$LABEOSV, 
                                      expr_sign_genes_all, 
                                      method = "pearson",
                                      adjust = "fdr")
sum(blood_eos_cor_all$p.adj < 0.05)
sum(blood_eos_cor_all$p < 0.05)
blood_eos_cor_df_all <- data.frame(cor_p_bec_all = signif(t(blood_eos_cor_all$p), 2),
                               cor_p_adj_bec_all = signif(t(blood_eos_cor_all$p.adj), 2))


feno_cor_all <- psych::corr.test(master.Table$FENRES, 
                                 expr_sign_genes_all, 
                                 method = "pearson",
                                 adjust = "fdr")
sum(feno_cor_all$p.adj < 0.05)
sum(feno_cor_all$p < 0.05)
feno_cor_df_all <- data.frame(cor_p_feno_all = signif(t(feno_cor_all$p), 2),
                          cor_p_adj_feno_all = signif(t(feno_cor_all$p.adj), 2))

# Combine results: 
genes_vs_bec_feno <- genes_vs_bec_feno %>%
  left_join(blood_eos_cor_df_all %>%
              tibble::rownames_to_column("Gene"), by = c("Gene" = "Gene")) %>%
  left_join(feno_cor_df_all %>%
              tibble::rownames_to_column("Gene"), by = c("Gene" = "Gene")) 

# Add classification :
genes_vs_bec_feno <- genes_vs_bec_feno %>%
  mutate(group_all = case_when(
    (cor_p_adj_feno_all < 0.05 & cor_p_adj_bec_all < 0.05) ~ "eos_feno_all",
    (cor_p_adj_feno_all < 0.05 & cor_p_adj_bec_all > 0.05) ~ "feno_all",
    (cor_p_adj_feno_all > 0.05 & cor_p_adj_bec_all < 0.05) ~ "eos_all", 
    .default = "none"
  ))

write.csv(genes_vs_bec_feno, "./th_high_th_low/genes_cor_with_eos_feno.csv", quote = F)

## 3. Run DE analysis for T2high genes for an intermediate group
# 3.1 create the subgroups              
master.Table.undet <- master.Table %>%
  mutate(group_feno_bec = case_when(
    (group_th == "undeterm" & FENRES > 25 & LABEOSV < 0.15) ~ "high_feno_low_bec",
    (group_th == "undeterm" & FENRES < 25 & LABEOSV > 0.3)  ~ "high_bec_low_feno",
    (group_th == "healthy" ~ "healthy"))) %>%
  filter(!is.na(group_feno_bec))

table(master.Table.undet$group_feno_bec)


## 3.2 Run the DGE 
library(edgeR)
  
# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table.undet$GenomeScan_ID)) %>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + group_feno_bec + age + gender + smoking.status, data = master.Table.undet)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")

DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgse 4.0 and above- to reproduce the results of edgeR 3.0

contrasts <- limma::makeContrasts(
  high_bec = group_feno_bechigh_bec_low_feno - group_feno_bechealthy,
  high_feno = group_feno_bechigh_feno_low_bec - group_feno_bechealthy,
  levels = design
)

qlf_high_bec   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_bec"])
qlf_high_feno <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_feno"])

summary(decideTests(qlf_high_bec)) 
summary(decideTests(qlf_high_feno)) 

de.results.both <- edgeR::topTags(
  qlf_high_bec,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  mutate (logFC_high_bec = round (logFC, 2),
          PValue_high_bec = signif(PValue, digits = 3),
          FDR_high_bec = signif(FDR, digits = 3)) %>%
  dplyr::select(c(Gene, logFC_high_bec, PValue_high_bec, FDR_high_bec)) %>%
  left_join(edgeR::topTags(
    qlf_high_feno,
    n=nrow(DGEL))$table %>%
      tibble::rownames_to_column("Gene") %>%
      mutate (logFC_high_feno = round (logFC, 2),
              PValue_high_feno = signif(PValue, digits = 3),
              FDR_high_feno = signif(FDR, digits = 3)) %>%
      dplyr::select(c(Gene, logFC_high_feno, PValue_high_feno, FDR_high_feno)), by = c("Gene" = "Gene")) 


# add p_val to the T2-high genes
de.results.high.undet <- de.results.high %>%
  left_join(de.results.both, by = c("Gene" = "Gene")) %>%
  mutate(adj_p_high_bec = signif(p.adjust(PValue_high_bec, method = 'fdr'), digits = 3),
         adj_p_high_feno = signif(p.adjust(PValue_high_feno, method = 'fdr'), digits = 3)) %>%
  mutate(gene_group = case_when(
    (adj_p_high_bec < 0.05 & adj_p_high_feno < 0.05) ~ "both_bec_feno",
    (adj_p_high_bec < 0.05 & adj_p_high_feno > 0.05) ~ "high_bec",
    (adj_p_high_bec > 0.05 & adj_p_high_feno < 0.05) ~ "high_feno")) %>%
  mutate(log_fc_direction_high_bec = if_else((logFC_high_bec * logFC) > 0, "same", "opposite"))
write.csv(de.results.high.undet, "./th_high_th_low/DE_genes_THhigh_healthy_high_bec_high_feno.csv",
          quote = F)

# table(de.results.high.undet$gene_group)
# 
# high_bec 
# 438 - same direction