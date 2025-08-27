library(dplyr)
library(psych)
library(fastDummies)
library(ggplot2)
# load module eigengenes, metadata, cell-types:
setwd("~/Work/RP2/ATLANTIS")

res.cibersort <- read.csv("./th_high_th_low/CIBERSORT.proportions.csv") %>%
  tibble::column_to_rownames ("X") %>%
  t() %>%
  as.data.frame() %>%

  tibble::rownames_to_column("Sample") 

master.Table <- read.csv("./th_high_th_low/master_table_th2.csv")

big_master_table <-  read.csv("./atlantis_patient_data.csv", header =TRUE, na.strings=c("","NA"))
additional_data <- read.csv("./Umi_dedup/Dif_expr/Asthma_groups/clinicaldata_bhr_feno_mld.csv", header =TRUE, na.strings=c("","NA"))%>%
  dplyr::select(-c(X,FENRES, MLD_ratio))
clinical_table <- big_master_table[,c('PT','PACKNO','BMI','PHADRES', 'DUR_DIS','AGE_DIAG','B_FEV1F','GINA',"MORE1EX",
                                      'ICS_MEAN_DOSE',
                                      'acq6_score','LABEOSV','LABNEUV','LABMACV', 'EOSP', 'NEUTROP',
                                      'B_RVTLCPNVF', 'B_FEV1PNVG','FENRES',
                                      'B_R520PNVR', 'B_SCONDPNVF','B_SACINPNVF', 'B_F50PNVR',
                                      'ICS', 'ICS_LABA', 'ICS_DDOSE_EQ', 'ICS_LABA_DDOSE_EQ')]
clinical_table <- clinical_table[!duplicated(clinical_table$PT), ] %>%
  mutate (GINA = as.integer(GINA))

clinical_cell_module_df <- read.csv('./th_high_th_low/WGCNA_th2low/WGCNA_modules_eigengenes.csv')%>%
  dplyr::rename(Sample = X) %>%
  dplyr::select(c(Sample, MEblack)) %>%
  left_join(res.cibersort, by = c("Sample" = "Sample")) %>%
  left_join(master.Table, by = c("Sample" = "GenomeScan_ID")) %>%
  left_join(clinical_table, by = c('PT'='PT')) %>%
  left_join(additional_data, by = c('PT'='PT')) %>%
  mutate(any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes")) %>%
  mutate(ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                                coalesce(ICS_DDOSE_EQ,0) + coalesce(ICS_LABA_DDOSE_EQ,0))) %>%
  dplyr::select(-c("ICS_LABA_DDOSE_EQ", "ICS_LABA_DDOSE_EQ", "PT", "X1", "filesize", "sample.id",
                   "ICS_MEAN_DOSE", "ICS_DDOSE_EQ")) %>%
  dummy_cols(select_columns = c("gender", "PHADRES", "MORE1EX", "smoking.status"),
             remove_first_dummy = TRUE,
             ignore_na = TRUE)

# Only numeric columns - no cell types
clinical_module_df <- clinical_cell_module_df %>%
  dplyr::select(-c("Basal.resting...Suprabasal", "Hillock.like", "Multiciliated.lineage",
                   "Club", "Goblet", "T.cell.lineage", "Dendritic.cells", "B_FEV1F", "asthma.status",
                   "MEblack"))
num_df <- clinical_module_df[sapply(clinical_module_df, is.numeric)]

cor_result <- psych::corr.test(clinical_cell_module_df$MEblack, num_df, method="pearson", adjust = "fdr")

cor_result$r
cor_result$p.adj


cor_clinical_module_df <- data.frame(cor_p = signif(t(cor_result$p), 2),
                                     cor_p_adj = signif(t(cor_result$p.adj), 2),
                                     cor_r = signif(t(cor_result$r), 2))

write.csv(cor_clinical_module_df, "./th_high_th_low/cor_blackME_vs_clinical.csv", quote = F)


## plot correlation between black module and cell types 
# Reshape to long format
df_long_cell_types <- clinical_cell_module_df %>%
  dplyr::select(c("Sample", "MEblack", "Basal.resting...Suprabasal", 
           "Hillock.like", "Multiciliated.lineage", "Club",
           "Goblet", "Dendritic.cells", "asthma.status")) %>%
  tidyr::pivot_longer(cols = all_of(c("Basal.resting...Suprabasal", 
                               "Hillock.like", "Multiciliated.lineage", "Club",
                               "Goblet", "Dendritic.cells")),
               names_to = "cell_types", values_to = "value") %>%
  mutate(cell_types = case_when(
    cell_types == "Basal.resting...Suprabasal" ~  "Basal resting",
    cell_types == "Hillock.like" ~  "Hillock like",
    cell_types == "Multiciliated.lineage" ~ "Multiciliated lineage",
    cell_types == "Dendritic.cells" ~ "Dendritic cells",
    .default = cell_types
  ))

# Compute correlation and p-values
stats <- df_long_cell_types %>%
  group_by(cell_types) %>%
  summarise(
    cor = cor(clinical_cell_module_df[["MEblack"]], value, method = "spearman"),
    pval = cor.test(clinical_cell_module_df[["MEblack"]], value, method = "spearman")$p.value,
    max_y = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# Merge stats into long data
df_long_cell_types <- df_long_cell_types %>%
  left_join(stats, by = "cell_types") %>%
  mutate(asthma.status = if_else(asthma.status == "H", "Healthy", "T2-low asthma"))

# Plot
plt_cor <- ggplot(df_long_cell_types, aes(x = .data[["MEblack"]], y = value)) +
  geom_point(aes(color = asthma.status)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ cell_types, scales = "free_y") +
  labs(x = "Black module eigengene", y = "Cell proportion") +
  geom_text(
    data = stats,
    aes(
      x = Inf, y = max_y + 0.02,
      label = paste0("r = ", round(cor, 2), "\n p = ", signif(pval, 3))
    ),
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "Healthy" = "#66C2A5",
    "T2-low asthma" = "#8DA0CB"
  )) +
  theme_minimal()
png('./th_high_th_low/plots/Cell_types_vs_MEBlack.png', res = 300,  width = 1880, height = 1400)
print(plt_cor)
dev.off()



