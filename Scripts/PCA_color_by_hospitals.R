#This script will perform PCA on the genes with top 2000 average logcpm, and color it by hospital
library(dplyr)
library(ggrepel)
library(edgeR)
library(rpca)
library(ggpubr)

setwd("/Users/tatiana/Work/RP2/ATLANTIS/")
master.Table <- read.csv("./th_high_th_low/master_table_th2.csv", header = TRUE) %>%
  filter (QC_check == 'YES') 

hospital.Table <- read.csv('./patient_data_final.csv',header = TRUE) %>%
  as.data.frame()%>%
  dplyr::select(c(PT,hospital))

master.Table <- master.Table%>%
  left_join(hospital.Table, by = c('PT' = 'PT'))%>%
  mutate(hospital = as.factor(hospital))

expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

# normalize data 
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
ATL_logCPM <- cpm(DGEL, log = T)

# calculate average logcpm per gene
#select genes with top 2000 average logcpm
Atlantis_aveLogCPM <- data.frame(Gene = rownames(ATL_logCPM), average_log_cpm = rowMeans(ATL_logCPM))
Atlantis_aveLogCPM  <- Atlantis_aveLogCPM [order(-Atlantis_aveLogCPM$average_log_cpm),]
Top_aveLogCPM <- Atlantis_aveLogCPM$Gene[1:2000]

# extract row counts for top 2000 genes
top_expression <- expression.data[rownames(expression.data) %in% Top_aveLogCPM,]

# voom Normalization
norm.expr.data.top <- top_expression[rowSums(top_expression) >= 10,] %>%
  limma::voom() %>%
  as.matrix()

norm.expr.data.top <- norm.expr.data.top[,master.Table$GenomeScan_ID]

norm.expr.data.pca.top <- norm.expr.data.top %>%
  t() %>% #transpose the matrix to calculate the components by sample 
  stats::prcomp(
    center = TRUE,
    scale. = FALSE)

## 1. Check hospitals
pc12_hosp <- ggbiplot::ggbiplot(norm.expr.data.pca.top, 
                                choices = 1:2, 
                                obs.scale = 1, 
                                var.scale = 1, 
                                ellipse = TRUE, 
                                groups = master.Table$hospital,
                                point.size = 1,
                                labels = NULL, 
                                var.axes = FALSE,
                                ellipse.prob = 0.95,
                                ellipse.fill = FALSE, 
                                ellipse.linewidth = 0.5) +
  theme_minimal() +
  guides(color=guide_legend("Hospitals"))


pc23_hosp <- ggbiplot::ggbiplot(norm.expr.data.pca.top, 
                                choices = 2:3, 
                                obs.scale = 1, 
                                var.scale = 1, 
                                ellipse = TRUE, 
                                groups = master.Table$hospital,
                                point.size = 1,
                                labels = NULL, 
                                var.axes = FALSE,
                                ellipse.prob = 0.95,
                                ellipse.fill = FALSE, 
                                ellipse.linewidth = 0.5) +
  theme_minimal()+
  guides(color=guide_legend("Hospitals"))


pc34_hosp <- ggbiplot::ggbiplot(norm.expr.data.pca.top, 
                                choices = 3:4, 
                                obs.scale = 1, 
                                var.scale = 1, 
                                ellipse = TRUE, 
                                groups = master.Table$hospital,
                                point.size = 1,
                                labels = NULL, 
                                var.axes = FALSE,
                                ellipse.prob = 0.95,
                                ellipse.fill = FALSE, 
                                ellipse.linewidth = 0.5) +
  theme_minimal()+
  guides(color=guide_legend("Hospitals"))


pc45_hosp <- ggbiplot::ggbiplot(norm.expr.data.pca.top, 
                                choices = 4:5, 
                                obs.scale = 1, 
                                var.scale = 1, 
                                ellipse = TRUE, 
                                groups = master.Table$hospital,
                                point.size = 1,
                                labels = NULL, 
                                var.axes = FALSE,
                                ellipse.prob = 0.95,
                                ellipse.fill = FALSE, 
                                ellipse.linewidth = 0.5) +
  theme_minimal()+
  guides(color=guide_legend("Hospitals"))


figure <- ggarrange(pc12_hosp, pc23_hosp,
                    pc34_hosp, pc45_hosp,
                    common.legend = TRUE,
                    legend="right",
                    labels = c("A", "B", "C", "D"),
                    nrow = 2, ncol = 2)

png("./th_high_th_low/plots/PC12345_hospitals.png",
    width = 2000, 
    height = 2000,
    res = 300)
print(figure)
dev.off()

