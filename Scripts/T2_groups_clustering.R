# This script will perform clustering of T2-high and T2-low asthma group and check which group is more heterogeneous 
# gap statistics is used to asses the optimal N of clusters
library(dplyr)
library(NbClust) 
library(cluster)
library(ggplot2)

setwd("Work/RP2/ATLANTIS")

# metadata
master.Table <- read.csv('./th_high_th_low/master_table_th2.csv')

# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

# normalize cpm
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
norm_expression <- as.data.frame(edgeR::cpm(DGEL,normalized.lib.sizes=TRUE, log=TRUE))

# function to create k clusters 
mycluster <- function(x, k) {list(cluster = cutree(hclust(dist(x, method = 'euclidian'), method = "complete"), k=k))}


# 1. Cluster t2low 
master.Table.low <- master.Table %>%
  filter(group_th == "low")
expression.low <- norm_expression[, master.Table.low$GenomeScan_ID]

# 1.1 Select Most Variable Genes:
# calculate coefficient of variation
row_stdev <- apply(expression.low, 1, sd)
row_mean <- apply(expression.low, 1, mean)

coef_of_var_low <- data.frame(Gene = rownames(expression.low),
                          stdev = row_stdev,
                          mean = row_mean, 
                          coef_of_variation = row_stdev/row_mean)%>%
  arrange(desc(coef_of_variation))

# select top 500 and 2000 most variable genes:
n_genes <- c(500, 2000)
gap_plot_list_low <- list()

for (n in n_genes) {
  top_n_variable_low <- as.matrix(expression.low[coef_of_var_low$Gene[1:n], master.Table.low$GenomeScan_ID])
  # Scale the Data
  top_n_variable_low_scaled <- scale(t(top_n_variable_low)) # scale by gene
  # 1.2 calculate gap statistics
  # Gap stat with factoextra package (just to check )
  dist_matrix_gap_low <- factoextra::fviz_nbclust(as.matrix(top_n_variable_low_scaled),
                                                  mycluster, 
                                                  method = c("gap")) +
    ggtitle(paste0("T2-low, top ", as.character(n)," variable genes"))
  png(paste0("./th_high_th_low/plots/gap_stat_t2low_", as.character(n), "_variable.png"),
      width = 600, height = 500, res = 150)
  print(dist_matrix_gap_low)
  dev.off()
  gap_plot_list_low[[paste0("t2low", as.character(n))]] <- dist_matrix_gap_low
}


# 2. Cluster t2high
master.Table.high <- master.Table %>%
  filter(group_th == "high")
expression.high <- norm_expression[, master.Table.high$GenomeScan_ID]

# 2.1 Select Most Variable Genes:
# calculate coefficient of variation
row_stdev <- apply(expression.high, 1, sd)
row_mean <- apply(expression.high, 1, mean)

coef_of_var_high <- data.frame(Gene = rownames(expression.high),
                          stdev = row_stdev,
                          mean = row_mean, 
                          coef_of_variation = row_stdev/row_mean)%>%
  arrange(desc(coef_of_variation))


# select top500 most variable genes:
n_genes <- c(500, 2000)
gap_plot_list_high <- list()

for (n in n_genes) {
  top_n_variable_high <- as.matrix(expression.high[coef_of_var_high$Gene[1:n], master.Table.high$GenomeScan_ID])
  # Scale the Data
  top_n_variable_high_scaled <- scale(t(top_n_variable_high)) # scale by gene
  # 1.2 calculate gap statistics
  # Gap stat with factoextra package (just to check )
  dist_matrix_gap_high <- factoextra::fviz_nbclust(as.matrix(top_n_variable_high_scaled),
                                                  mycluster, 
                                                  method = c("gap")) +
    ggtitle(paste0("T2-high, top ", as.character(n)," variable genes"))
  
  png(paste0("./th_high_th_low/plots/gap_stat_t2high_", as.character(n), "_variable.png"),
      width = 600, height = 500, res = 150)
  print(dist_matrix_gap_high)
  dev.off()
  gap_plot_list_high[[paste0("t2high", as.character(n))]] <- dist_matrix_gap_high
}

library(ggpubr)
combined_plot <- ggarrange(gap_plot_list_low[["t2low500"]], gap_plot_list_high[["t2high500"]],
            gap_plot_list_low[["t2low2000"]], gap_plot_list_high[["t2high2000"]],
            nrow = 2, ncol = 2,
            labels = c('A', 'B', "C", "D")) 

png(paste0("./th_high_th_low/plots/gap_all.png"),
    width = 1800, height = 1500, res = 200)
print(combined_plot)
dev.off()

# 3.0 Cluster T2-high and T2-low together

master.Table.low.high <- master.Table %>%
  filter(group_th %in% c("low", "high"))
expression.low.high <- norm_expression[, master.Table.low.high$GenomeScan_ID]

# 3.1 Select Most Variable Genes:
# calculate coefficient of variation
row_stdev <- apply(expression.low.high, 1, sd)
row_mean <- apply(expression.low.high, 1, mean)

coef_of_var <- data.frame(Gene = rownames(expression.low.high),
                          stdev = row_stdev,
                          mean = row_mean, 
                          coef_of_variation = row_stdev/row_mean) %>%
  arrange(desc(coef_of_variation))

# select top500 most variable genes:
top_500_variable_low_high <- as.matrix(expression.low.high[coef_of_var$Gene[1:500], master.Table.low.high$GenomeScan_ID])

heatmap <- pheatmap::pheatmap(top_500_variable_low_high,
                              annotation_col = master.Table.low.high %>%
                                dplyr::select(c(GenomeScan_ID, group_th, gender)) %>%
                                tibble::column_to_rownames("GenomeScan_ID"),
                              scale = "row",
                              clustering_method = "complete",  
                              show_rownames = FALSE,
                              show_colnames = FALSE)



# 4. Cluster all together (t2low, t2-high, healthy, undetermined)
# 4.1 Select Most Variable Genes:
# calculate coefficient of variation
row_stdev <- apply(norm_expression, 1, sd)
row_mean <- apply(norm_expression, 1, mean)

coef_of_var <- data.frame(Gene = rownames(norm_expression),
                          stdev = row_stdev,
                          mean = row_mean, 
                          coef_of_variation = row_stdev/row_mean)%>%
  arrange(desc(coef_of_variation))

# select top500 most variable genes:
top_500_variable_all <- as.matrix(norm_expression[coef_of_var$Gene[1:500], master.Table$GenomeScan_ID])

# scale the Data
top_500_variable_all_scaled <- scale(t(top_500_variable_all)) # scale by gene

# 4.2 calculate gap statistics
# calculate gap statistics using this function: as input use scaled matrix
#gap_stat_all <- clusGap(as.matrix(top_500_variable_all_scaled), FUN = mycluster,  K.max = 15, B = 100, d.power = 1)

# plot - repeat gap stat calculation
dist_matrix_gap_all <- factoextra::fviz_nbclust(as.matrix(top_500_variable_all_scaled), mycluster, 
                                                 method = c("gap"),
                                                 print.summary = T) +
  ggtitle("All hclust, euclidian, complete")

png("./th_high_th_low/plots/gap_stat_all.png",
    width = 600, height = 500, res = 150)
print(dist_matrix_gap_all)
dev.off()

## plot heatmap 
heatmap_all <- pheatmap::pheatmap(top_500_variable_all,
                              annotation_col = master.Table %>%
                                dplyr::select(c(GenomeScan_ID, group_th, gender)) %>%
                                tibble::column_to_rownames("GenomeScan_ID"),
                              scale = "row",
                              clustering_method = "complete",  # or "complete", "ward.D2", etc.
                              show_rownames = FALSE,
                              show_colnames = FALSE)
png("./th_high_th_low/plots/hclust_heatmap_all.png",
    width = 1000, height = 900, res = 150)
print(heatmap_all)
dev.off()

# 4.3 cut the tree at n = 4 
# Compute distance matrix (Euclidean by default)
dist_matrix <- dist(top_500_variable_all_scaled)
# Perform hierarchical clustering (default: complete linkage)
hc <- hclust(dist_matrix, method = "complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "", cex = 0.8)
# Cut the tree into 4 clusters
clusters <- cutree(hc, k = 4)

# Add cluster labels to the master table
master.Table.clusters <- master.Table %>%
  left_join(data.frame(ID = names(clusters), 
                       Clusters_4 = clusters), by = c("GenomeScan_ID" = "ID"))

table(master.Table.clusters$Clusters_4, master.Table.clusters$group_th)

chisq.test(table(master.Table.clusters$Clusters_4, master.Table.clusters$group_th))
chisq.test(table(master.Table.clusters$group_th, master.Table.clusters$Clusters_4))

# Cut the tree into 2 clusters
clusters_2 <- cutree(hc, k = 2)

# Add cluster labels to the master table
master.Table.clusters <- master.Table.clusters %>%
  left_join(data.frame(ID = names(clusters_2), 
                       Clusters_2 = clusters_2), by = c("GenomeScan_ID" = "ID"))

table(master.Table.clusters$Clusters_2, master.Table.clusters$group_th)

chisq.test(table(master.Table.clusters$Clusters_2, master.Table.clusters$group_th))




