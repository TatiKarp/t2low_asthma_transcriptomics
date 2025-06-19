library(dplyr)
library(ggrepel)
library(edgeR)
library(rpca)
library(ggpubr)
library("factoextra")


setwd("Work/RP2/ATLANTIS/")

# 1. prepare data
master.Table <- read.csv("./th_high_th_low/master_table_th2.csv", header = TRUE) %>%
  filter (QC_check == 'YES') 

hospital.Table <- read.csv('./patient_data_final.csv',header = TRUE) %>%
  as.data.frame()%>%
  dplyr::select(c(PT,hospital))


expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

# 1.1 Normalize expression data 
# Create an edgeR object, normalize

DGEL <- edgeR::DGEList(expression.data)
# keep <- edgeR::filterByExpr(DGEL) 
# DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
norm_expression <- edgeR::cpm(DGEL, normalized.lib.sizes=TRUE, log = FALSE)

# check filtering: 
keep_genes <- rowSums(norm_expression >= 50) >= 20 #  at least 50 CPM in at least 20 samples
filtered_expression <- norm_expression[keep_genes, ]

# 1. 2 Select most variable genes (top 2000 ) : 
# calculate coefficient of variation
row_stdev <- apply(filtered_expression, 1, sd)
row_mean <- apply(filtered_expression, 1, mean)

coef_of_var <- data.frame(Gene = rownames(filtered_expression),
                          stdev = row_stdev,
                          mean = row_mean, 
                          coef_of_variation = row_stdev/row_mean) %>%
  arrange(desc(coef_of_variation))

# select top500 most variable genes:
top_2000_variable <- as.matrix(norm_expression[coef_of_var$Gene[1:2000], master.Table$GenomeScan_ID])

# voom Normalization
norm.expr.data.top <- top_2000_variable %>%
  limma::voom() %>%
  as.matrix()

# 2. Calculate PC
norm.expr.data.pca.top <- norm.expr.data.top %>%
  t() %>% #transpose the matrix to calculate the components by sample 
  stats::prcomp(
    center = TRUE,
    scale. = FALSE)

# 2.1 Scree plot 
#calculate total variance explained by each principal component
var_explained_df <- data.frame(var_explained = 100*(norm.expr.data.pca.top$sdev)^2/sum((norm.expr.data.pca.top$sdev)^2)) %>%
  tibble::rownames_to_column("PC")

var_explained_df$PC <- factor(var_explained_df$PC, levels = paste0(1:nrow(var_explained_df))) # order PCs

scree_plt <- ggplot(var_explained_df[1:15,], 
       aes(x = PC,
           y = var_explained, 
           group = 1)) +
  geom_col() +
  geom_point() +
  geom_line() +
  labs(y = "% variance explained",
       x = "Principal Component") +
  geom_hline(yintercept = 1, linetype='dashed') +
  annotate("text", 
           x = 12,     # Adjust x position to where you want the label
           y = 2,     # Slightly above the dashed line
           label = "1% variance", 
           size = 4,
           hjust = 0) +
  theme_minimal() +
  theme(axis.text=element_text(size = 15),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "plain"),
        legend.position = "none",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6))

    

png("./th_high_th_low/plots/scree_pca_most_variable_mendian_dist.png",
    width=1200, height=850, res = 200)

print(scree_plt)
dev.off()

pc_sum_th2 <- c("2", "3", "4", "6", "8", "9", "11")
var_explained_th2 <- var_explained_df %>%
  filter(PC %in% pc_sum_th2) 
sum(var_explained_th2$var_explained)

# 3. Calculate distances: 

# Initialize result table
results <- data.frame(
  PC = character(),
  Group = character(),
  Mean_Distance = numeric(),
  Median_Distance = numeric(),
  stringsAsFactors = FALSE
)

# Loop over PC1 to PC11
 # calculate p-values: 

# Initialize comparison results table
comparison_results <- data.frame(
  PC = character(),
  Test = character(),
  P_Value = numeric(),
  Statistic = numeric(),
  stringsAsFactors = FALSE
)

for (pc_index in 1:11) { #PCs that explain at least 1% variability
  pc_name <- paste0("PC", pc_index)
  PC_values <- norm.expr.data.pca.top$x[, pc_name]
  
  # Extract values for "low" and "high" groups
  low_values <- PC_values[master.Table$group_th == "low"]
  high_values <- PC_values[master.Table$group_th == "high"]
  
  # Skip comparison if either group has fewer than 2 samples
  if (length(low_values) < 2 | length(high_values) < 2) next
  
  # Compute within-group distances
  low_distances <- as.vector(dist(low_values))
  high_distances <- as.vector(dist(high_values))
  
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(low_distances, high_distances, exact = FALSE)
  
  # Store result
  comparison_results <- rbind(
    comparison_results,
    data.frame(
      PC = pc_name,
      Test = "low_vs_high",
      P_Value = test_result$p.value,
      Statistic = test_result$statistic,
      stringsAsFactors = FALSE
    )
  )
}

# save all the distances - for boxplots
distance_data <- data.frame(
  PC = character(),
  Group = character(),
  Distance = numeric(),
  stringsAsFactors = FALSE
)

for (pc_index in 1:11) {
  pc_name <- paste0("PC", pc_index)
  PC_values <- norm.expr.data.pca.top$x[, pc_name]
  
  # Everyone
  distance_matrix_all <- as.matrix(dist(PC_values, method = "euclidean"))
  distances_all <- distance_matrix_all[lower.tri(distance_matrix_all)]
  distance_data <- rbind(
    distance_data,
    data.frame(
      PC = pc_name,
      Group = "everyone",
      Distance = distances_all,
      stringsAsFactors = FALSE
    )
  )
  
  # Groups
  groups <- c("low", "high", "healthy", "undeterm")
  for (group in groups) {
    group_values <- PC_values[master.Table$group_th == group]
    if (length(group_values) < 2) next
    
    distance_matrix <- as.matrix(dist(group_values, method = "euclidean"))
    distances <- distance_matrix[lower.tri(distance_matrix)]
    
    distance_data <- rbind(
      distance_data,
      data.frame(
        PC = pc_name,
        Group = group,
        Distance = distances,
        stringsAsFactors = FALSE
      )
    )
  }
}

# 3.2 Plot distances:
# 3.2.1 boxplot:

distance_data <- distance_data %>%
  filter(!(Group %in% c("everyone", "undeterm"))) %>%
  mutate(Group = recode(Group, 
                        high = "T2-high",
                        healthy = "Healthy", 
                        low = "T2-low")) 

distance_data$Group <- factor(distance_data$Group, levels = c("Healthy", "T2-low", "T2-high"))

# calculate p_value:
stat.test <- distance_data %>%
  group_by(PC) %>%
  rstatix::wilcox_test(Distance ~ Group) %>%
  mutate(p.adj = signif(p.adjust(p, method = 'fdr'), digits = 2)) %>%
  rstatix::add_y_position(step.increase = 0, fun = "median") %>%
  mutate(y.position = y.position + (y.position*3.5)) %>%
  mutate(y.position = case_when(
    (group1 == "Healthy" & group2 == "T2-low") ~ y.position,
    (group1 == "Healthy" & group2 == "T2-high") ~ y.position + (y.position*0.15),
    (group1 == "T2-low" & group2 == "T2-high") ~ y.position + (y.position*0.3)
  ))
  # mutate(y.position + as.numeric(factor(PC)) * 0.1))

distance_data$PC <- factor(distance_data$PC, levels = paste0("PC", 1:11)) # order PCs

bxplt_per_pc <- ggplot(distance_data %>%
                         filter(!(Group %in% c("everyone", "undeterm"))), 
                         aes(x = Group, y = Distance)) +
  geom_boxplot(aes(fill = Group), outliers = FALSE) +
  facet_wrap(~factor(PC, levels = paste0("PC", 1:11)), scale = "free_y") +
  scale_fill_manual(values = c(
    "Healthy" = "#66C2A5",
    "T2-low" = "#8DA0CB",
    "T2-high" = "#FFD92F"
  ))+
  theme_minimal() +
  labs(y = expression(atop(bold("Samples dissimilarities"), "Pairwise distance over each principal component"))) +
  ggpubr::stat_pvalue_manual(
    data = stat.test,
    y.position = "y.position",
    label = "p adj = {p.adj}",
    size = 3,
    remove.bracket = FALSE
  ) +
  theme(axis.text=element_text(size = 15),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "plain"),
        legend.position = "none",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size = 12),
        legend.key = element_rect(fill = "white"))
png("./th_high_th_low/plots/bxplt_mendian_dist_on_pca_per_group_most_variable.png",
    width = 2000, height = 2000, res = 200)
print(bxplt_per_pc)
dev.off()


# median
plt_mendian_dist <- ggplot(results, aes(x = PC, y = Median_Distance, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Median Euclidean Distances by PC, most variable",
    x = "Principal Component",
    y = "Mendian Euclidean Distance"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

png("./th_high_th_low/plots/mendian_dist_on_pca_per_group_most_variable.png",
    width=1200, height=850, res = 200)

print(plt_mendian_dist)
dev.off()



# combine scree plot and distances 
library(patchwork)

combined_plot <- 
  ((scree_plt | plot_spacer()) /
  bxplt_per_pc) +
  plot_layout(heights = c(2, 5), widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
        plot.subtitle = element_text(size = 14),
        plot.tag = element_text(face = "bold"))

png("./th_high_th_low/plots/Supplementary_combined_bxplt_mendian_dist_on_pca_per_group_most_variable.png",
    width = 3000, height = 3500,
    res = 300)
print(combined_plot)
dev.off()


## 3.2 Check T2 group 
pc12_t2_group <- ggbiplot::ggbiplot(norm.expr.data.pca.top, 
                                    choices = 1:2, 
                                    obs.scale = 1, 
                                    var.scale = 1, 
                                    ellipse = TRUE, 
                                    groups = master.Table$group_th,
                                    point.size = 1,
                                    labels = NULL, 
                                    var.axes = FALSE,
                                    ellipse.prob = 0.95,
                                    ellipse.fill = FALSE, 
                                    ellipse.linewidth = 0.5) +
  theme_minimal() +
  guides(color=guide_legend("group_th"))

