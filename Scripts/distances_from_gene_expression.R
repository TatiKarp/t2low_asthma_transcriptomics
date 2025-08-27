# this script will calculate distances on gene expression level between t2-high/low/healthy subjects
#library(ggrepel)
library(edgeR)
library(dplyr)
library(ggplot2)
#libraryggplot2#library(ggpubr)

setwd("~/Work/RP2/ATLANTIS")

# 1. prepare data
master.Table <- read.csv("./th_high_th_low/master_table_th2.csv", header = TRUE) %>%
  filter (QC_check == 'YES') 

expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

# 1.1 normalize data 
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
cpm_norm <- cpm(DGEL, normalized.lib.sizes = TRUE, log = F) # dim 12134


# 2. save all the distances - for boxplots
distance_data <- data.frame(
  PC = character(),
  Group = character(),
  Distance = numeric(),
  stringsAsFactors = FALSE
)

# 2.1 Everyone vs everyone (pairwise), subjects in the rows!
distance_matrix_all <- as.matrix(dist(t(cpm_norm), method = "euclidean"))
distances_all <- distance_matrix_all[lower.tri(distance_matrix_all)]
distance_data <- rbind(
  distance_data,
  data.frame(
    Group = "everyone",
    Distance = distances_all,
    stringsAsFactors = FALSE
  )
)
  

# 2.2 Loop over groups
groups <- c("low", "high", "healthy", "undeterm")
for (group in groups) {
  group_values <- cpm_norm[,master.Table$group_th == group]

  # calculate pairwise distances, subjects in the rows 
  distance_matrix <- as.matrix(dist(t(group_values), method = "euclidean"))
  distances <- distance_matrix[lower.tri(distance_matrix)]
  
  distance_data <- rbind(
    distance_data,
    data.frame(
      Group = group,
      Distance = distances,
      stringsAsFactors = FALSE
    )
  )
}

distance_data <- distance_data %>%
  filter(!(Group %in% c("everyone", "undeterm"))) %>%
  mutate(Group = recode(Group, 
                        high = "T2-high",
                        healthy = "Healthy", 
                        low = "T2-low")) 

distance_data$Group <- factor(distance_data$Group, levels = c("Healthy", "T2-low", "T2-high"))
  
# 3 Plot distances:
# calculate p_value:
stat.test <- distance_data %>%
  #filter(Group %in% c("high", "low")) %>%
  rstatix::wilcox_test(Distance ~ Group) %>%
  rstatix::add_y_position(step.increase = 0.5, fun = "median") %>%
  mutate(y.position = y.position + 45000)

bxplt <- ggplot(distance_data, 
                       aes(x = Group, y = Distance)) +
  geom_boxplot(aes(fill = Group), outliers = FALSE) +
  scale_fill_manual(values = c(
    "Healthy" = "#66C2A5",
    "T2-low" = "#8DA0CB",
    "T2-high" = "#FFD92F"
  )) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(y = expression(bold("Inter-sample gene expression dissimilarity"))) +
  #labs(title = "Median Euclidean Distances, all highly exressed gene", y = "Distance", x = "Group") +
  ggpubr::stat_pvalue_manual(
    data = stat.test,
    y.position = "y.position",
    label = "p = {p}",
    size = 3.5,
    remove.bracket = FALSE
  ) +
  theme(axis.text=element_text(size = 15),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        legend.position = "none",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size = 12),
        legend.key = element_rect(fill = "white"))

saveRDS(bxplt, "./th_high_th_low/plots/bxplt_median_dist_all_genes.rds")

png("./th_high_th_low/plots/bxplt_median_dist_all_genes.png",
    width = 1400, height = 1200, res = 200)
print(bxplt)
dev.off()


