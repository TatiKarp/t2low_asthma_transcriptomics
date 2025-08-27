# This script will create figure 1 for the paper - combining figures together 
library(patchwork)
library(ggplot2)
library(gridExtra)


setwd("~/Work/RP2/ATLANTIS")

sankey_plt <- readRDS("./th_high_th_low/plots/t2stable_sankey_plt_only_RNA.rds")
volcano_list <- readRDS("./th_high_th_low/plots/ATLANTIS_volcano_list.rds")
cell_types <- readRDS("./th_high_th_low/plots/Cell_types_th2_high_low_CIBERSORT.rds")
wgcna_modules <- readRDS("./th_high_th_low/plots/Eigengene_values_all_modules.rds")
distances_bxplt <- readRDS("./th_high_th_low/plots/bxplt_median_dist_all_genes.rds")


# Arrange first 3 figures together 
combined_plot_123 <- 
  sankey_plt /
  ((volcano_list[[1]]) | (volcano_list[[2]])) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
        plot.subtitle = element_text(size = 14))

png("./th_high_th_low/plots/combined_figure_1ABC.png",
    width = 3000, height = 2500,
    res = 300)
print(combined_plot_123)
dev.off()

# Arrange figure 4 and 5 together +distances 
combined_plot_45 <- 
  ((distances_bxplt | cell_types) +
     plot_layout(widths = c(1, 3))) /
  wgcna_modules +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
        plot.subtitle = element_text(size = 14),
        plot.tag = element_text(face = "bold"))

png("./th_high_th_low/plots/combined_figure_2ABC.png",
    width = 3000, height = 3500,
    res = 300)
print(combined_plot_45)
dev.off()



### all together, old version 
combined_plot <-  # First row
  ((volcano_list[[1]]) | (volcano_list[[2]])) /  # Second row (side by side)
  (sankey_plt | cell_types) / wgcna_modules +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
        plot.subtitle = element_text(size = 14))

png("./th_high_th_low/plots/combined_figure_1.png",
    width = 3500, height = 4700,
    res = 300)
print(combined_plot)
dev.off()
