# this plot will combine figure to create the final replication plot
library(ggplot2)
library(ggpubr)

## combine with U-biopred
setwd("./")
volcano_ubiopred_list <- readRDS("./th_high_th_low/plots/UBIOPRED_volcano_list.rds")
volcano_bamse_list <- readRDS("./th_high_th_low/plots/BAMSE_volcano_list.rds")

gsea_Ubio_list<- readRDS("./th_high_th_low/plots/GSEA_UBIOPRED_backgrouns_ATL_set.rds")
gsea_BAMSE_list <- readRDS("./th_high_th_low/plots/GSEA_BAMSE_backgrouns_ATL_set.rds")

volcano_bamse <-  ggarrange(volcano_bamse_list[[1]] , volcano_bamse_list[[2]],
                            nrow = 1, ncol = 2)  %>%
  annotate_figure(top = text_grob("BAMSE cohort", face = "bold", size = 15),
                  left = text_grob((expression(-log[10]~FDR)), rot = 90, size = 15))        
                            

volcano_ubio <-  ggarrange(volcano_ubiopred_list[[1]] , volcano_ubiopred_list[[2]],
                            nrow = 1, ncol = 2)  %>%
  annotate_figure(top = text_grob("U-BIOPRED cohort", face = "bold", size = 15),
                  left = text_grob((expression(-log[10]~FDR)), rot = 90, size = 15))        


all_volcano_BAMSE_Ubiopred <- ggarrange(volcano_bamse, volcano_ubio,
                                    nrow = 2, ncol = 1,
                                    labels = c('A', 'B')) 

gsea_BAMSE <- ggarrange(gsea_BAMSE_list[[1]], gsea_BAMSE_list[[2]],
                        nrow = 1, ncol = 2)  %>%
  annotate_figure(top = text_grob("BAMSE cohort", face = "bold", size = 15)) 

gsea_ubio <- ggarrange(gsea_Ubio_list[[1]], gsea_Ubio_list[[2]],
                        nrow = 1, ncol = 2)  %>%
  annotate_figure(top = text_grob("U-BIOPRED cohort", face = "bold", size = 15)) 


gsea_all <- ggarrange(gsea_BAMSE, gsea_ubio,
                      nrow = 2, ncol = 1) %>%
  annotate_figure(left = text_grob("enrichment score", rot = 90, size = 15))  

all_volcano_gsea_BAMSE_Ubiopred <- ggarrange(volcano_bamse, volcano_ubio,
                                        gsea_all,
                                        nrow = 3, ncol = 1,
                                        labels = c("A", "B", "C")) 

png("./th_high_th_low/plots/volcano_GSEA_BAMSE_Ubiopred.png",
    width = 2800, height = 4700,
    res = 300)
print(all_volcano_gsea_BAMSE_Ubiopred)
dev.off()

## Combine plots using patchwork

library(patchwork)

combined_plot <- (volcano_bamse / volcano_ubio / gsea_all) +                   # Arrange plots in a row
  plot_annotation(tag_levels = 'A'# Automatically labels plots as A, B, C
  ) & 
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
        plot.subtitle = element_text(size = 14))  


png("./th_high_th_low/plots/patchwork_volcano_GSEA_BAMSE_Ubiopred.png",
    width = 2800, height = 4700,
    res = 300)
print(combined_plot)
dev.off() # the same!