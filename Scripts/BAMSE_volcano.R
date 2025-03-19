### plot volcano plots for th2 high and low DE genes from BAMSE cohort ###
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")

de.results.high.BAMSE <- read.csv('./th_high_th_low/Bamse_replication/DE_genes_THhigh_healthy.legacyT.csv')
de.results.low.BAMSE <- read.csv('./th_high_th_low/Bamse_replication/DE_genes_THlow_healthy.legacyT.csv')

sum(de.results.high.BAMSE[de.results.high.BAMSE$FDR < 0.05,]$logFC > 0)
sum(de.results.high.BAMSE[de.results.high.BAMSE$FDR < 0.05,]$logFC < 0)

# plot high 
png("./th_high_th_low/plots/BAMSE_TH2_high_genes_volcano.png", width=1000, height=900,res = 150)

volcano_high_bamse <- ggplot(data=de.results.high.BAMSE, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color= ifelse((FDR < 0.05)&(logFC>0),'firebrick', ifelse((FDR < 0.05)&(logFC<0),'blue','gray')))) + 
  theme_classic()+
  geom_text_repel(aes(label=ifelse(((FDR < 0.01 )& (abs(logFC)>2)) | (FDR < 0.000001),
                                   (ifelse((!is.na(external_gene_name)), external_gene_name, '')), ''),
                      lineheight=0.5, hjust= 0.5, vjust= 0.4),
                  max.overlaps =50,
                  direction = 'both',
                  alpha = 0.6,
                  box.padding=0.3,
                  point.padding=0.5)+
  xlab (expression (log[2]~fold~change))+
  ylab (expression(-log[10]~FDR))+
  scale_color_identity(name = '', breaks= c('firebrick', 'blue','') ,
                       labels = c('higher in T2-high',
                                  'lower in T2-high', 
                                  ''), guide = "legend")+
  xlim(-5, 5) + 
  ylim(-0.1, 15.5) +
  ggtitle("T2-high vs controls") +
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 15))
print(volcano_high_bamse)
dev.off()

# plot low
png("./th_high_th_low/plots/BAMSE_TH2_low_genes_volcano.png", width=1000, height=900,res = 150)

volcano_low_bamse<- ggplot(data=de.results.low.BAMSE, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color= ifelse((FDR < 0.05)&(logFC>0),'firebrick', ifelse((FDR < 0.05)&(logFC<0),'blue','gray')))) + 
  theme_classic()+
  geom_text_repel(aes(label=ifelse(((FDR < 0.01 )& (abs(logFC)>2)) | (FDR < 0.000001),
                                   (ifelse((!is.na(external_gene_name)), external_gene_name, '')), ''),
                      lineheight=0.5, hjust= 0.5, vjust= 0.4),
                  max.overlaps =50,
                  direction = 'both',
                  alpha = 0.6,
                  box.padding=0.3,
                  point.padding=0.5)+
  xlab (expression (log[2]~fold~change))+
  ylab (expression(-log[10]~FDR))+
  scale_color_identity(name = '', breaks= c('firebrick', 'blue','gray') ,
                       labels = c('higher in T2-low',
                                  'lower in T2-low', 
                                  'not significant (FDR>0.05)'), guide = "legend")+
  xlim(-5, 5) + 
  ylim(-0.1, 15.5) +
  ggtitle("T2-low vs controls") +
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 15))
print(volcano_low_bamse)
dev.off()
saveRDS(list(volcano_high_bamse, volcano_low_bamse), "./th_high_th_low/plots/BAMSE_volcano_list.rds")


figure <- ggarrange(volcano_high_bamse, volcano_low_bamse,
                    legend="bottom",
                    labels = c("A", "B"),
                    ncol = 2)
png("./th_high_th_low/plots/Bamse_TH2_high_low_genes_volcano.png",
    width = 3800, 
    height = 2000,
    res = 300)
print(figure)
dev.off()

