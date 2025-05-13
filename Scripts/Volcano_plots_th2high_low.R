### plot volcano plots for th2 high and low DE genes ###
library(ggplot2)
library(dplyr)
library(ggrepel)


setwd("/Users/tatiana/Work/RP2/ATLANTIS")
source ("/Users/tatiana/Work/RP2/ATLANTIS_project/source_scripts/DE_scripts.R")

# upload the result table
de.results.high <- read.csv('./th_high_th_low/DE_genes_THhigh_healthy.csv')

de.results.low <- read.csv('./th_high_th_low/DE_genes_THlow_healthy.csv')

de.results.high.low <- read.csv('./th_high_th_low/DE_genes_THhigh_THlow.csv')

master.Table <-read.csv('./th_high_th_low/master_table_th2.csv')

#------------------
# plot volcano plot
#------------------

# th2 high 
png("./th_high_th_low/plots/TH2_high_genes_volcano.png",width=900, height=900,res = 150)
# volcano_high <- volcano.plot(de.results.high, de.results.high$logFC, de.results.high$FDR, "Th2 high vs healthy")
# volcano_high_p_val <- volcano.plot.p.val.y(de.results.high, de.results.high$logFC, de.results.high$FDR, 
volcano_high <- ggplot(data=de.results.high, aes(x=logFC, y=-log10(FDR))) +
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
                       labels = c('higher in T2-high',
                                  'lower in T2-high', 
                                  'not significant (FDR>0.05)'), guide = "legend")+
  xlim(-4, 4) + 
  ggtitle("T2-high vs controls") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="plain"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.y  = element_blank())
print(volcano_high)
dev.off()

# th2 low

png("./th_high_th_low/plots/TH2_low_genes_volcano.png",width=900, height=900,res = 150)
volcano_low <-  ggplot(data=de.results.low, aes(x=logFC, y=-log10(FDR))) +
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
  #labs(title = "Th2 low vs healthy")+
  xlab (expression (log[2]~fold~change))+
  ylab (expression(-log[10]~FDR))+
  scale_color_identity(name = '', breaks= c('firebrick', 'blue','gray') ,
                       labels = c('higher in T2-low',
                                  'lower in T2-low',
                                  'not significant (FDR>0.05)'), guide = "legend")+
  xlim(-4, 4) + 
  ylim(0, 11.5) +
  ggtitle("T2-low vs controls") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, face = "plain"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.6),
        legend.text=element_text(size = 12),
        legend.key = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.y  = element_blank())


print(volcano_low)
dev.off()

saveRDS(list(volcano_high, volcano_low), "./th_high_th_low/plots/ATLANTIS_volcano_list.rds")


figure <- ggarrange(volcano_high, volcano_low,
                    legend="bottom",
                    labels = c("A", "B"),
                    ncol = 2)

png("./th_high_th_low/plots/ATLANTIS_TH2_high_low_genes_volcano.png",
    width = 3800, 
    height = 2000,
    res = 300)
print(figure)
dev.off()
