### This script will check overlapped genes identified in ATLANTIS and BAMSE cohort when comparing T2-high vs healthy, T2-low vs healthy ###
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(fgsea)



setwd("/Users/tatiana/Work/RP2/ATLANTIS")

# Bamse result table:
de.results.high.BAMSE <- read.csv('./th_high_th_low/Bamse_replication/DE_genes_THhigh_healthy.legacyT.csv')
de.results.low.BAMSE <- read.csv('./th_high_th_low/Bamse_replication/DE_genes_THlow_healthy.legacyT.csv')

# ATLANTIS result table:
de.results.high.ATL <- read.csv('./th_high_th_low/DE_genes_THhigh_healthy.csv')
de.results.low.ATL <- read.csv('./th_high_th_low/DE_genes_THlow_healthy.csv')

######### GSEA ##########
## BAMSE as background

# define geneSets: 
geneSets_high <- list(Up_genes = de.results.high.ATL[de.results.high.ATL$logFC>0 & de.results.high.ATL$FDR<0.05,]$Gene,
                      Down_genes = de.results.high.ATL[de.results.high.ATL$logFC<0 & de.results.high.ATL$FDR<0.05,]$Gene)

# define background

ordered.de.results.high.BAMSE <- de.results.high.BAMSE%>%
  mutate(`-10logFDR` = ifelse(logFC>0, 1,-1)*-log10(FDR)) %>%
  mutate(`-10logPVal` = ifelse(logFC>0, 1,-1)*-log10(PValue)) %>%
  dplyr::arrange(desc(`-10logFDR`), PValue) 

BAMSE_high_ranks <- c(ordered.de.results.high.BAMSE$`-10logPVal`)

names(BAMSE_high_ranks) <- ordered.de.results.high.BAMSE$Gene

# GSEA two tailed
fgseaRes_both <-fgsea(pathways = geneSets_high,
                      stats = BAMSE_high_ranks,
                      minSize=15,
                      maxSize=2000,
                      eps = 0,
                      nPermSimple = 10000)

# plot 
source("/Users/tatiana/Work/RP2/ATLANTIS_project/source_scripts/GSEA_enrichment_plot.R")

up_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", NES], 2)
up_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", padj], 2)
up_bamse <- plotEnr(geneSets_high[['Up_genes']], BAMSE_high_ranks, 
              linecol =  "black", bincol =  "#76807E")+ #B4251A
  ggtitle(paste0('NES=', up_NES, ', padj=', up_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.high.BAMSE$`-10logFDR` > 0)) +
  ylab("T2-high higher genes")

down_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", NES], 2)
down_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", padj], 2)
down_bamse <- plotEnr(geneSets_high[['Down_genes']], BAMSE_high_ranks, 
                linecol =  "black", bincol =  "#76807E")+ #163A7D
  ggtitle(paste0('NES=', down_NES, ', padj=', down_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.high.BAMSE$`-10logFDR` < 0)) +
  ylab("T2-high lower genes")

all_plt <- ggarrange(up_bamse, down_bamse,
                     nrow = 2, 
                     labels = c("A",  "B"))
png("./th_high_th_low/plots/GSEA_BAMSE_background_ATL_set.png",
    width=1400, height=850, res = 200)

print(all_plt)
dev.off()

# save as r object 
saveRDS(list(up_bamse, down_bamse), "./th_high_th_low/plots/GSEA_BAMSE_backgrouns_ATL_set.rds")

