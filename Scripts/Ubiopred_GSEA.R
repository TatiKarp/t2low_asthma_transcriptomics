## This script will check overlapped genes identified in ATLANTIS and UBIOPRED cohort when comparing T2-high vs healthy, T2-low vs healthy ###
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(fgsea)



setwd("/Users/tatiana/Work/RP2/ATLANTIS")


# UBIOPRED result table:
de.results.high.UBiopred <- read.csv('./th_high_th_low/Ubiopred_replication/All_T2_High_vs_Healthy_Limma_DGE.csv')
de.results.low.UBiopred <- read.csv('./th_high_th_low/Ubiopred_replication/All_T2_Low_vs_Healthy_Limma_DGE.csv')

# ATLANTIS result table:
de.results.high.ATL <- read.csv('./th_high_th_low/DE_genes_THhigh_healthy.csv') %>%
  filter(external_gene_name != "") ## keep only gene with gene symbol names 
de.results.low.ATL <- read.csv('./th_high_th_low/DE_genes_THlow_healthy.csv')

## 1. ATLANTIS as background
# determine gene sets (up and down)
geneSets_high <- list(Up_genes = de.results.high.UBiopred[de.results.high.UBiopred$logFC>0 & de.results.high.UBiopred$adj.P.Val<0.05,]$Gene,
                      Down_genes = de.results.high.UBiopred[de.results.high.UBiopred$logFC<0 & de.results.high.UBiopred$adj.P.Val<0.05,]$Gene)

high_geneSet_up <- list(Up_genes = de.results.high.UBiopred[de.results.high.UBiopred$logFC>0 & de.results.high.UBiopred$adj.P.Val<0.05,]$Gene)
high_geneSet_down <- list(Down_genes = de.results.high.UBiopred[de.results.high.UBiopred$logFC<0 & de.results.high.UBiopred$adj.P.Val<0.05,]$Gene)

# Up_genes  : 57
# Down_genes: 33

de.results.high.ATL[de.results.high.ATL$FDR < 0.05, "external_gene_name"]
de.results.high.UBiopred[de.results.high.UBiopred$adj.P.Val < 0.05, "GeneSymbol"]

intersect(de.results.high.UBiopred[de.results.high.UBiopred$adj.P.Val < 0.05, "GeneSymbol"],de.results.high.ATL[de.results.high.ATL$FDR < 0.05, "external_gene_name"])
# 11 genes

########## GSEA  #############
## U-BIOPRED as background:

# define geneSets: 
geneSets_high <- list(Up_genes = de.results.high.ATL[de.results.high.ATL$logFC>0 & de.results.high.ATL$FDR<0.05,]$external_gene_name,
                      Down_genes = de.results.high.ATL[de.results.high.ATL$logFC<0 & de.results.high.ATL$FDR<0.05,]$external_gene_name)

# define background
ordered.de.results.high.UBiopred <- de.results.high.UBiopred %>%
  mutate(`-10logFDR` = ifelse(logFC>0, 1,-1)*-log10(adj.P.Val)) %>%
  mutate(`-10logPVal` = ifelse(logFC>0, 1,-1)*-log10(P.Value)) %>%
  dplyr::arrange(desc(`-10logFDR`), P.Value) 

Ubio_high_ranks <- c(ordered.de.results.high.UBiopred$`-10logPVal`)

names(Ubio_high_ranks) <- ordered.de.results.high.UBiopred$GeneSymbol

# two tailed GSEA
fgseaRes_both <-fgsea(pathways = geneSets_high,
                      stats = Ubio_high_ranks,
                      minSize=15,
                      maxSize=2000,
                      nPermSimple = 10000)

# plot
source("/Users/tatiana/Work/RP2/ATLANTIS_project/source_scripts/GSEA_enrichment_plot.R")

up_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", NES], 2)
up_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", padj], 2)
up_ubio <- plotEnr(geneSets_high[['Up_genes']], Ubio_high_ranks, 
              linecol =  "black", bincol =  "#8A9795")+ ##B4251A"
  ggtitle(paste0('NES=', up_NES, ', padj=', up_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.high.UBiopred$`-10logFDR` > 0)) +
  ylab("T2-high higher genes")

down_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", NES], 2)
down_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", padj], 2)
down_ubio <- plotEnr(geneSets_high[['Down_genes']], Ubio_high_ranks, 
                linecol =  "black", bincol =  "#8A9795")+ ##163A7D
  ggtitle(paste0('NES=', down_NES, ', padj=', down_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.high.UBiopred$`-10logFDR` < 0)) +
  ylab("T2-high lower genes")

all_plt <- ggarrange(up_ubio, down_ubio,
                     nrow = 2, 
                     labels = c("A",  "B"))
png("./th_high_th_low/plots/GSEA_UBIOPRED_backgrouns_ATL_set.png",
    width=1400, height=850, res = 200)

print(all_plt)
dev.off()


saveRDS(list(up_ubio, down_ubio), "./th_high_th_low/plots/GSEA_UBIOPRED_backgrouns_ATL_set.rds")





