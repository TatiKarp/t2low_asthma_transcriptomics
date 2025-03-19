# This script will perform pathways analysis for the gene modules that were found to be associated to T2-low asthma by WGCNA
library(clusterProfiler)
library(dplyr)
library(biomaRt)
library('org.Hs.eg.db')
library(enrichplot)


setwd('/Users/tatiana/Work/RP2/ATLANTIS')

genes_in_modules <- read.csv('./th_high_th_low/WGCNA_th2low/Genes_in_sign_modules.csv')
all_genes <- read.csv('./th_high_th_low/WGCNA_th2low/wgcna.modules.csv')

# create fucntion for GO
GO_all <- function(names_genes_in_modules, universe_genes, module) {
  pathways <- enrichGO(gene = as.character(names_genes_in_modules[names_genes_in_modules$`bwnet.colors` == module, ]$X),
                       universe = as.character(universe_genes),
                       OrgDb = 'org.Hs.eg.db',
                       keyType = "ENSEMBL",
                       ont = 'ALL',
                       pvalueCutoff = 0.05, 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05,
                       minGSSize = 10, 
                       maxGSSize = 10000)
  return(pathways)
}  
 
                  universe_genes = universe_genes)

GO_all <- lapply(c('greenyellow', 'black', 'purple'), GO_all, names_genes_in_modules = genes_in_modules,
                 universe_genes = all_genes$X)

# save datasets
modules <- c('greenyellow', 'black', 'purple')
# for( i in 1:length(KEGG_pathways)){
#   write.csv(KEGG_pathways[[i]]@result, paste0('./th_high_th_low/WGCNA_th2low/KEGG_pathways_in_module_',modules[i],'.csv'),
#             row.names = F)
# }

for( i in 1:length(GO_all)){
  write.table(GO_all[[i]]@result %>%
                mutate (pvalue = signif (pvalue, digits = 2),
                        p.adjust = signif(p.adjust, digits = 2),
                        qvalue = signif(qvalue, digits = 2)), paste0('./th_high_th_low/WGCNA_th2low/GO_pathways_in_module_', modules[i], '.csv'),
            row.names = F, sep = "\t", quote = F)
}



plot_GO <- function(GO_object){
  GO_BP2_up <- clusterProfiler::simplify(GO_object, cutoff=0.7, by="p.adjust", select_fun=min)
  #View(head(GO_BP2_up@result[, -8], 50))
    ## visualization
  ego2_BP_up <- pairwise_termsim(GO_BP2_up)
  p2_PB_up <- emapplot(ego2_BP_up, cex_label_category=.8, cex_line=.5) + coord_cartesian()
  p2_PB_up <- p2_PB_up + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                                               guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
  png(paste0('./th_high_th_low/WGCNA_th2low/plot_GO_pathways_in_module_',i,'.png'))
  a <- cowplot::plot_grid(p2_PB_up, rel_widths=c(1, 1.2))                     
  print(a)
  dev.off()
}

for( i in 1:length(GO_all)){
  plot_GO(GO_all[[i]])
}

## check 
pathways <- enrichGO(gene = as.character(genes_in_modules[genes_in_modules$`bwnet.colors` == 'black', ]$X),
                     universe = as.character(all_genes$X),
                     OrgDb = 'org.Hs.eg.db',
                     keyType = "ENSEMBL",
                     ont = 'ALL',
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 10000)
df <- pathways@result

