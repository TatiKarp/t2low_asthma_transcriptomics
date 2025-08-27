# This script will perform pathways analysis for the gene modules that were found to be associated to T2-low asthma by WGCNA
library(clusterProfiler)
library(dplyr)
library(biomaRt)
library('org.Hs.eg.db')
library(enrichplot)
library(ggplot2)


setwd("~/Work/RP2/ATLANTIS")

genes_in_modules <- read.csv('./th_high_th_low/WGCNA_th2low/Genes_in_sign_modules.csv')
all_genes <- read.csv('./th_high_th_low/WGCNA_th2low/wgcna.modules.csv')
## Add gene names
#ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# KEGG_all <- function(names_genes_in_modules, universe_genes, module) {
#   pathways <- enrichKEGG(gene = as.character(names_genes_in_modules[names_genes_in_modules$`bwnet$colors` == module, ]$entrezgene_id),
#                          organism = 'hsa',
#                          universe = as.character(universe_genes$entrezgene_id),
#                          pvalueCutoff = 0.05,
#                          pAdjustMethod = "BH",
#                          minGSSize = 5,
#                          qvalueCutoff = 1)
#   return(pathways)
# }
# create fucntion for GO
GO_all <- function(names_genes_in_modules, universe_genes, module) {
  pathways <- clusterProfiler::enrichGO(gene = as.character(names_genes_in_modules[names_genes_in_modules$`bwnet.colors` == module, ]$X),
                       universe = as.character(universe_genes),
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = 'ALL',
                       pvalueCutoff = 0.05, 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05,
                       minGSSize = 10, 
                       maxGSSize = 10000)
  return(pathways)
}  


# universe_genes <- getBM(
#   filters="ensembl_gene_id",
#   attributes=c("ensembl_gene_id", "entrezgene_id"),
#   values=all_genes$X,
#   mart=ensembl)
# 
# names_genes_in_modules <- getBM(
#   filters="ensembl_gene_id",
#   attributes=c("ensembl_gene_id", "entrezgene_id"),
#   values=genes_in_modules$X,
#   mart=ensembl) %>%
#   left_join(genes_in_modules%>%
#               tibble::rownames_to_column('ensembl_gene'), by = c('ensembl_gene_id'='X'))


### Calculate OR analysis for GO and KEGG pathways using lapply
# KEGG_pathways <- lapply(c('greenyellow', 'black', 'purple'), KEGG_all, names_genes_in_modules =names_genes_in_modules,
#                         universe_genes = universe_genes)

# GO_all <- lapply(c('greenyellow', 'black', 'purple'), GO_all, names_genes_in_modules =names_genes_in_modules,
#                   universe_genes = universe_genes)

GO_all_result <- lapply(c('greenyellow', 'black', 'purple'), GO_all, names_genes_in_modules = genes_in_modules,
                 universe_genes = all_genes$X)

# save datasets
modules <- c('greenyellow', 'black', 'purple')
# for( i in 1:length(KEGG_pathways)){
#   write.csv(KEGG_pathways[[i]]@result, paste0('./th_high_th_low/WGCNA_th2low/KEGG_pathways_in_module_',modules[i],'.csv'),
#             row.names = F)
# }

# save tables:

for( i in 1:length(GO_all_result)){
  write.table(GO_all_result[[i]]@result %>%
                mutate (pvalue = signif (pvalue, digits = 2),
                        p.adjust = signif(p.adjust, digits = 2),
                        qvalue = signif(qvalue, digits = 2)), paste0('./th_high_th_low/WGCNA_th2low/GO_pathways_in_module_', modules[i], '.csv'),
            row.names = F, sep = "\t", quote = F)
}


## Plot dotplots: 
plot_GO <- function(c){
  GO_one_module_plt <- enrichplot::dotplot(c)
  png(paste0('./th_high_th_low/WGCNA_th2low/dotplot_GO_pathways_in_module_', modules[i],'.png'),
             res = 300,  width = 1880, height = 1880)
  print(GO_one_module_plt)
  dev.off()
}

for( i in 1:length(GO_all_result)){
  if (dim(GO_all_result[[i]])[1] != 0) {
    
    plot_GO(GO_all_result[[i]]) 
  }
}


## plot pathways - this is not working
## use simplify to remove redundant terms

# GO_BP2_up <- clusterProfiler::simplify(GO_all[[1]], cutoff=0.7, by="p.adjust", select_fun=min)
# #View(head(GO_BP2_up@result[, -8], 50))
# 
# ## visualization
# ego2_BP_up <- pairwise_termsim(GO_BP2_up)
# p2_PB_up <- emapplot(ego2_BP_up, cex_label_category=.8, cex_line=.5) + coord_cartesian()
# p2_PB_up <- p2_PB_up + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
#                                              guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
# a <- cowplot::plot_grid(p2_PB_up, rel_widths=c(1, 1.2))                     
# 
# 
# plot_GO <- function(c){
#   GO_BP2_up <- clusterProfiler::simplify(GO_object, cutoff=0.7, by="p.adjust", select_fun=min)
#   #View(head(GO_BP2_up@result[, -8], 50))
#     ## visualization
#   ego2_BP_up <- pairwise_termsim(GO_BP2_up)
#   plt <- enrichplot::goplot(ego2_BP_up)
#   
#   p2_PB_up <- emapplot(ego2_BP_up, cex_label_category=.8, cex_line=.5) + coord_cartesian()
#   p2_PB_up <- p2_PB_up + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
#                                                guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
#   a <- cowplot::plot_grid(p2_PB_up, rel_widths=c(1, 1.2)) 
#   png(paste0('./th_high_th_low/WGCNA_th2low/plot_GO_pathways_in_module_', modules[i],'.png',
#              res = 300,  width = 1880, height = 1880))
#   print(a)
#   dev.off()
# }
# 
# for( i in 1:length(GO_all)){
#   plot_GO(GO_all[[i]])
# }


## add reactome:

library(ReactomePA)


all_gene_to_entrtzid <- clusterProfiler::bitr(all_genes$X, fromType = "ENSEMBL" , toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')

genes_in_modules_entrez <- genes_in_modules %>%
  left_join(all_gene_to_entrtzid, by = c("X" = "ENSEMBL"))

reactome_ORA_black <- enrichPathway(
  gene = genes_in_modules_entrez[genes_in_modules_entrez$`bwnet.colors` == 'black', ]$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe = all_gene_to_entrtzid$ENTREZID,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE)

reactome_ORA_black_plt <- enrichplot::dotplot(reactome_ORA_black)
png('./th_high_th_low/WGCNA_th2low/dotplot_reactome_pathways_in_module_black.png',
    res = 300,  width = 1880, height = 1880)
print(reactome_ORA_black_plt)
dev.off()

df_reactome_ORA_black <- reactome_ORA_black@result %>%
  filter(p.adjust < 0.05)
write.csv(df_reactome_ORA_black %>%
            mutate (pvalue = signif (pvalue, digits = 2),
                    p.adjust = signif(p.adjust, digits = 2),
                    qvalue = signif(qvalue, digits = 2)), 
          paste0('./th_high_th_low/WGCNA_th2low/REACTOME_pathways_in_module_',"black",'.csv'), row.names = F)

reactome_ORA_purple <- enrichPathway(
  gene = genes_in_modules_entrez[genes_in_modules_entrez$`bwnet.colors` == 'purple', ]$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe = all_gene_to_entrtzid$ENTREZID,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE)

reactome_ORA_purple_plt <- enrichplot::dotplot(reactome_ORA_purple)
png('./th_high_th_low/WGCNA_th2low/dotplot_reactome_pathways_in_module_purple.png',
    res = 300,  width = 1880, height = 1880)
print(reactome_ORA_purple_plt)
dev.off()

df_reactome_ORA_purple <- reactome_ORA_purple@result %>%
  filter(p.adjust < 0.05)
write.csv(df_reactome_ORA_purple %>%
            mutate (pvalue = signif (pvalue, digits = 2),
                    p.adjust = signif(p.adjust, digits = 2),
                    qvalue = signif(qvalue, digits = 2)),
          paste0('./th_high_th_low/WGCNA_th2low/REACTOME_pathways_in_module_',"purple",'.csv'), row.names = F)

reactome_ORA_greenyellow <- enrichPathway(
  gene = genes_in_modules_entrez[genes_in_modules_entrez$`bwnet.colors` == 'greenyellow', ]$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe = all_gene_to_entrtzid$ENTREZID,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE)

reactome_ORA_yellow_plt <- enrichplot::dotplot(reactome_ORA_greenyellow)
png('./th_high_th_low/WGCNA_th2low/dotplot_reactome_pathways_in_module_greenyellow.png',
    res = 300,  width = 1880, height = 1880)
print(reactome_ORA_yellow_plt)
dev.off()

df_reactome_ORA_greenyellow <- reactome_ORA_greenyellow@result %>%
  filter(p.adjust < 0.05)
write.csv(df_reactome_ORA_greenyellow %>%
            mutate (pvalue = signif (pvalue, digits = 2),
                    p.adjust = signif(p.adjust, digits = 2),
                    qvalue = signif(qvalue, digits = 2)),
          paste0('./th_high_th_low/WGCNA_th2low/REACTOME_pathways_in_module_',"greenyellow",'.csv'), row.names = F)

