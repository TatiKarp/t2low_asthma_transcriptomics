## check GO pathways using GSEA for T2-low non significant genes: 

library(clusterProfiler)
library(ReactomePA)


setwd("~/Work/RP2/ATLANTIS")

# load T2-low genes 
de.results.low <- read.csv("./th_high_th_low/DE_genes_THlow_healthy.csv") 

# create rank
de.results.low <- de.results.low %>%
  mutate(`-10logPval` = ifelse(logFC>0, 1,-1)*-log10(PValue)) %>%
  arrange(desc (`-10logPval`))

ranks <- de.results.low$`-10logPval`
names(ranks) <- de.results.low$Gene

## 1. Gene ontology 
gsea_GO <- clusterProfiler::gseGO(ranks,
                                  ont = "ALL",
                                  keyType = "ENSEMBL",
                                  OrgDb = 'org.Hs.eg.db',
                                  pvalueCutoff = 0.05,
                                  minGSSize = 10,
                                  maxGSSize = 10000)

gsea_GO_simp <- clusterProfiler::simplify(gsea_GO, cutoff=0.7, by="p.adjust", select_fun=min)
gsea_GO_df <- gsea_GO_simp@result %>%
  arrange(p.adjust) %>%
  dplyr::mutate_at(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), ~signif(., 3))

write.csv(gsea_GO_df, "./th_high_th_low/gsea_go_t2low.csv", quote = F)

# plot
gsea_GO_simp_plt <- enrichplot::dotplot(gsea_GO_simp)
png('./th_high_th_low/dotplot_gsea_go_t2low.png',
    res = 300,  width = 1880, height = 1880)
print(gsea_GO_simp_plt)
dev.off()

## 2. reactome
all_gene_to_entrtzid <- clusterProfiler::bitr(de.results.low$Gene, fromType = "ENSEMBL" , toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')

# 13.73% of input gene IDs are fail to map
 
de.results.low.etrez <- de.results.low %>%
  left_join(all_gene_to_entrtzid, by = c("Gene" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID))

ranks <- de.results.low.etrez$`-10logPval`
names(ranks) <- de.results.low.etrez$ENTREZID


gsea_reactome <- ReactomePA::gsePathway(geneList = ranks,
                       organism = "human",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       by = 'fgsea')
gsea_reactome_df <- gsea_reactome@result %>%
  arrange(p.adjust) %>%
  dplyr::mutate_at(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), ~signif(., 3))

write.csv(gsea_reactome_df, "./th_high_th_low/gsea_reactome_t2low.csv", quote = F)

# plot
gsea_reactome_simp_plt <- enrichplot::dotplot(gsea_reactome)
png('./th_high_th_low/dotplot_gsea_reactome_t2low.png',
    res = 300,  width = 1880, height = 1880)
print(gsea_reactome_simp_plt)
dev.off()