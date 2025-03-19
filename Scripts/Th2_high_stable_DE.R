## perform DE analysis of TH2-high stable asthma vs healthy
library(dplyr)
library(edgeR)
library(biomaRt)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")
master.Table <- read.csv("./th_high_th_low/master_table_th2_stable.csv") %>%
  filter(!is.na(stable_high))

# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

# Create a model matrix
design <- model.matrix(~0 + stable_high + age + gender + smoking.status, data = master.Table)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# Fit the model
DGEL <- edgeR::estimateDisp(DGEL, design)
#Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design)
contrasts <- limma::makeContrasts(
  highstable_healthy = stable_highYES - stable_highhealthy,
  levels = design
)

qlf  <- edgeR::glmQLFTest(fit, contrast = contrasts[,"highstable_healthy"])

## add gene names
ensembl38 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 108)

dataprobes <- getBM(filters = "ensembl_gene_id", attributes = c(
  "ensembl_gene_id",
  "external_gene_name"),
  values = as.character(rownames(expression.data)),
  mart = ensembl38
)

de.results.high <- edgeR::topTags(
  qlf,
  n=nrow(DGEL))$table %>% 
  tibble::rownames_to_column("Gene") %>%
  left_join(dataprobes, by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
          PValue = signif(PValue, digits = 3),
          FDR = signif(FDR, digits = 3))

write.csv(de.results.high,"./th_high_th_low/DE_genes_THhighstable_healthy.csv", row.names = F, quote = F)
#------------------
# plot volcano plot
#------------------

library(ggplot2)
library(ggrepel)

source ("/Users/tatiana/Work/RP2/ATLANTIS_project/source_scripts/DE_scripts.R")

de.results.high <- de.results.high %>%
  dplyr::rename(hgnc_symbol = external_gene_name)
# th2 high
png("./th_high_th_low/plots/TH2_high_stable_genes_volcano.png",width=900, height=900,res = 150)
volcano_low <- volcano.plot(de.results.high, de.results.high$logFC, de.results.high$FDR, "Th2 high stable vs healthy")
print(volcano_low)
dev.off()

