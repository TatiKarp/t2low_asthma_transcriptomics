library(dplyr)
library(edgeR)
library(biomaRt)
library(fgsea)


setwd("Work/RP2")
# 1. data load

master.Table <- read.csv('./ATLANTIS/th_high_th_low/master_table_th2.csv')

# add additional data 
big_master_table <-  read.csv('./ATLANTIS/atlantis_patient_data.csv', header = TRUE)%>%
  dplyr::select(c(PT, ICS, ICS_LABA, 
                  ICS_DDOSE_EQ, ICS_LABA_DDOSE_EQ))

big_master_table <- big_master_table[!duplicated(big_master_table$PT), ]

master.Table.TH.cor <- master.Table %>%
  filter(asthma.status == 'A') %>%
  left_join(big_master_table, by = c('PT'='PT')) %>%
    mutate(any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes")) %>%
  mutate(ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                                coalesce(ICS_DDOSE_EQ,0) + coalesce(ICS_LABA_DDOSE_EQ, 0)))
master.Table <- master.Table%>%
  left_join(master.Table.TH.cor%>%
              dplyr::select(c(GenomeScan_ID, any_ICS, ICS_dose_sum))) %>%
  mutate(ICS_dose_sum = if_else(is.na(ICS_dose_sum), 0, ICS_dose_sum))

# add expression data
expression.data <- read.csv('./ATLANTIS/Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()


# add gene names
ensembl38 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 108)

dataprobes <- getBM(filters = "ensembl_gene_id", attributes = c(
  "ensembl_gene_id",
  "external_gene_name"),
  values = as.character(rownames(expression.data)),
  mart = ensembl38
)


# 2. Run analysis by correcting for ICS
design <- model.matrix(~0 + group_th + age + gender + smoking.status + ICS_dose_sum, data = master.Table)
# Create an edgeR object, filter low expressed genes, normalize
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# Fit the model
DGEL <- edgeR::estimateDisp(DGEL, design)
# Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE)
# legacy = TRUE: for edgse 4.0 and above- to reproduce the results of edgeR 3.0

contrasts <- limma::makeContrasts(
  high_healthy = group_thhigh - group_thhealthy,
  low_healthy = group_thlow - group_thhealthy,
  high_low = group_thhigh - group_thlow,
  levels = design
)

qlf_high   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_healthy"])
qlf_low <- edgeR::glmQLFTest(fit, contrast = contrasts[,"low_healthy"])
qlf_high_low <- edgeR::glmQLFTest(fit, contrast = contrasts[,"high_low"])

summary(decideTests(qlf_high)) 
# -1*group_thhealthy 1*group_thhigh
# Down                                  63
# NotSig                             18851
# Up                                   153
summary(decideTests(qlf_low)) 
# -1*group_thhealthy 1*group_thlow
# Down                                  0
# NotSig                            19067
# Up                                    0


de.results.t2high.ICScorrected <- edgeR::topTags(
  qlf_high,
  n=nrow(DGEL))$table %>% 
  tibble::rownames_to_column("Gene") %>%
  left_join(dataprobes, by = c("Gene" = "ensembl_gene_id"))

## 3. Compare with non corrected results: 
# not corrected result table:
de.results.t2high<- read.csv("./ATLANTIS/th_high_th_low/DE_genes_THhigh_healthy.csv")

# define geneSets: 
geneSets_high_corrected <- list(Up_genes = de.results.t2high.ICScorrected[de.results.t2high.ICScorrected$logFC>0 & de.results.t2high.ICScorrected$FDR<0.05,]$Gene,
                      Down_genes = de.results.t2high.ICScorrected[de.results.t2high.ICScorrected$logFC<0 & de.results.t2high.ICScorrected$FDR<0.05,]$Gene)

## define background - uncorrected

ordered.de.results.t2high <- de.results.t2high %>%
  mutate(`-10logFDR` = ifelse(logFC>0, 1,-1)*-log10(FDR)) %>%
  mutate(`-10logPVal` = ifelse(logFC>0, 1,-1)*-log10(PValue)) %>%
  dplyr::arrange(desc(`-10logFDR`), PValue) 

t2high_ranks <- c(ordered.de.results.t2high$`-10logPVal`)

names(t2high_ranks) <- ordered.de.results.t2high$Gene

# GSEA two tailed
fgseaRes_both <-fgsea(pathways = geneSets_high_corrected,
                      stats = t2high_ranks,
                      minSize=15,
                      maxSize=2000,
                      eps = 0,
                      nPermSimple = 10000)

# plot 
source("./ATLANTIS_project/source_scripts/GSEA_enrichment_plot.R")

up_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", NES], 2)
up_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Up_genes", padj], 2)
up_corrected <- plotEnr(geneSets_high_corrected[['Up_genes']], t2high_ranks, 
                    linecol =  "black", bincol =  "#76807E")+ #B4251A
  ggtitle(paste0('NES=', up_NES, ', padj=', up_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.t2high$`-10logFDR` > 0)) +
  ylab("")

down_NES = round(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", NES], 2)
down_padj = signif(fgseaRes_both[fgseaRes_both$pathway == "Down_genes", padj], 2)
down_corrected <- plotEnr(geneSets_high[['Down_genes']], t2high_ranks, 
                      linecol =  "black", bincol =  "#76807E")+ #163A7D
  ggtitle(paste0('NES=', down_NES, ', padj=', down_padj)) +
  geom_vline(xintercept = sum(ordered.de.results.t2high$`-10logFDR` < 0)) +
  ylab("")

all_plt <- ggarrange(up_corrected, down_corrected,
                     nrow = 2, 
                     labels = c("A",  "B"))
png("./ATLANTIS/th_high_th_low/plots/GSEA_ICS_corrected_background_NOT_corrected.png",
    width=1400, height=850, res = 200)

print(all_plt)
dev.off()

