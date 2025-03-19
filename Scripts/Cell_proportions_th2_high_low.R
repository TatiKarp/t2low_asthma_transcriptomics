### This script will check the cell type proportions in TH2-high/ Th2-low and healthy in ATLATNIS

library(dplyr)
library(ggplot2)
library(ggpubr)
## Cibersort with matrix from Tessa (also used in Zaid's eQTL project) - with ENSEBLE IDs

setwd('/Users/tatiana/Work/RP2/ATLANTIS/')
## upload source Cibersort
source('/Users/tatiana/Work/RP2/ATLANTIS/Deconvolution/Deconvolution_Jos/CIBERSORT.R', verbose=TRUE)

master.Table <- read.csv('./th_high_th_low/master_table_th2.csv')

C <- read.csv("./Season/cell_types/CIBERSORTx_nose_subsampled_matrix_max200cells_ENSG_inferred_phenoclasses.CIBERSORTx_nose_subsampled_matrix_max200cells_ENSG_inferred_refsample.bm.K999.txt", sep='\t')
rownames(C)<-C$NAME
C <- C%>%
  dplyr::select(-NAME)

expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))

## cpm normalization
bulk = edgeR::cpm(expression.data)
bulk <- bulk %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Gene')

## overlap genes from signature and from bulk data
keep = intersect(rownames(C), bulk$Gene) ## keep 2034 out of the 2035

bulk <- bulk %>%
  filter(Gene %in% keep) %>%
  tibble::column_to_rownames("Gene")

Ref = C[keep,]

##############################################################
##################### run Cibersort ##########################

xf <- "./th_high_th_low/reference.tsv"
Ref %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames")%>%
  readr::write_tsv(
    file = xf
  )

yf <-  "./th_high_th_low//mixture.tsv"
bulk %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames") %>%
  readr::write_tsv(
    file = yf
  )

RESULTS <- CIBERSORT(sig_matrix = xf, mixture_file = yf, QN = FALSE, perm=100)
res.cibersort <- t(RESULTS[,1:(ncol(RESULTS)-3)]) %>%
  as.data.frame()

write.csv(res.cibersort, "./th_high_th_low/CIBERSORT.proportions.csv", row.names = T, quote = F)

### plot cell type proportions 
# prepare table for the plot
res.cibersort <- read.csv("./th_high_th_low/CIBERSORT.proportions.csv") %>%
  tibble::column_to_rownames ("X")

res.cibersort.plot <- res.cibersort %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer( -Sample,
    names_to = "cell_type",
    values_to = "proportion")

master.Table <- master.Table%>%
  left_join(res.cibersort.plot, by = c('GenomeScan_ID'='Sample'))

# add p-values
master.Table.plt <- master.Table %>%
  filter(group_th != 'undeterm') %>%
  mutate(cell_type = gsub("\\.", " ", cell_type)) %>%
  mutate(cell_type = if_else (cell_type =="Basal resting   Suprabasal", "Basal resting", cell_type)) %>%
  filter(cell_type != 'T cell lineage')

master.Table.plt$group_th <- factor(master.Table.plt$group_th, 
                                    levels = c("healthy", "low", "high"))

stat.test <- master.Table.plt%>%
  group_by(cell_type)%>%
  rstatix::wilcox_test(proportion ~ group_th, 
                       ref.group = 'healthy', 
                       p.adjust.method = "fdr")%>%
  rstatix::add_y_position(step.increase = 0.03)%>%
  rstatix::add_x_position(x = "cell_type", dodge = 0.7) 

# plot for the paper 
png('./th_high_th_low/plots/Cell_types_th2_high_low_CIBERSORT.png', res = 300,  width = 1880, height = 1400)
plt <- ggplot(master.Table.plt, aes(x = cell_type,y = proportion)) +
  geom_boxplot(aes(fill = group_th), outlier.shape = NA)+
  geom_point(aes(fill = group_th), size = 1, shape = 21, position=position_jitterdodge(0.1), alpha= 0.5)+
  ylab('Proportion')+
  xlab('Cell type')+
  scale_fill_manual(breaks = levels(factor(master.Table.plt$group_th)),
                    values = c("healthy"="#66C2A5","low"="#8DA0CB", "high"="#FFD92F"),
                    labels = c("Healthy", "T2-low", "T2-high"),
                    name = "")+
  # theme_classic()+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
  #       plot.margin = unit(c(1, 1, 1, 2), "lines"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"))+
  stat_pvalue_manual(stat.test, label = "p.adj = {p.adj}", size = 2.5, linetype = 1,  tip.length = 0.01, remove.bracket = FALSE)
#annotate("text",x = Inf, y = 1, hjust = 1, vjust = 0, label = "Wilcox test, p adjusted")
print(plt)
dev.off()

saveRDS(plt, "./th_high_th_low/plots/Cell_types_th2_high_low_CIBERSORT.rds")

library(colorBlindness)
cvdPlot(plt)
# Check median cell counts:
master.Table %>%
  group_by(cell_type, group_th) %>%
  summarise(Median = median(proportion)) 

summary (master.Table %>%
           filter(group_th == "low") %>%
           filter(cell_type == "Club"))

summary (master.Table %>%
           filter(group_th == "healthy") %>%
           filter(cell_type == "Club"))




