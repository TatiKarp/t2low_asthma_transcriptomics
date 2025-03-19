## construct WGCNA - identify gene modules - correlate with asthma status
library(dplyr)
library(DESeq2)
library(WGCNA)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)
library(biomaRt)
library(rstatix)


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
setwd('/Users/tatiana/Work/RP2/ATLANTIS')

# metadata
master.Table.low.healthy <- read.csv('./th_high_th_low/master_table_th2.csv')%>%
  filter(group_th == 'low' | group_th == 'healthy') #only Th2-low asthmatics


# add expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table.low.healthy$GenomeScan_ID))%>%
  as.matrix()

# QC - outlier detection ------------------------------------------------

gsg <- goodSamplesGenes(t(expression.data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes) ## exclude 4629 genes
table(gsg$goodSamples)

# remove genes that are detectd as outliers
expression.data <- expression.data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(expression.data)), method = "average")
plot(htree)

# pca - method 2

pca <- prcomp(t(expression.data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# exclude outlier samples if needed !!

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# fixing column names in colData
colData <- master.Table.low.healthy%>%
  tibble::column_to_rownames('GenomeScan_ID')

# making the rownames and column names identical
all(rownames(colData) %in% colnames(expression.data))
all(rownames(colData) == colnames(expression.data))

# create dds
dds <- DESeqDataSetFromMatrix(countData = expression.data,
                              colData = colData,
                              design = ~ 1) # not specifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ - https://edo98811.github.io/WGCNA_official_documentation/faq.html

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) # 14999 genes

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 7
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 15000,
                          TOMType = "signed", # the same direction
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

write.csv(module_eigengenes%>%
            tibble::rownames_to_column('Sample'), './th_high_th_low/wgcna_modiles_eigengenes.csv',
          row.names = FALSE)
# Print out a preview
head(module_eigengenes)

#write.csv(module_eigengenes, './th_high_th_low/WGCNA_th2low/WGCNA_modules_eigengenes.csv')

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# 6. Associate modules with Th2-low asthma/ healthy status 
# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(group_th == 'low', 1, 0)) %>% 
  dplyr::select(c(disease_state_bin))

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  tibble::column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[13],
             y = names(heatmap.data)[1:12],
             col = c("blue1", "skyblue", "white", "pink", "red"))


module.gene.mapping <- as.data.frame(bwnet$colors)
sign.modules <- module.gene.mapping %>% 
  filter(`bwnet$colors` %in% c('black', 'purple', 'greenyellow'))

all_new_gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                      filters = 'ensembl_gene_id', values = rownames(module.gene.mapping), mart = ensembl)

write.csv(sign.modules%>%
            tibble::rownames_to_column('gene')%>%
            left_join(all_new_gene, by = c('gene'='ensembl_gene_id'))%>%
            tibble::column_to_rownames('gene')%>%
            arrange(`bwnet$colors`),'./th_high_th_low/WGCNA_th2low/Genes_in_sign_modules.csv')
write.csv(module.gene.mapping%>%
            tibble::rownames_to_column('gene')%>%
            left_join(all_new_gene, by = c('gene'='ensembl_gene_id'))%>%
            tibble::column_to_rownames('gene')%>%
            arrange(`bwnet$colors`), './th_high_th_low/WGCNA_th2low/wgcna.modules.csv')

# plot eigengene data per group for each module

module_eigengenes_plt <- read.csv('./th_high_th_low/WGCNA_th2low/WGCNA_modules_eigengenes.csv')%>%
  tidyr::pivot_longer(-X, names_to= c('module'))%>%
  dplyr::rename(Sample = X)%>%
  left_join(colData%>%
              tibble::rownames_to_column('Sample'), by = c('Sample'='Sample'))
  
png('./th_high_th_low/WGCNA_th2low/Eigengene_values_sign_modules.png', res = 300,  width = 1880, height = 1400)
plt <- ggplot(module_eigengenes_plt%>%
         filter(module %in% paste0('ME',unique(sign.modules$`bwnet$colors`))), aes(module, value))+
  geom_boxplot(aes(fill = group_th),outlier.shape = NA)+
  geom_point(aes(fill = group_th),size = 1, shape = 21, position=position_jitterdodge(0.1), alpha= 0.5)+
  ylab('Eigengene values')+
  xlab('Module')+
  scale_fill_manual(breaks = module_eigengenes_plt$group_th,
                    values = c('healthy'='#66C2A5','low'='#8DA0CB'))+
  theme_classic()
print(plt)
dev.off()

stat.test <- module_eigengenes_plt%>%
  group_by(module)%>%
  rstatix::t_test(value~group_th)%>%
  rstatix::add_y_position(step.increase = 0.08) %>%
  rstatix::add_x_position(x = "module", dodge = 0) %>%
  mutate(p = round(p, digits = 3))%>%
  mutate(p.adj = round(p.adjust(p, method = 'fdr'), digits = 3))
# create a table only nominal p-value significant comparisons
stat.test.sign <- stat.test%>%
  filter(p < 0.05)
  

png('./th_high_th_low/WGCNA_th2low/Eigengene_values_all_modules.png', res = 300,  width = 2300, height = 1600)
plt <- ggplot(module_eigengenes_plt, aes(module, value))+
  geom_boxplot(aes(fill = group_th),outlier.shape = NA)+
  geom_point(aes(fill = group_th),size = 1, shape = 21, position=position_jitterdodge(0.1), alpha= 0.5)+
  ylab('Module eigengene value')+
  xlab('Module')+
  scale_fill_manual(breaks = levels(factor(module_eigengenes_plt$group_th)),
                    values = c('healthy'='#66C2A5','low'='#8DA0CB'), 
                    labels = c("Healthy", "T2-low"),
                    name = "") +
  theme_classic()+
  ggpubr::stat_pvalue_manual(stat.test, label = "p adj = {p.adj}", xmin = 'xmin', xmax ='xmax',size = 2.5, remove.bracket = TRUE)+
  ggpubr::stat_pvalue_manual(stat.test, label = "p = {p}", xmin = 'xmin', xmax ='xmax',size = 2.5, remove.bracket = TRUE,
                             y.position = stat.test$y.position + 0.03)+
  ylim(-0.4, 0.6) +
  #annotate("text",x = Inf, y = -0.4, hjust = 1, vjust = 0, label = "t-test, FDR adjusted") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"))

print(plt)
dev.off()
saveRDS(plt, "./th_high_th_low/plots/Eigengene_values_all_modules.rds")


# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile of each gene in the module. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# only interesting modules 
interesting.module.membership.measure <- module.membership.measure[c('MEblack','MEpurple','MEgreenyellow'),]
interesting.module.membership.measure.pvals <- module.membership.measure.pvals[c('MEblack','MEpurple','MEgreenyellow'),]

# create a data frame with all the correlations and p-values for interesting modules 
mod_names <- rownames(interesting.module.membership.measure)
names = colnames(interesting.module.membership.measure)
all_new_gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                      filters = 'ensembl_gene_id', values = names, mart = ensembl)

df.module.memberships <- as.data.frame(interesting.module.membership.measure)%>%
  tibble::rownames_to_column('rowname')%>%
  mutate(rowname = paste0('cor_', mod_names))%>%
  tibble::column_to_rownames('rowname')%>%
  t()%>%
  as.data.frame()%>%
  tibble::rownames_to_column('gene')%>%
  left_join(as.data.frame(t(as.data.frame(interesting.module.membership.measure.pvals)%>%
              tibble::rownames_to_column('rowname')%>%
              mutate(rowname = paste0('p_val_', mod_names))%>%
              tibble::column_to_rownames('rowname')))%>%
              tibble::rownames_to_column('gene'), by = c('gene'='gene'))%>%
  mutate(across(starts_with("p_val_"), ~p.adjust(.x, method = "BH"), .names = "{col}_adjusted"))%>%
  left_join(all_new_gene, by = c('gene' = 'ensembl_gene_id'))

## create separate table per module with module membership and p-value:
create_df_per_module <- function(module){
  df.one.module <- as.data.frame(module.membership.measure[module,, drop= FALSE])%>%
    tibble::rownames_to_column('rowname')%>%
    mutate(rowname = paste0('cor_', rowname))%>%
    tibble::column_to_rownames('rowname')%>%
    t()%>%
    as.data.frame()%>%
    tibble::rownames_to_column('gene')%>%
    filter(gene %in% c(rownames(module.gene.mapping)[module.gene.mapping$`bwnet$colors` == substring(module,3,nchar(module))]))%>%
    left_join(as.data.frame(interesting.module.membership.measure.pvals[module,, drop= FALSE])%>%
                              tibble::rownames_to_column('rowname')%>%
                              mutate(rowname = paste0('p_val_', rowname))%>%
                              tibble::column_to_rownames('rowname')%>%
                              t()%>%
                              as.data.frame()%>%
                              tibble::rownames_to_column('gene'),
              by = c('gene'='gene'))%>%
    mutate(across(starts_with("p_val_"), ~p.adjust(.x, method = "BH"), .names = "{col}_adjusted"))%>%
    left_join(all_new_gene, by = c('gene' = 'ensembl_gene_id'))%>%
    write.csv2(file = paste0('./th_high_th_low/WGCNA_th2low/', module, '.genes.module.memberships.csv'), row.names = FALSE)
}
sapply(c('MEblack','MEpurple','MEgreenyellow'), create_df_per_module)


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$disease_state_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

