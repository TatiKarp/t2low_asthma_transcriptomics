#This script will define T2-high/low stable and create a sankeyplot
library(dplyr)
library(readr)
library(ggplot2)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")
# master.Table <- read_csv("./th_high_th_low/master_table_th2.csv")

big_master_table <-  read_csv("./atlantis_patient_data.csv") %>%
  dplyr::select (c(PT, ASTHEA, VISIT, SYS_COR, BIO, LABEOSV, FENRES)) %>%
  filter(ASTHEA == "A") %>%
  mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES')) %>%
  filter(include == 'YES') %>%
  mutate(group_th = case_when(
    ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
    ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
    ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
    ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low',
    #is.na(LABEOSV) ~ 'labeos_NA',
    TRUE ~ 'undeterm'))

# delete sample with NA for blood eos at least at one visit 
na_visits_wide <- big_master_table %>%
  mutate(PT = as.character(PT)) %>%
  dplyr::select(c(PT, VISIT, LABEOSV)) %>%
  tidyr::pivot_wider (names_from = VISIT, values_from= LABEOSV) %>%
  mutate(Has_NA = if_any(c(`VISIT 1`, `VISIT 2`, `VISIT 3`), is.na))


sample_to_include <- na_visits_wide[na_visits_wide$Has_NA == FALSE, ]$PT


big_master_table <- big_master_table %>%
  filter(PT %in% sample_to_include)

stacked_plt <- ggplot(big_master_table, aes(x = VISIT, fill = group_th))+
  geom_bar(position = "fill")

big_master_table_wide <- big_master_table %>%
  mutate(PT = as.character(PT)) %>%
  dplyr::select(c(PT, VISIT, group_th)) %>%
  tidyr::pivot_wider (names_from = VISIT, values_from= group_th)


## sankey plot

library(ggalluvial)

## plot only for samples with RNA-seq data

master.Table <- read.csv("./th_high_th_low/master_table_th2.csv") #samples that have RNA data
plt_groups_RNA <- big_master_table %>%
  filter(PT %in% master.Table$PT) %>%
  mutate(group_th = as.factor(case_when(
    group_th == "high" ~ "T2-high",
    group_th == "low" ~ "T2-low",
    group_th == "undeterm" ~ "Undetermined"
  )))

#order levels 
plt_groups_RNA$group_th <- factor(plt_groups_RNA$group_th, levels = c("T2-high", "Undetermined", "T2-low"))

sankey_plt <- ggplot(plt_groups_RNA, aes(x = VISIT, stratum = group_th, 
                                      alluvium = PT, fill = group_th, label = group_th)) +
  scale_fill_manual(breaks = as.factor(plt_groups_RNA$group_th),
                    values = c('T2-low'= "#8DA0CB" , 'Undetermined' = "#B8C4C2",'T2-high'="#FFD92F"), 
                    name = "") +
  #scale_fill_brewer(type = "qual", palette = "Set2", name = "group") +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  ylab ("Number of participants") + 
  theme(legend.position = "bottom") +
  theme(axis.text=element_text(size=15),
        axis.title.x=element_blank(),
        axis.title.y =element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=12),
        legend.key = element_rect(fill = "white"))

library(colorBlindness)
cvdPlot(sankey_plt)

saveRDS(sankey_plt, "./th_high_th_low/plots/t2stable_sankey_plt_only_RNA.rds")

png("./th_high_th_low/plots/sankey_plt_only_RNA.png", width = 1400, height = 700, res = 150)
print(sankey_plt)
dev.off()



