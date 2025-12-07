#Aim 2 - Longitudinal evaluation of diet-induced gut microbiome alterations and their implications for gut health

#install lme4 packages (for Shirley's older R version)
install.packages("lme4", type = "source")
install.packages("lmerTest")

#load library
library(phyloseq)
library(tidyverse)
library(microbiome)
library(vegan)
library(ggplot2)
library(DESeq2)
library(lme4)
library(lmerTest)
library(emmeans)
library(pheatmap)
library(ape)
library(picante)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(patchwork)
library(tidyr)
library(tibble)
library(patchwork)
library(grid)


#### Load Data ####
#load filtered and rarefied data
load("ferm_rare.RData")

#check object loaded
ls()
ferm_rarefac
sample_data(ferm_rarefac) %>% head()



#### Alpha Diversity Linear Mixed Model (LMM) ####
##### Prepare Dataset #####
#prune zero-sum taxa and samples in the full dataset
ferm_rarefac <- prune_taxa(taxa_sums(ferm_rarefac) > 0, ferm_rarefac)
ferm_rarefac <- prune_samples(sample_sums(ferm_rarefac) > 0, ferm_rarefac)

#calculate Shannon diversity
alpha_shannon <- estimate_richness(ferm_rarefac, measures = "Shannon")

#calculate Faith's PD
#OTU table as matrix, samples in rows
otu_mat <- as(otu_table(ferm_rarefac), "matrix")
if(taxa_are_rows(ferm_rarefac)) otu_mat <- t(otu_mat)

#remove zero-sum taxa and samples
otu_mat <- otu_mat[, colSums(otu_mat) > 0, drop = FALSE]
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

#prune phylogenetic tree to match OTUs
#tree <- prune_taxa(colnames(otu_mat), phy_tree(ferm_rarefac))
tree <- phy_tree(ferm_rarefac)
common_otus <- intersect(colnames(otu_mat), tree$tip.label)
tree <- drop.tip(tree, setdiff(tree$tip.label, common_otus))
if(!is.rooted(tree)) tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)


#calculate Faith's PD
alpha_faith_pd <- pd(otu_mat, tree, include.root = TRUE)

#add alpha metrics to metadata
meta_full <- as.data.frame(sample_data(ferm_rarefac))
meta_full$Shannon <- alpha_shannon$Shannon[match(rownames(meta_full), rownames(alpha_shannon))]
meta_full$Faith_PD <- alpha_faith_pd$PD[match(rownames(meta_full), rownames(alpha_faith_pd))]

#update phyloseq object
sample_data(ferm_rarefac) <- sample_data(meta_full)


#subset phyloseq object based on intervention periods
#fresh vegetables: P1 -> P2 -> P3
ps_fresh <- subset_samples(ferm_rarefac, period %in% c("Base","VEG","WO1"))

#fermented vegetables: P1 -> P4 -> P5
ps_ferm <- subset_samples(ferm_rarefac, period %in% c("Base","FERM","WO2"))

#microbiome washouts: P1 -> P3 -> P5
ps_microb <- subset_samples(ferm_rarefac, period %in% c("Base","WO1","WO2"))



#prepare and check data for LMM
#ps_fresh
meta_fresh <- as(sample_data(ps_fresh), "data.frame")
#make participant_id and period as factor 
meta_fresh$participant_id <- factor(meta_fresh$participant_id)
meta_fresh$period <- factor(meta_fresh$period, levels = c("Base", "VEG", "WO1"))
str(meta_fresh[, c("participant_id", "period", "Shannon", "Faith_PD")])

#check for NAs
sum(is.na(meta_fresh$Shannon))
sum(is.na(meta_fresh$Faith_PD))
sum(is.na(meta_fresh$period))
sum(is.na(meta_fresh$participant_id))

#ps_ferm
meta_ferm <- as(sample_data(ps_ferm), "data.frame")
#make participant_id and period as factor 
meta_ferm$participant_id <- factor(meta_ferm$participant_id)
meta_ferm$period <- factor(meta_ferm$period, levels = c("Base", "FERM", "WO2"))
str(meta_ferm[, c("participant_id", "period", "Shannon", "Faith_PD")])

#check for NAs
sum(is.na(meta_ferm$Shannon))
sum(is.na(meta_ferm$Faith_PD))
sum(is.na(meta_ferm$period))
sum(is.na(meta_ferm$participant_id))

#ps_microb
meta_microb <- as(sample_data(ps_microb), "data.frame")
#make participant_id and period as factor 
meta_microb$participant_id <- factor(meta_microb$participant_id)
meta_microb$period <- factor(meta_microb$period, levels = c("Base", "WO1", "WO2"))
str(meta_microb[, c("participant_id", "period", "Shannon", "Faith_PD")])

#check for NAs
sum(is.na(meta_microb$Shannon))
sum(is.na(meta_microb$Faith_PD))
sum(is.na(meta_microb$period))
sum(is.na(meta_microb$participant_id))

##### Fresh LMM #####
###### Shannon  ######
#fresh shannon LMM (period as the fixed effects and participant_id as the random effect)
model_shannon_fresh <- lmer(Shannon ~ period + (1 | participant_id), data = meta_fresh)

#perform ANOVA to test and post-hoc test
anova(model_shannon_fresh)
emmeans(model_shannon_fresh, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_fresh_s <- emmeans(model_shannon_fresh, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_fresh_s <- emmeans(model_shannon_fresh, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_FRS <- contrasts_fresh_s %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_fresh_s$emmean[emm_df_fresh_s$period == g]) + 0.25
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_fresh_veg_shannon <- ggplot(emm_df_fresh_s, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  #add significance labels
  stat_pvalue_manual(
    pval_table_FRS,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.15
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


#box-like plot
LMM_fresh_veg_shannon_box <- ggplot(emm_df_fresh_s, aes(x = period, y = emmean)) +
  #box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_FRS,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.1
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "VEG" = "#E69F00", "WO1" = "#009E73"))

#save
ggsave("A2_LMM_Alpha_Shannon_Fresh.png", plot = LMM_fresh_veg_shannon, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_Shannon_Fresh_box.png", plot = LMM_fresh_veg_shannon_box, width = 6, height = 4, dpi = 300)


###### Faith's PD  ######
#fresh Faith LMM
model_faith_fresh <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_fresh)

#ANOVA and post-hoc contrasts
anova(model_faith_fresh)
emmeans(model_faith_fresh, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_fresh_f <- emmeans(model_faith_fresh, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_fresh_f <- emmeans(model_faith_fresh, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_FRF <- contrasts_fresh_f %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_fresh_f$emmean[emm_df_fresh_f$period == g]) + 2.5
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_fresh_veg_faith <- ggplot(emm_df_fresh_f, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FRF,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.15
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 


#box-like plot
LMM_fresh_veg_faith_box <- ggplot(emm_df_fresh_f, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FRF,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.1
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "VEG" = "#E69F00", "WO1" = "#009E73"))


#save plot
ggsave("A2_LMM_Alpha_FaithPD_Fresh.png", plot = LMM_fresh_veg_faith, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_FaithPD_Fresh_box.png", plot = LMM_fresh_veg_faith_box, width = 6, height = 4, dpi = 300)

##### Fermentation LMM  #####
###### Shannon ######

#ferm Shannon LMM
model_shannon_ferm <- lmer(Shannon ~ period + (1 | participant_id), data = meta_ferm)

#ANOVA and post-hoc contrasts
anova(model_shannon_ferm)
emmeans(model_shannon_ferm, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_ferm_s <- emmeans(model_shannon_ferm, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_ferm_s <- emmeans(model_shannon_ferm, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_FES <- contrasts_ferm_s %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_ferm_s$emmean[emm_df_ferm_s$period == g]) + 1
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_ferm_veg_shannon <- ggplot(emm_df_ferm_s, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FES,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.15
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 


#box-like plot
LMM_ferm_veg_shannon_box <- ggplot(emm_df_ferm_s, aes(x = period, y = emmean)) +
  #box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_FES,
    label = "label",
    tip.length = 0.03,
    step.increase = .5
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "FERM" = "#E69F00", "WO2" = "#009E73"))

#save plot
ggsave("A2_LMM_Alpha_Shannon_Ferm.png", plot = LMM_ferm_veg_shannon, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_Shannon_Ferm_box.png", plot = LMM_ferm_veg_shannon_box, width = 6, height = 4, dpi = 300)


###### Faith's PD ######
#ferm Faith LMM
model_faith_ferm <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_ferm)
#ANOVA and post-hoc contrasts
anova(model_faith_ferm)
emmeans(model_faith_ferm, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_ferm_f <- emmeans(model_faith_ferm, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_ferm_f <- emmeans(model_faith_ferm, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_FEF <- contrasts_ferm_f %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_ferm_f$emmean[emm_df_ferm_f$period == g]) + 3.25
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_ferm_veg_faith <- ggplot(emm_df_ferm_f, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Ferm Intervention – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FEF,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.1
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) 


#box-like plot
LMM_ferm_veg_faith_box <- ggplot(emm_df_ferm_f, aes(x = period, y = emmean)) +
  #box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_FEF,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.15
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "FERM" = "#E69F00", "WO2" = "#009E73"))

#save plot
ggsave("A2_LMM_Alpha_FaithPD_Ferm.png", plot = LMM_ferm_veg_faith, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_FaithPD_Ferm_box.png", plot = LMM_ferm_veg_faith_box, width = 6, height = 4, dpi = 300)

##### Microbiome LMM #####
###### Shannon ######
#microbiome Shannon LMM
model_shannon_microb <- lmer(Shannon ~ period + (1 | participant_id), data = meta_microb)

#ANOVA and post-hoc contrasts
anova(model_shannon_microb)
emmeans(model_shannon_microb, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_microb_s <- emmeans(model_shannon_microb, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_microb_s <- emmeans(model_shannon_microb, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_MS <- contrasts_microb_s %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_microb_s$emmean[emm_df_microb_s$period == g]) + 0.3
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_microbiome_shannon <- ggplot(emm_df_microb_s, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Microbiome – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add significance labels
  stat_pvalue_manual(
    pval_table_MS,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.05
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 


#box-like plot
LMM_microbiome_shannon_box <- ggplot(emm_df_microb_s, aes(x = period, y = emmean)) +
  #box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_MS,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.1
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Microbiome – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "WO1" = "#E69F00", "WO2" = "#009E73"))


#save plot
ggsave("A2_LMM_Alpha_Shannon_Microbiome.png", plot = LMM_microbiome_shannon, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_Shannon_Microbiome_box.png", plot = LMM_microbiome_shannon_box, width = 6, height = 4, dpi = 300)


###### Faith's PD ######
##microbiome Faith LMM
model_faith_microb <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_microb)
#ANOVA and post-hoc contrasts
anova(model_faith_microb)
emmeans(model_faith_microb, pairwise ~ period)

#extract the estimated marginal means from each period and convert to data frame
emm_df_microb_f <- emmeans(model_faith_microb, ~ period) %>% 
  as.data.frame()

#extract pairwise contrasts
contrasts_microb_f <- emmeans(model_faith_microb, pairwise ~ period)$contrasts %>%
  as.data.frame()

#prepare table for significance
pval_table_MF <- contrasts_microb_f %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_microb_f$emmean[emm_df_microb_f$period == g]) + 1
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 

#plot
LMM_microbiome_faith <- ggplot(emm_df_microb_f, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Microbiome – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # Add significance labels
  stat_pvalue_manual(
    pval_table_MF,
    label = "label",
    tip.length = 0.03,
    step.increase =  2.5
  ) +
  #increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) 


#box-like plot
LMM_microbiome_faith_box <- ggplot(emm_df_microb_f, aes(x = period, y = emmean)) +
  #box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_MF,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.15
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Microbiome – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "WO1" = "#E69F00", "WO2" = "#009E73"))


#save plot
ggsave("A2_LMM_Alpha_FaithPD_Microbiome.png", plot = LMM_microbiome_faith, width = 6, height = 4, dpi = 300)
ggsave("A2_LMM_Alpha_FaithPD_Microbiome_box.png", plot = LMM_microbiome_faith_box, width = 6, height = 4, dpi = 300)


#### Beta Diversity Linear Mixed Model (LMM) ####
##### Prepare Dataset #####
#prepare metadata
meta_fresh_df <- as.data.frame(sample_data(ps_fresh))
meta_fresh_df$SampleID <- rownames(meta_fresh_df)

meta_ferm_df <- as.data.frame(sample_data(ps_ferm))
meta_ferm_df$SampleID <- rownames(meta_ferm_df)

meta_microb_df <- as.data.frame(sample_data(ps_microb))
meta_microb_df$SampleID <- rownames(meta_microb_df)

##### Calculate Bray-Curtis and Weighted UniFrac Metrics #####
#beta diversity matrices calculations (bray and weighted unifrac)
#fresh
beta_bc_fresh <- phyloseq::distance(ps_fresh, method = "bray")
beta_wu_fresh <- phyloseq::distance(ps_fresh, method = "wunifrac")

#fermented
beta_bc_ferm <- phyloseq::distance(ps_ferm, method = "bray")
beta_wu_ferm <- phyloseq::distance(ps_ferm, method = "wunifrac")

#microbiome washouts
beta_bc_microb <- phyloseq::distance(ps_microb, method = "bray")
beta_wu_microb <- phyloseq::distance(ps_microb, method = "wunifrac")

##### Calculate Distance-to-baseline #####
#distance-to-baseline calculation function
distance_to_baseline <- function(dist_matrix, meta) {
  # Ensure metadata has SampleID
  if(!"SampleID" %in% colnames(meta)) meta$SampleID <- rownames(meta)
  
  # Convert dist object to matrix if needed
  if(inherits(dist_matrix, "dist")) dist_matrix <- as.matrix(dist_matrix)
  
  # Ensure rownames and colnames match SampleID
  if(is.null(rownames(dist_matrix))) rownames(dist_matrix) <- colnames(dist_matrix) <- meta$SampleID
  if(is.null(colnames(dist_matrix))) colnames(dist_matrix) <- rownames(dist_matrix) <- meta$SampleID
  
  dtb <- data.frame()
  
  for(pid in unique(meta$participant_id)) {
    samp_ids <- meta$SampleID[meta$participant_id == pid]
    
    # Identify baseline sample
    base_samp <- samp_ids[meta$period[meta$SampleID %in% samp_ids] == "Base"]
    
    if(length(base_samp) == 0) next
    
    # Compute distance to baseline for each sample
    for(s in samp_ids) {
      dtb <- rbind(dtb, data.frame(
        participant_id = pid,
        SampleID = s,
        period = meta$period[meta$SampleID == s],
        distance = dist_matrix[s, base_samp]
      ))
    }
  }
  
  return(dtb)
}


# Convert dist object to matrix first
beta_bc_fresh_mat <- as.matrix(beta_bc_fresh)
beta_wu_fresh_mat <- as.matrix(beta_wu_fresh)

beta_bc_ferm_mat <- as.matrix(beta_bc_ferm)
beta_wu_ferm_mat <- as.matrix(beta_wu_ferm)

beta_bc_microb_mat <- as.matrix(beta_bc_microb)
beta_wu_microb_mat <- as.matrix(beta_wu_microb)

#distance-to-baseline calculations
#fresh
dist_bc_fresh <- distance_to_baseline(beta_bc_fresh_mat, meta_fresh_df)
dist_wu_fresh <- distance_to_baseline(beta_wu_fresh_mat, meta_fresh_df)
#ferm
dist_bc_ferm <- distance_to_baseline(beta_bc_ferm_mat, meta_ferm_df)
dist_wu_ferm <- distance_to_baseline(beta_wu_ferm_mat, meta_ferm_df)
#microbiome washouts
dist_bc_microb <- distance_to_baseline(beta_bc_microb_mat, meta_microb_df)
dist_wu_microb <- distance_to_baseline(beta_wu_microb_mat, meta_microb_df)

##### LMM on Distance-to-baseline #####
#function fits LMM (period as fixed effect and participant_id as random effects) and returns ANOVA and post-hoc contract
run_lmm <- function(dist_df, response_label = "Distance") {
  # Ensure factors
  dist_df$participant_id <- factor(dist_df$participant_id)
  dist_df$period <- factor(dist_df$period, levels = unique(dist_df$period))
  
  # Fit linear mixed model
  lmm_model <- lmer(distance ~ period + (1 | participant_id), data = dist_df)
  
  # ANOVA table
  anova_res <- anova(lmm_model)
  
  # Post-hoc pairwise comparisons
  emmeans_res <- emmeans(lmm_model, pairwise ~ period)
  
  # Prepare estimated marginal means for plotting
  emm_df <- as.data.frame(emmeans(lmm_model, ~ period))
  emm_df$response <- response_label
  
  # Return as a list
  return(list(
    model = lmm_model,
    anova = anova_res,
    emmeans = emmeans_res,
    emm_df = emm_df
  ))
}

#Fresh Bray-Curtis
lmm_bc_fresh <- run_lmm(dist_bc_fresh, response_label = "Bray-Curtis")
lmm_wu_fresh <- run_lmm(dist_wu_fresh, response_label = "Weighted UniFrac")

# Fermented
lmm_bc_ferm <- run_lmm(dist_bc_ferm, response_label = "Bray-Curtis")
lmm_wu_ferm <- run_lmm(dist_wu_ferm, response_label = "Weighted UniFrac")

# Microbiome washouts
lmm_bc_microb <- run_lmm(dist_bc_microb, response_label = "Bray-Curtis")
lmm_wu_microb <- run_lmm(dist_wu_microb, response_label = "Weighted UniFrac")


#quick check
head(dist_bc_fresh)
str(dist_bc_fresh)
table(dist_bc_fresh$period)
table(dist_bc_fresh$participant_id)

##### Plots #####
#function to create combined line + boxplot figure
plot_distance_lmm <- function(dist_df, lmm_res, title_prefix) {
  #extract post-hoc contrasts for significance
  contrasts_df <- as.data.frame(lmm_res$emmeans$contrasts) %>%
    mutate(
      group1 = gsub(" -.*","",contrast),
      group2 = gsub(".*- ","",contrast),
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*"
      )
    ) %>%
    #stagger y positions to avoid overlap
    mutate(
      y.position = max(dist_df$distance) + seq(0.05, 0.05*nrow(.), by = 0.05)*max(dist_df$distance)
    ) %>%
    dplyr::select(group1, group2, sig, y.position)
  
  #line/trajectory plot
  line_plot <- ggplot(dist_df, aes(x = period, y = distance, group = participant_id)) +
    geom_line(alpha = 0.3, color = "grey") +
    stat_summary(aes(group=1), fun = mean, geom = "line", color = "blue", size = 1.2) +
    stat_summary(aes(group=1), fun.data = mean_se, geom = "ribbon", fill = "blue", alpha = 0.2) +
    labs(title = paste0(title_prefix, " - Trajectories"), x = "Period", y = "Distance to Baseline") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  #boxplot/violin with significance
  box_plot <- ggplot(dist_df, aes(x = period, y = distance, fill = period)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 3, color = "black", fill = "yellow") +
    ggpubr::stat_pvalue_manual(contrasts_df, label = "sig", tip.length = 0.01) +
    labs(title = paste0(title_prefix, " - Summary"), x = "Period", y = "Distance to Baseline") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  #combine side by side
  combined_plot <- line_plot | box_plot
  return(combined_plot)
}

#generate plots for all datasets
#fresh
beta_fresh_BC <- plot_distance_lmm(dist_bc_fresh, lmm_bc_fresh, "Fresh Bray-Curtis")
ggsave("A2_LMM_Beta_BC_Fresh.png", plot = beta_fresh_BC, width = 6, height = 4, dpi = 300)

beta_fresh_WUF <- plot_distance_lmm(dist_wu_fresh, lmm_wu_fresh, "Fresh Weighted UniFrac")
ggsave("A2_LMM_Beta_WUF_Fresh.png", plot = beta_fresh_WUF, width = 6, height = 4, dpi = 300)

#fermented
beta_ferm_BC <- plot_distance_lmm(dist_bc_ferm, lmm_bc_ferm, "Fermented Bray-Curtis")
ggsave("A2_LMM_Beta_BC_Ferm.png", plot = beta_ferm_BC, width = 6, height = 4, dpi = 300)

beta_ferm_WUF <- plot_distance_lmm(dist_wu_ferm, lmm_wu_ferm, "Fermented Weighted UniFrac")
ggsave("A2_LMM_Beta_WUF_Ferm.png", plot = beta_ferm_WUF, width = 6, height = 4, dpi = 300)


#microbiome washout
beta_microbiome_BC <- plot_distance_lmm(dist_bc_microb, lmm_bc_microb, "Washout Bray-Curtis")
ggsave("A2_LMM_Beta_BC_Microbiome.png", plot = beta_microbiome_BC, width = 6, height = 4, dpi = 300)

beta_microbiome_WUF <- plot_distance_lmm(dist_wu_microb, lmm_wu_microb, "Washout Weighted UniFrac")
ggsave("A2_LMM_Beta_WUF_Microbiome.png", plot = beta_microbiome_WUF, width = 6, height = 4, dpi = 300)



#function to extract post-hoc contrasts as a table
get_posthoc_table <- function(lmm_res, response_label = "Distance") {
  contrasts_df <- as.data.frame(lmm_res$emmeans$contrasts) %>%
    dplyr::mutate(
      response = response_label
    )
  return(contrasts_df)
}

#extract post-hoc contrasts 
posthoc_bc_fresh <- get_posthoc_table(lmm_bc_fresh, "Bray-Curtis - Fresh")
posthoc_wu_fresh <- get_posthoc_table(lmm_wu_fresh, "Weighted UniFrac - Fresh")

posthoc_bc_ferm <- get_posthoc_table(lmm_bc_ferm, "Bray-Curtis - Fermented")
posthoc_wu_ferm <- get_posthoc_table(lmm_wu_ferm, "Weighted UniFrac - Fermented")

posthoc_bc_microb <- get_posthoc_table(lmm_bc_microb, "Bray-Curtis - Washout")
posthoc_wu_microb <- get_posthoc_table(lmm_wu_microb, "Weighted UniFrac - Washout")

#combine all into a single table
posthoc_all <- bind_rows(
  posthoc_bc_fresh, posthoc_wu_fresh,
  posthoc_bc_ferm, posthoc_wu_ferm,
  posthoc_bc_microb, posthoc_wu_microb
)

#preview and export
head(posthoc_all)
write.csv(posthoc_all, "BetaDiversity_LMM_posthoc_contrasts.csv", row.names = FALSE)


#### DESeq2 LMM ####
#attempted to run the DESeq analysis by adding 1 to avoid all zeros, but ran into error
#code below for DESeq2 does not include +1 step

##### DESeq2 ##### 
#controls for repeated measures
design = ~ participant_id + period

##FRESH
dds_fresh <- phyloseq_to_deseq2(ps_fresh, ~ participant_id + period)

#ensure period is a factor with proper reference level
dds_fresh$period <- relevel(dds_fresh$period, "Base")

#DESeq2
dds_fresh <- DESeq(dds_fresh)

#contrasts
res_fresh_VEG_vs_Base <- results(dds_fresh, contrast = c("period", "VEG", "Base"))
res_fresh_WO1_vs_VEG  <- results(dds_fresh, contrast = c("period", "WO1", "VEG"))


##Ferm
dds_ferm<- phyloseq_to_deseq2(ps_ferm, ~ participant_id + period)

#ensure period is a factor with proper reference level
dds_ferm$period <- relevel(dds_ferm$period, "Base")

#DESeq2
dds_ferm <- DESeq(dds_ferm)

#contrasts
res_ferm_FERM_vs_Base <- results(dds_ferm, contrast = c("period", "FERM", "Base"))
res_ferm_WO2_vs_FERM  <- results(dds_ferm, contrast = c("period", "WO2", "FERM"))

##Microbiome
dds_microb<- phyloseq_to_deseq2(ps_microb, ~ participant_id + period)

#ensure period is a factor with proper reference level
dds_microb$period <- relevel(dds_microb$period, "Base")

#DESeq2
dds_microb <- DESeq(dds_microb)

#contrasts
res_microb_WO1_vs_Base <- results(dds_microb, contrast = c("period", "WO1", "Base"))
res_microb_WO2_vs_WO1  <- results(dds_microb, contrast = c("period", "WO2", "WO1"))

##### Variance-stabilizing Transformation (VST) #####
#prepare each subset for LMM analysis using VST
vst_fresh  <- varianceStabilizingTransformation(dds_fresh, blind = FALSE)
vst_df_fresh  <- as.data.frame(assay(vst_fresh))

vst_ferm   <- varianceStabilizingTransformation(dds_ferm, blind = FALSE)
vst_df_ferm   <- as.data.frame(assay(vst_ferm))

vst_microb <- varianceStabilizingTransformation(dds_microb, blind = FALSE)
vst_df_microb <- as.data.frame(assay(vst_microb))

##### LMM on Top Taxa #####
#convert to long format for merging with metadata
vst_fresh_long <- vst_df_fresh %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_fresh_df, by = "SampleID")

vst_ferm_long <- vst_df_ferm %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_ferm_df, by = "SampleID")

vst_microb_long <- vst_df_microb %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_microb_df, by = "SampleID")


#identify top significant taxa (n=10)
top_n <- 10

#fresh
res_fresh_df <- as.data.frame(res_fresh_VEG_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%  #move rownames to a column
  select(taxon, padj) %>%                #keep only the columns we need
  mutate(padj = as.numeric(padj)) %>%    #ensure padj is numeric
  filter(!is.na(padj)) %>%               #remove NA padj values
  arrange(padj)                          #sort by smallest padj

top_taxa_fresh <- res_fresh_df %>%
  slice_head(n = top_n) %>%              #select top n taxa
  pull(taxon)

#ferm
res_ferm_df <- as.data.frame(res_ferm_FERM_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%
  select(taxon, padj) %>%
  mutate(padj = as.numeric(padj)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top_taxa_ferm <- res_ferm_df %>%
  slice_head(n = top_n) %>%
  pull(taxon)

#microbiome
res_microb_df <- as.data.frame(res_microb_WO1_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%
  select(taxon, padj) %>%
  mutate(padj = as.numeric(padj)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top_taxa_microb <- res_microb_df %>%
  slice_head(n = top_n) %>%
  pull(taxon)


#filter the VST-long datasets for top taxa
vst_fresh_long_top <- vst_fresh_long %>%
  filter(taxon %in% top_taxa_fresh)
vst_ferm_long_top <- vst_ferm_long %>%
  filter(taxon %in% top_taxa_ferm)
vst_microb_long_top <- vst_microb_long %>%
  filter(taxon %in% top_taxa_microb)


#run LMM per taxon
run_lmm <- function(df) {
  df %>%
    group_by(taxon) %>%
    group_modify(~{
      tryCatch({
        m <- lmer(vst_abundance ~ period + (1 | participant_id), data = .x)
        emm <- emmeans(m, pairwise ~ period)
        tibble(
          model = list(m),
          emm = list(emm),
          contrasts = list(as.data.frame(emm$contrasts))
        )
      }, error = function(e) {
        tibble(model = NA, emm = NA, contrasts = NA)
      })
    })
}

lmm_results_fresh  <- run_lmm(vst_fresh_long_top)
lmm_results_ferm   <- run_lmm(vst_ferm_long_top)
lmm_results_microb <- run_lmm(vst_microb_long_top)

#double-check the periods used in LMM
unique(vst_fresh_long_top$period)
unique(vst_ferm_long_top$period)
unique(vst_microb_long_top$period)


#double-check if LMM was sucessful 
lmm_results_fresh %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )

lmm_results_ferm %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )

lmm_results_microb %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )


##### Plots #####
#function to tidy EMMs for plotting
emm_fresh <- lmm_results_fresh %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))

emm_ferm <- lmm_results_ferm %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))

emm_microb <- lmm_results_microb %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))


#create function to merge taxa
add_taxon_labels_safe <- function(lmm_results, physeq_obj) {
  #extract taxonomy table
  tax_df <- as.data.frame(tax_table(physeq_obj))
  tax_df$taxon <- rownames(tax_df)
  
  #extract emmeans from LMM
  emm_df <- lmm_results %>%
    dplyr::filter(!is.na(model)) %>%
    dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
    dplyr::select(taxon, emm_df) %>%
    tidyr::unnest(cols = c(emm_df))
  
  #merge taxonomy
  emm_named <- emm_df %>%
    left_join(tax_df %>% select(taxon, Family, Genus), by = "taxon") %>%
    #create label with fallback and uniqueness
    mutate(
      taxon_label = case_when(
        !is.na(Genus) ~ paste(Family, Genus, sep = "_"),
        !is.na(Family) ~ Family,
        TRUE ~ taxon
      ),
      #append OTU ID to make labels unique
      taxon_label = paste0(taxon_label, "_", taxon)
    )
  
  return(emm_named)
}


#apply function to subset
emm_fresh_named  <- add_taxon_labels_safe(lmm_results_fresh, ps_fresh)
emm_ferm_named   <- add_taxon_labels_safe(lmm_results_ferm, ps_ferm)
emm_microb_named <- add_taxon_labels_safe(lmm_results_microb, ps_microb)


#line plot
emm_line_plot_taxa <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean, group = taxon_label, color = taxon_label)) +
    geom_point(size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    theme_bw() +
    labs(
      y = "Estimated VST Abundance",
      x = "Period",
      title = title,
      color = "Taxon"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

#generate and save 
deseq_fresh_taxa <- emm_line_plot_taxa(emm_fresh_named, title = "Top 10 Significant Taxa: Fresh")
deseq_ferm_taxa   <- emm_line_plot_taxa(emm_ferm_named, title = "Top 10 Significant Taxa: Fermented")
deseq_microb_taxa <- emm_line_plot_taxa(emm_microb_named, title = "Top 10 Significant Taxa: Microbiome Washout")

ggsave("A2_LMM_DESeq2_Fresh.png", plot = deseq_fresh_taxa, width = 8, height = 6, dpi = 300)
ggsave("A2_LMM_DESeq2_Ferm.png", plot = deseq_ferm_taxa, width = 8, height = 6, dpi = 300)
ggsave("A2_LMM_DESeq2_Microbiome.png", plot = deseq_microb_taxa, width = 8, height = 6, dpi = 300)


#facet plots
emm_line_plot_taxa_facet <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean)) +
    geom_point(size = 3, color = "steelblue") +
    geom_line(aes(group = 1), color = "steelblue") +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    facet_wrap(~ taxon_label, scales = "free_y") +
    theme_bw() +
    labs(
      y = "Estimated VST Abundance",
      x = "Period",
      title = title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#generate plots and save
deseq_fresh_facet  <- emm_line_plot_taxa_facet(emm_fresh_named, title = "Top 10 Significant Taxa: Fresh")
deseq_ferm_facet   <- emm_line_plot_taxa_facet(emm_ferm_named, title = "Top 10 Significant Taxa: Fermented")
deseq_microb_facet <- emm_line_plot_taxa_facet(emm_microb_named, title = "Top 10 Significant Taxa: Microbiome Washout")

ggsave("A2_LMM_DESeq2_Fresh_facet.png", plot = deseq_fresh_facet, width = 8, height = 6, dpi = 300)
ggsave("A2_LMM_DESeq2_Ferm_facet.png", plot = deseq_ferm_facet, width = 8, height = 6, dpi = 300)
ggsave("A2_LMM_DESeq2_Microbiome_facet.png", plot = deseq_microb_facet, width = 8, height = 6, dpi = 300)



#test
#save as CSV
emm_fresh_named <- emm_fresh_named %>% mutate(Subset = "Fresh")
emm_ferm_named <- emm_ferm_named %>% mutate(Subset = "Fermented")
emm_microb_named <- emm_microb_named %>% mutate(Subset = "Microbiome")

all_top_taxa_points <- bind_rows(emm_fresh_named, emm_ferm_named, emm_microb_named)

write.csv(all_top_taxa_points, "Top10_Taxa_All_Subsets_Points.csv", row.names = FALSE)


###combine graphs
#add the subset info and combine

#fresh
emm_df_fresh_s$metric <- "Shannon"
emm_df_fresh_s$subset <- "Fresh"

emm_df_fresh_f$metric <- "Faith_PD"
emm_df_fresh_f$subset <- "Fresh"

#ferm
emm_df_ferm_s$metric <- "Shannon"
emm_df_ferm_s$subset <- "Fermented"

emm_df_ferm_f$metric <- "Faith_PD"
emm_df_ferm_f$subset <- "Fermented"

#microbiome washouts
emm_df_microb_s$metric <- "Shannon"
emm_df_microb_s$subset <- "Microbiome_Washout"

emm_df_microb_f$metric <- "Faith_PD"
emm_df_microb_f$subset <- "Microbiome_Washout"

#combine all data
plot_data <- bind_rows(
  emm_df_fresh_s, emm_df_fresh_f,
  emm_df_ferm_s, emm_df_ferm_f,
  emm_df_microb_s, emm_df_microb_f
)

plot_data

#reorder and change label names
period_labels <- c(
  "Base" = "Base",
  "VEG"  = "VEG",
  "FERM" = "FERM",
  "WO1"  = "WO1",
  "WO2"  = "WO2"
)

metric_labels <- c(
  "Shannon" = "Shannon",
  "Faith's PD" = "Faith's PD"
)

combined_plot <- ggplot(plot_data, aes(x = period, y = emmean, color = period)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  facet_grid(metric ~ subset, scales = "free_y",
             labeller = labeller(metric = metric_labels)) +
  scale_x_discrete(labels = period_labels) +
  ylab("Alpha Diversity") +
  xlab("Period") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey90", color = "black", size = 0.5)
  )


ggsave("LMM_Alpha_combined.png", plot = combined_plot, width = 6, height = 4, dpi = 300)

A2_LMM_DESeq2_Top10_Fresh.png


#### Figure Plots ####
##Legend 
#create a legend (couldn't extract all periods color as one legend)
period_colors <- c(
  "Base" = "#4E79A7",
  "VEG" = "#B699C6",
  "WO1" = "#F28E2B",
  "FERM" = "#59A14F",
  "WO2" = "#E15759")

#create a data frame to generate a legend, y is arbitary here because we are only interested in the legend
legend_data <- data.frame(
  period = factor(names(period_colors), levels = c("Base", "VEG", "WO1", "FERM", "WO2")),
  y = 1
)

#generate a ggplot object solely to create the legend
legend_plot <- ggplot(legend_data, aes(x = period, y = y, fill = period)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = period_colors, name = "Period") +  # colors remain the same
  theme_void() +
  theme(legend.position = "right")  
#extract the legend from the ggplot object
shared_legend <- get_legend(legend_plot)
#wrapextracted legend so it can be treated as a patchwork element
legend_wrap <- wrap_elements(full = shared_legend)

#save legend
ggsave("Figure_3_LMM_Alpha_Faith_combined_legend.png", plot = legend_wrap, width = 6, height = 4, dpi = 300)



##LMM alpha diversity (only Faith's PD)
"Base" = "#4E79A7",
"VEG" = "#B699C6",
"WO1" = "#F28E2B",
"FERM" = "#59A14F",
"WO2" = "#E15759")


LMM_fresh_veg_faith_box_plot <- ggplot(emm_df_fresh_f, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FRF,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.1
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_fill_manual(values = c("Base" = "#4E79A7", "VEG" = "#B699C6", "WO1" = "#F28E2B"))



LMM_ferm_veg_faith_box_plot <- ggplot(emm_df_ferm_f, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FEF,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.15
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )  +
  scale_fill_manual(values = c("Base" = "#4E79A7", "FERM" = "#59A14F", "WO2" = "#E15759"))


LMM_microbiome_faith_box_plot <- ggplot(emm_df_microb_f, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_MF,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.15
  ) +
  
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Microbiome – Faith's PD") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )  +
  scale_fill_manual(values = c("Base" = "#4E79A7", "WO1" = "#F28E2B", "WO2" = "#E15759"))

#plot all three Faith LMM 
F3_combined_plot <- LMM_microbiome_faith_box_plot + plot_spacer() + LMM_fresh_veg_faith_box_plot + plot_spacer() + LMM_ferm_veg_faith_box_plot +
  plot_layout(ncol = 5, widths = c(1, 0.1, 1, 0.1, 1))
F3_combined_plot

#save plot
ggsave("Figure_3_LMM_Alpha_Faith_combined.png", plot = F3_combined_plot, width = 6, height = 4, dpi = 300)

##LMM alpha diversity (Shannon for S1)

LMM_fresh_veg_shannon_box_plot <- ggplot(emm_df_fresh_s, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  #EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  #add significance labels
  stat_pvalue_manual(
    pval_table_FRS,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.1
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_fill_manual(values = c("Base" = "#4E79A7", "VEG" = "#B699C6", "WO1" = "#F28E2B"))

LMM_ferm_veg_shannon_box_plot <- ggplot(emm_df_ferm_s, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_FES,
    label = "label",
    tip.length = 0.03,
    step.increase = .5
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_fill_manual(values = c("Base" = "#4E79A7", "FERM" = "#59A14F", "WO2" = "#E15759"))



LMM_microbiome_shannon_box_plot <- ggplot(emm_df_microb_s, aes(x = period, y = emmean)) +
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.7, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table_MS,
    label = "label",
    tip.length = 0.03,
    step.increase = -0.1
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Microbiome – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_fill_manual(values = c("Base" = "#4E79A7", "WO1" = "#F28E2B", "WO2" = "#E15759"))

S1_LMM_Shannon_combined_plot <- LMM_microbiome_shannon_box_plot + 
  plot_spacer() + 
  LMM_fresh_veg_shannon_box_plot + 
  plot_spacer() + 
  LMM_ferm_veg_shannon_box_plot+
  plot_layout(ncol = 5, widths = c(1, 0.1, 1, 0.1, 1))
S1_LMM_Shannon_combined_plot

#save plot
ggsave("Sup_Figure_1_LMM_Alpha_Shannon_combined.png", plot = S1_LMM_Shannon_combined_plot, width = 6, height = 4, dpi = 300)


##LMM beta diversity (only Bray-Curtis)
# Combine plots
beta_combined_BC <- (beta_fresh_BC) /
  (beta_ferm_BC)  /
  (beta_microbiome_BC) + 
  plot_annotation(title = "Bray-Curtis Beta Diversity LMM Analysis")

# Save
ggsave("Beta_BrayCurtis_Combined.png", plot = beta_combined_BC, width = 6, height = 12, dpi = 300)


beta_combined_BC <- (beta_fresh_BC) /
  (beta_ferm_BC)  /
  (beta_microbiome_BC) + 
  plot_annotation(title = "Bray-Curtis Beta Diversity LMM Analysis")

# Save the combined figure
ggsave("Beta_Diversity_Combined.png", plot = beta_combined, width = 12, height = 12, dpi = 300)




wu_posthoc <- posthoc_all %>%
  filter(grepl("Weighted UniFrac", response)) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    direction = ifelse(estimate > 0, "away from baseline", "toward baseline")
  ) %>%
  select(Metric = response,
         Contrast = contrast,
         Estimate = estimate,
         SE = SE,
         t_ratio = t.ratio,
         p_value = p.value,
         Significance = significance,
         Direction = direction)

head(wu_posthoc)

