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


#load filtered and rarefied data
load("ferm_rare.RData")


#check object loaded
ls()
ferm_rarefac
sample_data(ferm_rarefac) %>% head()




#### Alpha Diversity ####
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

#### Linear Mixed Model ####

##Fresh Shannon
model_shannon_fresh <- lmer(Shannon ~ period + (1 | participant_id), data = meta_fresh)

#ANOVA for fixed effects
anova(model_shannon_fresh)

#post-hoc tests (period comparisons)
emmeans(model_shannon_fresh, pairwise ~ period)

emm_fresh_shannon <- emmeans(model_shannon_fresh, ~ period)
emm_df_fresh <- as.data.frame(emm_fresh_shannon)
ggplot(emm_df_fresh, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  theme_minimal()

##Fresh Faith
model_faith_fresh <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_fresh)

# ANOVA
anova(model_faith_fresh)

# Post-hoc contrasts
emmeans(model_faith_fresh, pairwise ~ period)

# Plot
emm_fresh_faith <- emmeans(model_faith_fresh, ~ period)
emm_df_faith_fresh <- as.data.frame(emm_fresh_faith)
ggplot(emm_df_faith_fresh, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  theme_minimal()


##Ferm Shannon
model_shannon_ferm <- lmer(Shannon ~ period + (1 | participant_id), data = meta_ferm)
anova(model_shannon_ferm)
emmeans(model_shannon_ferm, pairwise ~ period)

# Plot
emm_shannon_ferm <- emmeans(model_shannon_ferm, ~ period)
emm_df_shannon_ferm <- as.data.frame(emm_shannon_ferm)
ggplot(emm_df_shannon_ferm, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  theme_minimal()

##Ferm Faith
model_faith_ferm <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_ferm)
anova(model_faith_ferm)
emmeans(model_faith_ferm, pairwise ~ period)

emm_faith_ferm <- emmeans(model_faith_ferm, ~ period)
emm_df_faith_ferm <- as.data.frame(emm_faith_ferm)
ggplot(emm_df_faith_ferm, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  theme_minimal()


##Microbiome Shannon
model_shannon_microb <- lmer(Shannon ~ period + (1 | participant_id), data = meta_microb)
anova(model_shannon_microb)
emmeans(model_shannon_microb, pairwise ~ period)

emm_shannon_microb <- emmeans(model_shannon_microb, ~ period)
emm_df_shannon_microb <- as.data.frame(emm_shannon_microb)
ggplot(emm_df_shannon_microb, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  theme_minimal()

##Microbiome Faith
model_faith_microb <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_microb)
anova(model_faith_microb)
emmeans(model_faith_microb, pairwise ~ period)

emm_faith_microb <- emmeans(model_faith_microb, ~ period)
emm_df_faith_microb <- as.data.frame(emm_faith_microb)
ggplot(emm_df_faith_microb, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  theme_minimal()

###Combine Graph
#add the subset info and combine

#fresh
emm_df_fresh$metric <- "Shannon"
emm_df_fresh$subset <- "Fresh"

emm_df_faith_fresh$metric <- "Faith_PD"
emm_df_faith_fresh$subset <- "Fresh"

#ferm
emm_df_shannon_ferm$metric <- "Shannon"
emm_df_shannon_ferm$subset <- "Fermented"

emm_df_faith_ferm$metric <- "Faith_PD"
emm_df_faith_ferm$subset <- "Fermented"

#microbiome washouts
emm_df_shannon_microb$metric <- "Shannon"
emm_df_shannon_microb$subset <- "Microbiome Washout"

emm_df_faith_microb$metric <- "Faith_PD"
emm_df_faith_microb$subset <- "Microbiome Washout"

#combine all data
plot_data <- bind_rows(
  emm_df_fresh, emm_df_faith_fresh,
  emm_df_shannon_ferm, emm_df_faith_ferm,
  emm_df_shannon_microb, emm_df_faith_microb
)

#plot all
ggplot(plot_data, aes(x=period, y=emmean, color=period)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  facet_grid(metric ~ subset, scales = "free_y") +
  ylab("Alpha Diversity") +
  xlab("Period") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, size=1), # add border
    strip.background = element_rect(fill="grey90", color="black", size=0.5) # optional: facet label background
  )