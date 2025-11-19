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
model_shannon <- lmer(Shannon ~ period + (1 | participant_id), data = meta_fresh)

#ANOVA for fixed effects
anova(model_shannon)

#post-hoc tests (period comparisons)
emmeans(model_shannon, pairwise ~ period)





