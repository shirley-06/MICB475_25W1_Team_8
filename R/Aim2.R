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

#load filtered and rarefied data
load("ferm_rare.RData")

#check object loaded
ls()

#check metadata columns
colnames(sample_data(ferm_rarefac))


#subset groups
#fresh vegetables: P1 -> P2 -> P3
ps_fresh <- subset_samples(ferm_rarefac, period %in% c("Base","VEG","WO1"))

#fermented vegetables: P1 -> P4 -> P5
ps_ferm <- subset_samples(ferm_rarefac, period %in% c("Base","FERM","WO2"))

#fermented vegetables: P1 -> P4 -> P5
ps_microb <- subset_samples(ferm_rarefac, period %in% c("Base","FERM","WO2"))


#remove taxa with zero counts after subsetting 
ps_fresh <- prune_taxa(taxa_sums(ps_fresh) > 0, ps_fresh)
ps_ferm  <- prune_taxa(taxa_sums(ps_ferm)  > 0, ps_ferm)
ps_microb  <- prune_taxa(taxa_sums(ps_microb)  > 0, ps_microb)

#### Alpha Diversity ####