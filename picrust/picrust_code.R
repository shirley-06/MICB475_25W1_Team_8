#Install packages
#BiocManager::install('ALDEx2')

library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggpicrust2)
library(ape)
library(picante)
library(vegan)

#Loading files
meta <- read_tsv(file="fermentation_metadata.tsv")
ko_data <- read_tsv(file="pred_metagenome_unstrat.tsv")
pathway_data <- read_tsv(file="path_abun_unstrat.tsv")

#Formatting imported tables
meta <- meta %>%
  rename('sample_id'='sample-id')

pathway_data <- as.data.frame(pathway_data)

#####MetaCyc Analysis#####

###Pulling sample IDs from tables###
#Filtering Metadata for P2 and P4#
metadata_p2p4 <- meta %>%
  filter(period %in% c("FERM","VEG"))

#Metadata sample IDs#
meta_samp_ids_p2p4 <- select(metadata_p2p4, sample_id) %>%
  as.data.frame()

#Pathway/MetaCyc sample IDs#
pathway_samp_ids <- colnames(pathway_data) %>%
  as.data.frame()
pathway_samp_ids <- pathway_samp_ids[-1,] %>%
  as.data.frame()
colnames(pathway_samp_ids)[1] <- "sample_id"


##Creating list of shared sample IDs betweeen metadata and metacyc, for P2/P4##
common_ids_pwy <- intersect(pathway_samp_ids, meta_samp_ids_p2p4)
colnames(common_ids_pwy)[1] <- "sample_name"
print(common_ids_pwy)

###Removing missing sample IDs from metadata###
meta_removed_p2p4 <- subset(metadata_p2p4, sample_id %in% common_ids_pwy$sample_name)

###Removing missing sample IDs from pathway###
pathway_data <- column_to_rownames(pathway_data, var ="pathway")
common_ids_pwy_list <- as.character(common_ids_pwy$sample_name)
pathway_removed_p2p4 <- pathway_data[, common_ids_pwy_list]

#Running differential analysis

daa_results_pwy = pathway_daa(abundance = pathway_removed_p2p4,
                          metadata = meta_removed_p2p4,
                          group = "period",
                          daa_method = "ALDEx2")

####Running same protocol - filtered tables for FERM and VEG####

#Using same metadata objects from MetaCyc analysis












