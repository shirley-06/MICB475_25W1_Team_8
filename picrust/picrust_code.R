#Install packages
#BiocManager::install('ALDEx2')

library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggpicrust2)
library(ape)
library(picante)
library(vegan)

#loading files
meta <- read_tsv(file="fermentation_metadata.tsv")
ko_data <- read_tsv(file="combined_KO_predicted.tsv")
ec_data <- read_tsv(file="combined_EC_predicted.tsv")
pathway_data <- read_tsv(file="path_abun_unstrat.tsv")

#Transposing ko_data / swapping rows and columns
transposed_ko_matrix <- t(ko_data)
ko_transposed <- as.data.frame(transposed_ko_matrix)
colnames(ko_transposed) <- ko_transposed[1, ]
ko_transposed <- ko_transposed[-1, ]

#Formatting tables
meta <- meta %>%
  rename('sample_name'='sample-id')

pathway_data <- as.data.frame(pathway_data)

#Formatting tables to remove missing sample IDs
meta_cols <- select(meta, sample_name)
meta_cols <- as.data.frame(meta_cols)
colnames(meta_cols)[1] <- "sample_id"

pathway_cols <- colnames(pathway_data)
pathway_cols <- as.data.frame(pathway_cols)
pathway_cols <- pathway_cols[-1,]
pathway_cols <- as.data.frame(pathway_cols)
colnames(pathway_cols)[1] <- "sample_id"

missing_ids <- setdiff(meta_cols, pathway_cols)
colnames(missing_ids)[1] <- "sample_name"
print(missing_ids)

meta_removed <- anti_join(meta, missing_ids, by = "sample_name")

pathway_data <- column_to_rownames(pathway_data, var = "pathway")

#Performing differential analysis for pathways
daa_results = pathway_daa(abundance = pathway_data,
                          metadata = meta_removed,
                          group = "period",
                          daa_method = "ALDEx2")


