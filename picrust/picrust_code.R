
library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggpicrust2)
library(ape)
library(picante)
library(vegan)

meta <- read_tsv(file="fermentation_metadata.tsv")
ko_data <- read_tsv(file="combined_KO_predicted.tsv")
ec_data <- read_tsv(file="combined_EC_predicted.tsv")
pathway_data <- read_tsv(file="path_abun_unstrat.tsv")

transposed_ko_matrix <- t(ko_data)
ko_transposed <- as.data.frame(transposed_ko_matrix)
colnames(ko_transposed) <- ko_transposed[1, ]
ko_transposed <- ko_transposed[-1, ]

meta <- meta %>%
  rename('sample_name'='sample-id')

pathway_data <- as.data.frame(pathway_data)

meta_cols <- select(meta, sample_name)
meta_cols <- as.data.frame(meta_cols)
colnames(meta_cols)[1] <- "sample_id"

pathway_cols <- colnames(pathway_data)
pathway_cols <- as.data.frame(pathway_cols)
pathway_cols <- pathway_cols[-1,]
pathway_cols <- as.data.frame(pathway_cols)
colnames(pathway_cols)[1] <- "sample_id"

missing_ids <- setdiff(meta_cols, pathway_cols)
colnames(missing_ids)[1] <- "sample_id"
print(missing_ids)