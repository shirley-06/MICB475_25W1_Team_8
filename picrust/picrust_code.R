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

###Removing missing sample IDs###
meta_removed <- anti_join(meta, missing_ids, by = "sample_name")

pathway_data <- column_to_rownames(pathway_data, var = "pathway")

###Performing differential analysis for pathways###

daa_results_pwy = pathway_daa(abundance = pathway_data,
                          metadata = meta_removed,
                          group = "period",
                          daa_method = "ALDEx2")

saveRDS(daa_results_pwy, file = "daa_results_pwy.rds")

annotated_daa_pwy <- pathway_annotation(pathway = "MetaCyc",
                   daa_results_df = daa_results,
                   ko_to_kegg = FALSE)

pathway_graph <- pathway_errorbar(abundance = daa_results,
                  daa_results_df = pathway_ann_daa,
                  Group = meta_removed$period,
                  p_values_threshold = 5e-11,
                  order = "pathway class",
                  ko_to_kegg = FALSE,
                  p_value_bar = TRUE,
                  x_lab = "pathway_name")


####Trying with filtered tables for FERM and VEG####

##Creating filtered metadata with non-matching sample IDs + non-relevant period removed##
metadata_removed_filt <- meta_removed %>%
  filter(period %in% c("FERM","VEG"))

#Extracting list of sample IDs - FERM/VEG period#
meta_cols_filt <- metadata_removed_filt %>%
  select(sample_name)
meta_cols_filt <- as.data.frame(meta_cols_filt)
colnames(meta_cols_filt)[1] <- "sample_id"

pathway_cols_filt <- colnames(pathway_data)
pathway_cols_filt <- as.data.frame(pathway_cols_filt)
colnames(pathway_cols_filt)[1] <- "sample_id"

#Creating list of intersecting samples between period-filtered metadata and pathway_data#
meta_removed_filt <- intersect(pathway_cols_filt, meta_cols_filt)
meta_removed_filt <- unlist(meta_removed_filt)

#Selecting sample IDs in pathway_data that match filtered metadata
pathway_data_filt <- pathway_data[,meta_removed_filt]

daa_results_pwy_filtered = pathway_daa(abundance = pathway_data_filt,
                              metadata = metadata_removed_filt,
                              group = "period",
                              daa_method = "ALDEx2")







