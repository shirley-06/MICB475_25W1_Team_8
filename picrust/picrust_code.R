#Install packages
#install.packages("pheatmap")

library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggpicrust2)
library(ape)
library(picante)
library(vegan)
library(GGally)
library(pheatmap)
library(ggplot2)
library(scales)
library(stringr)

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

###This is not needed for every dataset##
###Some sample ID discrepancies between our tables need to be resolved before using ggpicrust###

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

###Running differential analysis with LinDA###

daa_results_pwy = pathway_daa(abundance = pathway_removed_p2p4,
                          metadata = meta_removed_p2p4,
                          group = "period",
                          daa_method = "LinDA")

pwy_annotated <- pathway_annotation(pathway = "MetaCyc",
                   daa_results_df = daa_results_pwy,
                   ko_to_kegg = FALSE)

#Checking number of significant hits#
significant_pwy <- pwy_annotated %>% 
  filter(p_adjust < 0.05)
#Did not include log2 FC - would only have 1 result for heatmap

##Plotting log2 FC with ggplot

L2FC_metacyc <- ggplot(data = significant_pwy, aes(x=description, 
                                             y=log2FoldChange)) +
  geom_col(fill="steelblue") +
  labs(x="Pathway Name", y= "Log2 Fold Change") +
  coord_flip()

ggsave("L2FC_plot.png", plot = L2FC_metacyc, width=9, height=5)

##Generating MetaCyc PCA plot##
colnames(meta_removed_p2p4)[1] <- "sample_name"

pca_pwy <- pathway_pca(abundance = pathway_removed_p2p4,
                       metadata = meta_removed_p2p4,
                       group = "period")

ggsave("pca_pwy.png", plot=pca_pwy, width=6, height=6, units="in")

##Generating heatmap##
pwy_sig_features <- significant_pwy %>%
  pull("feature") %>%
  unique()

pwy_rel_abun <- pathway_data %>%
  apply(2,function(x) x/sum(x)) %>%
  as.data.frame()

pwy_stats = pwy_rel_abun %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(pwy_sig_features)) %>%
  cor(method = "spearman")

color_palette <- colorRampPalette(c("blue","white","red"))(40)
breaks <- seq(-1, 1, length.out=41)

pwy_heatmap <- pheatmap(pwy_stats,
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",
                       clustering_method = "complete",
                       color = color_palette, breaks = breaks,
                       main = "",
                       fontsize_row = 10, fontsize_col = 10,
                       filename = "pwy_heatmap.png",height= 9, width = 9)


####Running same protocol - filtered tables for FERM and VEG####
#Using the same metadata objects from MetaCyc analysis#

#KO sample IDs#
ko_data <- column_to_rownames(ko_data, var = "function")
ko_samp_ids <- colnames(ko_data) %>%
  as.data.frame()
colnames(ko_samp_ids)[1] <- "sample_id"

##Creating list of shared sample IDs betweeen metadata and ko, for P2/P4##
common_ids_ko <- intersect(ko_samp_ids, meta_samp_ids_p2p4)
colnames(common_ids_ko)[1] <- "sample_name"
print(common_ids_ko)

###Removing missing sample IDs from metadata###
meta_removed_p2p4_ko <- subset(metadata_p2p4, sample_id %in% common_ids_ko$sample_name)

###Removing missing sample IDs from pathway###
common_ids_ko_list <- as.character(common_ids_ko$sample_name)
ko_removed_p2p4 <- ko_data[, common_ids_ko_list]

###Running differential analysis with LinDA###

daa_results_ko = pathway_daa(abundance = ko_removed_p2p4,
                                   metadata = meta_removed_p2p4_ko,
                                   group = "period",
                                   daa_method = "LinDA",
                                   select=NULL, reference=NULL)

ko_annotated <- pathway_annotation(pathway = "KO",
                                   daa_results_df = daa_results_ko,
                                   ko_to_kegg = TRUE)

#saveRDS(ko_annotated, file = "ko_annotated.rds")

ko_annotated = readRDS("ko_annotated.rds")

#Checking number of significant hits#
significant_ko <- ko_annotated %>% 
  filter(p_adjust < 5e-2, abs(log2FoldChange)>2)

#Manually separating duplicate uncharacterized proteins for visualization#
significant_ko[10, "pathway_name"] <- "uncharacterized protein "

##Generating log2 FC plot with ggplot##

L2FC_ko <- ggplot(data = significant_ko, aes(x=pathway_name, 
                                             y=log2FoldChange, 
                                             fill=str_wrap(pathway_description, width=40))) +
  geom_col() +
  labs(x="Pathway Name", y= "Log2 Fold Change", fill="Pathway Description") +
  guides(fill = guide_legend(byrow = TRUE)) + 
  theme(legend.key.spacing.y = unit(0.2,"cm")) + #Adjusting spacing between legend keys#
  coord_flip()

ggsave("L2FC_ko.png", plot=L2FC_ko, width=12, height=6, units="in")

##Generating KEGG PCA plot##
colnames(meta_removed_p2p4_ko)[1] <- "sample_name"

pca_ko <- pathway_pca(abundance = ko_removed_p2p4,
                      metadata = meta_removed_p2p4_ko,
                      group = "period")

ggsave("pca_ko.png", plot=pca_ko, width=6, height=6, units="in")

##Generating heatmap##
ko_sig_features <- significant_ko %>%
  pull("feature") %>%
  unique()

ko_rel_abun <- ko_data %>%
  apply(2,function(x) x/sum(x)) %>%
  as.data.frame()

ko_stats = ko_rel_abun %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(ko_sig_features)) %>%
  cor(method = "spearman")

color_palette <- colorRampPalette(c("blue","white","red"))(40)
breaks <- seq(-1, 1, length.out=41)

ko_heatmap <- pheatmap(ko_stats,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = color_palette, breaks = breaks,
         main = "",
         fontsize_row = 10, fontsize_col = 10,
         filename = "ko_heatmap.png",height= 9, width = 9)


