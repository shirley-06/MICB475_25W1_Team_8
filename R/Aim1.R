#Aim 1 - Assess impact of fresh versus fermented vegetable consumption on gut microbial diversity and composition

#!/usr/bin/env Rscript
# NOTE: to install phyloseq, please use the following code instead of the usual "install.packages" function:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#load library
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library (picante)
library (ggpubr)
library(DESeq2)
library(forcats)
library(patchwork)

#### Generate Phyloseq Object ####
##### Import Data #####
#change file paths as necessary (import as tibble)
#import metadata
metafp <- "Export/fermentation_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

#import OTU table
otufp <- "Export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

#import taxonomy table
taxfp <- "Export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

#import phylogenetic tree
phylotreefp <- "Export/tree.nwk"
phylotree <- read.tree(phylotreefp)

##### Format Data #####
###### OTU Table ###### 
# OTU tables should be a matrix with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

#save everything except first column (OTU ID) into a matrix 
#take the otu data and index it to remove the first column
otu_mat <- as.matrix(otu[,-1])
#make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
#use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

###### Sample Metadata ###### 
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

###### Taxonomy ###### 
#convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
#save everything except feature IDs 
tax_mat <- tax_mat[,-1]
#make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
#make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

##### Create phyloseq Object #####
#merge all objects into a phyloseq object
fermentation_ps_raw <- phyloseq(OTU, SAMP, TAX, phylotree)

#view components of phyloseq object with the following commands
otu_table(fermentation_ps_raw)
sample_data(fermentation_ps_raw)
tax_table(fermentation_ps_raw)
phy_tree(fermentation_ps_raw)

##### Filter phyloseq Object #####
#remove unwanted taxa: mitochondria, chloroplasts, Eukaryota, and Archaea
ps_taxa_filter <- subset_taxa(
  fermentation_ps_raw, 
  !Order %in% c("o__Chloroplast") &   # remove chloroplasts
    !Family %in% c("f__Mitochondria") & # remove mitochondria
    !Domain %in% c("d__Eukaryota", "d__Archaea","Unassigned") # remove Eukaryota and Archaea
)


#subset to only 'female' samples from host_sex category and remove 'NA' samples in Group category
ps_metadata_filter <- subset_samples (ps_taxa_filter, 
                             host_sex == "female" &
                            Group != "not applicable")

#create a new column and group 'Control' and 'Antibiotics' group as one group
sample_data(ps_metadata_filter)$Group_new <- ifelse(
  sample_data(ps_metadata_filter)$Group %in% c("Control", "Antibiotics"),
  "CTRL_AB",
  sample_data(ps_metadata_filter)$Group
)

ps_filter <- ps_metadata_filter


#view data set after filtering
#look at the taxa table
table(tax_table(fermentation_ps_raw)[, "Domain"])
table(tax_table(ps_filter)[, "Domain"])

#look at metadata table
table(sample_data(ps_filter)$Group_new)
table(sample_data(fermentation_ps_raw)$host_sex)
table(sample_data(ps_filter)$host_sex)

#look at number of samples before and after filtering
nsamples(fermentation_ps_raw)
nsamples(ps_filter)
#% of samples kept --> ~80% were kept after filtering
nsamples(ps_filter) / nsamples(fermentation_ps_raw)

#look at number of ASVs before and after filtering
ntaxa(fermentation_ps_raw)
ntaxa(ps_filter)
#% of ASVs kept --> ~99.3% were kept after filtering
ntaxa (ps_filter) / ntaxa(fermentation_ps_raw) 

#make helper function to look at the top 5 OTUs from any phyloseq object before and after filtering
get_top_otus <- function(ps_obj, n = 5) {
  otu_df <- as.data.frame(as.table(as.matrix(otu_table(ps_obj))))
  colnames(otu_df) <- c("OTU", "Sample", "Count")
  otu_df[order(-otu_df$Count), ][1:n, ]
}

#look at top 5 OTUs before and after filtering
top5_raw <- get_top_otus(fermentation_ps_raw)
top5_filtered <- get_top_otus(ps_filter)

#print results
top5_raw
top5_filtered


#### Rarefaction ####
#summarize sequencing depth per group
#convert sample_data to a data frame and add total reads per sample
sample_df <- as.data.frame(sample_data(ps_filter))
sample_df$TotalReads <- sample_sums(ps_filter)

#summarize min and max sequencing depth per Group_new to determine rarefaction depth
sample_df %>%
  group_by(Group_new) %>%
  summarise(
    MinReads = min(TotalReads),
    MaxReads = max(TotalReads),
    MeanReads = mean(TotalReads),
    .groups = "drop"
  )
#rarefaction depth (minimum reads = 31218)

#clean OTU table to remove zero-count OTUs and samples
ps_clean <- prune_taxa(taxa_sums(ps_filter) > 0, ps_filter)
ps_clean <- prune_samples(sample_sums(ps_clean) > 0, ps_clean)

#set rngseed the same number each time to reproduce this exact analysis
ferm_rarefac <- rarefy_even_depth(ps_clean,
                               rngseed = 1,
                               sample.size = 31218,
                               replace = FALSE,
                               trimOTUs = TRUE,
                               verbose = TRUE)

#check total reads per sample — all should equal 31218
table(sample_sums(ferm_rarefac))

#check for how many OTUs remain
ntaxa(ps_filter) #after filtering
ntaxa(ps_clean) #after cleaning (before rarefaction)
ntaxa(ferm_rarefac) #after rarefaction

#check how many samples remain
nsamples(ps_filter) #after filtering
nsamples(ps_clean) #after cleaning (before rarefaction)
nsamples(ferm_rarefac) #after rarefaction

#confirm your grouping variable is still intact
table(sample_data(ferm_rarefac)$Group_new)


##### Save Raw and Rarefaction files #####
save(ps_filter, file="ferm_final.RData")
save(ferm_rarefac, file="ferm_rare.RData")


#### Alpha Diversity ####
##### Shannon and Faith's PD Calculations #####
#Shannon
#compute Shannon
alpha_shannon <- estimate_richness(ferm_rarefac, measures = "Shannon")

#Faith's PD
#extract OTU table from phyloseq object and convert to matrix
#if taxa are stored as rows in the phyloseq object, transpose the matrix
otu_mat <- as(otu_table(ferm_rarefac), "matrix")
if (taxa_are_rows(ferm_rarefac)) {
  otu_mat <- t(otu_mat)
}

#compute Faith's PD
alpha_faith_pd <- pd(otu_mat, phy_tree(ferm_rarefac), include.root = TRUE)


#add Shannon and Faith PD metrics to metadata as new columns
sd <- as.data.frame(sample_data(ferm_rarefac))
sd$Shannon <- alpha_shannon$Shannon[match(rownames(sd), rownames(alpha_shannon))]
sd$Faith_PD <- alpha_faith_pd$PD[match(rownames(sd), rownames(alpha_faith_pd))]

#convert sample_data back into the phyloseq object
sample_data(ferm_rarefac) <- sample_data(sd)

#verify
head(sample_data(ferm_rarefac))
sample_variables(ferm_rarefac)

##### Shannon and Faith's PD Plots #####
#extract and convert the sample data from phyloseq to a true data.frame
sd_df <- data.frame(sample_data(ferm_rarefac), check.names = TRUE, stringsAsFactors = FALSE)
# Check class
class(sd_df)

###### FRESH ######
#filter for groups and periods for the intervention period
sd_plot_veg <- sd_df %>%
  filter(Group_new %in% c("Constipation", "CTRL_AB") &
           period %in% c("Base", "VEG")) %>%
  select(Subject = participant_id, Group_new, period, Shannon, Faith_PD)

#make period a factor for plotting
sd_plot_veg$period <- factor(sd_plot_veg$period,
                         levels = c("Base", "VEG"),
                         labels = c("Base", "Veg"))
#rename groups in 'Group_new' and reorder group
sd_plot_veg$Group_new <- recode(sd_plot_veg$Group_new,
                            "CTRL_AB" = "Healthy",
                            "Constipation" = "Constipation")
sd_plot_veg$Group_new <- factor(sd_plot_veg$Group_new, levels = c("Healthy", "Constipation"))

#preview
head(sd_plot_veg)


#boxplot of Shannon diversity with Wilcox test annotations
Veg_Shannon <- ggplot(sd_plot_veg, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Veg" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE,
                     label.y = max(sd_plot_veg$Shannon) * 1.05,  
                     label.x = 1.5                            
  ) +
  theme_minimal() +
  labs(title = "Shannon Alpha Diversity Before and After Fresh Vegetable",
       x = "Period",
       y = "Shannon Index",
       fill = "Period") +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


#boxplot of Faith's PD with Wilcoxon test annotations
Veg_FaithPD <- ggplot(sd_plot_veg, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Veg" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE,
                     label.y = max(sd_plot_veg$Faith_PD) * 1.05,  
                     label.x = 1.5                            
  ) +
  theme_minimal() +
  labs(title = "Faith's Phylogenetic Diversity Before and After Fresh Vegetable",
       x = "Period",
       y = "Faith's PD",
       fill = "Period") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))


#saving
ggsave("A1_Alpha_Shannon_Fresh.png", plot = Veg_Shannon, width = 6, height = 4, dpi = 300)
ggsave("A1_Alpha_FaithPD_Fresh.png", plot = Veg_FaithPD, width = 6, height = 4, dpi = 300)


###### FERM ######
#filter for groups and periods for the intervention period
sd_plot_ferm <- sd_df %>%
  filter(Group_new %in% c("Constipation", "CTRL_AB") &
           period %in% c("WO1", "FERM")) %>%
  select(Subject = participant_id, Group_new, period, Shannon, Faith_PD)

#make period a factor for plotting
sd_plot_ferm$period <- factor(sd_plot_ferm$period,
                         levels = c("WO1", "FERM"),
                         labels = c("WO1", "Ferm"))
#rename groups in 'Group_new' and reorder group
sd_plot_ferm$Group_new <- recode(sd_plot_ferm$Group_new,
                            "CTRL_AB" = "Healthy",
                            "Constipation" = "Constipation")
sd_plot_ferm$Group_new <- factor(sd_plot_ferm$Group_new, levels = c("Healthy", "Constipation"))

#preview
head(sd_plot_ferm)


# Boxplot of Shannon diversity with Wilcox test annotations
Ferm_Shannon <- ggplot(sd_plot_ferm, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("WO1" = "#56B4E9", "Ferm" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE,
                     label.y = max(sd_plot_ferm$Shannon) * 1.05,  
                     label.x = 1.5                            
  ) +
  theme_minimal() +
  labs(title = "Shannon Alpha Diversity Before and After Fermented Vegetable",
       x = "Period",
       y = "Shannon Index",
       fill = "Period") +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Boxplot of Faith's PD with Wilcoxon test annotations
Ferm_FaithPD <- ggplot(sd_plot_ferm, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("WO1" = "#56B4E9", "Ferm" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE,
                     label.y = max(sd_plot_ferm$Faith_PD) * 1.05,  
                     label.x = 1.5                            
  ) +
  theme_minimal() +
  labs(title = "Faith's Phylogenetic Diversity Before and After Fermented Vegetable",
       x = "Period",
       y = "Faith's PD",
       fill = "Period") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))


#saving
ggsave("A1_Alpha_Shannon_Ferm.png", plot = Ferm_Shannon, width = 6, height = 4, dpi = 300)
ggsave("A1_Alpha_FaithPD_Ferm.png", plot = Ferm_FaithPD, width = 6, height = 4, dpi = 300)


#### Beta Diversity ####
#subset to VEG and FERM periods
ps_beta <- subset_samples(ferm_rarefac, period %in% c("VEG", "FERM"))
ps_beta <- prune_taxa(taxa_sums(ps_beta) > 0, ps_beta)  # remove zero-count OTUs

#prepare metadata
meta_df <- data.frame(sample_data(ps_beta), check.names = TRUE, stringsAsFactors = FALSE)

#make period a factor
meta_df$period <- factor(meta_df$period, levels = c("VEG", "FERM"))

#rename groups and make period a factor for plotting
meta_df$Group_new <- ifelse(meta_df$Group_new == "CTRL_AB", "Healthy", meta_df$Group_new)
meta_df$Group_new <- factor(meta_df$Group_new, levels = c("Healthy", "Constipation"))

##### Bray-Curtis Calculations #####
#extract OTU table as matrix and ensure taxa are in columns
otu_mat <- as(otu_table(ps_beta), "matrix")
#calculate bray curtis distance
if(taxa_are_rows(ps_beta)) otu_mat <- t(otu_mat)
bray_dist <- vegdist(otu_mat, method = "bray")

#PERMANOVA
adonis_bc <- adonis2(bray_dist ~ period, data = meta_df)
print(adonis_bc)

#perform classical multidimensional scaling on the Bray–Curtis distance matrix
#k = 2 extracts the first two principal coordinate axes 
pcoa_bc <- cmdscale(bray_dist, eig = TRUE, k = 2)

#calculate percent variance explained by each axis
eig_vals <- pcoa_bc$eig
var_exp <- round((eig_vals / sum(eig_vals)) * 100, 2)

#build a tidy data frame for plotting
pcoa_df <- data.frame(
  Sample = rownames(meta_df),
  Axis1 = pcoa_bc$points[,1],
  Axis2 = pcoa_bc$points[,2],
  period = meta_df$period,
  Group_new = meta_df$Group_new
)

##### Bray-Curtis PCoA Plot #####
bray_PCoA <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = period, shape = Group_new)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#F8766D", "FERM" = "#619CFF")) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances (VEG vs FERM)",
    x = paste0("PCoA 1 (", var_exp[1], "%)"),
    y = paste0("PCoA 2 (", var_exp[2], "%)"),
    color = "Period",
    shape = "Group"
  )


#plot PCoA (facet by group)
bray_PCoA_facet <- bray_PCoA +
  facet_wrap(~Group_new) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

#save plots
ggsave("A1_Beta_BC_FreshvsFerm.png", bray_PCoA, width = 6, height = 4, dpi = 300)
ggsave("A1_Beta_BC_FreshvsFerm_facet.png", bray_PCoA_facet, width = 6, height = 4, dpi = 300)


##### Weighted UniFrac Calculations #####
wunifrac_dist <- phyloseq::distance(ps_beta, method = "wunifrac")

#PERMANOVA
adonis_wu <- adonis2(wunifrac_dist ~ period, data = meta_df)
print(adonis_wu)

#perform classical multidimensional scaling on the Weighted UniFrac distance matrix
#k = 2 extracts the first two principal coordinate axes 
pcoa_wu <- cmdscale(as.matrix(wunifrac_dist), eig = TRUE, k = 2)

#build a tidy data frame for plotting
pcoa_df_wu <- data.frame(
  Sample = rownames(meta_df),
  Axis1 = pcoa_wu$points[,1],
  Axis2 = pcoa_wu$points[,2],
  period = meta_df$period,
  Group_new = meta_df$Group_new
)

##### Weighted UniFrac PCoA Plot #####
WUniFrac_PCoA <- ggplot(pcoa_df_wu, aes(x = Axis1, y = Axis2, color = period, shape = Group_new)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#F8766D", "FERM" = "#619CFF")) +
  theme_minimal() +
  labs(title = "PCoA of Weighted UniFrac Distances (VEG vs FERM)",
       x = paste0("PCoA 1 (", var_exp[1], "%)"),
       y = paste0("PCoA 2 (", var_exp[2], "%)"),
       color = "Period",
       shape = "Group")

#plot PCoA (facet by group)
WUniFrac_PCoA_facet <- WUniFrac_PCoA +
  facet_wrap(~Group_new) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # <-- adds facet border
  )


#save plots
ggsave("A1_Beta_WUF_FreshvsFerm.png", WUniFrac_PCoA, width = 6, height = 4, dpi = 300)
ggsave("A1_Beta_WUF_FreshvsFerm_facet.png", WUniFrac_PCoA_facet, width = 6, height = 4, dpi = 300)


#### DESeq2 ####
#add 1 to OTU counts to handle zeros
ps_plus1 <- transform_sample_counts(ferm_rarefac, function(x) x + 1)

#subset phyloseq object for Base vs VEG -- 'base' is the reference group
ps_veg <- subset_samples(ps_plus1, period %in% c("Base", "VEG"))
#convert to DESeq2 object and run DESeq
dds_VEG <- phyloseq_to_deseq2(ps_veg, ~ period)
dds_VEG <- DESeq(dds_VEG)
#extract results
res_VEG <- results(dds_VEG, contrast = c("period", "VEG", "Base"), tidy = TRUE)


#subset phyloseq object for WO1 vs FERM -- 'WO1' is the reference group
ps_ferm <- subset_samples(ps_plus1, period %in% c("WO1", "FERM"))
#convert to DESeq2 object and run DESeq
dds_FERM <- phyloseq_to_deseq2(ps_ferm, ~ period)
dds_FERM <- DESeq(dds_FERM)
#extract results
res_FERM <- results(dds_FERM, contrast = c("period", "FERM", "WO1"), tidy = TRUE)

#subset phyloseq object for VEG vs FERM -- 'VEG' is the reference group
ps_deseq <- subset_samples(ps_plus1, period %in% c("VEG", "FERM"))
#convert to DESeq2 object and run DESeq
dds_VEG_FERM <- phyloseq_to_deseq2(ps_deseq, ~ period)
dds_VEG_FERM <- DESeq(dds_VEG_FERM)
#extract results
res_VEG_FERM <- results(dds_VEG_FERM, contrast = c("period", "FERM", "VEG"), tidy = TRUE)


#function to create volcano plot from DESeq2 results
plot_volcano <- function(res_df, title = "Volcano Plot") {
  #add a significance column
  res_df <- res_df %>%
    mutate(Significant = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    ))
  #create volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9", "NS" = "grey70")) +
    theme_minimal() +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10 Adjusted P-value",
         color = "Significance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  
  return(p)
}

#volcano plots
volcano_VEG <- plot_volcano(res_VEG, title = "DESeq2: Base vs VEG")
volcano_FERM <- plot_volcano(res_FERM, title = "DESeq2: WO1 vs FERM")
volcano_VEG_FERM <- plot_volcano(res_VEG_FERM, title = "DESeq2: VEG vs FERM")

#preview plots
volcano_VEG
volcano_FERM
volcano_VEG_FERM

#save plots
ggsave("A1_DESeq2_Volcano_BasevsVEG.png", plot = volcano_VEG, width = 6, height = 4, dpi = 300)
ggsave("A1_DESeq2_Volcano_WO1vsFERM.png", plot = volcano_FERM, width = 6, height = 4, dpi = 300)
ggsave("A1_DESeq2_Volcano_VEGvsFERM.png", plot = volcano_VEG_FERM, width = 6, height = 4, dpi = 300)



#function to get significant taxa with enrichment direction
get_sig_table <- function(deseq_res, group1, group2, alpha = 0.05) {
  deseq_res %>%
    filter(!is.na(padj)) %>%               # remove NA padj
    filter(padj < alpha) %>%               # significance threshold
    mutate(
      Comparison = paste(group1, "vs", group2),
      Enriched_in = ifelse(log2FoldChange > 0, group1, group2)
    ) %>%
    select(Taxon = row, Comparison, log2FoldChange, padj, Enriched_in)
}

#apply to each comparison
sig_VEG <- get_sig_table(res_VEG, "VEG", "Base")
sig_FERM <- get_sig_table(res_FERM, "FERM", "WO1")
sig_VEG_FERM <- get_sig_table(res_VEG_FERM, "FERM", "VEG")

#combine all into one table and preview
sig_all <- bind_rows(sig_VEG, sig_FERM, sig_VEG_FERM)
head(sig_all)

#save DESeq2 significant taxa as CSV
write.csv(sig_all, "A1_DESeq2_significant_taxa_all_comparisons.csv", row.names = FALSE)

#merge DESeq2 analysis with taxonomy table
#extract taxonomy table
taxa_df <- as.data.frame(tax_table(ps_plus1))
taxa_df$ASV <- rownames(taxa_df)

#merge taxa into significant result table
sig_VEG_merged <- left_join(
  dplyr::rename(sig_VEG, ASV = Taxon),
  taxa_df,
  by = "ASV")

sig_FERM_merged <- left_join(
  dplyr::rename(sig_FERM, ASV = Taxon),
  taxa_df,
  by = "ASV")

sig_VEG_FERM_merged <- left_join(
  dplyr::rename(sig_VEG_FERM, ASV = Taxon),
  taxa_df,
  by = "ASV")

#combine all significant merged tables  
sig_all_merged <- bind_rows(sig_VEG_merged, sig_FERM_merged, sig_VEG_FERM_merged)

#export
write.csv(sig_all_merged, "A1_DESeq2_significant_taxa_with_taxonomy.csv", row.names = FALSE)

#get top 10 enriched/depleted taxa per group
top10_taxa <- sig_all_merged %>%
  group_by(Enriched_in) %>%
  arrange(desc(abs(log2FoldChange))) %>%  # use absolute fold change to capture strongest changes
  slice_head(n = 10) %>%
  ungroup()

#preview and save
top10_taxa
write.csv(top10_taxa, "A1_DESeq2_top10_taxa_per_group.csv", row.names = FALSE)

#list for top 10 enriched/depleted taxa per comparison and per group
top10_taxa_by_comparison <- sig_all_merged %>%
  group_by(Comparison, Enriched_in) %>%
  arrange(desc(abs(log2FoldChange))) %>%  #strongest fold changes first
  slice_head(n = 10) %>%
  ungroup()

#preview and save
top10_taxa_by_comparison
write.csv(top10_taxa_by_comparison, "A1_DESeq2_top10_taxa_per_comparison_group.csv", row.names = FALSE)

#make a label for taxa (use Genus if available, otherwise Family)
heatmap_data <- top10_taxa_by_comparison %>%
  mutate(
    #remove 'g__' and 'f__' prefixes
    Genus = sub("^g__", "", Genus),
    Family = sub("^f__", "", Family),
    TaxonLabel = ifelse(!is.na(Genus) & Genus != "", Genus, Family) #make TaxonLabel
  ) %>%
  select(Comparison, TaxonLabel, log2FoldChange) %>%
  mutate(TaxonLabel = fct_reorder(TaxonLabel, log2FoldChange))

#heatmap plot of top 10 taxa
heatmap_data <- heatmap_data %>%
  mutate(
    Comparison = factor(
      recode(Comparison,
             "VEG vs Base" = "Base vs VEG",
             "FERM vs WO1" = "WO1 vs FERM",
             "FERM vs VEG" = "VEG vs FERM"),
      levels = c("Base vs VEG", "WO1 vs FERM", "VEG vs FERM")
    )
  )

Top10_heatmap <- ggplot(heatmap_data, aes(x = Comparison, y = TaxonLabel, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    name = "log2FC"
  ) +
  labs(
    x = "Comparison",
    y = "Taxon (Genus/Family)",
    title = "Top 10 Differentially Abundant Taxa Across Comparisons"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("A1_DESeq2_Heatmap_Top10Taxa.png", Top10_heatmap, width = 15, height = 10, dpi = 300)



###Group-stratified DESeq2 analyis

#Subset phyloseq object by group
#use +1 data
#CTRL_AB
ps_ctrl_ab <- subset_samples(ps_plus1, Group_new == "CTRL_AB")
ps_ctrl_ab <- prune_taxa(taxa_sums(ps_ctrl_ab) > 0, ps_ctrl_ab)
#Constipation
ps_constipation <- subset_samples(ps_plus1, Group_new == "Constipation")
ps_constipation <- prune_taxa(taxa_sums(ps_constipation) > 0, ps_constipation)


#run DESeq2 for HEALTHY group
#subset: Base vs VEG
ps_ctrl_ab_veg <- subset_samples(ps_ctrl_ab, period %in% c("Base", "VEG"))

dds_ctrl_ab_VEG <- phyloseq_to_deseq2(ps_ctrl_ab_veg, ~ period)
dds_ctrl_ab_VEG <- DESeq(dds_ctrl_ab_VEG)

res_ctrl_ab_VEG <- results(dds_ctrl_ab_VEG, contrast = c("period", "VEG", "Base"), tidy = TRUE)


#subset: WO1 vs FERM
ps_ctrl_ab_ferm <- subset_samples(ps_ctrl_ab, period %in% c("WO1", "FERM"))

dds_ctrl_ab_FERM <- phyloseq_to_deseq2(ps_ctrl_ab_ferm, ~ period)
dds_ctrl_ab_FERM <- DESeq(dds_ctrl_ab_FERM)

res_ctrl_ab_FERM <- results(dds_ctrl_ab_FERM, contrast = c("period", "FERM", "WO1"), tidy = TRUE)

# subset: VEG vs FERM
ps_ctrl_ab_veg_ferm <- subset_samples(ps_ctrl_ab, period %in% c("VEG", "FERM"))

dds_ctrl_ab_VEG_FERM <- phyloseq_to_deseq2(ps_ctrl_ab_veg_ferm, ~ period)
dds_ctrl_ab_VEG_FERM <- DESeq(dds_ctrl_ab_VEG_FERM)

res_ctrl_ab_VEG_FERM <- results(dds_ctrl_ab_VEG_FERM, contrast = c("period", "FERM", "VEG"), tidy = TRUE)


#run DESEQ2 for constipation
#subset: Base vs VEG
ps_const_veg <- subset_samples(ps_constipation, period %in% c("Base", "VEG"))

dds_const_VEG <- phyloseq_to_deseq2(ps_const_veg, ~ period)
dds_const_VEG <- DESeq(dds_const_VEG)

res_const_VEG <- results(dds_const_VEG, contrast = c("period", "VEG", "Base"), tidy = TRUE)

#subset: WO1 vs FERM
ps_const_ferm <- subset_samples(ps_constipation, period %in% c("WO1", "FERM"))

dds_const_FERM <- phyloseq_to_deseq2(ps_const_ferm, ~ period)
dds_const_FERM <- DESeq(dds_const_FERM)

res_const_FERM <- results(dds_const_FERM, contrast = c("period", "FERM", "WO1"), tidy = TRUE)

#subset: VEG vs FERM
ps_const_veg_ferm <- subset_samples(ps_constipation, period %in% c("VEG", "FERM"))

dds_const_VEG_FERM <- phyloseq_to_deseq2(ps_const_veg_ferm, ~ period)
dds_const_VEG_FERM <- DESeq(dds_const_VEG_FERM)

res_const_VEG_FERM <- results(dds_const_VEG_FERM, contrast = c("period", "FERM", "VEG"), tidy = TRUE)

#plot
v1 <- plot_volcano(res_ctrl_ab_VEG,  "Healthy: VEG vs Base")
v2 <- plot_volcano(res_ctrl_ab_FERM,  "Healthy: FERM vs WO1")
v3 <- plot_volcano(res_ctrl_ab_VEG_FERM,  "Healthy: FERM vs VEG")

v4 <- plot_volcano(res_const_VEG, "Constipation: VEG vs Base")
v5 <- plot_volcano(res_const_FERM, "Constipation: FERM vs WO1")
v6 <- plot_volcano(res_const_VEG_FERM, "Constipation: FERM vs VEG")

#combine into 6 panel
combined_volcano <- (v1 | v2 | v3) /
  (v4 | v5 | v6)

combined_volcano

ggsave("A1_DESeq2_STRAT_volcano_combined.png", combined_volcano, width = 15, height = 10, dpi = 300)


# Helper function to clean plots
clean_volcano <- function(p) {
  p + 
    theme(
      panel.grid = element_blank(),        # remove grid
      legend.position = "none",            # remove legend
      panel.border = element_rect(color = "black", fill = NA, size = 1)  # add border
    )
}

# Apply to all plots
volcano_ctrl_base_veg  <- clean_volcano(v1)
volcano_ctrl_wo1_ferm  <- clean_volcano(v2)
volcano_ctrl_veg_ferm  <- clean_volcano(v3)

volcano_const_base_veg <- clean_volcano(v4)
volcano_const_wo1_ferm <- clean_volcano(v5)
volcano_const_veg_ferm <- clean_volcano(v5)

# Combine plots with patchwork
combined_volcano <- (volcano_ctrl_base_veg | volcano_ctrl_wo1_ferm | volcano_ctrl_veg_ferm) /
  (volcano_const_base_veg | volcano_const_wo1_ferm | volcano_const_veg_ferm)



#function to get significant taxa with direction
get_sig_taxa <- function(res_df, group1, group2, alpha = 0.05) {
  res_df %>%
    filter(!is.na(padj)) %>%               # remove NA padj
    filter(padj < alpha) %>%               # significant taxa
    mutate(
      Comparison = paste(group1, "vs", group2),
      Enriched_in = ifelse(log2FoldChange > 0, group1, group2)
    ) %>%
    select(Taxon = row, Comparison, log2FoldChange, padj, Enriched_in)
}

#healthy group
sig_ctrl_veg  <- get_sig_taxa(res_ctrl_ab_VEG,  "VEG", "Base")
sig_ctrl_ferm  <- get_sig_taxa(res_ctrl_ab_FERM,  "FERM", "WO1")
sig_ctrl_veg_ferm  <- get_sig_taxa(res_ctrl_ab_VEG_FERM,  "FERM", "VEG")

#constipation group
sig_const_veg <- get_sig_taxa(res_const_VEG, "VEG", "Base")
sig_const_ferm <- get_sig_taxa(res_const_FERM, "FERM", "WO1")
sig_const_veg_ferm <- get_sig_taxa(res_const_VEG_FERM, "FERM", "VEG")

#combine all into one table
sig_all <- bind_rows(
  sig_ctrl_base_veg, sig_ctrl_wo1_ferm, sig_ctrl_veg_ferm,
  sig_const_base_veg, sig_const_wo1_ferm, sig_const_veg_ferm
)

#extract taxonomy table as a data.frame
taxa_df <- as.data.frame(tax_table(ps_plus1))
taxa_df$ASV <- rownames(taxa_df)  # add ASV/OTU IDs as a column

#merge taxonomy with significant taxa
sig_all_merged <- sig_all %>%
  left_join(dplyr::rename(taxa_df, Taxon = ASV), by = "Taxon")

#create a clean taxon label using Genus if available, otherwise Family
sig_all_merged <- sig_all_merged %>%
  mutate(
    Genus = sub("^g__", "", Genus),
    Family = sub("^f__", "", Family),
    TaxonLabel = ifelse(!is.na(Genus) & Genus != "", Genus, Family)
  )

#preview
head(sig_all_merged)

#save as CSV
write.csv(sig_all_merged, "A1_DESeq2_STRAT_significant_taxa_with_taxonomy.csv", row.names = FALSE)









#### Manuscript Plots ####
##### Alpha Diversity #####
#custom colors
my_colors <- c(
  "Base" = "#4E79A7",
  "Veg" = "#B699C6",
  "WO1" = "#F28E2B",
  "Ferm" = "#59A14F",
  "WO2" = "#E15759"
)

# Adjust Shannon and Faith PD plots
Veg_Shannon <- ggplot(sd_plot_veg, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(title = NULL, x = "Period", y = "Shannon Index") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

Veg_FaithPD <- ggplot(sd_plot_veg, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(title = NULL, x = "Period", y = "Faith's PD") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

Ferm_Shannon <- ggplot(sd_plot_ferm, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(title = NULL, x = "Period", y = "Shannon Index") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

Ferm_FaithPD <- ggplot(sd_plot_ferm, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(title = NULL, x = "Period", y = "Faith's PD") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

#combine alpha plots
Alpha_Shannon_FaithPD_Combined <- 
  (Veg_Shannon + Veg_FaithPD) /
  (Ferm_Shannon + Ferm_FaithPD) +
  plot_layout(widths = c(1, 2)) +  # slightly wider right column
  plot_annotation(
    title = "Alpha Diversity (Fresh vs Fermented Vegetable Interventions)"
  )
Alpha_Shannon_FaithPD_Combined

#save
ggsave("A1_Figure_Alpha_Shannon_FaithPD_Combined.png", Alpha_Shannon_FaithPD_Combined, width = 15, height = 10, dpi = 300)


##### Beta Diversity #####
#adjust BC and WUF PCoA plots
#BC PCoA
bray_PCoA_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = period, shape = Group_new)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#F8766D", "FERM" = "#619CFF")) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances (VEG vs FERM)",
    x = paste0("PCoA 1 (", var_exp[1], "%)"),
    y = paste0("PCoA 2 (", var_exp[2], "%)"),
    color = "Period",
    shape = "Group"
  )+
  theme(
    panel.grid = element_blank(),        
    panel.border = element_rect(      
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    legend.position = "none"
  )

#WUF PCoA
WUniFrac_PCoA_plot <- ggplot(pcoa_df_wu, aes(x = Axis1, y = Axis2, color = period, shape = Group_new)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#F8766D", "FERM" = "#619CFF")) +
  theme_minimal() +
  labs(title = "PCoA of Weighted UniFrac Distances (VEG vs FERM)",
       x = paste0("PCoA 1 (", var_exp[1], "%)"),
       y = paste0("PCoA 2 (", var_exp[2], "%)"),
       color = "Period",
       shape = "Group") +
  theme(
    panel.grid = element_blank(),        
    panel.border = element_rect(      
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

#combine beta plots
Beta_BC_WUF_PCoA_Combined <- bray_PCoA_plot + plot_spacer() + WUniFrac_PCoA_plot +
  plot_layout(ncol = 3, widths = c(1, 0.1, 1))
Beta_BC_WUF_PCoA_Combined

#save
ggsave("A1_Figure_Beta_BC_WUF_PCoA_Combined.png", Beta_BC_WUF_PCoA_Combined, width = 15, height = 8, dpi = 300)

##### DESeq2 #####
#adjust DESeq2 plots
#remove legend for plotting
volcano_VEG_plot <- volcano_VEG + theme(legend.position = "none")
volcano_FERM_plot <- volcano_FERM + theme(legend.position = "none")

#combine volcano plots
combined_volcano <- volcano_VEG_plot + plot_spacer() + volcano_FERM_plot + plot_spacer() + volcano_VEG_FERM +
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1)) &
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
combined_volcano

#save
ggsave("A1_Figure_DESeq2_Volcano.png", combined_volcano, width = 15, height = 8, dpi = 300)

#DESeq2 volcano plot to be combined with the Top10_heatmap plot manually (could not format scale of both plots to fit as one figure )
Top10_heatmap
