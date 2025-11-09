#Alpha-Beta Diversity Plots

#!/usr/bin/env Rscript
# NOTE: to install phyloseq, please use the following code instead of the usual "install.packages" function:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#load library
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library (picante)
library(dplyr)
library(tidyr)
library(ggplot2)
library (ggpubr)

#### Load data ####
# Change file paths as necessary, below we are importing them as tibble
metafp <- "Export/fermentation_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "Export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "Export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "Export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix 
##take the otu data and index it to remove the first column
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
fermentation_ps_raw <- phyloseq(OTU, SAMP, TAX, phylotree)

# View components of phyloseq object with the following commands
otu_table(fermentation_ps_raw)
sample_data(fermentation_ps_raw)
tax_table(fermentation_ps_raw)
phy_tree(fermentation_ps_raw)



#### Filter phyloseq object ####
# Remove unwanted taxa: mitochondria, chloroplasts, Eukaryota, and Archaea
ps_taxa_filter <- subset_taxa(
  fermentation_ps_raw, 
  !Order %in% c("o__Chloroplast") &   # remove chloroplasts
    !Family %in% c("f__Mitochondria") & # remove mitochondria
    !Domain %in% c("d__Eukaryota", "d__Archaea","Unassigned") # remove Eukaryota and Archaea
)


#subset to only female samples and remove NA samples in Group samples
ps_metadata_filter <- subset_samples (ps_taxa_filter, 
                             host_sex == "female" &
                            Group != "not applicable")

#Create a new column to group 'Control' and 'Antibiotics' group as one group
sample_data(ps_metadata_filter)$Group_new <- ifelse(
  sample_data(ps_metadata_filter)$Group %in% c("Control", "Antibiotics"),
  "CTRL_AB",
  sample_data(ps_metadata_filter)$Group
)

ps_filter <- ps_metadata_filter


#feature-based filtering
#ps_abund_filter <- filter_taxa(
#  ps_filter,
#  function(x) sum(x) > (0.00005 * sum(sample_sums(ps_filter))),
#  prune = TRUE
#)

#number of taxa removed after feature-based filtering
#ntaxa(ps_filter) - ntaxa(ps_abund_filter)

#number of reads removed after feature-based filtering
#sum(otu_table(ps_filter)) - sum(otu_table(ps_abund_filter))

#number of taxa remaining
#ntaxa(ps_abund_filter)
#not using feature-based filtering, it removed a lot of samples


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



#look at the top 5 OTUs before and after filtering
# Define a helper function to get top 5 OTUs from any phyloseq object
get_top_otus <- function(ps_obj, n = 5) {
  otu_df <- as.data.frame(as.table(as.matrix(otu_table(ps_obj))))
  colnames(otu_df) <- c("OTU", "Sample", "Count")
  otu_df[order(-otu_df$Count), ][1:n, ]
}

# Use it before and after filtering
top5_raw <- get_top_otus(fermentation_ps_raw)
top5_filtered <- get_top_otus(ps_filter)

# Print results
top5_raw
top5_filtered



#### Rarefaction ####
##Summarize sequencing depth per group
# Convert sample_data to a data frame and add total reads per sample
sample_df <- as.data.frame(sample_data(ps_filter))
sample_df$TotalReads <- sample_sums(ps_filter)

# Summarize min and max sequencing depth per Group_new to determine rarefaction depth
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


# Rarefaction
# set rngseed the same number each time to reproduce this exact analysis
ferm_rarefac <- rarefy_even_depth(ps_clean,
                               rngseed = 1,
                               sample.size = 31218,
                               replace = FALSE,
                               trimOTUs = TRUE,
                               verbose = TRUE)


#check rarefied object
# Check total reads per sample — all should equal 31218
table(sample_sums(ferm_rarefac))

# Check how many OTUs remain
ntaxa(ps_filter) #after filtering
ntaxa(ps_clean) #after cleaning (before rarefaction)
ntaxa(ferm_rarefac) #after rarefaction

# Check how many samples remain
nsamples(ps_filter) #after filtering
nsamples(ps_clean) #after cleaning (before rarefaction)
nsamples(ferm_rarefac) #after rarefaction

# Confirm your grouping variable is still intact
table(sample_data(ferm_rarefac)$Group_new)


##### Saving #####
save(ps_filter, file="ferm_final.RData")
save(ferm_rarefac, file="ferm_rare.RData")


#### Alpha Diversity ####
#Shannon
alpha_shannon <- estimate_richness(ferm_rarefac, measures = "Shannon")

#Faith's PD
otu_mat <- as(otu_table(ferm_rarefac), "matrix")
if (taxa_are_rows(ferm_rarefac)) {
  otu_mat <- t(otu_mat)
}

alpha_faith_pd <- pd(otu_mat, phy_tree(ferm_rarefac), include.root = TRUE)


# Add Shannon and Faith PD metrics to metadata as new columns
sd <- as.data.frame(sample_data(ferm_rarefac))
sd$Shannon <- alpha_shannon$Shannon[match(rownames(sd), rownames(alpha_shannon))]
sd$Faith_PD <- alpha_faith_pd$PD[match(rownames(sd), rownames(alpha_faith_pd))]

# Put the updated sample_data back into the phyloseq object
sample_data(ferm_rarefac) <- sample_data(sd)

# Verify
head(sample_data(ferm_rarefac))
sample_variables(ferm_rarefac)


#Alpha Diversity Plots
# extract and convert the sample data from phyloseq to a true data.frame
sd_df <- data.frame(sample_data(ferm_rarefac), check.names = TRUE, stringsAsFactors = FALSE)
# Check class
class(sd_df)

#Alpha plots for FRESH vegetables
# Filter only the groups and periods you care about
sd_plot_veg <- sd_df %>%
  filter(Group_new %in% c("Constipation", "CTRL_AB") &
           period %in% c("Base", "VEG")) %>%
  select(Subject = participant_id, Group_new, period, Shannon, Faith_PD)

# Make period a factor for plotting
sd_plot_veg$period <- factor(sd_plot_veg$period,
                         levels = c("Base", "VEG"),
                         labels = c("Base", "Veg"))
# Recode Group_new and reorder group
sd_plot_veg$Group_new <- recode(sd_plot_veg$Group_new,
                            "CTRL_AB" = "Healthy",
                            "Constipation" = "Constipation")
sd_plot_veg$Group_new <- factor(sd_plot_veg$Group_new, levels = c("Healthy", "Constipation"))

# Preview
head(sd_plot_veg)


# Boxplot of Shannon diversity with Wilcox test annotations
Veg_Shannon <- ggplot(sd_plot_veg, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Veg" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE) +
  theme_minimal() +
  labs(title = "Shannon Alpha Diversity Before and After Fresh Vegetable",
       x = "Period",
       y = "Shannon Index",
       fill = "Period") +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

# Boxplot of Faith's PD with Wilcoxon test annotations
Veg_FaithPD <- ggplot(sd_plot_veg, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Veg" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE) +
  theme_minimal() +
  labs(title = "Faith's Phylogenetic Diversity Before and After Fresh Vegetable",
       x = "Period",
       y = "Faith's PD",
       fill = "Period") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"))


#saving
ggsave("Fresh_Veg_Shannon.png", plot = Veg_Shannon, width = 6, height = 4, dpi = 300)
ggsave("Fresh_Veg_FaithPD.png", plot = Veg_FaithPD, width = 6, height = 4, dpi = 300)


#Alpha plots for FERM vegetables
# Filter only the groups and periods you care about
sd_plot_ferm <- sd_df %>%
  filter(Group_new %in% c("Constipation", "CTRL_AB") &
           period %in% c("WO1", "FERM")) %>%
  select(Subject = participant_id, Group_new, period, Shannon, Faith_PD)

# Make period a factor for plotting
sd_plot_ferm$period <- factor(sd_plot_ferm$period,
                         levels = c("WO1", "FERM"),
                         labels = c("Base", "Ferm"))
# Recode Group_new and reorder group
sd_plot_ferm$Group_new <- recode(sd_plot_ferm$Group_new,
                            "CTRL_AB" = "Healthy",
                            "Constipation" = "Constipation")
sd_plot_ferm$Group_new <- factor(sd_plot_ferm$Group_new, levels = c("Healthy", "Constipation"))

# Preview
head(sd_plot_ferm)


# Boxplot of Shannon diversity with Wilcox test annotations
Ferm_Shannon <- ggplot(sd_plot_ferm, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Ferm" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE) +
  theme_minimal() +
  labs(title = "Shannon Alpha Diversity Before and After Fermented Vegetable",
       x = "Period",
       y = "Shannon Index",
       fill = "Period") +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

# Boxplot of Faith's PD with Wilcoxon test annotations
Ferm_FaithPD <- ggplot(sd_plot_ferm, aes(x = period, y = Faith_PD, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Group_new) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "Ferm" = "#E69F00")) +
  stat_compare_means(aes(group = period),
                     method = "wilcox.test",
                     label = "p.format",
                     paired = FALSE) +
  theme_minimal() +
  labs(title = "Faith's Phylogenetic Diversity Before and After Fermented Vegetable",
       x = "Period",
       y = "Faith's PD",
       fill = "Period") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"))


#saving
ggsave("Ferm_Veg_Shannon.png", plot = Ferm_Shannon, width = 6, height = 4, dpi = 300)
ggsave("Ferm_Veg_FaithPD.png", plot = Ferm_FaithPD, width = 6, height = 4, dpi = 300)


#### Beta Diversity ####
# Subset phyloseq object to periods VEG and FERM
ps_veg_ferm <- subset_samples(ferm_rarefac, period %in% c("VEG", "FERM"))

# Remove OTUs with zero counts after subsetting
ps_veg_ferm <- prune_taxa(taxa_sums(ps_veg_ferm) > 0, ps_veg_ferm)

#Bray-Curtis
# Extract OTU table and make sure taxa are in columns
otu_mat <- as(otu_table(ps_veg_ferm), "matrix")
if(taxa_are_rows(ps_veg_ferm)) {
  otu_mat <- t(otu_mat)
}
# Bray-Curtis distance
bray_dist <- vegdist(otu_mat, method = "bray")

#covert sample data to data frame for PERMAOVA
meta_df <- data.frame(sample_data(ps_veg_ferm), check.names = TRUE, stringsAsFactors = FALSE)
class(meta_df)
#Make period a factor for PERMAOVA
meta_df$period <- factor(meta_df$period, levels = c("VEG", "FERM"))

# PERMANOVA
adonis_bc <- adonis2(bray_dist ~ period, data = meta_df)
adonis_bc

# Perform PCoA (classical multidimensional scaling)
pcoa_bc <- cmdscale(bray_dist, eig = TRUE, k = 2)  # 2 dimensions

# Create a data frame for plotting
pcoa_df <- data.frame(
  Sample = rownames(meta_df),
  Axis1 = pcoa_bc$points[,1],
  Axis2 = pcoa_bc$points[,2],
  period = meta_df$period,
  Group_new = meta_df$Group_new
)

# Make period a factor for consistent colors
pcoa_df$period <- factor(pcoa_df$period, levels = c("VEG", "FERM"))

# Rename groups so CTRL_AB -> Healthy
pcoa_df$Group_new <- ifelse(pcoa_df$Group_new == "CTRL_AB", "Healthy", pcoa_df$Group_new)

# Make Group_new a factor to control facet order: Healthy left, Constipation right
pcoa_df$Group_new <- factor(pcoa_df$Group_new, levels = c("Healthy", "Constipation"))


# Plot PCoA
bray_PCoA <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = period, shape = Group_new)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#16A7A1", "FERM" = "#A09BC2")) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances (VEG vs FERM)",
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Period",
    shape = "Group"
  )

# Faceted PCoA plot
bray_PCoA_facet <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = period)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = period), type = "norm", linetype = 2, alpha = 0.3) +
  scale_color_manual(values = c("VEG" = "#16A7A1", "FERM" = "#A09BC2")) +
  facet_wrap(~Group_new) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances by Group",
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Period"
  ) +
  theme(strip.text = element_text(face = "bold"))


ggsave("FreshvsFerm_Bray_PCoA.png", plot = bray_PCoA, width = 6, height = 4, dpi = 300)
ggsave("FreshvsFerm_Bray_PCoA_facet.png", plot = bray_PCoA_facet, width = 6, height = 4, dpi = 300)


#Weighted UniFrac distance
# ---------------------------
# Subset phyloseq object to periods VEG and FERM
# ---------------------------
ps_wunifrac <- subset_samples(ferm_rarefac, period %in% c("VEG", "FERM"))

# Remove OTUs with zero counts
ps_wunifrac <- prune_taxa(taxa_sums(ps_wunifrac) > 0, ps_wunifrac)

# ---------------------------
# Weighted UniFrac distance
# ---------------------------
wunifrac_dist <- phyloseq::distance(ps_wunifrac, method = "wunifrac")

# ---------------------------
# PERMANOVA
# ---------------------------
meta_df <- data.frame(sample_data(ps_wunifrac), check.names = TRUE, stringsAsFactors = FALSE)
meta_df$period <- factor(meta_df$period, levels = c("VEG", "FERM"))

adonis_wu <- adonis2(wunifrac_dist ~ period, data = meta_df)
adonis_wu






# extract and convert the sample data from phyloseq to a true data.frame
sd_df <- data.frame(sample_data(ferm_rarefac), check.names = TRUE, stringsAsFactors = FALSE)
# Check class
class(sd_df)

#Beta Diversity 
# Filter only the groups and periods you care about
sd_plot_veg <- sd_df %>%
  filter(Group_new %in% c("Constipation", "CTRL_AB") &
           period %in% c("Base", "VEG")) %>%
  select(Subject = participant_id, Group_new, period, Shannon, Faith_PD)





