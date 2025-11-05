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
    !Domain %in% c("d__Eukaryota", "d__Archaea") # remove Eukaryota and Archaea
)


#subset to only female samples and remove NA samples in Group samples
ps_filter <- subset_samples (ps_taxa_filter, 
                             host_sex == "female" ) &
                            Group != "not applicable")

#Create a new column to group 'Control' and 'Antibiotics' group as one group
sample_data(ps_filter)$Group_new <- ifelse(
  sample_data(ps_filter)$Group %in% c("Control", "Antibiotics"),
  "CTRL_AB",
  sample_data(ps_filter)$Group
)

table(sample_data(ps_filter)$Group_new)

###viewing OTU table after filtering
otu_df <- as.data.frame(otu_table(ps_filter))
head(otu_df) 

# Rarefy samples
# set rngseed the same number each time to reproduce this exact analysis
# t transposes the table to use rarecurve function
# cex decreases font size
#rarefaction curve using final filtered table
rarecurve(t(as.data.frame(otu_table(ps_filter))), cex=0.1)
ferm_rare <- rarefy_even_depth(ps_filter, rngseed = 1, sample.size = 35798)



##### Saving #####
save(ps_filter, file="ferm_final.RData")
save(ferm_rare, file="ferm_rare.RData")