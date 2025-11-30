#install lme4 packages (for Shirley's older R version)
install.packages("lme4", type = "source")
install.packages("lmerTest")

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
library(ggsignif)
library(ggpubr)
library(dplyr)
library(patchwork)
library(tidyr)
library(tibble)



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

##### Linear Mixed Model (LMM) #####

##fresh shannon LMM (period as the fixed effects and participant_id as the random effect)
model_shannon_fresh <- lmer(Shannon ~ period + (1 | participant_id), data = meta_fresh)

#perform ANOVA to test 
anova(model_shannon_fresh)

#post-hoc tests (period comparisons to test which periods differ from each other)
emmeans(model_shannon_fresh, pairwise ~ period)

#extract the estimated marginal means from each period
emm_fresh_shannon <- emmeans(model_shannon_fresh, ~ period)
#convert estimated marginal means output into a data frame for plotting
emm_df_fresh <- as.data.frame(emm_fresh_shannon)

#plot the estimated means 
LMM_fresh_veg_shannon <- ggplot(emm_df_fresh, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("LMM_Fresh_Veg_Shannon.png", plot = LMM_fresh_veg_shannon, width = 6, height = 4, dpi = 300)


# Pairwise contrasts from the LMM
contrasts_fresh <- emmeans(model_shannon_fresh, pairwise ~ period)$contrasts %>%
  as.data.frame()

# Correctly prepare table for plotting significance
pval_table <- contrasts_fresh %>%
  mutate(
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".*- ", "", contrast),
    y.position = sapply(group1, function(g) {
      max(emm_df_fresh$emmean[emm_df_fresh$period == g]) + 0.25
    }),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*"
    )
  ) %>%
  filter(p.value < 0.05) 



LMM_fresh_veg_shannon_sig <- ggplot(emm_df_fresh, aes(x = period, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.15
  ) +
  
  # Increase y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))  # add 5% below, 20% above



LMM_fresh_veg_shannon_box <- ggplot(emm_df_fresh, aes(x = period, y = emmean)) +
  
  # box-like rectangles for CI, filled by period
  geom_rect(
    aes(
      xmin = as.numeric(period) - 0.3,
      xmax = as.numeric(period) + 0.3,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = period
    ),
    alpha = 0.5, color = "black",
    inherit.aes = FALSE
  ) +
  
  # EMM points in the middle
  geom_point(size = 3, color = "black") +
  
  # Add significance labels
  stat_pvalue_manual(
    pval_table,
    label = "label",
    tip.length = 0.03,
    step.increase = 0.1
  ) +
  
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Shannon") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c("Base" = "#56B4E9", "VEG" = "#E69F00", "WO1" = "#009E73"))


















##Fresh Faith LMM
model_faith_fresh <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_fresh)

#ANOVA
anova(model_faith_fresh)

#post-hoc contrasts
emmeans(model_faith_fresh, pairwise ~ period)

#extract the estimated marginal means from each period
emm_fresh_faith <- emmeans(model_faith_fresh, ~ period)
#convert estimated marginal means output into a data frame for plotting
emm_df_faith_fresh <- as.data.frame(emm_fresh_faith)


#plot
LMM_fresh_veg_faith <- ggplot(emm_df_faith_fresh, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fresh Intervention – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#save plot
ggsave("LMM_Fresh_Veg_Faith.png", plot = LMM_fresh_veg_faith, width = 6, height = 4, dpi = 300)




##Ferm Shannon LMM
model_shannon_ferm <- lmer(Shannon ~ period + (1 | participant_id), data = meta_ferm)
#ANOVA and post-hoc contrasts
anova(model_shannon_ferm)
emmeans(model_shannon_ferm, pairwise ~ period)

#extract the estimated marginal means from each period and convert output into a data frame for plotting
emm_shannon_ferm <- emmeans(model_shannon_ferm, ~ period)
emm_df_shannon_ferm <- as.data.frame(emm_shannon_ferm)

#plot
LMM_ferm_veg_shannon <- ggplot(emm_df_shannon_ferm, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#save plot
ggsave("LMM_Ferm_Veg_Shannon.png", plot = LMM_ferm_veg_shannon, width = 6, height = 4, dpi = 300)




##Ferm Faith LMM
model_faith_ferm <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_ferm)
#ANOVA and post-hoc contrasts
anova(model_faith_ferm)
emmeans(model_faith_ferm, pairwise ~ period)

#extract the estimated marginal means from each period and convert output into a data frame for plotting
emm_faith_ferm <- emmeans(model_faith_ferm, ~ period)
emm_df_faith_ferm <- as.data.frame(emm_faith_ferm)

#plot
LMM_ferm_veg_faith <- ggplot(emm_df_faith_ferm, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Fermentation Intervention – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#save plot
ggsave("LMM_Ferm_Veg_Faith.png", plot = LMM_ferm_veg_faith, width = 6, height = 4, dpi = 300)


##Microbiome Shannon LMM
model_shannon_microb <- lmer(Shannon ~ period + (1 | participant_id), data = meta_microb)
#ANOVA and post-hoc contrasts
anova(model_shannon_microb)
emmeans(model_shannon_microb, pairwise ~ period)

#extract the estimated marginal means from each period and convert output into a data frame for plotting
emm_shannon_microb <- emmeans(model_shannon_microb, ~ period)
emm_df_shannon_microb <- as.data.frame(emm_shannon_microb)
#plot
LMM_microbiome_shannon <- ggplot(emm_df_shannon_microb, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Shannon Diversity") +
  xlab("Period") +
  ggtitle("Microbiome – Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#save plot
ggsave("LMM_Microbiome_Shannon.png", plot = LMM_ferm_veg_faith, width = 6, height = 4, dpi = 300)



##Microbiome Faith
model_faith_microb <- lmer(Faith_PD ~ period + (1 | participant_id), data = meta_microb)
#ANOVA and post-hoc contrasts
anova(model_faith_microb)
emmeans(model_faith_microb, pairwise ~ period)

#extract the estimated marginal means from each period and convert output into a data frame for plotting
emm_faith_microb <- emmeans(model_faith_microb, ~ period)
emm_df_faith_microb <- as.data.frame(emm_faith_microb)

#plot
LMM_microbiome_faith <- ggplot(emm_df_faith_microb, aes(x=period, y=emmean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  ylab("Faith's PD") +
  xlab("Period") +
  ggtitle("Microbiome – Faith's PD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#save plot
ggsave("LMM_Microbiome_Faith.png", plot = LMM_microbiome_faith, width = 6, height = 4, dpi = 300)


###combine graphs
#add the subset info and combine

#fresh
emm_df_fresh$metric <- "Shannon"
emm_df_fresh$subset <- "Fresh"

emm_df_faith_fresh$metric <- "Faith_PD"
emm_df_faith_fresh$subset <- "Fresh"

#ferm
emm_df_shannon_ferm$metric <- "Shannon"
emm_df_shannon_ferm$subset <- "Fermented"

emm_df_faith_ferm$metric <- "Faith_PD"
emm_df_faith_ferm$subset <- "Fermented"

#microbiome washouts
emm_df_shannon_microb$metric <- "Shannon"
emm_df_shannon_microb$subset <- "Microbiome Washout"

emm_df_faith_microb$metric <- "Faith_PD"
emm_df_faith_microb$subset <- "Microbiome Washout"

#combine all data
plot_data <- bind_rows(
  emm_df_fresh, emm_df_faith_fresh,
  emm_df_shannon_ferm, emm_df_faith_ferm,
  emm_df_shannon_microb, emm_df_faith_microb
)

#plot all
combined_plot <- ggplot(plot_data, aes(x=period, y=emmean, color=period)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.2) +
  facet_grid(metric ~ subset, scales = "free_y") +
  ylab("Alpha Diversity") +
  xlab("Period") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, size=1), # add border
    strip.background = element_rect(fill="grey90", color="black", size=0.5) # optional: facet label background
  )

ggsave("LMM_Alpha_combined.png", plot = combined_plot, width = 6, height = 4, dpi = 300)


#### Beta Diversity ####
#prepare metadata
meta_fresh_df <- as.data.frame(sample_data(ps_fresh))
meta_fresh_df$SampleID <- rownames(meta_fresh_df)

meta_ferm_df <- as.data.frame(sample_data(ps_ferm))
meta_ferm_df$SampleID <- rownames(meta_ferm_df)

meta_microb_df <- as.data.frame(sample_data(ps_microb))
meta_microb_df$SampleID <- rownames(meta_microb_df)


#beta diversity matrices calculations (bray and weighted unifrac)
#fresh
beta_bc_fresh <- phyloseq::distance(ps_fresh, method = "bray")
beta_wu_fresh <- phyloseq::distance(ps_fresh, method = "wunifrac")

#fermented
beta_bc_ferm <- phyloseq::distance(ps_ferm, method = "bray")
beta_wu_ferm <- phyloseq::distance(ps_ferm, method = "wunifrac")

#microbiome washouts
beta_bc_microb <- phyloseq::distance(ps_microb, method = "bray")
beta_wu_microb <- phyloseq::distance(ps_microb, method = "wunifrac")


#distance-to-baseline calculation function
distance_to_baseline <- function(dist_matrix, meta) {
  # Ensure metadata has SampleID
  if(!"SampleID" %in% colnames(meta)) meta$SampleID <- rownames(meta)
  
  # Convert dist object to matrix if needed
  if(inherits(dist_matrix, "dist")) dist_matrix <- as.matrix(dist_matrix)
  
  # Ensure rownames and colnames match SampleID
  if(is.null(rownames(dist_matrix))) rownames(dist_matrix) <- colnames(dist_matrix) <- meta$SampleID
  if(is.null(colnames(dist_matrix))) colnames(dist_matrix) <- rownames(dist_matrix) <- meta$SampleID
  
  dtb <- data.frame()
  
  for(pid in unique(meta$participant_id)) {
    samp_ids <- meta$SampleID[meta$participant_id == pid]
    
    # Identify baseline sample
    base_samp <- samp_ids[meta$period[meta$SampleID %in% samp_ids] == "Base"]
    
    if(length(base_samp) == 0) next
    
    # Compute distance to baseline for each sample
    for(s in samp_ids) {
      dtb <- rbind(dtb, data.frame(
        participant_id = pid,
        SampleID = s,
        period = meta$period[meta$SampleID == s],
        distance = dist_matrix[s, base_samp]
      ))
    }
  }
  
  return(dtb)
}



# Convert dist object to matrix first
beta_bc_fresh_mat <- as.matrix(beta_bc_fresh)
beta_wu_fresh_mat <- as.matrix(beta_wu_fresh)

beta_bc_ferm_mat <- as.matrix(beta_bc_ferm)
beta_wu_ferm_mat <- as.matrix(beta_wu_ferm)

beta_bc_microb_mat <- as.matrix(beta_bc_microb)
beta_wu_microb_mat <- as.matrix(beta_wu_microb)

#distance-to-baseline calculations
#fresh
dist_bc_fresh <- distance_to_baseline(beta_bc_fresh_mat, meta_fresh_df)
dist_wu_fresh <- distance_to_baseline(beta_wu_fresh_mat, meta_fresh_df)
#ferm
dist_bc_ferm <- distance_to_baseline(beta_bc_ferm_mat, meta_ferm_df)
dist_wu_ferm <- distance_to_baseline(beta_wu_ferm_mat, meta_ferm_df)
#microbiome washouts
dist_bc_microb <- distance_to_baseline(beta_bc_microb_mat, meta_microb_df)
dist_wu_microb <- distance_to_baseline(beta_wu_microb_mat, meta_microb_df)


#Linear Mixed Models (Distance-to-Baseline)
#function fits LMM (period as fixed effect and participant_id as random effects) and returns ANOVA and post-hoc contract
run_lmm <- function(dist_df, response_label = "Distance") {
  # Ensure factors
  dist_df$participant_id <- factor(dist_df$participant_id)
  dist_df$period <- factor(dist_df$period, levels = unique(dist_df$period))
  
  # Fit linear mixed model
  lmm_model <- lmer(distance ~ period + (1 | participant_id), data = dist_df)
  
  # ANOVA table
  anova_res <- anova(lmm_model)
  
  # Post-hoc pairwise comparisons
  emmeans_res <- emmeans(lmm_model, pairwise ~ period)
  
  # Prepare estimated marginal means for plotting
  emm_df <- as.data.frame(emmeans(lmm_model, ~ period))
  emm_df$response <- response_label
  
  # Return as a list
  return(list(
    model = lmm_model,
    anova = anova_res,
    emmeans = emmeans_res,
    emm_df = emm_df
  ))
}

#Fresh Bray-Curtis
lmm_bc_fresh <- run_lmm(dist_bc_fresh, response_label = "Bray-Curtis")
lmm_wu_fresh <- run_lmm(dist_wu_fresh, response_label = "Weighted UniFrac")

# Fermented
lmm_bc_ferm <- run_lmm(dist_bc_ferm, response_label = "Bray-Curtis")
lmm_wu_ferm <- run_lmm(dist_wu_ferm, response_label = "Weighted UniFrac")

# Microbiome washouts
lmm_bc_microb <- run_lmm(dist_bc_microb, response_label = "Bray-Curtis")
lmm_wu_microb <- run_lmm(dist_wu_microb, response_label = "Weighted UniFrac")


#quick check
head(dist_bc_fresh)
str(dist_bc_fresh)
table(dist_bc_fresh$period)
table(dist_bc_fresh$participant_id)

#Plots for beta diversity 
# Function to create combined line + boxplot figure
plot_distance_lmm <- function(dist_df, lmm_res, title_prefix) {
  
  # Extract post-hoc contrasts for significance
  contrasts_df <- as.data.frame(lmm_res$emmeans$contrasts) %>%
    mutate(
      group1 = gsub(" -.*","",contrast),
      group2 = gsub(".*- ","",contrast),
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    # staggered y positions to avoid overlap
    mutate(
      y.position = max(dist_df$distance) + seq(0.05, 0.05*nrow(.), by = 0.05)*max(dist_df$distance)
    ) %>%
    dplyr::select(group1, group2, sig, y.position)
  
  # Line/trajectory plot
  line_plot <- ggplot(dist_df, aes(x = period, y = distance, group = participant_id)) +
    geom_line(alpha = 0.3, color = "grey") +
    stat_summary(aes(group=1), fun = mean, geom = "line", color = "blue", size = 1.2) +
    stat_summary(aes(group=1), fun.data = mean_se, geom = "ribbon", fill = "blue", alpha = 0.2) +
    labs(title = paste0(title_prefix, " - Trajectories"), x = "Period", y = "Distance to Baseline") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  #boxplot/violin with significance
  box_plot <- ggplot(dist_df, aes(x = period, y = distance, fill = period)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 3, color = "black", fill = "yellow") +
    ggpubr::stat_pvalue_manual(contrasts_df, label = "sig", tip.length = 0.01) +
    labs(title = paste0(title_prefix, " - Summary"), x = "Period", y = "Distance to Baseline") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  #combine side by side
  combined_plot <- line_plot | box_plot
  return(combined_plot)
}

#generate plots for all datasets
#fresh
beta_fresh_BC <- plot_distance_lmm(dist_bc_fresh, lmm_bc_fresh, "Fresh Bray-Curtis")
ggsave("D2B_Fresh_Bray-Curtis.png", plot = beta_fresh_BC, width = 6, height = 4, dpi = 300)

beta_fresh_WUF <- plot_distance_lmm(dist_wu_fresh, lmm_wu_fresh, "Fresh Weighted UniFrac")
ggsave("D2B_Fresh_WeightedUF.png", plot = beta_fresh_WUF, width = 6, height = 4, dpi = 300)

#fermented
beta_ferm_BC <- plot_distance_lmm(dist_bc_ferm, lmm_bc_ferm, "Fermented Bray-Curtis")
ggsave("D2B_Ferm_Bray-Curtis.png", plot = beta_ferm_BC, width = 6, height = 4, dpi = 300)

beta_ferm_WUF <- plot_distance_lmm(dist_wu_ferm, lmm_wu_ferm, "Fermented Weighted UniFrac")
ggsave("D2B_Ferm_WeightedUF.png", plot = beta_ferm_WUF, width = 6, height = 4, dpi = 300)


#microbiome washout
beta_microbiome_BC <- plot_distance_lmm(dist_bc_microb, lmm_bc_microb, "Washout Bray-Curtis")
ggsave("D2B_microbiome_Bray-Curtis.png", plot = beta_microbiome_BC, width = 6, height = 4, dpi = 300)

beta_microbiome_WUF <- plot_distance_lmm(dist_wu_microb, lmm_wu_microb, "Washout Weighted UniFrac")
ggsave("D2B_microbiome_WeightedUF.png", plot = beta_microbiome_WUF, width = 6, height = 4, dpi = 300)



#function to extract post-hoc contrasts as a table
get_posthoc_table <- function(lmm_res, response_label = "Distance") {
  contrasts_df <- as.data.frame(lmm_res$emmeans$contrasts) %>%
    dplyr::mutate(
      response = response_label
    )
  return(contrasts_df)
}

#extract post-hoc contrasts 
posthoc_bc_fresh <- get_posthoc_table(lmm_bc_fresh, "Bray-Curtis - Fresh")
posthoc_wu_fresh <- get_posthoc_table(lmm_wu_fresh, "Weighted UniFrac - Fresh")

posthoc_bc_ferm <- get_posthoc_table(lmm_bc_ferm, "Bray-Curtis - Fermented")
posthoc_wu_ferm <- get_posthoc_table(lmm_wu_ferm, "Weighted UniFrac - Fermented")

posthoc_bc_microb <- get_posthoc_table(lmm_bc_microb, "Bray-Curtis - Washout")
posthoc_wu_microb <- get_posthoc_table(lmm_wu_microb, "Weighted UniFrac - Washout")

#combine all into a single table
posthoc_all <- bind_rows(
  posthoc_bc_fresh, posthoc_wu_fresh,
  posthoc_bc_ferm, posthoc_wu_ferm,
  posthoc_bc_microb, posthoc_wu_microb
)

#preview and export
head(posthoc_all)
write.csv(posthoc_all, "BetaDiversity_LMM_posthoc_contrasts.csv", row.names = FALSE)


#### DESeq ####
#controls for repeated measures
design = ~ participant_id + period

#run DESeq2 for each subset
#I also tried to run the DESeq analysis by adding 1 to avoid all zeros, but ran into error
#code below does not include that +1 step
#fresh
dds_fresh <- phyloseq_to_deseq2(ps_fresh, ~ participant_id + period)

#ensure period is a factor with proper reference level
dds_fresh$period <- relevel(dds_fresh$period, "Base")

# Run DESeq2
dds_fresh <- DESeq(dds_fresh)

# Contrasts
res_fresh_VEG_vs_Base <- results(dds_fresh, contrast = c("period", "VEG", "Base"))
res_fresh_WO1_vs_VEG  <- results(dds_fresh, contrast = c("period", "WO1", "VEG"))


#repeat DESEQ analysis for ferm subset
dds_ferm<- phyloseq_to_deseq2(ps_ferm, ~ participant_id + period)
dds_ferm$period <- relevel(dds_ferm$period, "Base")
dds_ferm <- DESeq(dds_ferm)
res_ferm_FERM_vs_Base <- results(dds_ferm, contrast = c("period", "FERM", "Base"))
res_ferm_WO2_vs_FERM  <- results(dds_ferm, contrast = c("period", "WO2", "FERM"))

#repeat DESEQ analysis for microbiome subset
dds_microb<- phyloseq_to_deseq2(ps_microb, ~ participant_id + period)
dds_microb$period <- relevel(dds_microb$period, "Base")
dds_microb <- DESeq(dds_microb)
res_microb_WO1_vs_Base <- results(dds_microb, contrast = c("period", "WO1", "Base"))
res_microb_WO2_vs_WO1  <- results(dds_microb, contrast = c("period", "WO2", "WO1"))


#prepare each subset for LMM analysis using variance-stabilizing transformation (VST)
vst_fresh  <- varianceStabilizingTransformation(dds_fresh, blind = FALSE)
vst_df_fresh  <- as.data.frame(assay(vst_fresh))

vst_ferm   <- varianceStabilizingTransformation(dds_ferm, blind = FALSE)
vst_df_ferm   <- as.data.frame(assay(vst_ferm))

vst_microb <- varianceStabilizingTransformation(dds_microb, blind = FALSE)
vst_df_microb <- as.data.frame(assay(vst_microb))

#convert to long format for merging with metadata
vst_fresh_long <- vst_df_fresh %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_fresh_df, by = "SampleID")

vst_ferm_long <- vst_df_ferm %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_ferm_df, by = "SampleID")

vst_microb_long <- vst_df_microb %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(
    cols = -taxon,
    names_to = "SampleID",
    values_to = "vst_abundance"
  ) %>%
  left_join(meta_microb_df, by = "SampleID")


#identify top significant taxa (n=20)
top_n <- 10

#fresh
res_fresh_df <- as.data.frame(res_fresh_VEG_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%  # move rownames to a column
  select(taxon, padj) %>%                # keep only the columns we need
  mutate(padj = as.numeric(padj)) %>%    # ensure padj is numeric
  filter(!is.na(padj)) %>%               # remove NA padj values
  arrange(padj)                          # sort by smallest padj

top_taxa_fresh <- res_fresh_df %>%
  slice_head(n = top_n) %>%              # select top n taxa
  pull(taxon)

#ferm
res_ferm_df <- as.data.frame(res_ferm_FERM_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%
  select(taxon, padj) %>%
  mutate(padj = as.numeric(padj)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top_taxa_ferm <- res_ferm_df %>%
  slice_head(n = top_n) %>%
  pull(taxon)

#microbiome
res_microb_df <- as.data.frame(res_microb_WO1_vs_Base) %>%
  rownames_to_column(var = "taxon") %>%
  select(taxon, padj) %>%
  mutate(padj = as.numeric(padj)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top_taxa_microb <- res_microb_df %>%
  slice_head(n = top_n) %>%
  pull(taxon)


#filter the VST-long datasets for top taxa
vst_fresh_long_top <- vst_fresh_long %>%
  filter(taxon %in% top_taxa_fresh)
vst_ferm_long_top <- vst_ferm_long %>%
  filter(taxon %in% top_taxa_ferm)
vst_microb_long_top <- vst_microb_long %>%
  filter(taxon %in% top_taxa_microb)


#run LMM per taxon
run_lmm <- function(df) {
  df %>%
    group_by(taxon) %>%
    group_modify(~{
      tryCatch({
        m <- lmer(vst_abundance ~ period + (1 | participant_id), data = .x)
        emm <- emmeans(m, pairwise ~ period)
        tibble(
          model = list(m),
          emm = list(emm),
          contrasts = list(as.data.frame(emm$contrasts))
        )
      }, error = function(e) {
        tibble(model = NA, emm = NA, contrasts = NA)
      })
    })
}

lmm_results_fresh  <- run_lmm(vst_fresh_long_top)
lmm_results_ferm   <- run_lmm(vst_ferm_long_top)
lmm_results_microb <- run_lmm(vst_microb_long_top)



#check the periods used in LMM
unique(vst_fresh_long_top$period)
unique(vst_ferm_long_top$period)
unique(vst_microb_long_top$period)


#check if LMM was sucessful 
lmm_results_fresh %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )

lmm_results_ferm %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )

lmm_results_microb %>%
  mutate(success = !is.na(model)) %>%
  summarise(
    total = n(),
    succeeded = sum(success),
    failed = sum(!success)
  )


#DESeq2 plots
#function to tidy EMMs for plotting
emm_fresh <- lmm_results_fresh %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))

emm_ferm <- lmm_results_ferm %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))

emm_microb <- lmm_results_microb %>%
  dplyr::filter(!is.na(model)) %>%
  dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
  dplyr::select(taxon, emm_df) %>%
  tidyr::unnest(cols = c(emm_df))

#combind plots
plot_emm_overall <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean, group = taxon, color = taxon)) +
    geom_point(size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    theme_bw() +
    labs(
      y = "VST abundance (estimated)",
      x = "Period",
      title = title,
      color = "Taxon"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

#apply function to all subset and save plot
deseq_fresh_combined <- plot_emm_overall(emm_fresh, title = "Top 10 Significant Taxa: Fresh")
ggsave("DESeq2_Fresh_combined.png", plot = deseq_fresh_combined, width = 6, height = 4, dpi = 300)

deseq_ferm_combined <- plot_emm_overall(emm_ferm, title = "Top 10 Significant Taxa: Fermented")
ggsave("DESeq2_Ferm_combined.png", plot = deseq_ferm_combined, width = 6, height = 4, dpi = 300)

deseq_microb_combined <- plot_emm_overall(emm_microb, title = "Top 10 Significant Taxa: Microbiome Washout")
ggsave("DESeq2_Microb_combined.png", plot = deseq_microb_combined, width = 6, height = 4, dpi = 300)



#facet plots
plot_emm_facets <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean)) +
    geom_point(size = 3, color = "steelblue") +
    geom_line(aes(group = 1), color = "steelblue") +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    facet_wrap(~ taxon, scales = "free_y") +
    theme_bw() +
    labs(
      y = "VST abundance (estimated)",
      x = "Period",
      title = title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#apply function to all subset and save plot
deseq_fresh_facet <- plot_emm_facets(emm_fresh, title = "Top 10 Significant Taxa: Fresh")
ggsave("DESeq2_Fresh_facet.png", plot = deseq_fresh_facet, width = 6, height = 4, dpi = 300)

deseq_ferm_facet <- plot_emm_facets(emm_ferm, title = "Top 10 Significant Taxa: Fermented")
ggsave("DESeq2_Ferm_facet.png", plot = deseq_ferm_facet, width = 6, height = 4, dpi = 300)

deseq_microb_facet <- plot_emm_facets(emm_microb, title = "Top 10 Significant Taxa: Microbiome Washout")
ggsave("DESeq2_Microb_facet.png", plot = deseq_microb_facet, width = 6, height = 4, dpi = 300)


#create function to merge taxa
add_taxon_labels_safe <- function(lmm_results, physeq_obj) {
  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(physeq_obj))
  tax_df$taxon <- rownames(tax_df)
  
  # Extract emmeans from LMM
  emm_df <- lmm_results %>%
    dplyr::filter(!is.na(model)) %>%
    dplyr::mutate(emm_df = purrr::map(emm, ~ as.data.frame(.x$emmeans))) %>%
    dplyr::select(taxon, emm_df) %>%
    tidyr::unnest(cols = c(emm_df))
  
  # Merge taxonomy
  emm_named <- emm_df %>%
    left_join(tax_df %>% select(taxon, Family, Genus), by = "taxon") %>%
    # Create label with fallback and uniqueness
    mutate(
      taxon_label = case_when(
        !is.na(Genus) ~ paste(Family, Genus, sep = "_"),
        !is.na(Family) ~ Family,
        TRUE ~ taxon
      ),
      # Append OTU ID to make labels unique
      taxon_label = paste0(taxon_label, "_", taxon)
    )
  
  return(emm_named)
}


#apply function to subset
emm_fresh_named  <- add_taxon_labels_safe(lmm_results_fresh, ps_fresh)
emm_ferm_named   <- add_taxon_labels_safe(lmm_results_ferm, ps_ferm)
emm_microb_named <- add_taxon_labels_safe(lmm_results_microb, ps_microb)


#update line plots
plot_emm_overall_safe <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean, group = taxon_label, color = taxon_label)) +
    geom_point(size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    theme_bw() +
    labs(
      y = "VST abundance (estimated)",
      x = "Period",
      title = title,
      color = "Taxon"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

#generate and save 
deseq_fresh_combined  <- plot_emm_overall_safe(emm_fresh_named, title = "Top 10 Significant Taxa: Fresh")
deseq_ferm_combined   <- plot_emm_overall_safe(emm_ferm_named, title = "Top 10 Significant Taxa: Fermented")
deseq_microb_combined <- plot_emm_overall_safe(emm_microb_named, title = "Top 10 Significant Taxa: Microbiome Washout")

ggsave("DESeq2_Fresh_combined_safe.png", plot = deseq_fresh_combined, width = 8, height = 6, dpi = 300)
ggsave("DESeq2_Ferm_combined_safe.png", plot = deseq_ferm_combined, width = 8, height = 6, dpi = 300)
ggsave("DESeq2_Microb_combined_safe.png", plot = deseq_microb_combined, width = 8, height = 6, dpi = 300)


#update facet plots
plot_emm_facets_safe <- function(emm_df, title = "LMM Top Taxa") {
  ggplot(emm_df, aes(x = period, y = emmean)) +
    geom_point(size = 3, color = "steelblue") +
    geom_line(aes(group = 1), color = "steelblue") +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    facet_wrap(~ taxon_label, scales = "free_y") +
    theme_bw() +
    labs(
      y = "VST abundance (estimated)",
      x = "Period",
      title = title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#generate plots and save
deseq_fresh_facet  <- plot_emm_facets_safe(emm_fresh_named, title = "Top 10 Significant Taxa: Fresh")
deseq_ferm_facet   <- plot_emm_facets_safe(emm_ferm_named, title = "Top 10 Significant Taxa: Fermented")
deseq_microb_facet <- plot_emm_facets_safe(emm_microb_named, title = "Top 10 Significant Taxa: Microbiome Washout")

ggsave("DESeq2_Fresh_facet_safe.png", plot = deseq_fresh_facet, width = 8, height = 6, dpi = 300)
ggsave("DESeq2_Ferm_facet_safe.png", plot = deseq_ferm_facet, width = 8, height = 6, dpi = 300)
ggsave("DESeq2_Microb_facet_safe.png", plot = deseq_microb_facet, width = 8, height = 6, dpi = 300)
