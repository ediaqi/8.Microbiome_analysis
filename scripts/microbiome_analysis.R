#################################
# Microbiome analysis in R: alpha diversity, beta diversity, and relative abundance
# Author: Iva Sunic
# Date: 20.5.2025.
#################################

# Loading and installing required packages
packages <- c("phyloseq", "ggplot2", "vegan", "ggpubr", "dplyr", "tidyr", "RColorBrewer", "readr")
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

# Loading data (ASV, taxonomy, metadata)
otu_table <- read.delim("example_asv.tsv", row.names = 1, check.names = FALSE)
taxonomy <- read.delim("example_taxonomy.tsv", row.names = 1)
metadata <- read.csv("example_metadata.csv", row.names = 1)

# Creating phyloseq object
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(taxonomy))
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, SAM)


# Alpha Diversity
#################################

# Calculations
alpha_df <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_df$SampleID <- rownames(alpha_df)
alpha_df <- cbind(sample_data(ps)[rownames(alpha_df), ], alpha_df)
alpha_df

# Plot
alpha_plot <- plot_richness(ps, x = "Healthy", measures = c("Shannon", "Simpson")) +
  theme_minimal() +
  ggtitle("Alpha Diversity")
alpha_plot


# Beta Diversity
#################################

# Calculation of Bray-Curtis distances
bray_dist <- phyloseq::distance(ps, method = "bray")
bray_matrix <- as.matrix(bray_dist)
bray_matrix

# PCoA Ordination Plot
ord <- ordinate(ps, method = "PCoA", distance = "bray")
beta_plot <- plot_ordination(ps, ord, color = "Healthy") +
  geom_point(size = 3) +
  ggtitle("Beta Diversity (PCoA - Bray-Curtis)") +
  theme_minimal()
beta_plot

# Relative Abundance: Barplot
#################################

# Transformation to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Top 10 taxa
top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:10])
ps_top <- prune_taxa(top_taxa, ps_rel)

# Barplot
bar_plot <- plot_bar(ps_top, fill = "Phylum") +
  facet_wrap(~Healthy) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Top 10 Taxa - Relative Abundance")
bar_plot


# Relative Abundance: Violin Plot
#################################

# Violin plot of the most abundant taxons
top_taxon <- top_taxa[1]
df <- psmelt(ps_rel)
violin_df <- df %>% filter(OTU == top_taxon)

violin_plot <- ggplot(violin_df, aes(x = Healthy, y = Abundance)) +
  geom_violin(fill = "skyblue") +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  ggtitle(paste("Violin Plot for", top_taxon))
violin_plot

