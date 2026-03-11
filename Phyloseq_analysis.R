
# ============================================================================== #
# ITS2 Symbiodiniaceae Diversity Analysis
# Ningaloo 2023 — Phyloseq Pipeline
#
# Description:
#   Loads pre-processed DADA2/LULU output data, calls phyloseq object created before,
#   performs basic quality filtering based on sequencing depth, and prepares
#   the object for downstream diversity analyses.
#
# Input files:
#   - Metadata_Ningaloo_2023.csv      : Sample metadata
#   - ASVtablepostLULU.csv            : ASV table (post-LULU curation)
#   - BlastTable.2023ITS2.PRELULU.csv : BLASTn taxonomy assignments
#   - LULU_phyloseq_object_2023.rds   : Saved phyloseq object (pre-built)
#
# ============================================================================== #

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggh4x) 
library(vegan)

mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2), "cm"),
                 #plot.title = element_text(size = 11, face = "bold"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 10, hjust = 0.5),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 10),
                 legend.text.align = 0,
                 legend.title = element_text(size = 10),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 10, hjust = 0),
                 panel.spacing = unit(.5, "cm"),
                 text = element_text(family = "Helvetica Neue"))

# ============================================================================== #
# 1. READ IN SEQUENCING DATA ####
# ============================================================================== #

seqtab.nochim <- read.csv(
  "/Users/teresa/Desktop/WAsymbionts/ITS2_data/OutputDADA_AllASVs.csv",
  sep     = ",",
  row.names = 1
)

# ============================================================================== #
# 2. CONSTRUCT PHYLOSEQ OBJECT (see archived code below) ####
# ============================================================================== #

# NOTE: The phyloseq object was constructed in a previous session and saved
# as an .rds file. The code used to build it is preserved below for reference.

# -- 2a. Sample Metadata -------------------------------------------------------
# meta <- read.csv(
#   "/Users/teresa/Desktop/WAsymbionts/ITS2_data/Metadata_Ningaloo_2023.csv",
#   sep       = ",",
#   row.names = 1
# )
#
# # Move row names to an explicit SampleID column and reorder to front
# meta$SampleID <- rownames(meta)
# rownames(meta) <- NULL
# meta <- meta[, c("SampleID", setdiff(colnames(meta), "SampleID"))]
# head(meta)

# -- 2b. OTU Table (post-LULU) -------------------------------------------------
# otu <- as.matrix(
#   read.csv(
#     "/Users/teresa/Desktop/WAsymbionts/ITS2_data/ASVtablepostLULU.csv",
#     sep       = ",",
#     row.names = 1
#   )
# )
#
# # Ensure all values are numeric and restore row/column names
# otu <- apply(otu, 2, as.numeric)
# rownames(otu) <- rownames(as.matrix(read.csv("ASVtablepostLULU.csv", sep = ",", row.names = 1)))
# colnames(otu) <- colnames(as.matrix(read.csv("ASVtablepostLULU.csv", sep = ",", row.names = 1)))

# -- 2c. Taxonomy Table --------------------------------------------------------
# taxa <- as.matrix(
#   read.csv(
#     "/Users/teresa/Desktop/WAsymbionts/ITS2_data/BlastTable.2023ITS2.PRELULU.csv",
#     sep       = ",",
#     row.names = 1
#   )
# )
# head(taxa)

# -- 2d. Assemble Phyloseq Object ----------------------------------------------
# PS <- phyloseq(
#   otu_table(otu,  taxa_are_rows = FALSE),
#   sample_data(meta),
#   tax_table(taxa)
# )

# ============================================================================== #
# 3. LOAD PRE-BUILT PHYLOSEQ OBJECT ####
# ============================================================================== #


ps <- readRDS("/Users/teresa/Desktop/WAsymbionts/ITS2_data/LULU_phyloseq_object_2023.rds")

# ============================================================================== #
# 4. QUALITY FILTERING — SEQUENCING DEPTH ####
# ============================================================================== #

# Inspect the distribution of per-sample sequencing depth
hist(sample_sums(ps),
     main = "Distribution of Sequencing Depth",
     xlab = "Total Reads per Sample",
     col  = "steelblue")

# Remove samples below the minimum read-depth threshold
# (typically a small number of ITS2 samples fall below 1,000 reads)
ps2 <- prune_samples(sample_sums(ps) >= 10000, ps)

# Verify filtering result
hist(sample_sums(ps2),
     main = "Sequencing Depth After Filtering (≥10,000 reads)",
     xlab = "Total Reads per Sample",
     col  = "steelblue")

# NOTE: Taxonomic filtering is not required here — all ASVs returned matches
# to Symbiodiniaceae and no non-target sequences were detected.
# The line below is retained for reference only:
# ps3 <- subset_taxa(ps2, (Class != "NA"))

# Rename filtered object for downstream use
ps <- ps2

sample_data(ps)
tax_table(ps)
otu_table(ps)

# ============================================================================== #
# 5. RELATIVE ABUNDANCE ####
# ============================================================================== #

#change to relative abundance within each sample, then melt to use with ggplot
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df <- psmelt(ps_rel)

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Type_label)) 
  
df <- df %>%
  mutate(Type = recode(Type,
                            "C1232.."    = "Cladocopium C1232",
                            "C21.12."    = "Cladocopium C21.12",
                            "C21.3d.C3k" = "Cladocopium C21.3d.C3k",
                            "C3.."       = "Cladocopium C3" ))
df <- df %>%
   mutate(Species_label = recode(Species,
                                      "A.millepora" = "italic('A. millepora')",
                                      "A.tenuis"    = "italic('A. tenuis')" ))

legend_labels <- list(
    "Cladocopium C1232"     = expression(italic("Cladocopium") ~ "C1232"),
    "Cladocopium C21.12"    = expression(italic("Cladocopium") ~ "C21.12"),
    "Cladocopium C21.3d.C3k"= expression(italic("Cladocopium") ~ "C21.3d.C3k"),
    "Cladocopium C3"        = expression(italic("Cladocopium") ~ "C3"))

rel_abundance <- ggplot(df, aes(x = Sample, y = Abundance, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Species_label * Site, scales = "free_x",
             labeller = labeller(Species_label = label_parsed)) +
  scale_fill_manual(
    values = c(
      "Cladocopium C1232"      = "#8FA8BF",
      "Cladocopium C21.12"     = "#C4876A",
      "Cladocopium C21.3d.C3k" = "#6B8C6B",
      "Cladocopium C3"         = "#C9A96E"),
    name   = "Symbiodiniaceae species",
    labels = legend_labels) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)) +
  ylab("Relative Abundance") +
  theme(axis.title.y = element_text(size = 10, hjust = 0.5)) +
  xlab("") +
  xlab("") +
  scale_x_discrete(labels = NULL) + # remove this to bring back sample labels
# scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(axis.ticks.x = element_blank()) +
  theme(
        panel.spacing.y = unit(0.1, "cm"),
        strip.text      = element_text(size = 10, hjust = 0, lineheight = 0.85),
        axis.ticks.x    = element_blank(),
        axis.title.y    = element_text(size = 10, hjust = 0.5)) +
  mytheme

# Same plot but facets in a row

rel_abundance_row <- ggplot(df, aes(x = Sample, y = Abundance, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Species_label * Site, scales = "free_x",
             labeller = labeller(Species_label = label_parsed),
             nrow = 1) +
  scale_fill_manual(
    values = c(
      "Cladocopium C1232"      = "#8FA8BF",
      "Cladocopium C21.12"     = "#C4876A",
      "Cladocopium C21.3d.C3k" = "#6B8C6B",
      "Cladocopium C3"         = "#C9A96E"),
    name   = "Symbiodiniaceae species",
    labels = legend_labels) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)) +
  ylab("Relative Abundance") +
  theme(axis.title.y = element_text(size = 10, hjust = 0.5)) +
  xlab("") +
  xlab("") +
  scale_x_discrete(labels = NULL) + # remove this to bring back sample labels
  # scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(axis.ticks.x = element_blank()) +
  theme(
    panel.spacing.y = unit(0.1, "cm"),
    strip.text      = element_text(size = 10, hjust = 0, lineheight = 0.85),
    axis.ticks.x    = element_blank(),
    axis.title.y    = element_text(size = 10, hjust = 0.5)) +
  mytheme

print(rel_abundance_row)

rel_abundance %>%
  ggsave(
    filename = "rel_abundance.pdf",
    path     = "figures",
    device   = cairo_pdf,
    units    = "cm",
    width    = 16,
    height   = 11
  )

# ============================================================================== #
# 6. ALPHA DIVERSITY ####
# ============================================================================== #

# plot_richness() is native to phyloseq, 
# works directly on my ps object and calculates multiple indices at once. 

sample_data(ps)

alpha_div <- plot_richness(
  ps,
  x       = "Species",
  measures = c("Shannon", "Simpson"),
  color   = "Site"
)$data

# Rename index labels for clean facet titles (optional — matches reference fig)
alpha_div <- alpha_div %>%
  mutate(variable = recode(variable,
                         "Shannon" = "Shannon Diversity",
                         "Simpson" = "Simpson Diversity"))

# Site colours matching the relative abundance palette
site_colours <- c(
  "Oyster Stacks" = "#C4876A",
  "Pelican Point"  = "#8FA8BF"
)

alpha_plot <- ggplot(alpha_div, aes(x = Species, y = value, colour = Site)) +
  geom_boxplot(
    aes(colour = Site),
    outlier.shape = NA,
    width         = 0.5,
    position      = position_dodge(width = 0.7)) +
  geom_jitter(
    aes(colour = Site),
    size     = 1.5,
    alpha    = 0.8,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_x_discrete(
    labels = c(
      "A.millepora" = expression(italic("A. millepora")),
      "A.tenuis"    = expression(italic("A. tenuis")))) +
  scale_colour_manual(values = site_colours, name = "Site") +
  labs(
    title = "Alpha Diversity",
    x     = NULL,
    y     = "Diversity Index") +
  mytheme +
  theme(plot.title = element_text(face = "bold")) +
  theme(
    legend.position    = c(1, 0.05),   # x=right, y=bottom in normalised coords
    legend.justification = c("right", "bottom"))

alpha_plot <- alpha_plot +
    facetted_pos_scales(y = list(
        variable == "Shannon Diversity" ~ scale_y_continuous(
            limits = c(0, 2),
            expand = c(0, 0),
            labels = c("", "0.5", "1", "1.5", "2")),
        variable == "Simpson Diversity" ~ scale_y_continuous(
            limits = c(0, 1),
            expand = c(0, 0),
            labels = c("", "0.25", "0.5", "0.75", "1")) ))
print(alpha_plot)

# STATS
erich <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Extract sample IDs and diversity estimates
Sample   <- row.names(erich)
alpha    <- data.frame(Sample, erich)

# Extract sample metadata and ensure Sample column name matches
s <- data.frame(sample_data(ps))
s$Sample <- row.names(s)

# Merge diversity estimates with metadata
alphadiv <- merge(alpha, s, by = "Sample")


# SUMMARY MEANS BY GROUP

cat("\n--- Mean Shannon by Site ---\n")
tapply(alphadiv$Shannon, alphadiv$Site, mean)

cat("\n--- Mean Simpson by Site ---\n")
tapply(alphadiv$Simpson, alphadiv$Site, mean)

cat("\n--- Mean Shannon by Species ---\n")
tapply(alphadiv$Shannon, alphadiv$Species, mean)

cat("\n--- Mean Simpson by Species ---\n")
tapply(alphadiv$Simpson, alphadiv$Species, mean)

cat("\n--- Mean Shannon by Site x Species ---\n")
tapply(alphadiv$Shannon, list(alphadiv$Site, alphadiv$Species), mean)

cat("\n--- Mean Simpson by Site x Species ---\n")
tapply(alphadiv$Simpson, list(alphadiv$Site, alphadiv$Species), mean)


# TEST FOR NORMALITY (Shapiro-Wilk)

cat("\n--- Shapiro-Wilk Normality Tests ---\n")

shapiro_shannon <- shapiro.test(alphadiv$Shannon)
shapiro_simpson <- shapiro.test(alphadiv$Simpson)

cat("Shannon: W =", shapiro_shannon$statistic,
    " p =", shapiro_shannon$p.value, "\n")
cat("Simpson: W =", shapiro_simpson$statistic,
    " p =", shapiro_simpson$p.value, "\n")

# Flag normality (p > 0.05 = normal)
shannon_normal <- shapiro_shannon$p.value > 0.05
simpson_normal <- shapiro_simpson$p.value > 0.05

cat("\nShannon normal:", shannon_normal,
    "| Simpson normal:", simpson_normal, "\n")

# Both indices are non-normal (Shapiro p < 0.05), so Wilcoxon rank-sum tests
# are used for two-group comparisons (2 sites, 2 species).

cat("\n--- Wilcoxon Tests: Site ---\n")
print(wilcox.test(Shannon ~ Site, data = alphadiv))
print(wilcox.test(Simpson ~ Site, data = alphadiv))

cat("\n--- Wilcoxon Tests: Species ---\n")
print(wilcox.test(Shannon ~ Species, data = alphadiv))
print(wilcox.test(Simpson ~ Species, data = alphadiv))


# INTERACTION: SITE x SPECIES (4 GROUPS)
# Kruskal-Wallis if non-normal (>2 groups); two-way ANOVA if normal
# Create an interaction factor for Kruskal-Wallis.

alphadiv$Site_Species <- interaction(alphadiv$Site, alphadiv$Species)

cat("\n--- Interaction Tests: Site x Species ---\n")

if (shannon_normal) {
  cat("Shannon ~ Site * Species: two-way ANOVA\n")
  print(summary(aov(Shannon ~ Site * Species, data = alphadiv)))
} else {
  cat("Shannon ~ Site x Species: Kruskal-Wallis\n")
  print(kruskal.test(Shannon ~ Site_Species, data = alphadiv))
}

if (simpson_normal) {
  cat("Simpson ~ Site * Species: two-way ANOVA\n")
  print(summary(aov(Simpson ~ Site * Species, data = alphadiv)))
} else {
  cat("Simpson ~ Site x Species: Kruskal-Wallis\n")
  print(kruskal.test(Simpson ~ Site_Species, data = alphadiv))
}

# With more than two groups, Kruskal-Wallis is used.
# Creating interaction factor to test all four site x species combinations

alphadiv$Site_Species <- interaction(alphadiv$Site, alphadiv$Species)

cat("\n--- Kruskal-Wallis Tests: Site x Species ---\n")
print(kruskal.test(Shannon ~ Site_Species, data = alphadiv))
print(kruskal.test(Simpson ~ Site_Species, data = alphadiv))

# ============================================================================== #
# 7. BETA DIVERSITY ####
# ============================================================================== #

# Extract OTU table (samples as rows) and metadata
otu  <- as.data.frame(otu_table(ps))
if (!taxa_are_rows(ps)) otu <- otu  else otu <- t(otu)
meta <- as.data.frame(sample_data(ps))


# COMPUTE DISSIMILARITY MATRICES

# -- Bray-Curtis (abundance-weighted) 
dist_bray <- vegdist(otu, method = "bray")

# -- Jaccard (presence-absence) 
# Binary transform converts counts to 0/1 before computing Jaccard distance
otu_pa    <- decostand(otu, method = "pa")   # presence-absence transform
dist_jacc <- vegdist(otu_pa, method = "jaccard", binary = TRUE)

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

# PERMANOVA — do groups differ in centroid?
# Tests whether Species, Site, and their interaction explain significant
# variation in symbiont community composition.
# set.seed() ensures reproducibility of permutation results.

set.seed(123)

cat("\n--- PERMANOVA: Bray-Curtis ---\n")
perm_bray <- adonis2(
  dist_bray ~ Species * Site,
  data       = meta,
  permutations = 999
)
print(perm_bray)

cat("\n--- PERMANOVA: Jaccard ---\n")
perm_jacc <- adonis2(
  dist_jacc ~ Species * Site,
  data       = meta,
  permutations = 999
)
print(perm_jacc)

# PERMDISP — do groups differ in dispersion?

# PERMANOVA can be significant due to differences in dispersion rather than
# location. betadisper tests whether within-group spread differs between groups.
# If significant, interpret PERMANOVA results with caution.

# Create combined Site x Species grouping factor
group_interaction <- interaction(meta$Species, meta$Site)

cat("\n--- PERMDISP: Bray-Curtis ---\n")
disp_bray <- betadisper(dist_bray, group_interaction)
print(permutest(disp_bray, permutations = 999))

cat("\n--- PERMDISP: Jaccard ---\n")
disp_jacc <- betadisper(dist_jacc, group_interaction)
print(permutest(disp_jacc, permutations = 999))

# ============================================================================== #
# 8. ORDINATION ####
# ============================================================================== #

site_colours <- c(
  "Oyster Stacks" = "#C4876A",
  "Pelican Point"  = "#8FA8BF"
)
species_shapes <- c("A.millepora" = 16, "A.tenuis" = 17)  # circle, triangle

# PCoA
pcoa_bray <- ordinate(ps, method = "PCoA", distance = dist_bray)
pcoa_jacc <- ordinate(ps, method = "PCoA", distance = dist_jacc)

# Extract % variance explained for axis labels
bray_eig <- round(pcoa_bray$values$Relative_eig[1:2] * 100, 1)
jacc_eig <- round(pcoa_jacc$values$Relative_eig[1:2] * 100, 1)

# Build PCoA data frames
pcoa_bray_df <- data.frame(
  Axis1   = pcoa_bray$vectors[, 1],
  Axis2   = pcoa_bray$vectors[, 2],
  Site    = meta$Site,
  Species = meta$Species
)

pcoa_jacc_df <- data.frame(
  Axis1   = pcoa_jacc$vectors[, 1],
  Axis2   = pcoa_jacc$vectors[, 2],
  Site    = meta$Site,
  Species = meta$Species
)

## 8.1. PCoA plot — Bray-Curtis ####
plot_pcoa_bray <- ggplot(pcoa_bray_df, aes(x = Axis1, y = Axis2,
                                           colour = Site, shape = Species)) +
  stat_ellipse(aes(group = interaction(Species, Site), colour = Site),
               type = "t", level = 0.95, linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.6) +
  scale_colour_manual(values = site_colours) +
  scale_shape_manual(
    values = species_shapes,
    labels = c(
      "A.millepora" = expression(italic("A. millepora")),
      "A.tenuis"    = expression(italic("A. tenuis"))
    )
  ) +
  labs(
    title = "PCoA — Bray-Curtis",
    x     = paste0("PCoA1 (", bray_eig[1], "%)"),
    y     = paste0("PCoA2 (", bray_eig[2], "%)")
  ) +
  mytheme +
  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position = "none") # only keep it for plotting combined figure


## 8.2. PCoA plot — Jaccard ####

plot_pcoa_jacc <- ggplot(pcoa_jacc_df, aes(x = Axis1, y = Axis2,
                                           colour = Site, shape = Species)) +
  stat_ellipse(aes(group    = interaction(Species, Site),
                   colour   = Site,
                   linetype = Species),
               type = "norm", level = 0.95, linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.6) +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(
    name   = "Species",
    values = species_shapes,
    labels = c(
      "A.millepora" = expression(italic("A. millepora")),
      "A.tenuis"    = expression(italic("A. tenuis")))) +
  scale_linetype_manual(
    name   = "Species",
    values = c("A.millepora" = "solid", "A.tenuis" = "dashed"),
    labels = c(
      "A.millepora" = expression(italic("A. millepora")),
      "A.tenuis"    = expression(italic("A. tenuis")))) +
  labs(
    title = "PCoA — Jaccard",
    x     = paste0("PCoA1 (", jacc_eig[1], "%)"),
    y     = paste0("PCoA2 (", jacc_eig[2], "%)")) +
  mytheme

## 8.3. NMDS ####

set.seed(123)
nmds_bray <- metaMDS(otu,    distance = "bray",    k = 2, trymax = 100)
nmds_jacc <- metaMDS(otu_pa, distance = "jaccard", k = 2, trymax = 100,
                                             binary = TRUE)

cat("\nNMDS stress (Bray-Curtis):", nmds_bray$stress,
           "\nNMDS stress (Jaccard):    ", nmds_jacc$stress, "\n")

# Stress < 0.1 = good, < 0.2 = acceptable, > 0.2 = interpret with caution
# Not keeping Jaccard

nmds_bray_df <- data.frame(
    NMDS1   = nmds_bray$points[, 1],
    NMDS2   = nmds_bray$points[, 2],
    Site    = meta$Site,
    Species = meta$Species)

set.seed(123)
nmds_bray <- metaMDS(otu,    method = "bray", k = 2, trymax = 100)

cat("\nNMDS stress (Bray-Curtis):", nmds_bray$stress)
# Stress < 0.1 = good, < 0.2 = acceptable, > 0.2 = interpret with caution

nmds_bray$converged   # should be TRUE
nmds_bray$stress      # should be 0.05767

nmds_bray_df <- data.frame(
  NMDS1   = nmds_bray$points[, 1],
  NMDS2   = nmds_bray$points[, 2],
  Site    = meta$Site,
  Species = meta$Species
)

## 8.4. NMDS plot — Bray-Curtis ####

plot_nmds_bray <- ggplot(nmds_bray_df, aes(x = NMDS1, y = NMDS2,
                                           colour = Site, shape = Species)) +
  stat_ellipse(aes(group    = interaction(Species, Site),
                   colour   = Site,
                   linetype = Species),
               type = "norm", level = 0.95, linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.6) +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(
    name   = "Species",
    values = species_shapes,
    labels = c(
      "A.millepora" = expression(italic("A. millepora")),
      "A.tenuis"    = expression(italic("A. tenuis"))
    ),
    # Merge linetype into the shape legend via override.aes
    guide  = guide_legend(override.aes = list(
      linetype = c("solid", "dashed"),
      shape    = c(16, 17)
    ))
  ) +
  # Suppress separate linetype legend — merged into shape legend above
  scale_linetype_manual(
    values = c("A.millepora" = "solid", "A.tenuis" = "dashed"),
    guide  = "none"
  ) +
  annotate("text", x = Inf, y = -Inf,
           label = paste("Stress =", round(nmds_bray$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 3, colour = "grey40") +
  labs(title = "NMDS — Bray-Curtis", x = "NMDS1", y = "NMDS2") +
  theme(plot.title = element_text(face = "bold")) +
  mytheme

## 8.5. Together ####
# Plot PCoA Bray-Curtis | NMDS Bray-Curtis side by side.
# Jaccard PCoA omitted — Jaccard results reported in PERMANOVA table only.

library(patchwork)

final_figure <- rel_abundance_row / (alpha_plot | plot_pcoa_bray | plot_nmds_bray) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag   = element_text(size = 11, face = "bold", family = "Helvetica Neue"),
      plot.title = element_text(size = 11, face = "bold", family = "Helvetica Neue")))

final_figure <- final_figure & 
  theme(
    text        = element_text(size = 14, family = "Helvetica Neue"),
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 13),
    strip.text  = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    plot.tag    = element_text(size = 13, face = "bold")
  )

ggsave(
  filename = "Combined_RelAbund_AlphaBetaDiv.pdf",
  plot     = final_figure,
  width    = 18,
  height   = 12,
  units    = "in",
  device   = cairo_pdf  
)

# Font sizes edited in Affinity Designer


