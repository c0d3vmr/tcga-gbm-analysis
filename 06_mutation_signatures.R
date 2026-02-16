# ============================================
# Script 06: Mutation Signatures & Patterns
# Analyze specific mutation patterns and hotspots
# Date: 2024
# ============================================

library(tidyverse)
library(maftools)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  MUTATION SIGNATURES & PATTERNS            ║\n")
cat("║  Deep dive into specific mutations         ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD DATA
# ============================================

cat("Loading data...\n")

maf_data <- readRDS("data/raw/gbm_mutations.rds")
clinical <- read.csv("data/raw/clinical_data.csv")
pathway_genes <- read.csv("data/reference/pathway_genes.csv")

# Prepare clinical data
clinical_for_maf <- clinical %>%
  mutate(Tumor_Sample_Barcode = bcr_patient_barcode)

# Create MAF object
maf <- read.maf(
  maf = maf_data,
  clinicalData = clinical_for_maf,
  isTCGA = TRUE,
  verbose = FALSE
)

cat(paste("✓ Loaded", nrow(maf_data), "mutations\n"))
cat(paste("✓ Loaded", length(unique(maf_data$Tumor_Sample_Barcode)), "samples\n\n"))

# ============================================
# 2. IDENTIFY HOTSPOT MUTATIONS
# ============================================

cat("Identifying mutation hotspots...\n")

# Find recurrent mutations (same gene + position)
hotspots <- maf_data %>%
  filter(Hugo_Symbol %in% pathway_genes$Gene) %>%
  mutate(
    mutation_id = paste(Hugo_Symbol, 
                        HGVSp_Short, 
                        sep = "_")
  ) %>%
  count(Hugo_Symbol, HGVSp_Short, mutation_id, Variant_Classification) %>%
  filter(n >= 2) %>%  # Occurs in at least 2 patients
  arrange(desc(n))

cat(paste("✓ Found", nrow(hotspots), "recurrent mutations (≥2 patients)\n\n"))

if(nrow(hotspots) > 0) {
  cat("Top 10 hotspot mutations:\n")
  for(i in 1:min(10, nrow(hotspots))) {
    cat(sprintf("  %2d. %s %s: %d patients (%s)\n",
                i,
                hotspots$Hugo_Symbol[i],
                hotspots$HGVSp_Short[i],
                hotspots$n[i],
                hotspots$Variant_Classification[i]))
  }
  cat("\n")
}

write.csv(hotspots, "results/tables/Table10_Hotspot_Mutations.csv", row.names = FALSE)

# ============================================
# 3. FIGURE 9: TP53 LOLLIPOP PLOT
# ============================================

cat("Creating Figure 9: TP53 Mutation Lollipop Plot...\n")

# Check if TP53 has mutations
tp53_muts <- maf_data %>% filter(Hugo_Symbol == "TP53")

if(nrow(tp53_muts) > 0) {
  png("results/figures/Figure09_TP53_Lollipop.png", 
      width = 14, height = 8, res = 300, units = "in")
  
  lollipopPlot(
    maf = maf,
    gene = "TP53",
    AACol = "HGVSp_Short",
    showMutationRate = TRUE,
    labelPos = "all",
    refSeqID = "NM_000546"
  )
  
  dev.off()
  
  cat("✓ Figure 9 saved: TP53 lollipop plot\n\n")
} else {
  cat("⚠ No TP53 mutations found for lollipop plot\n\n")
}

# ============================================
# 4. FIGURE 10: PTEN MUTATION DISTRIBUTION
# ============================================

cat("Creating Figure 10: PTEN Mutation Distribution...\n")

# Get PTEN mutations
pten_muts <- maf_data %>% 
  filter(Hugo_Symbol == "PTEN") %>%
  mutate(
    mutation_type = Variant_Classification,
    protein_change = HGVSp_Short
  )

if(nrow(pten_muts) > 0) {
  
  # Create visualization
  png("results/figures/Figure10_PTEN_Distribution.png", 
      width = 14, height = 8, res = 300, units = "in")
  
  # Try lollipop plot for PTEN
  tryCatch({
    lollipopPlot(
      maf = maf,
      gene = "PTEN",
      AACol = "HGVSp_Short",
      showMutationRate = TRUE,
      labelPos = "all"
    )
  }, error = function(e) {
    # If lollipop fails, create a bar plot instead
    pten_summary <- pten_muts %>%
      count(Variant_Classification) %>%
      arrange(desc(n))
    
    par(mar = c(5, 10, 4, 2))
    barplot(pten_summary$n, 
            names.arg = pten_summary$Variant_Classification,
            las = 2, 
            main = "PTEN Mutation Types",
            xlab = "Number of Patients",
            col = brewer.pal(8, "Set2"),
            horiz = TRUE)
  })
  
  dev.off()
  
  cat("✓ Figure 10 saved: PTEN mutation distribution\n\n")
  
} else {
  cat("⚠ No PTEN mutations found\n\n")
}

# ============================================
# 5. MUTATION IMPACT CLASSIFICATION
# ============================================

cat("Classifying mutation functional impact...\n")

# Classify mutations by predicted impact
mutation_impact <- maf_data %>%
  filter(Hugo_Symbol %in% pathway_genes$Gene) %>%
  mutate(
    impact = case_when(
      Variant_Classification %in% c("Nonsense_Mutation", "Frame_Shift_Del", 
                                    "Frame_Shift_Ins") ~ "High (Loss of function)",
      Variant_Classification %in% c("Missense_Mutation", "In_Frame_Del", 
                                    "In_Frame_Ins") ~ "Moderate (Altered function)",
      Variant_Classification %in% c("Splice_Site", "Splice_Region") ~ "Moderate (Splicing)",
      Variant_Classification %in% c("Silent", "3'UTR", "5'UTR", 
                                    "Intron") ~ "Low (Silent/Non-coding)",
      TRUE ~ "Unknown"
    )
  )

impact_summary <- mutation_impact %>%
  count(impact, Hugo_Symbol) %>%
  group_by(Hugo_Symbol) %>%
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  ) %>%
  ungroup()

cat("\nMutation impact distribution:\n")
impact_overall <- mutation_impact %>%
  count(impact) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

for(i in 1:nrow(impact_overall)) {
  cat(sprintf("  %-30s: %4d (%.1f%%)\n",
              impact_overall$impact[i],
              impact_overall$n[i],
              impact_overall$percentage[i]))
}
cat("\n")

write.csv(impact_summary, "results/tables/Table11_Mutation_Impact.csv", row.names = FALSE)

# ============================================
# 6. FIGURE 11: MUTATION IMPACT BY GENE
# ============================================

cat("Creating Figure 11: Mutation Impact Classification...\n")

# Get top 10 mutated genes
top_genes <- impact_summary %>%
  group_by(Hugo_Symbol) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  head(10) %>%
  pull(Hugo_Symbol)

# Filter for top genes
impact_plot_data <- impact_summary %>%
  filter(Hugo_Symbol %in% top_genes) %>%
  mutate(
    Hugo_Symbol = factor(Hugo_Symbol, levels = rev(top_genes)),
    impact = factor(impact, levels = c(
      "High (Loss of function)",
      "Moderate (Altered function)",
      "Moderate (Splicing)",
      "Low (Silent/Non-coding)",
      "Unknown"
    ))
  )

p_impact <- ggplot(impact_plot_data, 
                   aes(x = Hugo_Symbol, y = n, fill = impact)) +
  geom_col(position = "stack", color = "black", size = 0.3) +
  coord_flip() +
  labs(
    title = "Functional Impact of Mutations by Gene",
    subtitle = "Top 10 mutated pathway genes",
    x = "Gene",
    y = "Number of Mutations",
    fill = "Predicted Impact"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_manual(values = c(
    "High (Loss of function)" = "#D73027",
    "Moderate (Altered function)" = "#FC8D59",
    "Moderate (Splicing)" = "#FEE090",
    "Low (Silent/Non-coding)" = "#E0F3F8",
    "Unknown" = "#CCCCCC"
  ))

ggsave("results/figures/Figure11_Mutation_Impact.png", 
       p_impact, width = 12, height = 8, dpi = 300)

cat("✓ Figure 11 saved\n\n")

# ============================================
# 7. VARIANT ALLELE FREQUENCY (VAF) ANALYSIS
# ============================================

cat("Analyzing Variant Allele Frequency (VAF)...\n")

# Check if VAF data is available
if("t_alt_count" %in% colnames(maf_data) && "t_depth" %in% colnames(maf_data)) {
  
  vaf_data <- maf_data %>%
    filter(Hugo_Symbol %in% pathway_genes$Gene) %>%
    filter(t_depth > 0) %>%
    mutate(
      vaf = t_alt_count / t_depth,
      vaf_category = case_when(
        vaf < 0.2 ~ "Low (<20%)",
        vaf < 0.4 ~ "Moderate (20-40%)",
        vaf >= 0.4 ~ "High (≥40%)"
      )
    )
  
  cat(paste("✓ VAF calculated for", nrow(vaf_data), "mutations\n"))
  cat(paste("  Mean VAF:", round(mean(vaf_data$vaf, na.rm = TRUE) * 100, 1), "%\n"))
  cat(paste("  Median VAF:", round(median(vaf_data$vaf, na.rm = TRUE) * 100, 1), "%\n\n"))
  
  # ============================================
  # 8. FIGURE 12: VAF DISTRIBUTION
  # ============================================
  
  cat("Creating Figure 12: VAF Distribution...\n")
  
  p_vaf <- ggplot(vaf_data, aes(x = vaf * 100, fill = vaf_category)) +
    geom_histogram(bins = 30, color = "black", alpha = 0.8) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 55, y = Inf, 
             label = "50% (expected clonal)", 
             vjust = 2, color = "red", fontface = "bold") +
    labs(
      title = "Variant Allele Frequency Distribution",
      subtitle = "Pathway gene mutations in GBM",
      x = "Variant Allele Frequency (%)",
      y = "Number of Mutations",
      fill = "VAF Category"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    ) +
    scale_fill_manual(values = c("#FEE5D9", "#FCAE91", "#FB6A4A"))
  
  ggsave("results/figures/Figure12_VAF_Distribution.png", 
         p_vaf, width = 10, height = 7, dpi = 300)
  
  cat("✓ Figure 12 saved\n\n")
  
  # Save VAF data
  vaf_summary <- vaf_data %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      n_mutations = n(),
      mean_vaf = round(mean(vaf, na.rm = TRUE), 3),
      median_vaf = round(median(vaf, na.rm = TRUE), 3),
      sd_vaf = round(sd(vaf, na.rm = TRUE), 3)
    ) %>%
    arrange(desc(n_mutations))
  
  write.csv(vaf_summary, "results/tables/Table12_VAF_Summary.csv", row.names = FALSE)
  
} else {
  cat("⚠ VAF data not available in MAF file\n")
  cat("Creating placeholder Figure 12...\n\n")
  
  # Create a placeholder plot
  p_placeholder <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "VAF data not available\nin this dataset", 
             size = 8, fontface = "bold") +
    theme_void() +
    labs(title = "Variant Allele Frequency Distribution",
         subtitle = "Data not available in TCGA MAF file")
  
  ggsave("results/figures/Figure12_VAF_Distribution.png", 
         p_placeholder, width = 10, height = 7, dpi = 300)
  
  cat("✓ Figure 12 (placeholder) saved\n\n")
}

# ============================================
# 9. GENE-SPECIFIC MUTATION PATTERNS
# ============================================

cat("Analyzing gene-specific mutation patterns...\n")

# Focus on key genes
key_genes <- c("TP53", "PTEN", "EGFR", "NF1", "PIK3CA", "PIK3R1", "RB1")

gene_patterns <- maf_data %>%
  filter(Hugo_Symbol %in% key_genes) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    total_mutations = n(),
    unique_patients = n_distinct(Tumor_Sample_Barcode),
    missense = sum(Variant_Classification == "Missense_Mutation"),
    nonsense = sum(Variant_Classification == "Nonsense_Mutation"),
    frameshift = sum(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins")),
    splice = sum(Variant_Classification %in% c("Splice_Site", "Splice_Region")),
    inframe = sum(Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins"))
  ) %>%
  mutate(
    pct_missense = round(100 * missense / total_mutations, 1),
    pct_lof = round(100 * (nonsense + frameshift) / total_mutations, 1)
  ) %>%
  arrange(desc(unique_patients))

cat("\nGene-specific mutation patterns:\n")
cat(sprintf("%-10s | Patients | Total Muts | Missense | LOF\n", "Gene"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for(i in 1:nrow(gene_patterns)) {
  cat(sprintf("%-10s | %8d | %10d | %7.1f%% | %6.1f%%\n",
              gene_patterns$Hugo_Symbol[i],
              gene_patterns$unique_patients[i],
              gene_patterns$total_mutations[i],
              gene_patterns$pct_missense[i],
              gene_patterns$pct_lof[i]))
}
cat("\n")

write.csv(gene_patterns, "results/tables/Table13_Gene_Patterns.csv", row.names = FALSE)

# ============================================
# 10. MUTATION CLUSTERING ANALYSIS
# ============================================

cat("Analyzing mutation clustering patterns...\n")

# For genes with position data, check for clustering
genes_with_positions <- maf_data %>%
  filter(Hugo_Symbol %in% key_genes) %>%
  filter(!is.na(Protein_position)) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    n_mutations = n(),
    n_unique_positions = n_distinct(Protein_position),
    clustering_ratio = n_mutations / n_unique_positions
  ) %>%
  arrange(desc(clustering_ratio))

cat("\nMutation clustering analysis:\n")
cat("(Ratio > 1 indicates mutations cluster at specific positions)\n\n")
for(i in 1:nrow(genes_with_positions)) {
  cat(sprintf("  %s: %.2f mutations per position\n",
              genes_with_positions$Hugo_Symbol[i],
              genes_with_positions$clustering_ratio[i]))
}
cat("\n")

write.csv(genes_with_positions, "results/tables/Table14_Mutation_Clustering.csv", 
          row.names = FALSE)

# ============================================
# 11. SUMMARY STATISTICS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  MUTATION PATTERN SUMMARY                  ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Hotspot Mutations:\n")
if(nrow(hotspots) > 0) {
  cat(sprintf("  • Total recurrent mutations: %d\n", nrow(hotspots)))
  cat(sprintf("  • Most frequent: %s %s (%d patients)\n",
              hotspots$Hugo_Symbol[1],
              hotspots$HGVSp_Short[1],
              hotspots$n[1]))
} else {
  cat("  • No recurrent mutations found\n")
}
cat("\n")

cat("Functional Impact:\n")
cat(sprintf("  • High impact (LOF): %.1f%%\n",
            impact_overall$percentage[impact_overall$impact == "High (Loss of function)"]))
cat(sprintf("  • Moderate impact: %.1f%%\n",
            sum(impact_overall$percentage[grepl("Moderate", impact_overall$impact)])))
cat("\n")

cat("Gene-Specific Patterns:\n")
cat(sprintf("  • TP53: %d patients, %.1f%% missense\n",
            gene_patterns$unique_patients[gene_patterns$Hugo_Symbol == "TP53"],
            gene_patterns$pct_missense[gene_patterns$Hugo_Symbol == "TP53"]))
cat(sprintf("  • PTEN: %d patients, %.1f%% LOF\n",
            gene_patterns$unique_patients[gene_patterns$Hugo_Symbol == "PTEN"],
            gene_patterns$pct_lof[gene_patterns$Hugo_Symbol == "PTEN"]))
cat("\n")

# ============================================
# 12. FINAL SUMMARY
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  ANALYSIS COMPLETE!                        ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Figures created:\n")
cat("  ✓ Figure 9: TP53 lollipop plot\n")
cat("  ✓ Figure 10: PTEN mutation distribution\n")
cat("  ✓ Figure 11: Mutation impact classification\n")
cat("  ✓ Figure 12: VAF distribution\n")
cat("\n")

cat("Tables created:\n")
cat("  ✓ Table 10: Hotspot mutations\n")
cat("  ✓ Table 11: Mutation impact summary\n")
cat("  ✓ Table 12: VAF summary\n")
cat("  ✓ Table 13: Gene-specific patterns\n")
cat("  ✓ Table 14: Mutation clustering\n")
cat("\n")

cat("═══════════════════════════════════════════════\n")
cat("Progress: 12/24 figures complete! (50%)\n")
cat("Next: Script 07 - Therapeutic Vulnerability Mapping\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")