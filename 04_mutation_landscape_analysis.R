# ============================================
# Script 04: Mutation Landscape Analysis
# Comprehensive analysis of mutation patterns in GBM
# ============================================

library(tidyverse)
library(maftools)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  MUTATION LANDSCAPE ANALYSIS               ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD DATA
# ============================================

cat("Loading data...\n")

# Load mutation data
maf_data <- readRDS("data/raw/gbm_mutations.rds")
clinical <- read.csv("data/raw/clinical_data.csv")
pathway_genes <- read.csv("data/reference/pathway_genes.csv")
mutation_matrix <- read.csv("data/processed/mutation_matrix.csv")

cat(paste("✓ Loaded", nrow(maf_data), "total mutations\n"))
cat(paste("✓ Loaded", nrow(clinical), "patient clinical records\n"))
cat(paste("✓ Loaded", nrow(pathway_genes), "pathway genes\n"))
cat(paste("✓ Loaded mutation matrix:", nrow(mutation_matrix), "patients\n\n"))

# ============================================
# 2. CREATE MAF OBJECT
# ============================================

cat("Creating MAF object for maftools...\n")

# Prepare clinical data with correct column name
clinical_for_maf <- clinical %>%
  mutate(Tumor_Sample_Barcode = bcr_patient_barcode)

# Create MAF object
maf <- read.maf(
  maf = maf_data,
  clinicalData = clinical_for_maf,
  isTCGA = TRUE,
  verbose = FALSE
)

cat("✓ MAF object created\n\n")

# ============================================
# 3. BASIC MUTATION STATISTICS
# ============================================

cat("╔════════════════════════════════════════════╗\n")
cat("║  BASIC MUTATION STATISTICS                 ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# Get summary statistics
maf_summary <- getSampleSummary(maf)
gene_summary <- getGeneSummary(maf)

cat(paste("Total samples:", nrow(maf_summary), "\n"))
cat(paste("Total genes with mutations:", nrow(gene_summary), "\n"))
cat(paste("Median mutations per sample:", 
          median(maf_summary$total), "\n"))
cat(paste("Mean mutations per sample:", 
          round(mean(maf_summary$total), 1), "\n\n"))

# ============================================
# 4. FIGURE 1: ONCOPLOT (TOP MUTATED GENES)
# ============================================

cat("Creating Figure 1: Oncoplot...\n")

# Get top 20 mutated genes overall (not just pathway genes)
top_genes <- gene_summary$Hugo_Symbol[1:20]

png("results/figures/Figure01_Oncoplot.png", 
    width = 14, height = 10, res = 300, units = "in")

oncoplot(
  maf = maf,
  top = 20,
  showTumorSampleBarcodes = FALSE,
  fontSize = 0.8,
  legendFontSize = 10,
  annotationFontSize = 10,
  titleFontSize = 14,
  sortByAnnotation = TRUE,
  keepGeneOrder = FALSE,
  removeNonMutated = TRUE
)

dev.off()

cat("✓ Figure 1 saved: results/figures/Figure01_Oncoplot.png\n\n")

# ============================================
# 5. FIGURE 2: MUTATION TYPE DISTRIBUTION
# ============================================

cat("Creating Figure 2: Mutation Type Distribution...\n")

# Get mutation type counts
mutation_types <- maf_data %>%
  count(Variant_Classification) %>%
  arrange(desc(n)) %>%
  mutate(
    percentage = round(100 * n / sum(n), 1),
    Variant_Classification = factor(Variant_Classification, 
                                    levels = Variant_Classification)
  )

# Create bar plot
p_mut_types <- ggplot(mutation_types, 
                      aes(x = Variant_Classification, y = n, fill = Variant_Classification)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n(", percentage, "%)")), 
            vjust = -0.3, size = 3) +
  labs(
    title = "Distribution of Mutation Types in GBM",
    subtitle = paste("Total mutations:", nrow(maf_data)),
    x = "Mutation Type",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave("results/figures/Figure02_Mutation_Types.png", 
       p_mut_types, width = 12, height = 8, dpi = 300)

cat("✓ Figure 2 saved: results/figures/Figure02_Mutation_Types.png\n\n")

# Save table
write.csv(mutation_types, "results/tables/Table01_Mutation_Types.csv", row.names = FALSE)

# ============================================
# 6. CALCULATE TUMOR MUTATION BURDEN (TMB)
# ============================================

cat("Calculating Tumor Mutation Burden...\n")

# TMB = number of mutations per patient
tmb_data <- maf_summary %>%
  as.data.frame() %>%
  select(Tumor_Sample_Barcode, total) %>%
  rename(patient_barcode = Tumor_Sample_Barcode, tmb = total) %>%
  mutate(
    patient_id = substr(patient_barcode, 1, 12),
    tmb_category = case_when(
      tmb < 25 ~ "Low (<25)",
      tmb < 50 ~ "Medium (25-50)",
      tmb < 100 ~ "High (50-100)",
      TRUE ~ "Very High (>100)"
    ),
    tmb_category = factor(tmb_category, 
                          levels = c("Low (<25)", "Medium (25-50)", 
                                     "High (50-100)", "Very High (>100)"))
  )

cat(paste("✓ TMB calculated for", nrow(tmb_data), "patients\n"))
cat(paste("  Mean TMB:", round(mean(tmb_data$tmb), 1), "\n"))
cat(paste("  Median TMB:", median(tmb_data$tmb), "\n"))
cat(paste("  Range:", min(tmb_data$tmb), "-", max(tmb_data$tmb), "\n\n"))

# Save TMB data
write.csv(tmb_data, "data/processed/tumor_mutation_burden.csv", row.names = FALSE)

# ============================================
# 7. FIGURE 3: CO-MUTATION ANALYSIS
# ============================================

cat("Creating Figure 3: Co-mutation Heatmap...\n")

# Focus on pathway genes for co-mutation analysis
pathway_mut_matrix <- mutation_matrix %>%
  select(-patient_id) %>%
  as.matrix()

# Calculate co-occurrence (genes mutated together)
# For each pair of genes, count patients with both mutated
n_patients <- nrow(mutation_matrix)
comutation_matrix <- matrix(0, 
                            nrow = ncol(pathway_mut_matrix), 
                            ncol = ncol(pathway_mut_matrix))
rownames(comutation_matrix) <- colnames(pathway_mut_matrix)
colnames(comutation_matrix) <- colnames(pathway_mut_matrix)

for(i in 1:ncol(pathway_mut_matrix)) {
  for(j in 1:ncol(pathway_mut_matrix)) {
    # Count patients with both genes mutated
    both_mutated <- sum(pathway_mut_matrix[, i] == 1 & 
                          pathway_mut_matrix[, j] == 1)
    comutation_matrix[i, j] <- both_mutated
  }
}

# Convert to percentage
comutation_pct <- 100 * comutation_matrix / n_patients

# Create heatmap
png("results/figures/Figure03_Comutation_Heatmap.png", 
    width = 12, height = 12, res = 300, units = "in")

pheatmap(
  comutation_pct,
  color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(50),
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Co-mutation Frequency in GBM Pathway Genes\n(% of patients with both genes mutated)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "grey90"
)

dev.off()

cat("✓ Figure 3 saved: results/figures/Figure03_Comutation_Heatmap.png\n\n")

# ============================================
# 8. STATISTICAL TEST FOR CO-OCCURRENCE
# ============================================

cat("Testing for significant co-occurrence patterns...\n")

# Test all pairs with Fisher's exact test
comutation_tests <- data.frame()

for(i in 1:(ncol(pathway_mut_matrix)-1)) {
  for(j in (i+1):ncol(pathway_mut_matrix)) {
    gene1 <- colnames(pathway_mut_matrix)[i]
    gene2 <- colnames(pathway_mut_matrix)[j]
    
    # Create contingency table
    both_mut <- sum(pathway_mut_matrix[, i] == 1 & pathway_mut_matrix[, j] == 1)
    gene1_only <- sum(pathway_mut_matrix[, i] == 1 & pathway_mut_matrix[, j] == 0)
    gene2_only <- sum(pathway_mut_matrix[, i] == 0 & pathway_mut_matrix[, j] == 1)
    neither <- sum(pathway_mut_matrix[, i] == 0 & pathway_mut_matrix[, j] == 0)
    
    contingency <- matrix(c(both_mut, gene1_only, gene2_only, neither), nrow = 2)
    
    # Fisher's exact test
    test <- fisher.test(contingency)
    
    comutation_tests <- rbind(comutation_tests, data.frame(
      Gene1 = gene1,
      Gene2 = gene2,
      Both_Mutated = both_mut,
      Odds_Ratio = test$estimate,
      P_Value = test$p.value
    ))
  }
}

# Adjust for multiple testing
comutation_tests$P_Adjusted <- p.adjust(comutation_tests$P_Value, method = "BH")

# Identify significant co-occurrences
significant_comutations <- comutation_tests %>%
  filter(P_Adjusted < 0.05) %>%
  arrange(P_Adjusted)

cat(paste("✓ Found", nrow(significant_comutations), 
          "significant co-mutation pairs (FDR < 0.05)\n\n"))

if(nrow(significant_comutations) > 0) {
  cat("Top significant co-mutations:\n")
  for(i in 1:min(5, nrow(significant_comutations))) {
    cat(sprintf("  %s + %s: OR=%.2f, p=%.2e\n",
                significant_comutations$Gene1[i],
                significant_comutations$Gene2[i],
                significant_comutations$Odds_Ratio[i],
                significant_comutations$P_Adjusted[i]))
  }
  cat("\n")
}

# Save results
write.csv(comutation_tests, "results/tables/Table02_Comutation_Tests.csv", row.names = FALSE)

# ============================================
# 9. FIGURE 4: TMB DISTRIBUTION
# ============================================

cat("Creating Figure 4: Tumor Mutation Burden Distribution...\n")

# Create histogram with categories
p_tmb <- ggplot(tmb_data, aes(x = tmb, fill = tmb_category)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8) +
  geom_vline(xintercept = median(tmb_data$tmb), 
             linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = median(tmb_data$tmb) + 10, y = Inf, 
           label = paste("Median =", median(tmb_data$tmb)), 
           vjust = 2, color = "red", fontface = "bold") +
  labs(
    title = "Tumor Mutation Burden Distribution in GBM",
    subtitle = paste0("n = ", nrow(tmb_data), " patients"),
    x = "Total Mutations per Patient",
    y = "Number of Patients",
    fill = "TMB Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "YlOrRd")

ggsave("results/figures/Figure04_TMB_Distribution.png", 
       p_tmb, width = 10, height = 7, dpi = 300)

cat("✓ Figure 4 saved: results/figures/Figure04_TMB_Distribution.png\n\n")

# TMB category summary
tmb_summary <- tmb_data %>%
  count(tmb_category) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

write.csv(tmb_summary, "results/tables/Table03_TMB_Summary.csv", row.names = FALSE)

# ============================================
# 10. PATHWAY-SPECIFIC MUTATION FREQUENCIES
# ============================================

cat("Calculating pathway-specific mutation frequencies...\n")

# Merge pathway annotation with mutation matrix
pathway_mutation_freq <- data.frame(
  Gene = colnames(mutation_matrix)[-1],  # Exclude patient_id column
  Frequency = colMeans(mutation_matrix[, -1])
) %>%
  left_join(pathway_genes, by = "Gene") %>%
  mutate(
    Frequency_Pct = round(100 * Frequency, 1),
    N_Patients = round(Frequency * nrow(mutation_matrix))
  ) %>%
  arrange(desc(Frequency))

cat("\nTop 10 mutated pathway genes:\n")
for(i in 1:min(10, nrow(pathway_mutation_freq))) {
  cat(sprintf("  %2d. %-10s: %5.1f%% (%s pathway)\n", 
              i,
              pathway_mutation_freq$Gene[i],
              pathway_mutation_freq$Frequency_Pct[i],
              pathway_mutation_freq$Pathway[i]))
}
cat("\n")

write.csv(pathway_mutation_freq, 
          "results/tables/Table04_Pathway_Gene_Frequencies.csv", 
          row.names = FALSE)

# ============================================
# 11. SUMMARY STATISTICS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  SUMMARY STATISTICS                        ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Mutation Landscape:\n")
cat(sprintf("  • Total patients analyzed: %d\n", nrow(mutation_matrix)))
cat(sprintf("  • Pathway genes with mutations: %d\n", ncol(mutation_matrix) - 1))
cat(sprintf("  • Median TMB: %d mutations/patient\n", median(tmb_data$tmb)))
cat(sprintf("  • Most mutated gene: %s (%.1f%%)\n", 
            pathway_mutation_freq$Gene[1],
            pathway_mutation_freq$Frequency_Pct[1]))
cat("\n")

cat("Mutation Types:\n")
cat(sprintf("  • Most common type: %s (%.1f%%)\n",
            mutation_types$Variant_Classification[1],
            mutation_types$percentage[1]))
cat("\n")

cat("Co-mutation Patterns:\n")
cat(sprintf("  • Significant co-mutations: %d pairs\n", 
            nrow(significant_comutations)))
if(nrow(significant_comutations) > 0) {
  cat(sprintf("  • Strongest association: %s + %s (OR=%.2f)\n",
              significant_comutations$Gene1[1],
              significant_comutations$Gene2[1],
              significant_comutations$Odds_Ratio[1]))
}
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
cat("Figure 1: Oncoplot (top 20 genes)\n")
cat("Figure 2: Mutation type distribution\n")
cat("Figure 3: Co-mutation heatmap\n")
cat("Figure 4: TMB distribution\n")
cat("\n")

cat("Tables created:\n")
cat("Table 1: Mutation type summary\n")
cat("Table 2: Co-mutation statistical tests\n")
cat("Table 3: TMB categories\n")
cat("Table 4: Pathway gene frequencies\n")
cat("\n")

cat("Data files created:\n")
cat("data/processed/tumor_mutation_burden.csv\n")
cat("\n")