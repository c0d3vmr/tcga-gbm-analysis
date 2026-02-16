# ============================================
# Assignment 2: Preprocess Expression Data
# Enhanced version with better error handling
# Date: 2024
# ============================================

library(tidyverse)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  RNA-seq Expression Preprocessing          ║\n")
cat("║  This takes 10-15 minutes                  ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD RNA-SEQ DATA
# ============================================

cat("Loading RNA-seq data...\n")

rnaseq <- readRDS("data/raw/gbm_rnaseq.rds")
pathway_genes <- read.csv("data/reference/pathway_genes.csv")

cat(paste("✓ RNA-seq dimensions:", nrow(rnaseq), "genes x", ncol(rnaseq), "samples\n"))
cat(paste("✓ Loaded", nrow(pathway_genes), "pathway genes\n\n"))

# ============================================
# 2. EXTRACT COUNT MATRIX
# ============================================

cat("Extracting count matrix...\n")

# Get unstranded counts
counts <- assay(rnaseq, "unstranded")

# Get gene information
gene_info <- rowData(rnaseq)

cat(paste("✓ Extracted counts for", nrow(counts), "genes\n"))

# Extract patient IDs from column names (first 12 characters)
cat("Standardizing patient IDs...\n")
colnames(counts) <- substr(colnames(counts), 1, 12)

# Remove duplicate samples (keep first occurrence)
n_before <- ncol(counts)
counts <- counts[, !duplicated(colnames(counts))]
n_after <- ncol(counts)

cat(paste("✓ Removed", n_before - n_after, "duplicate samples\n"))
cat(paste("✓ Final:", n_after, "unique patients\n\n"))

# ============================================
# 3. FILTER LOWLY EXPRESSED GENES
# ============================================

cat("Filtering lowly expressed genes...\n")

# Keep genes with at least 10 counts in at least 10% of samples
min_samples <- ceiling(0.1 * ncol(counts))
keep_genes <- rowSums(counts >= 10) >= min_samples

counts_filtered <- counts[keep_genes, ]

cat(paste("✓ Kept", nrow(counts_filtered), "out of", nrow(counts), "genes\n"))
cat(paste("  (", round(100 * nrow(counts_filtered) / nrow(counts), 1), "% retained)\n\n"))

# ============================================
# 4. NORMALIZE WITH DESeq2 VST
# ============================================

cat("Normalizing with DESeq2 VST...\n")
cat("(This is the slow part - be patient!)\n\n")

# Create DESeq2 object
# Create simple colData (metadata)
coldata <- data.frame(
  patient_id = colnames(counts_filtered),
  condition = rep("GBM", ncol(counts_filtered))
)
rownames(coldata) <- colnames(counts_filtered)

cat("Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ 1  # No design, just normalization
)

cat("Running variance stabilizing transformation...\n")
cat("(This takes ~5-10 minutes - don't worry if it seems stuck!)\n\n")

# Variance stabilizing transformation
vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

cat("✓ VST normalization complete!\n\n")

# ============================================
# 5. EXTRACT PATHWAY GENES
# ============================================

cat("Extracting pathway genes from normalized data...\n")

# Get gene symbols from gene_info
gene_info_filtered <- gene_info[keep_genes, ]

# Create mapping
gene_mapping <- data.frame(
  ensembl_id = rownames(vst_matrix),
  gene_symbol = gene_info_filtered$gene_name
)

cat(paste("✓ Created gene mapping for", nrow(gene_mapping), "genes\n"))

# Filter for pathway genes
pathway_indices <- which(gene_mapping$gene_symbol %in% pathway_genes$Gene)
pathway_expression <- vst_matrix[pathway_indices, ]

# Replace row names with gene symbols
rownames(pathway_expression) <- gene_mapping$gene_symbol[pathway_indices]

# Remove any genes that didn't map properly (NAs)
pathway_expression <- pathway_expression[!is.na(rownames(pathway_expression)), ]

cat(paste("✓ Extracted", nrow(pathway_expression), "pathway genes\n"))

# Check which pathway genes are missing
missing_genes <- setdiff(pathway_genes$Gene, rownames(pathway_expression))
if(length(missing_genes) > 0) {
  cat(paste("\n⚠ Note:", length(missing_genes), "pathway genes not found in expression data:\n"))
  cat(paste("  ", paste(head(missing_genes, 5), collapse=", ")))
  if(length(missing_genes) > 5) cat(", ...")
  cat("\n  (This is normal - some genes may not be expressed in GBM)\n\n")
} else {
  cat("✓ All pathway genes found!\n\n")
}

# ============================================
# 6. TRANSPOSE AND SAVE
# ============================================

cat("Saving processed expression data...\n")

# Transpose so patients are rows, genes are columns
expression_matrix <- t(pathway_expression) %>%
  as.data.frame() %>%
  rownames_to_column("patient_id")

write.csv(expression_matrix, "data/processed/expression_matrix.csv", row.names = FALSE)
cat("✓ Saved pathway gene expression: data/processed/expression_matrix.csv\n")

# Also save full normalized expression for later analyses
full_expression <- vst_matrix %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("patient_id")

write.csv(full_expression, "data/processed/expression_full.csv", row.names = FALSE)
cat("✓ Saved full expression data: data/processed/expression_full.csv\n\n")

# ============================================
# 7. VERIFY OVERLAP WITH MUTATION DATA
# ============================================

cat("Checking overlap with mutation data...\n")

if(file.exists("data/processed/mutation_matrix.csv")) {
  mut_matrix <- read.csv("data/processed/mutation_matrix.csv")
  
  common_patients <- intersect(expression_matrix$patient_id, mut_matrix$patient_id)
  
  cat(paste("✓ Patients in expression data:", nrow(expression_matrix), "\n"))
  cat(paste("✓ Patients in mutation data:", nrow(mut_matrix), "\n"))
  cat(paste("✓ Patients with BOTH datasets:", length(common_patients), "\n"))
  
  overlap_pct <- round(100 * length(common_patients) / nrow(expression_matrix), 1)
  cat(paste("✓ Overlap:", overlap_pct, "%\n\n"))
  
  if(overlap_pct < 70) {
    cat("⚠ Warning: Low overlap between datasets. Check patient ID formats.\n\n")
  }
} else {
  cat("⚠ Mutation data not found yet - run script 02 first\n\n")
}

# ============================================
# 8. QUALITY CONTROL PLOTS
# ============================================

cat("Creating quality control plots...\n")

# PCA
cat("  - Running PCA...\n")
pca <- prcomp(t(vst_matrix), scale. = FALSE)
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  patient_id = colnames(vst_matrix)
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5, size = 2, color = "steelblue") +
  labs(title = "PCA of VST-Normalized Expression",
       subtitle = paste(nrow(pca_data), "GBM samples"),
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "% variance)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

ggsave("results/figures/QC_PCA.png", p_pca, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved PCA plot: results/figures/QC_PCA.png\n")

# Sample correlation heatmap (for pathway genes only - smaller)
cat("  - Creating sample correlation heatmap...\n")

# Use pathway expression for faster correlation
sample_cor <- cor(pathway_expression)

png("results/figures/QC_Sample_Correlation.png", width = 10, height = 10, res = 300, units = "in")
heatmap(sample_cor, 
        main = "Sample-to-Sample Correlation (Pathway Genes)",
        col = colorRampPalette(c("blue", "white", "red"))(50),
        labRow = FALSE, labCol = FALSE)
dev.off()

cat("  ✓ Saved correlation heatmap: results/figures/QC_Sample_Correlation.png\n")

# Gene expression distribution
cat("  - Creating expression distribution plot...\n")

# Sample 1000 genes for faster plotting
set.seed(123)
sample_genes <- sample(1:nrow(vst_matrix), min(1000, nrow(vst_matrix)))
expr_long <- vst_matrix[sample_genes, ] %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "patient", values_to = "expression")

p_dist <- ggplot(expr_long, aes(x = expression)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Distribution of Normalized Expression Values",
       subtitle = "VST-transformed RNA-seq data",
       x = "VST Expression",
       y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

ggsave("results/figures/QC_Expression_Distribution.png", p_dist, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved distribution plot: results/figures/QC_Expression_Distribution.png\n\n")

# ============================================
# 9. SUMMARY STATISTICS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  EXPRESSION SUMMARY                        ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat(sprintf("Total samples: %d\n", ncol(vst_matrix)))
cat(sprintf("Total genes (after filtering): %d\n", nrow(vst_matrix)))
cat(sprintf("Pathway genes extracted: %d\n", nrow(pathway_expression)))
cat("\n")

cat("Expression range (pathway genes):\n")
cat(sprintf("  Min: %.2f\n", min(pathway_expression)))
cat(sprintf("  Mean: %.2f\n", mean(pathway_expression)))
cat(sprintf("  Max: %.2f\n", max(pathway_expression)))
cat("\n")

# Top expressed pathway genes
top_expr <- sort(rowMeans(pathway_expression), decreasing = TRUE)
cat("Top 10 expressed pathway genes:\n")
for(i in 1:min(10, length(top_expr))) {
  cat(sprintf("  %2d. %-10s (mean VST = %.2f)\n", 
              i, names(top_expr)[i], top_expr[i]))
}
cat("\n")

# ============================================
# FINAL SUMMARY
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  PROCESSING COMPLETE!                      ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Files created:\n")
cat("  ✓ data/processed/expression_matrix.csv (pathway genes)\n")
cat("  ✓ data/processed/expression_full.csv (all genes)\n")
cat("  ✓ results/figures/QC_PCA.png\n")
cat("  ✓ results/figures/QC_Sample_Correlation.png\n")
cat("  ✓ results/figures/QC_Expression_Distribution.png\n")
cat("\n")

cat("═══════════════════════════════════════════════\n")
cat("Ready for next step: Integrate multi-omic data!\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")