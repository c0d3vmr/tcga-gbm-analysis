# ============================================
# Assignment 2: Preprocess Mutations and CNV
# FIXED VERSION - No Errors
# Date: 2024
# ============================================

library(tidyverse)
library(maftools)

# Create processed data directory
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  Mutation & CNV Preprocessing              ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD DATA
# ============================================

cat("Loading data files...\n")

maf_data <- readRDS("data/raw/gbm_mutations.rds")
cnv_data <- readRDS("data/raw/gbm_cnv.rds")
clinical <- read.csv("data/raw/clinical_data.csv")
pathway_genes <- read.csv("data/reference/pathway_genes.csv")

cat(paste("✓ Loaded", nrow(maf_data), "mutations\n"))
cat(paste("✓ Loaded", nrow(pathway_genes), "pathway genes\n"))
cat(paste("✓ Loaded", nrow(clinical), "patient clinical records\n\n"))

# ============================================
# 2. PREPARE CLINICAL DATA FOR MAFTOOLS
# ============================================

cat("Preparing clinical data for maftools...\n")

# maftools expects a column called "Tumor_Sample_Barcode"
# But TCGA clinical data has "bcr_patient_barcode"
# We need to create the right column name

clinical_for_maf <- clinical %>%
  mutate(
    # Create Tumor_Sample_Barcode from patient barcode
    Tumor_Sample_Barcode = bcr_patient_barcode
  )

cat("✓ Clinical data prepared\n\n")

# ============================================
# 3. CREATE MAF OBJECT (for maftools)
# ============================================

cat("Creating MAF object...\n")

# Extract patient IDs (first 12 characters of barcode)
maf_data$patient_id <- substr(maf_data$Tumor_Sample_Barcode, 1, 12)

# Create maftools object with FIXED clinical data
maf <- read.maf(
  maf = maf_data,
  clinicalData = clinical_for_maf,  # Use the fixed version!
  isTCGA = TRUE
)

cat("✓ MAF object created successfully\n\n")

# Save summary
maf_summary <- getSampleSummary(maf)
write.csv(maf_summary, "data/processed/mutation_summary.csv", row.names = FALSE)

cat("✓ Mutation summary saved\n\n")

# ============================================
# 4. FOCUS ON PATHWAY GENES
# ============================================

cat("Extracting pathway gene mutations...\n")

# Extract mutations in pathway genes only
pathway_mutations <- maf_data %>%
  filter(Hugo_Symbol %in% pathway_genes$Gene) %>%
  filter(Variant_Classification %in% c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site"
  ))

cat(paste("✓ Found", nrow(pathway_mutations), "mutations in pathway genes\n\n"))

# Create binary mutation matrix (patient x gene)
cat("Creating mutation matrix...\n")

mutation_matrix <- pathway_mutations %>%
  distinct(patient_id, Hugo_Symbol) %>%
  mutate(mutated = 1) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = mutated,
    values_fill = 0
  )

write.csv(mutation_matrix, "data/processed/mutation_matrix.csv", row.names = FALSE)

cat(paste("✓ Mutation matrix saved:", nrow(mutation_matrix), "patients x", 
          ncol(mutation_matrix)-1, "genes\n\n"))

# ============================================
# 5. PROCESS CNV DATA (SKIP IF NULL)
# ============================================

cat("Processing CNV data...\n")

# Check if CNV data exists (we skipped it in download)
if(is.null(cnv_data)) {
  cat("⚠ CNV data is NULL (was skipped in download)\n")
  cat("Creating empty placeholder files...\n\n")
  
  # Create empty CNV matrices as placeholders
  cnv_matrix <- data.frame(patient_id = mutation_matrix$patient_id)
  cnv_binary <- data.frame(patient_id = mutation_matrix$patient_id)
  cnv_freq <- data.frame(Gene_Symbol = character(), 
                         cnv_status = character(),
                         n = integer(),
                         frequency = numeric())
  
  write.csv(cnv_matrix, "data/processed/cnv_matrix.csv", row.names = FALSE)
  write.csv(cnv_binary, "data/processed/cnv_binary.csv", row.names = FALSE)
  write.csv(cnv_freq, "data/processed/cnv_frequencies.csv", row.names = FALSE)
  
  cat("✓ Created placeholder CNV files\n\n")
  
} else {
  cat("Processing CNV data...\n")
  
  # Get sample columns (everything except first 3)
  sample_cols <- colnames(cnv_data)[4:ncol(cnv_data)]
  
  # Extract patient IDs from sample columns
  cnv_data$patient_ids <- sapply(sample_cols, function(x) substr(x, 1, 12))
  
  # Filter for pathway genes
  cnv_pathway <- cnv_data %>%
    filter(Gene_Symbol %in% pathway_genes$Gene)
  
  # Reshape to long format
  cnv_long <- cnv_pathway %>%
    select(Gene_Symbol, all_of(sample_cols)) %>%
    pivot_longer(
      cols = -Gene_Symbol,
      names_to = "sample_barcode",
      values_to = "cnv_score"
    ) %>%
    mutate(patient_id = substr(sample_barcode, 1, 12))
  
  # Classify CNV status
  cnv_long <- cnv_long %>%
    mutate(
      cnv_status = case_when(
        cnv_score > 0.3 ~ "Amplification",
        cnv_score < -0.3 ~ "Deletion",
        TRUE ~ "Neutral"
      )
    )
  
  # Create CNV matrix (patient x gene)
  cnv_matrix <- cnv_long %>%
    select(patient_id, Gene_Symbol, cnv_score) %>%
    pivot_wider(
      names_from = Gene_Symbol,
      values_from = cnv_score,
      values_fn = mean  # If multiple measurements, take mean
    )
  
  write.csv(cnv_matrix, "data/processed/cnv_matrix.csv", row.names = FALSE)
  
  # Create binary CNV matrix (amplified/deleted/neutral)
  cnv_binary <- cnv_long %>%
    select(patient_id, Gene_Symbol, cnv_status) %>%
    pivot_wider(
      names_from = Gene_Symbol,
      values_from = cnv_status,
      values_fn = function(x) x[1]  # Take first if multiple
    )
  
  write.csv(cnv_binary, "data/processed/cnv_binary.csv", row.names = FALSE)
  
  # CNV frequency per gene
  cnv_freq <- cnv_long %>%
    filter(cnv_status != "Neutral") %>%
    count(Gene_Symbol, cnv_status) %>%
    group_by(Gene_Symbol) %>%
    mutate(total_samples = length(unique(cnv_long$patient_id))) %>%
    mutate(frequency = n / total_samples) %>%
    arrange(desc(frequency))
  
  write.csv(cnv_freq, "data/processed/cnv_frequencies.csv", row.names = FALSE)
  
  cat("✓ CNV data processed successfully\n\n")
}

# ============================================
# 6. SUMMARY STATISTICS
# ============================================

cat("Calculating summary statistics...\n")

# Mutation frequency per gene
mut_freq <- pathway_mutations %>%
  count(Hugo_Symbol) %>%
  mutate(frequency = n / length(unique(maf_data$patient_id))) %>%
  arrange(desc(frequency))

write.csv(mut_freq, "data/processed/mutation_frequencies.csv", row.names = FALSE)

cat("✓ Summary statistics calculated\n\n")

# ============================================
# 7. DISPLAY TOP RESULTS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  TOP MUTATED GENES IN GBM                  ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

top_mutations <- head(mut_freq, 10)
for(i in 1:nrow(top_mutations)) {
  gene <- top_mutations$Hugo_Symbol[i]
  freq <- round(top_mutations$frequency[i] * 100, 1)
  n <- top_mutations$n[i]
  cat(sprintf("%2d. %-10s: %5.1f%% (%d patients)\n", i, gene, freq, n))
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

cat("Files created in data/processed/:\n")
cat("  ✓ mutation_summary.csv\n")
cat("  ✓ mutation_matrix.csv\n")
cat("  ✓ mutation_frequencies.csv\n")
cat("  ✓ cnv_matrix.csv\n")
cat("  ✓ cnv_binary.csv\n")
cat("  ✓ cnv_frequencies.csv\n")
cat("\n")

cat("Summary:\n")
cat(sprintf("  • Mutation matrix: %d patients × %d genes\n", 
            nrow(mutation_matrix), ncol(mutation_matrix)-1))
cat(sprintf("  • CNV matrix: %d patients\n", nrow(cnv_matrix)))
cat(sprintf("  • Top mutated gene: %s (%.1f%%)\n", 
            mut_freq$Hugo_Symbol[1], mut_freq$frequency[1]*100))
cat("\n")

cat("═══════════════════════════════════════════════\n")
cat("Ready for next step: Expression preprocessing!\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")