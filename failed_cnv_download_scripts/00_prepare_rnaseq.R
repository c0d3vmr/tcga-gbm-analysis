# ============================================
# Prepare Already-Downloaded RNA-seq Files
# This combines the 411 individual files into one dataset
# Date: 2024
# ============================================

library(TCGAbiolinks)
library(SummarizedExperiment)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  Preparing RNA-seq Data                    ║\n")
cat("║  (Combining 411 downloaded files)          ║\n")
cat("║  This takes 5-10 minutes                   ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. CHECK IF ALREADY PREPARED
# ============================================

if(file.exists("data/raw/gbm_rnaseq.rds")) {
  cat("✓ RNA-seq .rds file already exists!\n")
  
  # Verify it's valid
  tryCatch({
    rnaseq <- readRDS("data/raw/gbm_rnaseq.rds")
    cat(paste("  Dimensions:", nrow(rnaseq), "genes x", ncol(rnaseq), "samples\n"))
    cat("\n✓ File is valid. Nothing to do!\n\n")
    quit(save = "no")
  }, error = function(e) {
    cat("✗ File exists but is corrupted. Will re-create it.\n\n")
    file.remove("data/raw/gbm_rnaseq.rds")
  })
}

# ============================================
# 2. VERIFY RAW FILES EXIST
# ============================================

cat("Checking for downloaded RNA-seq files...\n")

rnaseq_files <- list.files(
  "data/raw/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification",
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = FALSE
)

cat(paste("✓ Found", length(rnaseq_files), "RNA-seq files\n\n"))

if(length(rnaseq_files) == 0) {
  cat("✗ ERROR: No .tsv files found!\n")
  cat("The files should be in:\n")
  cat("  data/raw/TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/\n\n")
  stop("No RNA-seq files to process")
}

# ============================================
# 3. CREATE QUERY OBJECT (for GDCprepare)
# ============================================

cat("Creating query object...\n")

# We need to create a query object even though files are already downloaded
# This tells GDCprepare where to find the files
query_rnaseq <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

cat("✓ Query created\n\n")

# ============================================
# 4. PREPARE (COMBINE) THE FILES
# ============================================

cat("╔════════════════════════════════════════════╗\n")
cat("║  COMBINING FILES INTO ONE DATASET          ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("This will:\n")
cat("  1. Read all 411 individual .tsv files\n")
cat("  2. Combine them into one matrix\n")
cat("  3. Add gene annotations\n")
cat("  4. Save as .rds file\n\n")

cat("Starting... (this takes 5-10 minutes)\n")
cat("You'll see progress updates below:\n\n")

# This is the key step - it combines all files
tryCatch({
  
  rnaseq_data <- GDCprepare(
    query = query_rnaseq,
    directory = "data/raw/",
    summarizedExperiment = TRUE
  )
  
  cat("\n✓ Files combined successfully!\n\n")
  
  # Show dimensions
  cat("Dataset information:\n")
  cat(paste("  Genes:", nrow(rnaseq_data), "\n"))
  cat(paste("  Samples:", ncol(rnaseq_data), "\n"))
  cat(paste("  Matrix size:", nrow(rnaseq_data) * ncol(rnaseq_data), "values\n\n"))
  
  # ============================================
  # 5. SAVE THE COMBINED DATA
  # ============================================
  
  cat("Saving combined dataset...\n")
  
  saveRDS(rnaseq_data, "data/raw/gbm_rnaseq.rds")
  
  file_size <- file.size("data/raw/gbm_rnaseq.rds") / 1024^2
  cat(paste("✓ Saved:", round(file_size, 1), "MB\n\n"))
  
  # ============================================
  # 6. VERIFY THE SAVED FILE
  # ============================================
  
  cat("Verifying saved file...\n")
  
  # Try to load it back
  test_load <- readRDS("data/raw/gbm_rnaseq.rds")
  
  cat("✓ File loads correctly\n")
  cat(paste("✓ Contains", nrow(test_load), "genes x", ncol(test_load), "samples\n\n"))
  
  # ============================================
  # SUCCESS!
  # ============================================
  
  cat("\n")
  cat("╔════════════════════════════════════════════╗\n")
  cat("║  SUCCESS! RNA-seq Data Ready               ║\n")
  cat("╚════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Created file: data/raw/gbm_rnaseq.rds\n")
  cat(paste("  Size:", round(file_size, 1), "MB\n"))
  cat(paste("  Genes:", nrow(rnaseq_data), "\n"))
  cat(paste("  Samples:", ncol(rnaseq_data), "\n\n"))
  
  cat("═══════════════════════════════════════════════\n")
  cat("NEXT STEP:\n")
  cat("Run script: 03_preprocess_expression.R\n")
  cat("═══════════════════════════════════════════════\n")
  cat("\n")
  
}, error = function(e) {
  
  cat("\n")
  cat("╔════════════════════════════════════════════╗\n")
  cat("║  ERROR DURING FILE PREPARATION             ║\n")
  cat("╚════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("Error message:", e$message, "\n\n")
  
  cat("TROUBLESHOOTING:\n")
  cat("─────────────────────────────────────────────\n")
  cat("1. Check that you have enough RAM (need ~4-8 GB)\n")
  cat("2. Close other programs to free up memory\n")
  cat("3. Restart RStudio and try again\n\n")
  
  cat("If error persists:\n")
  cat("  The individual files might be corrupted\n")
  cat("  You may need to re-download RNA-seq data\n\n")
  
})

cat("═══════════════════════════════════════════════\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")