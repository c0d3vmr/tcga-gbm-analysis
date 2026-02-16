# ============================================
# Alternative CNV Download: GISTIC2 Data
# This uses pre-processed CNV calls (easier!)
# Date: 2024
# ============================================

library(TCGAbiolinks)
library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  Downloading GISTIC2 CNV Data              ║\n")
cat("║  (Pre-processed, smaller files)            ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# Check if already exists
if(file.exists("data/raw/gbm_cnv_gistic.rds")) {
  cat("✓ GISTIC2 CNV data already exists!\n\n")
  cnv_data <- readRDS("data/raw/gbm_cnv_gistic.rds")
  cat(paste("Dimensions:", nrow(cnv_data), "genes x", ncol(cnv_data)-3, "samples\n\n"))
  quit(save = "no")
}

cat("Downloading GISTIC2 copy number data...\n")
cat("This is faster because files are pre-processed!\n\n")

tryCatch({
  
  # Query for GISTIC2 data
  query_gistic <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number Scores",
    workflow.type = "GISTIC - Copy Number Score"
  )
  
  cat("✓ Query successful\n")
  cat(paste("Found", nrow(getResults(query_gistic)), "files\n\n"))
  
  cat("Downloading... (should be faster than gene-level CNV)\n\n")
  
  GDCdownload(
    query = query_gistic,
    directory = "data/raw/"
  )
  
  cat("\n✓ Download complete! Preparing data...\n")
  
  # Prepare
  cnv_gistic <- GDCprepare(query_gistic, directory = "data/raw/")
  
  # Save
  saveRDS(cnv_gistic, "data/raw/gbm_cnv_gistic.rds")
  
  cat("\n")
  cat("╔════════════════════════════════════════════╗\n")
  cat("║  SUCCESS! GISTIC2 CNV Downloaded           ║\n")
  cat("╚════════════════════════════════════════════╝\n")
  cat("\n")
  cat(paste("Dimensions:", nrow(cnv_gistic), "genes x", 
            ncol(cnv_gistic)-3, "samples\n"))
  cat("\nSaved to: data/raw/gbm_cnv_gistic.rds\n\n")
  
  cat("NOTE: GISTIC scores are different from gene-level scores:\n")
  cat("  -2 = Deletion\n")
  cat("  -1 = Loss\n")
  cat("   0 = Neutral\n")
  cat("  +1 = Gain\n")
  cat("  +2 = Amplification\n\n")
  
  cat("NEXT STEP:\n")
  cat("Modify script 02 to use GISTIC scores instead\n\n")
  
}, error = function(e) {
  
  cat("\n✗ GISTIC2 download also failed\n")
  cat("Error:", e$message, "\n\n")
  
  cat("╔════════════════════════════════════════════╗\n")
  cat("║  ALL AUTOMATIC METHODS FAILED              ║\n")
  cat("╚════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("This is a known issue with TCGA CNV data.\n\n")
  
  cat("YOUR OPTIONS:\n")
  cat("─────────────────────────────────────────────\n")
  cat("1. Proceed WITHOUT CNV (recommended)\n")
  cat("   - Your mutation data is excellent\n")
  cat("   - Can still do great analysis with mutations + expression\n")
  cat("   - Skip CNV-related figures\n\n")
  
  cat("2. Manual download from GDC Portal\n")
  cat("   - Visit: https://portal.gdc.cancer.gov/\n")
  cat("   - Search: TCGA-GBM CNV\n")
  cat("   - Download manually\n")
  cat("   - Takes 30+ mins of manual work\n\n")
  
  cat("3. Use cBioPortal data (different source)\n")
  cat("   - I can provide code for this\n")
  cat("   - More reliable but different format\n\n")
  
  cat("RECOMMENDATION: Option 1 (proceed without CNV)\n")
  cat("Your analysis is already strong without it!\n\n")
})

cat("═══════════════════════════════════════════════\n")
cat("Script completed\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")