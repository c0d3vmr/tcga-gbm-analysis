# ============================================
# Script 07: Therapeutic Vulnerability Mapping
# Identify actionable mutations and drug targets
# Date: 2024
# ============================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  THERAPEUTIC VULNERABILITY MAPPING         ║\n")
cat("║  Identifying actionable drug targets       ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD DATA
# ============================================

cat("Loading data...\n")

mutation_matrix <- read.csv("data/processed/mutation_matrix.csv")
pathway_status <- read.csv("data/processed/pathway_status.csv")
maf_data <- readRDS("data/raw/gbm_mutations.rds")

cat(paste("✓ Loaded data for", nrow(mutation_matrix), "patients\n\n"))

# ============================================
# 2. DEFINE DRUG-GENE ASSOCIATIONS
# ============================================

cat("Defining drug-gene associations...\n")

# Create comprehensive drug database for GBM
drug_database <- data.frame(
  Drug = c(
    # EGFR inhibitors
    "Erlotinib", "Gefitinib", "Afatinib", "Lapatinib",
    
    # PI3K/AKT/mTOR inhibitors
    "Everolimus", "Temsirolimus", "Buparlisib", "Alpelisib",
    
    # MDM2 inhibitors
    "Idasanutlin", "Milademetan",
    
    # CDK4/6 inhibitors
    "Palbociclib", "Ribociclib", "Abemaciclib",
    
    # PARP inhibitors
    "Olaparib", "Veliparib",
    
    # Standard chemotherapy
    "Temozolomide", "Bevacizumab",
    
    # MEK inhibitors
    "Trametinib", "Cobimetinib",
    
    # PDGFR inhibitors
    "Imatinib", "Sunitinib"
  ),
  
  Target_Gene = c(
    # EGFR inhibitors
    "EGFR", "EGFR", "EGFR", "EGFR",
    
    # PI3K/mTOR
    "MTOR", "MTOR", "PIK3CA", "PIK3CA",
    
    # MDM2
    "MDM2", "MDM2",
    
    # CDK4/6
    "CDK4", "CDK4", "CDK4",
    
    # PARP
    "ATM", "ATM",
    
    # Chemotherapy
    "General", "General",
    
    # MEK
    "NF1", "NF1",
    
    # PDGFR
    "PDGFRA", "PDGFRA"
  ),
  
  Drug_Class = c(
    rep("EGFR inhibitor", 4),
    rep("PI3K/mTOR inhibitor", 4),
    rep("MDM2 inhibitor", 2),
    rep("CDK4/6 inhibitor", 3),
    rep("PARP inhibitor", 2),
    rep("Chemotherapy", 2),
    rep("MEK inhibitor", 2),
    rep("PDGFR inhibitor", 2)
  ),
  
  FDA_Status = c(
    # EGFR
    "Approved (other cancers)", "Approved (other cancers)", 
    "Approved (other cancers)", "Approved (other cancers)",
    
    # PI3K/mTOR
    "Approved (other cancers)", "Approved (other cancers)", 
    "Clinical trial", "Approved (breast cancer)",
    
    # MDM2
    "Clinical trial", "Clinical trial",
    
    # CDK4/6
    "Approved (breast cancer)", "Approved (breast cancer)", 
    "Approved (breast cancer)",
    
    # PARP
    "Approved (ovarian)", "Clinical trial",
    
    # Chemo
    "Approved (GBM)", "Approved (GBM)",
    
    # MEK
    "Approved (melanoma)", "Approved (melanoma)",
    
    # PDGFR
    "Approved (other cancers)", "Approved (other cancers)"
  ),
  
  Evidence_Level = c(
    # EGFR
    "2B", "2B", "3", "3",
    
    # PI3K/mTOR
    "2B", "2B", "3", "2B",
    
    # MDM2
    "3", "3",
    
    # CDK4/6
    "2B", "2B", "2B",
    
    # PARP
    "3", "3",
    
    # Chemo
    "1A", "1A",
    
    # MEK
    "2B", "3",
    
    # PDGFR
    "2B", "2B"
  ),
  
  stringsAsFactors = FALSE
)

cat(paste("✓ Defined", nrow(drug_database), "drug-gene associations\n\n"))

write.csv(drug_database, "results/tables/Table15_Drug_Database.csv", row.names = FALSE)

# ============================================
# 3. CALCULATE THERAPEUTIC VULNERABILITY
# ============================================

cat("Calculating therapeutic vulnerabilities for each patient...\n")

# Initialize vulnerability matrix
vulnerability_matrix <- mutation_matrix %>%
  select(patient_id)

# Helper function to safely check if gene is mutated
has_gene <- function(gene_name) {
  if(gene_name %in% colnames(mutation_matrix)) {
    return(mutation_matrix[[gene_name]] == 1)
  } else {
    return(rep(FALSE, nrow(mutation_matrix)))
  }
}

# For each drug class, determine vulnerability
# EGFR inhibitors
vulnerability_matrix$EGFR_inhibitor <- ifelse(
  has_gene("EGFR"), "High", "Low"
)

# PI3K/mTOR inhibitors (PTEN loss, PIK3CA/PIK3R1 mutations)
vulnerability_matrix$PI3K_mTOR_inhibitor <- case_when(
  (has_gene("PTEN") | has_gene("PIK3CA") | has_gene("PIK3R1")) ~ "High",
  TRUE ~ "Moderate"
)

# MDM2 inhibitors (MDM2 amp + TP53 WT)
# Note: MDM2 might not be in mutation data (it's usually amplified, not mutated)
vulnerability_matrix$MDM2_inhibitor <- case_when(
  has_gene("MDM2") & !has_gene("TP53") ~ "High",
  !has_gene("TP53") ~ "Moderate",
  TRUE ~ "Low"
)

# CDK4/6 inhibitors (CDKN2A/B loss, CDK4/6 alterations, RB1 intact)
vulnerability_matrix$CDK4_6_inhibitor <- case_when(
  (has_gene("CDKN2A") | has_gene("CDKN2B") | 
     has_gene("CDK4") | has_gene("CDK6")) & !has_gene("RB1") ~ "High",
  !has_gene("RB1") ~ "Moderate",
  TRUE ~ "Low"
)

# PARP inhibitors (ATM mutations, DNA repair defects)
vulnerability_matrix$PARP_inhibitor <- ifelse(
  has_gene("ATM"), "Moderate", "Low"
)

# MEK inhibitors (NF1 loss)
vulnerability_matrix$MEK_inhibitor <- ifelse(
  has_gene("NF1"), "Moderate", "Low"
)

# PDGFR inhibitors
vulnerability_matrix$PDGFR_inhibitor <- ifelse(
  has_gene("PDGFRA"), "Moderate", "Low"
)

# Temozolomide (standard)
vulnerability_matrix$Temozolomide <- "Moderate"

cat("✓ Vulnerability calculated for", nrow(vulnerability_matrix), "patients\n\n")

# Save vulnerability matrix
write.csv(vulnerability_matrix, "data/processed/therapeutic_vulnerabilities.csv", 
          row.names = FALSE)

# ============================================
# 4. FIGURE 13: THERAPEUTIC VULNERABILITY HEATMAP
# ============================================

cat("Creating Figure 13: Therapeutic Vulnerability Heatmap...\n")

# Convert to numeric for heatmap
vuln_numeric <- vulnerability_matrix %>%
  select(-patient_id) %>%
  mutate(across(everything(), ~case_when(
    . == "High" ~ 3,
    . == "Moderate" ~ 2,
    . == "Low" ~ 1,
    TRUE ~ 0
  ))) %>%
  as.matrix()

rownames(vuln_numeric) <- vulnerability_matrix$patient_id

# Transpose so drugs are rows, patients are columns
vuln_numeric <- t(vuln_numeric)

# Better row names
rownames(vuln_numeric) <- c(
  "EGFR inhibitors",
  "PI3K/mTOR inhibitors",
  "MDM2 inhibitors",
  "CDK4/6 inhibitors",
  "PARP inhibitors",
  "MEK inhibitors",
  "PDGFR inhibitors",
  "Temozolomide"
)

# Sort patients by total vulnerability score
col_order <- order(colSums(vuln_numeric), decreasing = TRUE)
vuln_numeric <- vuln_numeric[, col_order]

# Create heatmap
png("results/figures/Figure13_Therapeutic_Vulnerabilities.png", 
    width = 14, height = 8, res = 300, units = "in")

pheatmap(
  vuln_numeric,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  breaks = seq(0, 3, length.out = 51),
  main = "Therapeutic Vulnerability Landscape in GBM",
  legend_labels = c("Low", "Moderate", "High"),
  fontsize = 12,
  fontsize_row = 11,
  cellwidth = NA,
  cellheight = 20,
  legend_breaks = c(1, 2, 3),
  border_color = "grey90"
)

dev.off()

cat("✓ Figure 13 saved\n\n")

# ============================================
# 5. CALCULATE DRUGGABILITY SCORES
# ============================================

cat("Calculating patient druggability scores...\n")

# Get only the drug columns (exclude patient_id)
drug_columns <- setdiff(colnames(vulnerability_matrix), "patient_id")

# Count high-priority targets per patient
druggability_scores <- vulnerability_matrix %>%
  mutate(
    n_high = rowSums(select(., all_of(drug_columns)) == "High"),
    n_moderate = rowSums(select(., all_of(drug_columns)) == "Moderate"),
    n_low = rowSums(select(., all_of(drug_columns)) == "Low"),
    druggability_score = n_high * 3 + n_moderate * 2 + n_low * 1,
    druggability_category = case_when(
      n_high >= 2 ~ "Highly druggable",
      n_high == 1 ~ "Moderately druggable",
      n_moderate >= 3 ~ "Moderately druggable",
      TRUE ~ "Limited options"
    )
  )

cat("\nDruggability distribution:\n")
druggability_summary <- druggability_scores %>%
  count(druggability_category) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

for(i in 1:nrow(druggability_summary)) {
  cat(sprintf("  %-25s: %3d patients (%.1f%%)\n",
              druggability_summary$druggability_category[i],
              druggability_summary$n[i],
              druggability_summary$percentage[i]))
}
cat("\n")

write.csv(druggability_scores, "data/processed/druggability_scores.csv", 
          row.names = FALSE)

# ============================================
# 6. FIGURE 14: ACTIONABLE MUTATION FREQUENCIES
# ============================================

cat("Creating Figure 14: Actionable Mutation Frequencies...\n")

# Define actionable genes
actionable_genes <- c("EGFR", "PTEN", "PIK3CA", "PIK3R1", "NF1", 
                      "MDM2", "CDK4", "CDK6", "CDKN2A", "CDKN2B", 
                      "ATM", "PDGFRA")

# Calculate frequencies
actionable_freq <- mutation_matrix %>%
  select(all_of(c("patient_id", intersect(actionable_genes, colnames(mutation_matrix))))) %>%
  pivot_longer(-patient_id, names_to = "Gene", values_to = "Mutated") %>%
  group_by(Gene) %>%
  summarise(
    n_patients = sum(Mutated),
    percentage = 100 * mean(Mutated)
  ) %>%
  arrange(desc(percentage)) %>%
  mutate(Gene = factor(Gene, levels = Gene))

# Add drug class annotation
actionable_freq <- actionable_freq %>%
  mutate(
    drug_class = case_when(
      Gene == "EGFR" ~ "EGFR inhibitors",
      Gene %in% c("PTEN", "PIK3CA", "PIK3R1") ~ "PI3K/mTOR inhibitors",
      Gene == "MDM2" ~ "MDM2 inhibitors",
      Gene %in% c("CDK4", "CDK6", "CDKN2A", "CDKN2B") ~ "CDK4/6 inhibitors",
      Gene == "ATM" ~ "PARP inhibitors",
      Gene == "NF1" ~ "MEK inhibitors",
      Gene == "PDGFRA" ~ "PDGFR inhibitors",
      TRUE ~ "Other"
    )
  )

p_actionable <- ggplot(actionable_freq, 
                       aes(x = Gene, y = percentage, fill = drug_class)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(round(percentage, 1), "%\n(", n_patients, ")")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(
    title = "Frequency of Actionable Mutations in GBM",
    subtitle = "Genes with FDA-approved or clinical trial drug associations",
    x = "Gene",
    y = "Patients with Mutation (%)",
    fill = "Drug Class"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, max(actionable_freq$percentage) * 1.15))

ggsave("results/figures/Figure14_Actionable_Mutations.png", 
       p_actionable, width = 12, height = 8, dpi = 300)

cat("✓ Figure 14 saved\n\n")

write.csv(actionable_freq, "results/tables/Table16_Actionable_Frequencies.csv", 
          row.names = FALSE)

# ============================================
# 7. FIGURE 15: DRUG-PATHWAY NETWORK
# ============================================

cat("Creating Figure 15: Drug-Pathway Association Network...\n")

# Calculate pathway-drug associations
pathway_drug_assoc <- data.frame(
  Pathway = c(
    rep("RTK/PI3K", 4),
    rep("TP53", 2),
    rep("RB/Cell Cycle", 3)
  ),
  Drug_Class = c(
    "EGFR inhibitors",
    "PI3K/mTOR inhibitors",
    "MEK inhibitors",
    "PDGFR inhibitors",
    "MDM2 inhibitors",
    "Temozolomide",
    "CDK4/6 inhibitors",
    "CDK4/6 inhibitors",
    "Temozolomide"
  ),
  N_Patients = c(
    sum(vulnerability_matrix$EGFR_inhibitor == "High"),
    sum(vulnerability_matrix$PI3K_mTOR_inhibitor == "High"),
    sum(vulnerability_matrix$MEK_inhibitor == "Moderate"),
    sum(vulnerability_matrix$PDGFR_inhibitor == "Moderate"),
    sum(vulnerability_matrix$MDM2_inhibitor == "High"),
    nrow(vulnerability_matrix),
    sum(vulnerability_matrix$CDK4_6_inhibitor == "High"),
    sum(vulnerability_matrix$CDK4_6_inhibitor == "Moderate"),
    nrow(vulnerability_matrix)
  )
)

# Create bubble plot
p_network <- ggplot(pathway_drug_assoc, 
                    aes(x = Pathway, y = Drug_Class, size = N_Patients, 
                        color = Pathway)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = N_Patients), color = "white", 
            fontface = "bold", size = 3) +
  labs(
    title = "Drug-Pathway Association Network",
    subtitle = "Number of patients with therapeutic vulnerabilities",
    x = "Signaling Pathway",
    y = "Drug Class",
    size = "Number of\nPatients"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right"
  ) +
  scale_size_continuous(range = c(5, 20)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

ggsave("results/figures/Figure15_Drug_Pathway_Network.png", 
       p_network, width = 12, height = 8, dpi = 300)

cat("✓ Figure 15 saved\n\n")

# ============================================
# 8. FIGURE 16: PATIENT DRUGGABILITY SCORES
# ============================================

cat("Creating Figure 16: Patient Druggability Distribution...\n")

p_druggability <- ggplot(druggability_scores, 
                         aes(x = druggability_score, 
                             fill = druggability_category)) +
  geom_histogram(bins = 15, color = "black", alpha = 0.8) +
  geom_vline(xintercept = median(druggability_scores$druggability_score),
             linetype = "dashed", color = "red", size = 1) +
  annotate("text", 
           x = median(druggability_scores$druggability_score) + 1, 
           y = Inf, 
           label = paste("Median =", median(druggability_scores$druggability_score)),
           vjust = 2, color = "red", fontface = "bold") +
  labs(
    title = "Patient Druggability Score Distribution",
    subtitle = "Based on actionable mutations and therapeutic options",
    x = "Druggability Score\n(High = more treatment options)",
    y = "Number of Patients",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "RdYlGn", direction = 1)

ggsave("results/figures/Figure16_Druggability_Scores.png", 
       p_druggability, width = 10, height = 7, dpi = 300)

cat("✓ Figure 16 saved\n\n")

# ============================================
# 9. COMBINATION THERAPY OPPORTUNITIES
# ============================================

cat("Identifying combination therapy opportunities...\n")

# Find patients with multiple high-priority targets
combination_opportunities <- druggability_scores %>%
  filter(n_high >= 2) %>%
  select(patient_id, n_high, druggability_score)

# Get specific combinations
combo_details <- vulnerability_matrix %>%
  filter(patient_id %in% combination_opportunities$patient_id) %>%
  pivot_longer(-patient_id, names_to = "Drug", values_to = "Vulnerability") %>%
  filter(Vulnerability == "High") %>%
  group_by(patient_id) %>%
  summarise(
    drug_combinations = paste(Drug, collapse = " + "),
    n_drugs = n()
  ) %>%
  arrange(desc(n_drugs))

cat(paste("\n✓ Identified", nrow(combination_opportunities), 
          "patients with ≥2 high-priority targets\n"))

if(nrow(combo_details) > 0) {
  cat("\nTop combination therapy candidates:\n")
  for(i in 1:min(5, nrow(combo_details))) {
    cat(sprintf("  Patient %s: %s\n",
                combo_details$patient_id[i],
                combo_details$drug_combinations[i]))
  }
  cat("\n")
}

write.csv(combo_details, "results/tables/Table17_Combination_Therapy.csv", 
          row.names = FALSE)

# ============================================
# 10. SUMMARY STATISTICS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  THERAPEUTIC VULNERABILITY SUMMARY         ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Actionable Mutations:\n")
cat(sprintf("  • Most frequent: %s (%.1f%%)\n",
            actionable_freq$Gene[1],
            actionable_freq$percentage[1]))
cat(sprintf("  • Patients with ≥1 actionable mutation: %d (%.1f%%)\n",
            sum(druggability_scores$n_high > 0 | druggability_scores$n_moderate > 0),
            100 * mean(druggability_scores$n_high > 0 | druggability_scores$n_moderate > 0)))
cat("\n")

cat("Druggability:\n")
cat(sprintf("  • Highly druggable patients: %d (%.1f%%)\n",
            druggability_summary$n[druggability_summary$druggability_category == "Highly druggable"],
            druggability_summary$percentage[druggability_summary$druggability_category == "Highly druggable"]))
cat(sprintf("  • Median druggability score: %d\n",
            median(druggability_scores$druggability_score)))
cat("\n")

cat("Drug Classes:\n")
drug_class_summary <- vulnerability_matrix %>%
  summarise(across(-patient_id, ~sum(. == "High"))) %>%
  pivot_longer(everything(), names_to = "Drug", values_to = "N_High") %>%
  arrange(desc(N_High))

for(i in 1:min(5, nrow(drug_class_summary))) {
  cat(sprintf("  • %s: %d high-priority patients\n",
              drug_class_summary$Drug[i],
              drug_class_summary$N_High[i]))
}
cat("\n")

cat("Combination Therapy:\n")
cat(sprintf("  • Patients eligible for combination therapy: %d\n",
            nrow(combination_opportunities)))
cat(sprintf("  • Average drugs per combo patient: %.1f\n",
            mean(combo_details$n_drugs)))
cat("\n")

# ============================================
# 11. FINAL SUMMARY
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  ANALYSIS COMPLETE!                        ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Figures created:\n")
cat("  ✓ Figure 13: Therapeutic vulnerability heatmap\n")
cat("  ✓ Figure 14: Actionable mutation frequencies\n")
cat("  ✓ Figure 15: Drug-pathway network\n")
cat("  ✓ Figure 16: Patient druggability scores\n")
cat("\n")

cat("Tables created:\n")
cat("  ✓ Table 15: Drug database\n")
cat("  ✓ Table 16: Actionable mutation frequencies\n")
cat("  ✓ Table 17: Combination therapy opportunities\n")
cat("\n")

cat("Data files created:\n")
cat("  ✓ data/processed/therapeutic_vulnerabilities.csv\n")
cat("  ✓ data/processed/druggability_scores.csv\n")
cat("\n")

cat("═══════════════════════════════════════════════\n")
cat("Progress: 16/24 figures complete! (67%)\n")
cat("Next: Script 08 - Survival Analysis\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")