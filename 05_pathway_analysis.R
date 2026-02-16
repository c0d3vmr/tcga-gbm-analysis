# ============================================
# Script 05: Pathway Analysis
# Analyze mutation patterns at pathway level
# Date: 2024
# ============================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(gridExtra)
library(grid)

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  PATHWAY ANALYSIS                          ║\n")
cat("║  Analyzing mutation patterns by pathway    ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

# ============================================
# 1. LOAD DATA
# ============================================

cat("Loading data...\n")

mutation_matrix <- read.csv("data/processed/mutation_matrix.csv")
pathway_genes <- read.csv("data/reference/pathway_genes.csv")
clinical <- read.csv("data/raw/clinical_data.csv")
tmb_data <- read.csv("data/processed/tumor_mutation_burden.csv")

cat(paste("✓ Loaded mutation data for", nrow(mutation_matrix), "patients\n"))
cat(paste("✓ Loaded", nrow(pathway_genes), "pathway genes\n\n"))

# ============================================
# 2. CALCULATE PATHWAY ALTERATION STATUS
# ============================================

cat("Calculating pathway alteration status for each patient...\n")

# Get genes by pathway
rtk_genes <- pathway_genes %>% filter(Pathway == "RTK_PI3K") %>% pull(Gene)
tp53_genes <- pathway_genes %>% filter(Pathway == "TP53") %>% pull(Gene)
rb_genes <- pathway_genes %>% filter(Pathway == "RB_CellCycle") %>% pull(Gene)

# Find which genes are in mutation matrix
rtk_genes_avail <- intersect(rtk_genes, colnames(mutation_matrix))
tp53_genes_avail <- intersect(tp53_genes, colnames(mutation_matrix))
rb_genes_avail <- intersect(rb_genes, colnames(mutation_matrix))

cat(paste("  RTK/PI3K pathway:", length(rtk_genes_avail), "genes available\n"))
cat(paste("  TP53 pathway:", length(tp53_genes_avail), "genes available\n"))
cat(paste("  RB/Cell Cycle pathway:", length(rb_genes_avail), "genes available\n\n"))

# For each patient, determine if pathway is altered (≥1 mutation)
pathway_status <- mutation_matrix %>%
  mutate(
    RTK_altered = rowSums(select(., all_of(rtk_genes_avail))) > 0,
    TP53_altered = rowSums(select(., all_of(tp53_genes_avail))) > 0,
    RB_altered = rowSums(select(., all_of(rb_genes_avail))) > 0,
    
    # Count mutations per pathway
    RTK_count = rowSums(select(., all_of(rtk_genes_avail))),
    TP53_count = rowSums(select(., all_of(tp53_genes_avail))),
    RB_count = rowSums(select(., all_of(rb_genes_avail))),
    
    # Total pathway alterations
    total_pathways_altered = RTK_altered + TP53_altered + RB_altered
  )

cat("✓ Pathway status calculated\n\n")

# Save pathway status
write.csv(pathway_status, "data/processed/pathway_status.csv", row.names = FALSE)

# ============================================
# 3. FIGURE 5: PATHWAY ALTERATION FREQUENCIES
# ============================================

cat("Creating Figure 5: Pathway Alteration Frequencies...\n")

# Calculate frequencies
pathway_freq <- data.frame(
  Pathway = c("RTK/PI3K", "TP53", "RB/Cell Cycle"),
  N_Patients = c(
    sum(pathway_status$RTK_altered),
    sum(pathway_status$TP53_altered),
    sum(pathway_status$RB_altered)
  ),
  Percentage = c(
    100 * mean(pathway_status$RTK_altered),
    100 * mean(pathway_status$TP53_altered),
    100 * mean(pathway_status$RB_altered)
  )
) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway))

# Create bar plot
p_pathway_freq <- ggplot(pathway_freq, aes(x = Pathway, y = Percentage, fill = Pathway)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%\n(", N_Patients, " pts)")),
            vjust = -0.3, size = 5, fontface = "bold") +
  labs(
    title = "Pathway Alteration Frequencies in GBM",
    subtitle = paste("n =", nrow(mutation_matrix), "patients"),
    x = "Signaling Pathway",
    y = "Patients with ≥1 Mutation (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

ggsave("results/figures/Figure05_Pathway_Frequencies.png", 
       p_pathway_freq, width = 10, height = 8, dpi = 300)

cat("✓ Figure 5 saved\n\n")

# Save table
write.csv(pathway_freq, "results/tables/Table05_Pathway_Frequencies.csv", row.names = FALSE)

# ============================================
# 4. PATHWAY CO-OCCURRENCE ANALYSIS
# ============================================

cat("Analyzing pathway co-occurrence patterns...\n")

# Count combinations
pathway_combinations <- pathway_status %>%
  count(RTK_altered, TP53_altered, RB_altered) %>%
  mutate(
    pattern = case_when(
      !RTK_altered & !TP53_altered & !RB_altered ~ "None",
      RTK_altered & !TP53_altered & !RB_altered ~ "RTK only",
      !RTK_altered & TP53_altered & !RB_altered ~ "TP53 only",
      !RTK_altered & !TP53_altered & RB_altered ~ "RB only",
      RTK_altered & TP53_altered & !RB_altered ~ "RTK + TP53",
      RTK_altered & !TP53_altered & RB_altered ~ "RTK + RB",
      !RTK_altered & TP53_altered & RB_altered ~ "TP53 + RB",
      RTK_altered & TP53_altered & RB_altered ~ "All three"
    ),
    percentage = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

cat("\nPathway combination patterns:\n")
for(i in 1:nrow(pathway_combinations)) {
  cat(sprintf("  %-15s: %3d patients (%.1f%%)\n",
              pathway_combinations$pattern[i],
              pathway_combinations$n[i],
              pathway_combinations$percentage[i]))
}
cat("\n")

write.csv(pathway_combinations, "results/tables/Table06_Pathway_Combinations.csv", 
          row.names = FALSE)

# ============================================
# 5. STATISTICAL TESTS FOR CO-OCCURRENCE
# ============================================

cat("Testing for significant pathway co-occurrence...\n")

# Fisher's exact tests for pairwise co-occurrence
test_rtk_tp53 <- fisher.test(table(pathway_status$RTK_altered, pathway_status$TP53_altered))
test_rtk_rb <- fisher.test(table(pathway_status$RTK_altered, pathway_status$RB_altered))
test_tp53_rb <- fisher.test(table(pathway_status$TP53_altered, pathway_status$RB_altered))

cooccurrence_tests <- data.frame(
  Comparison = c("RTK vs TP53", "RTK vs RB", "TP53 vs RB"),
  Odds_Ratio = c(test_rtk_tp53$estimate, test_rtk_rb$estimate, test_tp53_rb$estimate),
  P_Value = c(test_rtk_tp53$p.value, test_rtk_rb$p.value, test_tp53_rb$p.value),
  Interpretation = c(
    ifelse(test_rtk_tp53$estimate > 1, "Co-occurrence", "Mutual exclusivity"),
    ifelse(test_rtk_rb$estimate > 1, "Co-occurrence", "Mutual exclusivity"),
    ifelse(test_tp53_rb$estimate > 1, "Co-occurrence", "Mutual exclusivity")
  )
)

cat("\nCo-occurrence test results:\n")
for(i in 1:nrow(cooccurrence_tests)) {
  cat(sprintf("  %s: OR=%.2f, p=%.4f (%s)\n",
              cooccurrence_tests$Comparison[i],
              cooccurrence_tests$Odds_Ratio[i],
              cooccurrence_tests$P_Value[i],
              cooccurrence_tests$Interpretation[i]))
}
cat("\n")

write.csv(cooccurrence_tests, "results/tables/Table07_Cooccurrence_Tests.csv", 
          row.names = FALSE)

# ============================================
# 6. FIGURE 6: PATHWAY CO-OCCURRENCE VENN DIAGRAM
# ============================================

cat("Creating Figure 6: Pathway Co-occurrence Venn Diagram...\n")

# Calculate overlap numbers
n_rtk_only <- sum(pathway_status$RTK_altered & !pathway_status$TP53_altered & !pathway_status$RB_altered)
n_tp53_only <- sum(!pathway_status$RTK_altered & pathway_status$TP53_altered & !pathway_status$RB_altered)
n_rb_only <- sum(!pathway_status$RTK_altered & !pathway_status$TP53_altered & pathway_status$RB_altered)
n_rtk_tp53 <- sum(pathway_status$RTK_altered & pathway_status$TP53_altered & !pathway_status$RB_altered)
n_rtk_rb <- sum(pathway_status$RTK_altered & !pathway_status$TP53_altered & pathway_status$RB_altered)
n_tp53_rb <- sum(!pathway_status$RTK_altered & pathway_status$TP53_altered & pathway_status$RB_altered)
n_all_three <- sum(pathway_status$RTK_altered & pathway_status$TP53_altered & pathway_status$RB_altered)

# Create Venn diagram
png("results/figures/Figure06_Pathway_Venn_Diagram.png", 
    width = 10, height = 10, res = 300, units = "in")

grid.newpage()
draw.triple.venn(
  area1 = sum(pathway_status$RTK_altered),
  area2 = sum(pathway_status$TP53_altered),
  area3 = sum(pathway_status$RB_altered),
  n12 = sum(pathway_status$RTK_altered & pathway_status$TP53_altered),
  n23 = sum(pathway_status$TP53_altered & pathway_status$RB_altered),
  n13 = sum(pathway_status$RTK_altered & pathway_status$RB_altered),
  n123 = n_all_three,
  category = c("RTK/PI3K", "TP53", "RB/Cell Cycle"),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  alpha = 0.3,
  cex = 2.5,
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05, 0.05),
  main = "Pathway Alteration Overlap in GBM"
)

dev.off()

cat("✓ Figure 6 saved\n\n")

# ============================================
# 7. FIGURE 7: PATIENT STRATIFICATION HEATMAP
# ============================================

cat("Creating Figure 7: Patient Stratification Heatmap...\n")

# Create matrix for heatmap (pathways x patients)
heatmap_data <- pathway_status %>%
  select(patient_id, RTK_altered, TP53_altered, RB_altered) %>%
  column_to_rownames("patient_id") %>%
  t() %>%
  as.matrix()

# Convert logical to numeric
heatmap_data <- heatmap_data * 1

# Row names
rownames(heatmap_data) <- c("RTK/PI3K", "TP53", "RB/Cell Cycle")

# Order columns by total pathways altered
col_order <- order(colSums(heatmap_data), decreasing = TRUE)
heatmap_data <- heatmap_data[, col_order]

png("results/figures/Figure07_Patient_Stratification.png", 
    width = 14, height = 6, res = 300, units = "in")

pheatmap(
  heatmap_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = c("white", "red"),
  border_color = NA,
  main = "Patient Stratification by Pathway Alterations",
  legend_labels = c("Wild-type", "Altered"),
  fontsize = 12,
  fontsize_row = 12
)

dev.off()

cat("✓ Figure 7 saved\n\n")

# ============================================
# 8. FIGURE 8: PATHWAY MUTATION BURDEN
# ============================================

cat("Creating Figure 8: Pathway Mutation Burden Distribution...\n")

# Reshape data for plotting
pathway_burden <- pathway_status %>%
  select(patient_id, RTK_count, TP53_count, RB_count) %>%
  pivot_longer(-patient_id, names_to = "Pathway", values_to = "Mutations") %>%
  mutate(
    Pathway = case_when(
      Pathway == "RTK_count" ~ "RTK/PI3K",
      Pathway == "TP53_count" ~ "TP53",
      Pathway == "RB_count" ~ "RB/Cell Cycle"
    ),
    Pathway = factor(Pathway, levels = c("RTK/PI3K", "TP53", "RB/Cell Cycle"))
  )

# Create violin plot
p_burden <- ggplot(pathway_burden, aes(x = Pathway, y = Mutations, fill = Pathway)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(alpha = 0.2, width = 0.2, size = 0.5) +
  labs(
    title = "Pathway Mutation Burden Distribution",
    subtitle = "Number of mutations per pathway per patient",
    x = "Pathway",
    y = "Number of Mutations"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

ggsave("results/figures/Figure08_Pathway_Burden.png", 
       p_burden, width = 10, height = 8, dpi = 300)

cat("✓ Figure 8 saved\n\n")

# ============================================
# 9. PATHWAY BURDEN STATISTICS
# ============================================

cat("Calculating pathway mutation burden statistics...\n\n")

burden_summary <- pathway_burden %>%
  group_by(Pathway) %>%
  summarise(
    Mean = round(mean(Mutations), 2),
    Median = median(Mutations),
    SD = round(sd(Mutations), 2),
    Min = min(Mutations),
    Max = max(Mutations),
    Patients_with_mutations = sum(Mutations > 0),
    Pct_with_mutations = round(100 * sum(Mutations > 0) / n(), 1)
  )

cat("Pathway mutation burden summary:\n")
print(burden_summary)
cat("\n")

write.csv(burden_summary, "results/tables/Table08_Pathway_Burden_Summary.csv", 
          row.names = FALSE)

# ============================================
# 10. PATIENT SUBGROUPS BY PATHWAY PATTERNS
# ============================================

cat("Defining patient subgroups based on pathway patterns...\n")

# Create clinically meaningful subgroups
patient_subgroups <- pathway_status %>%
  mutate(
    subgroup = case_when(
      total_pathways_altered == 0 ~ "No pathway alterations",
      total_pathways_altered == 1 ~ "Single pathway",
      total_pathways_altered == 2 ~ "Two pathways",
      total_pathways_altered == 3 ~ "All pathways",
      TRUE ~ "Other"
    ),
    subgroup = factor(subgroup, levels = c(
      "No pathway alterations",
      "Single pathway",
      "Two pathways",
      "All pathways"
    ))
  )

subgroup_summary <- patient_subgroups %>%
  count(subgroup) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

cat("\nPatient subgroup distribution:\n")
for(i in 1:nrow(subgroup_summary)) {
  cat(sprintf("  %-25s: %3d patients (%.1f%%)\n",
              subgroup_summary$subgroup[i],
              subgroup_summary$n[i],
              subgroup_summary$percentage[i]))
}
cat("\n")

# Save subgroups
write.csv(patient_subgroups, "data/processed/patient_subgroups.csv", row.names = FALSE)
write.csv(subgroup_summary, "results/tables/Table09_Subgroup_Summary.csv", row.names = FALSE)

# ============================================
# 11. SUMMARY STATISTICS
# ============================================

cat("\n")
cat("╔════════════════════════════════════════════╗\n")
cat("║  PATHWAY ANALYSIS SUMMARY                  ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat("\n")

cat("Pathway Alteration Frequencies:\n")
cat(sprintf("  • RTK/PI3K: %.1f%% (%d patients)\n", 
            pathway_freq$Percentage[1], pathway_freq$N_Patients[1]))
cat(sprintf("  • TP53: %.1f%% (%d patients)\n", 
            pathway_freq$Percentage[2], pathway_freq$N_Patients[2]))
cat(sprintf("  • RB/Cell Cycle: %.1f%% (%d patients)\n", 
            pathway_freq$Percentage[3], pathway_freq$N_Patients[3]))
cat("\n")

cat("Pathway Co-occurrence:\n")
cat(sprintf("  • All three pathways: %d patients (%.1f%%)\n",
            n_all_three,
            100 * n_all_three / nrow(pathway_status)))
cat(sprintf("  • Two pathways: %d patients\n",
            sum(pathway_status$total_pathways_altered == 2)))
cat(sprintf("  • Single pathway: %d patients\n",
            sum(pathway_status$total_pathways_altered == 1)))
cat(sprintf("  • No pathways: %d patients\n",
            sum(pathway_status$total_pathways_altered == 0)))
cat("\n")

cat("Average Mutations per Pathway:\n")
cat(sprintf("  • RTK/PI3K: %.1f mutations\n", 
            mean(pathway_status$RTK_count)))
cat(sprintf("  • TP53: %.1f mutations\n", 
            mean(pathway_status$TP53_count)))
cat(sprintf("  • RB/Cell Cycle: %.1f mutations\n", 
            mean(pathway_status$RB_count)))
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
cat("  ✓ Figure 5: Pathway alteration frequencies\n")
cat("  ✓ Figure 6: Pathway co-occurrence Venn diagram\n")
cat("  ✓ Figure 7: Patient stratification heatmap\n")
cat("  ✓ Figure 8: Pathway mutation burden distribution\n")
cat("\n")

cat("Tables created:\n")
cat("  ✓ Table 5: Pathway frequencies\n")
cat("  ✓ Table 6: Pathway combinations\n")
cat("  ✓ Table 7: Co-occurrence tests\n")
cat("  ✓ Table 8: Pathway burden summary\n")
cat("  ✓ Table 9: Patient subgroups\n")
cat("\n")

cat("Data files created:\n")
cat("  ✓ data/processed/pathway_status.csv\n")
cat("  ✓ data/processed/patient_subgroups.csv\n")
cat("\n")

cat("═══════════════════════════════════════════════\n")
cat("Progress: 8/24 figures complete! (33%)\n")
cat("Next: Script 06 - Mutation Signatures & Patterns\n")
cat("═══════════════════════════════════════════════\n")
cat("\n")