# TCGA Glioblastoma Genomic Analysis

**A comprehensive multi-omic analysis identifying therapeutic vulnerabilities in glioblastoma**

## Author

Vishakha Medha Rao,
The Awty International School,
December 2025

## Overview

This project performs integrated genomic analysis of 292 glioblastoma patients from The Cancer Genome Atlas (TCGA) to:
- Characterize the mutation landscape across key signaling pathways
- Identify patient-specific therapeutic vulnerabilities
- Correlate genomic alterations with clinical outcomes
- Generate personalized treatment recommendations

## Key Findings

1. **74.3%** of patients have actionable mutations for targeted therapy
2. **TP53** is the most frequently altered gene (46.4%)
3. Multiple signaling pathways are commonly disrupted:
   - RTK/PI3K pathway: 82.9%
   - TP53 pathway: 42.1%
   - RB/Cell Cycle: 14.7%
4. Age is the strongest independent predictor of survival

## Methods

### Data
- **Source**: TCGA-GBM (The Cancer Genome Atlas)
- **Patients**: 292
- **Data types**: Somatic mutations, clinical outcomes
- **Genes analyzed**: 28 across 3 core pathways

### Analysis Pipeline
1. **Mutation landscape** - MAF file analysis, oncoplot visualization
2. **Pathway analysis** - RTK/PI3K, TP53, RB/Cell Cycle pathway assessment
3. **Mutation patterns** - Hotspot identification, functional impact prediction
4. **Therapeutic mapping** - Drug-gene associations, druggability scoring
5. **Survival analysis** - Kaplan-Meier curves, Cox proportional hazards
6. **Case studies** - Patient-specific treatment recommendations

## Repository Structure

```
TCGA_GBM_Assignment2/
├── scripts/              # 9 R analysis scripts
├── data/
│   ├── raw/             # Downloaded TCGA data
│   ├── processed/       # Cleaned datasets
│   └── reference/       # Pathway gene lists
├── results/
│   ├── figures/         # 26 publication-quality figures
│   └── tables/          # 21 summary tables
└── README.md
```

## Key Figures

- **Figure 1**: Mutation landscape (oncoplot)
- **Figure 5**: Pathway alteration frequencies
- **Figure 13**: Therapeutic vulnerability heatmap
- **Figure 17**: Gene-based survival curves
- **Figure 22**: Patient case studies
- **Figure 23**: Literature validation

## Technologies Used

- **R** (primary analysis)
- **TCGAbiolinks** (data acquisition)
- **maftools** (mutation analysis)
- **survival/survminer** (survival analysis)
- **ggplot2** (visualization)
- **tidyverse** (data manipulation)

## Clinical Implications

This analysis demonstrates that:
1. Most GBM patients harbor actionable mutations
2. Patient-specific treatment strategies can be derived from genomic profiles
3. Combination therapies targeting multiple pathways may be beneficial
4. Despite identifying targets, outcomes remain poor - highlighting need for novel approaches

## Reproducibility

All scripts are documented and can be run sequentially (01 through 09). 
See individual script headers for specific requirements and dependencies.

## Acknowledgments

Data provided by The Cancer Genome Atlas Research Network.

## Citation

- TCGA Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. *Nature* 455, 1061-1068 (2008).
