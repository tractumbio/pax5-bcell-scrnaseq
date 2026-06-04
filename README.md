# Pax5 Single-Cell Analysis

Reproducible single-cell RNA-seq analysis accompanying the Pax5 manuscript.  
This repository contains all code to reproduce the data processing pipeline and Figures 4a–e.

---

## Overview

We characterise Pax5 function across B-cell development using single-cell RNA-seq.  
A published mouse B-cell atlas (Lee et al., *Nature Communications*) is used as a reference  
to annotate an experimental query dataset comparing **wild-type (WT)** and **mutant (MUT)** animals.  
Differential expression and gene set enrichment are then used to quantify how Pax5 loss  
perturbs its transcriptional programme across discrete developmental stages.

---

## Repository structure

```
.
├── data/                          # Input data (not tracked — see data/README.md)
├── plots/                         # Output directory for figures and tables
│
├── data_processing/
│   └── label_transfer_and_annotation.R   # Full preprocessing and annotation pipeline
│
└── figures/
    ├── 4a-c_umap_plots/
    │   └── fig4a-c_plots.R                # UMAP visualisations of B-cell labels
    ├── 4d_pax5_diffexp/
    │   └── pax5_DE_by_celltype.R          # Pax5 differential expression per cell type
    └── 4e_pax5_target_enrichment/
        └── pax5_target_enrichment.R       # Pax5 target gene set enrichment (fgsea)
```

---

## Reproducing the analysis

All scripts should be run **from the repository root** so that relative paths resolve correctly.

### Step 1 — Data processing and annotation

```r
source("data_processing/label_transfer_and_annotation.R")
```

This script:
- Preprocesses reference and query Seurat objects (cell cycle regression, SCTransform)
- Integrates them using SCT-aware anchor-based integration
- Transfers B-cell labels from the reference via majority-vote assignment
- Sub-classifies Mature B clusters using marker gene module scores
- Saves the fully annotated object to `data/anselm_labeled.rds`

> **Input:** `data/anselm_demux_filtered.rds`, `data/lee_etat_natcom_data.rds`, `data/marker_modules.rds`  
> **Output:** `data/anselm_labeled.rds`, `plots/Fig4_umap_data.csv`

---

### Step 2 — UMAP visualisations (Figures 4a–c)

```r
source("figures/4a-c_umap_plots/fig4a-c_plots.R")
```

Generates UMAP plots coloured by cluster identity, fine-grained labels (FinalLab),  
merged labels (FinalLab2), and mature B sub-classification (FinalLab3).  
Also exports a per-cell CSV with UMAP coordinates and all metadata.

> **Input:** `data/anselm_labeled.rds`  
> **Output:** `figures/4a-c_umap_plots/*.pdf`, `figures/4a-c_umap_plots/fig4_per_cell_metadata_umap.csv`

---

### Step 3 — Pax5 differential expression (Figure 4d)

```r
source("figures/4d_pax5_diffexp/pax5_DE_by_celltype.R")
```

Tests Pax5 expression differences between MUT and WT cells within each B-cell stage  
using two complementary methods:
- **SCT / Wilcoxon** — on PrepSCTFindMarkers-corrected normalised data
- **RNA / negative binomial** — on raw UMI counts

Results are visualised as lollipop plots with dot size encoding % cells expressing Pax5.

> **Input:** `data/anselm_labeled.rds`  
> **Output:** `figures/4d_pax5_diffexp/*.csv`, `figures/4d_pax5_diffexp/*.pdf`

---

### Step 4 — Pax5 target gene set enrichment (Figure 4e)

```r
source("figures/4e_pax5_target_enrichment/pax5_target_enrichment.R")
```

Runs genome-wide differential expression (negative binomial, RNA counts) per cell type,  
then uses **fgsea** to test whether Pax5-activated and Pax5-repressed gene sets are  
coordinately dysregulated in MUT cells. Genes are ranked by average log2 fold change.  
DE results are cached on first run — subsequent runs skip the DE step automatically.

> **Input:** `data/anselm_labeled.rds`, `data/pax5_targets.rds`  
> **Output:** `figures/4e_pax5_target_enrichment/*.csv`, `figures/4e_pax5_target_enrichment/*.pdf`

---

## Dependencies

All required packages are checked and installed automatically when each script is run.  
The analysis was developed with the following versions:

| Package | Source |
|---------|--------|
| `Seurat` | Bioconductor / CRAN |
| `fgsea` | Bioconductor |
| `dplyr` | CRAN |
| `ggplot2` | CRAN |
| `stringr` | CRAN |

R version ≥ 4.2 is recommended.

---

## Data availability

Raw and processed data are available from the authors upon request.  
See [`data/README.md`](data/README.md) for a description of each required file.

---

## Citation

> manuscript in preparation
