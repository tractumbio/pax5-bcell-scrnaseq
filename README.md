# Pax5 single-cell and BCR repertoire analysis

Reproducible analysis code accompanying the manuscript on Pax5 function across
B-cell development. This repository regenerates **Figures 4–6** and all
associated supporting data from a single command.

- **Reference atlas:** Lee et al., *Nature Communications* (mouse B-cell atlas)
- **Query data:** experimental scRNA-seq + 10x BCR (VDJ) of **wild-type (WT)**
  vs **Pax5-mutant (MUT)** mouse B cells
- **Data archive (Zenodo):** https://doi.org/10.5281/zenodo.20553303

---

## 1. What this repository produces

| Figure | Analysis | Script |
|--------|----------|--------|
| 4a–c | UMAP of B-cell label annotation | `output/fig4/scripts/fig4abc_umap_and_markers.R` |
| 4d (markers) | Canonical stage-marker dot plot | same script |
| 4d (DE) | Pax5 differential expression (MUT vs WT) | `output/fig4/scripts/fig4d_pax5_diffexp.R` |
| 4e | Pax5 target gene-set enrichment (fgsea) | `output/fig4/scripts/fig4e_pax5_target_enrichment.R` |
| 5a | ProB DEG heatmaps | `output/fig5/scripts/fig5a_proB_deg_heatmaps.R` |
| 5b | Per-stage gene-set enrichment (Hallmark/Reactome/GO/TF) | `output/fig5/scripts/fig5b_proB_gsea.R` |
| 6 | BCR VDJ repertoire analysis | `output/fig6/scripts/fig6_VDJ_analysis.R` |

---

## 2. Repository layout

```
.
├── README.md                 ← this file
├── METHODS.txt               ← full methods text (for the manuscript)
├── session_info.txt          ← exact R + package versions used
├── main_script.R             ← single entry point: runs the whole pipeline
├── setup.R                   ← installs + loads all required R packages
├── LICENSE                   ← MIT (code)
├── CITATION.cff              ← how to cite this repository
│
├── data/                     ← input data (NOT in git; downloaded from Zenodo)
│   └── README.md             ← data manifest + Zenodo link
│
└── output/                   ← all generated results, organised per figure
    ├── data_processing/
    │   ├── scripts/          ← label transfer & annotation pipeline
    │   └── supporting_data/  ← UMAP coordinates + label table
    ├── fig4/
    │   ├── scripts/          ← analysis code for Figure 4
    │   ├── plots/            ← figure PDFs
    │   └── supporting_data/  ← source-data CSVs / cached results
    ├── fig5/  {scripts, plots, supporting_data}
    └── fig6/  {scripts, plots, supporting_data}
```

Each figure folder is self-contained: the **code**, the **plots**, and the
**source data** behind every panel live together.

---

## 3. System requirements

- **Operating system:** macOS, Linux, or Windows
- **R:** version **4.5.x** (developed on 4.5.1; see `session_info.txt`)
- **RAM:** ≥ 16 GB recommended (integration + genome-wide DE are memory-heavy)
- **Disk:** ~5 GB for data + outputs
- **Internet:** required on first run (package install + Zenodo download)

---

## 4. Installing R

**macOS**
1. Install [R from CRAN](https://cran.r-project.org/bin/macosx/) (the `.pkg`).
2. (Optional) Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/).
3. Command-line tools may be needed for some packages: `xcode-select --install`.

**Linux (Ubuntu/Debian)**
```bash
sudo apt update
sudo apt install r-base r-base-dev libcurl4-openssl-dev libssl-dev \
                 libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
                 libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

**Windows**
1. Install [R from CRAN](https://cran.r-project.org/bin/windows/base/).
2. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (matching your R version).
3. (Optional) Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/).

Verify:
```bash
R --version      # should report 4.5.x
```

---

## 5. Reproducing the analysis (one command)

```bash
git clone https://github.com/tractumbio/pax5-bcell-scrnaseq.git
cd pax5-bcell-scrnaseq
Rscript main_script.R
```

`main_script.R` will, in order:
1. **Install** any missing CRAN/Bioconductor packages (via `setup.R`).
2. **Download** all input data from Zenodo into `data/`.
3. **Create** the `output/` directory tree.
4. **Run** every figure script in dependency order and report timing.

> **Runtime:** approximately **45–90 minutes** on a modern laptop. The label
> transfer (Stage 1) and genome-wide negative-binomial DE (Fig 4e) are the
> slowest steps. Fig 4e caches its DE result to
> `output/fig4/supporting_data/pax5_target_DE_cache.rds`; delete this file to
> force a re-run.

### Running a single figure
Each script can be run on its own **from the repository root** (so relative
paths resolve), provided the data are present and Stage 1 / Fig 4e have been
run first where required:
```bash
Rscript -e 'source("output/fig4/scripts/fig4d_pax5_diffexp.R")'
```

---

## 6. Reproducibility notes

- Exact package versions are listed in [`session_info.txt`](session_info.txt).
- For a fully pinned environment, initialise [`renv`](https://rstudio.github.io/renv/)
  before the first run:
  ```r
  install.packages("renv"); renv::init()
  ```
  then install the versions from `session_info.txt`.
- Random seeds: Seurat steps (UMAP, clustering, integration) use Seurat's
  internal defaults; minor numerical differences across platforms/BLAS builds
  are expected but do not change the biological conclusions.

---

## 7. Data availability

All input files are archived on Zenodo and downloaded automatically:

> **DOI: [10.5281/zenodo.20553303](https://doi.org/10.5281/zenodo.20553303)**

See [`data/README.md`](data/README.md) for the file manifest.

---

## 8. Citation

If you use this code, please cite the manuscript and this repository — see
[`CITATION.cff`](CITATION.cff).

## 9. License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). Data on
Zenodo are released under their own (CC-BY) terms.
