# setup.R
# Central package installation and loading for all analysis scripts.
# Source this file at the top of any script or via main_script.R.

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "ggrepel",
  "tidyr",
  "stringr",
  "pheatmap",
  "RColorBrewer",
  "jsonlite"
)

bioc_packages <- c(
  "Seurat",
  "SeuratObject",
  "fgsea"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing (CRAN): ", pkg)
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing (Bioconductor): ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(jsonlite)

message("All packages loaded.")
