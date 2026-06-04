# ─────────────────────────────────────────────────────────────────────────────
# main_script.R
# Pax5 single-cell analysis — master run script
#
# This is the single entry point for the full analysis. Running this script
# is sufficient to reproduce all figures and tables from scratch on a new
# machine. It proceeds through the following stages automatically:
#
# ── STAGE 0: Setup ────────────────────────────────────────────────────────────
#
#   a) Package installation
#      All required CRAN and Bioconductor packages are checked and installed
#      if missing. This includes Seurat, fgsea, ggplot2, dplyr, tidyr, and
#      stringr. No manual installation is needed.
#
#   b) Data download
#      All input .rds files are downloaded automatically from a public Google
#      Drive link as a single zip archive and extracted into data/. If the
#      files are already present (e.g. on a second run), the download is
#      skipped. The following files are downloaded:
#        - anselm_demux_filtered.rds   raw query Seurat object
#        - lee_etat_natcom_data.rds     published B-cell reference atlas
#        - marker_modules.rds           mature B-cell marker gene modules
#        - pax5_targets.rds             Pax5-activated and Pax5-repressed gene sets
#      Note: anselm_labeled.rds is produced by Stage 1 and does not need
#      to be downloaded.
#
#   c) Directory creation
#      All output directories are created if they do not already exist.
#
# ── STAGE 1: Data processing and annotation ───────────────────────────────────
#
#   Script: data_processing/label_transfer_and_annotation.R
#
#   Both the reference (Lee et al., Nature Communications) and query datasets
#   are normalised with SCTransform while regressing out cell cycle variation.
#   The two objects are then integrated using Seurat's SCT-aware anchor-based
#   integration to build a shared embedding. Cell type labels from the reference
#   are transferred to query cells via FindTransferAnchors/TransferData, and a
#   majority-vote rule assigns each Seurat cluster its dominant label (FinalLab).
#   Sub-types are collapsed into a cleaner set (FinalLab2), and Mature B clusters
#   are further resolved into functional sub-types (FoB, MzB, GerminalCenter etc.)
#   based on marker gene module scores (FinalLab3). The fully annotated object is
#   saved to data/anselm_labeled.rds.
#
# ── STAGE 2: UMAP visualisations — Figures 4a-c ──────────────────────────────
#
#   Script: figures/4a-c_umap_plots/fig4a-c_plots.R
#
#   Four UMAP plots are generated and saved as PDFs: raw Seurat clusters,
#   fine-grained labels (FinalLab), merged labels (FinalLab2), and mature B
#   sub-classification (FinalLab3). Cell type ordering in legends follows the
#   developmental sequence: PreproB → ProB → PreB → Immature B → FoB → MzB.
#   A per-cell CSV is exported with UMAP coordinates and all metadata columns.
#
# ── STAGE 3: Pax5 differential expression — Figure 4d ────────────────────────
#
#   Script: figures/4d_pax5_diffexp/pax5_DE_by_celltype.R
#
#   Pax5 expression is compared between MUT and WT cells within each B-cell
#   stage using two complementary methods:
#     - SCT/Wilcoxon: run on PrepSCTFindMarkers-corrected normalised data,
#       the Seurat-recommended approach for SCT-normalised objects.
#     - RNA/negative binomial: run on raw UMI counts, modelling count
#       overdispersion explicitly.
#   Results are visualised as lollipop plots (log2FC on x, cell type on y,
#   dot size = % cells expressing Pax5) with adjusted p-values labelled.
#
# ── STAGE 4: Pax5 target gene set enrichment — Figure 4e ─────────────────────
#
#   Script: figures/4e_pax5_target_enrichment/pax5_target_enrichment.R
#
#   Genome-wide differential expression (negative binomial, raw RNA counts)
#   is run for every expressed gene in each B-cell stage (MUT vs WT). Genes
#   are ranked by average log2 fold change and passed to fgsea, which tests
#   whether Pax5-activated and Pax5-repressed gene sets (from pax5_targets.rds)
#   are coordinately enriched at the top or bottom of the ranked list. DE
#   results are cached after the first run so subsequent runs skip the slow
#   testing step. Outputs include NES summary plots, leading edge gene tables,
#   and per-cell-type barplots of leading edge gene fold changes.
#
# ── USAGE ─────────────────────────────────────────────────────────────────────
#
#   From the terminal:
#     Rscript main_script.R
#
#   From an R session:
#     source("main_script.R")
#
#   The script must be run from the repository root directory.
#   Runtime is approximately 45-90 minutes depending on hardware, with the
#   label transfer and genome-wide DE steps being the most computationally
#   intensive.
#
# ─────────────────────────────────────────────────────────────────────────────

# Detect repo root from script location (works with both Rscript and source())
args     <- commandArgs(trailingOnly = FALSE)
file_arg <- args[startsWith(args, "--file=")]
if (length(file_arg)) {
  root_dir <- normalizePath(dirname(sub("--file=", "", file_arg)))
} else {
  root_dir <- getwd()
}
setwd(root_dir)
message("Working directory set to: ", root_dir)

# ── Package installation ───────────────────────────────────────────────────────

cran_packages <- c("dplyr", "ggplot2", "ggrepel", "tidyr", "stringr")
bioc_packages <- c("Seurat", "SeuratObject", "fgsea")

message(strrep("─", 70))
message("Checking required packages ...")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager ...")
  install.packages("BiocManager")
}

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

message("All packages present.")
message(strrep("─", 70))

# ── Download input data from Zenodo ───────────────────────────────────────────

zenodo_doi    <- "10.5281/zenodo.20545876"
zenodo_id     <- sub(".*zenodo\\.", "", zenodo_doi)
zenodo_api    <- paste0("https://zenodo.org/api/records/", zenodo_id)

if (!dir.exists("data")) dir.create("data", recursive = TRUE)

message("Fetching file list from Zenodo (", zenodo_doi, ") ...")

tryCatch({
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

  meta      <- jsonlite::fromJSON(zenodo_api)
  files_df  <- meta$files
  rds_files <- files_df[grepl("\\.rds$", files_df$key), ]

  if (nrow(rds_files) == 0) stop("No .rds files found in Zenodo record.")

  message("Found ", nrow(rds_files), " files — downloading to data/ ...")

  for (i in seq_len(nrow(rds_files))) {
    fname  <- rds_files$key[i]
    furl   <- rds_files$links$self[i]
    fdest  <- file.path("data", fname)
    message("  Downloading: ", fname)
    download.file(furl, destfile = fdest, method = "curl", extra = "-L", quiet = FALSE)
  }

  message("All data files downloaded to data/")

}, error = function(e) {
  stop("Failed to download data from Zenodo: ", e$message,
       "\nPlease download manually from: https://doi.org/", zenodo_doi,
       " and place .rds files in data/")
})

message(strrep("─", 70))

# ── Output directory structure ────────────────────────────────────────────────

dirs <- c(
  "data",
  "plots",
  "figures/4a-c_umap_plots",
  "figures/4d_pax5_diffexp",
  "figures/4e_pax5_target_enrichment"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created directory: ", d)
  }
}

# ── Helper: run a script with timing and error handling ──────────────────────

run_script <- function(path) {
  message("\n", strrep("─", 70))
  message("RUNNING: ", path)
  message(strrep("─", 70))
  t0 <- proc.time()
  tryCatch(
    source(path, chdir = FALSE),
    error = function(e) {
      message("\nERROR in ", path, ":\n", e$message)
      stop(e)
    }
  )
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  message("\nCompleted: ", path, " (", elapsed, " min)")
}

# ── Pipeline ──────────────────────────────────────────────────────────────────

t_start <- proc.time()

run_script("data_processing/label_transfer_and_annotation.R")
run_script("figures/4a-c_umap_plots/fig4a-c_plots.R")
run_script("figures/4d_pax5_diffexp/pax5_DE_by_celltype.R")
run_script("figures/4e_pax5_target_enrichment/pax5_target_enrichment.R")

total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)

message("\n", strrep("═", 70))
message("PIPELINE COMPLETE — total runtime: ", total, " min")
message(strrep("═", 70))
message("\nOutput locations:")
message("  Annotated object   →  data/anselm_labeled.rds")
message("  Figure 4a-c        →  figures/4a-c_umap_plots/")
message("  Figure 4d          →  figures/4d_pax5_diffexp/")
message("  Figure 4e          →  figures/4e_pax5_target_enrichment/")
