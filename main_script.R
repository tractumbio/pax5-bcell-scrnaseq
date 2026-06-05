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
#      All input data files are downloaded automatically from the Zenodo
#      archive (DOI 10.5281/zenodo.20553303) into data/. Includes the raw and
#      annotated Seurat objects, the reference atlas, marker modules, Pax5
#      target sets, combined VDJ contigs, and the MSigDB gene-set zip.
#      See data/README.md for the full manifest.
#
#   c) Directory creation
#      All output directories under output/ are created if they do not exist.
#
# ── STAGE 1: Data processing and annotation ───────────────────────────────────
#
#   Script: output/data_processing/scripts/01_label_transfer_and_annotation.R
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
#   saved to data/bcells_annotated.rds.
#
# ── STAGE 2: UMAP visualisations — Figures 4a-c ──────────────────────────────
#
#   Script: output/fig4/scripts/fig4abc_umap_and_markers.R
#
#   Four UMAP plots are generated and saved as PDFs: raw Seurat clusters,
#   fine-grained labels (FinalLab), merged labels (FinalLab2), and mature B
#   sub-classification (FinalLab3). Cell type ordering in legends follows the
#   developmental sequence: PreproB → ProB → PreB → Immature B → FoB → MzB.
#   A per-cell CSV is exported with UMAP coordinates and all metadata columns.
#
# ── STAGE 3: Pax5 differential expression — Figure 4d ────────────────────────
#
#   Script: output/fig4/scripts/fig4d_pax5_diffexp.R
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
#   Script: output/fig4/scripts/fig4e_pax5_target_enrichment.R
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

# ── Package installation and loading ──────────────────────────────────────────

message(strrep("─", 70))
message("Checking and loading packages ...")
source("setup.R")
message(strrep("─", 70))

# ── Download input data from Zenodo ───────────────────────────────────────────

zenodo_doi    <- "10.5281/zenodo.20553303"
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

  # Also download the MSigDB zip if not already present
  msigdb_zip <- "data/msigdb_v2026.1.Mm_files_to_download_locally.zip"
  if (!file.exists(msigdb_zip)) {
    msigdb_files <- files_df[grepl("\\.zip$", files_df$key), ]
    if (nrow(msigdb_files) > 0) {
      message("  Downloading: ", msigdb_files$key[1])
      download.file(msigdb_files$links$self[1], destfile = msigdb_zip,
                    method = "curl", extra = "-L", quiet = FALSE)
    } else {
      message("  Note: MSigDB zip not found in Zenodo record — download manually if needed.")
    }
  } else {
    message("  MSigDB zip already present — skipping.")
  }

}, error = function(e) {
  stop("Failed to download data from Zenodo: ", e$message,
       "\nPlease download manually from: https://doi.org/", zenodo_doi,
       " and place files in data/")
})

message(strrep("─", 70))

# ── Output directory structure ────────────────────────────────────────────────

dirs <- c(
  "data",
  "output/data_processing/supporting_data",
  "output/fig4/plots", "output/fig4/supporting_data",
  "output/fig5/plots", "output/fig5/supporting_data",
  "output/fig6/plots", "output/fig6/supporting_data"
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

run_script("output/data_processing/scripts/01_label_transfer_and_annotation.R")
run_script("output/fig4/scripts/fig4abc_umap_and_markers.R")
run_script("output/fig4/scripts/fig4d_pax5_diffexp.R")
run_script("output/fig4/scripts/fig4e_pax5_target_enrichment.R")
run_script("output/fig5/scripts/fig5a_proB_deg_heatmaps.R")
run_script("output/fig5/scripts/fig5b_proB_gsea.R")
run_script("output/fig6/scripts/fig6_VDJ_analysis.R")

total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)

message("\n", strrep("═", 70))
message("PIPELINE COMPLETE — total runtime: ", total, " min")
message(strrep("═", 70))
message("\nOutput locations:")
message("  Annotated object   →  data/bcells_annotated.rds")
message("  Data processing    →  output/data_processing/")
message("  Figure 4 (a-e)     →  output/fig4/")
message("  Figure 5 (a-b)     →  output/fig5/")
message("  Figure 6 (VDJ)     →  output/fig6/")
