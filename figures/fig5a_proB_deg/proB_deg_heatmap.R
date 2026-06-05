# Heatmap of top 50 differentially expressed genes in ProB cells (MUT vs WT)
#
# This script reads pre-computed differential expression results from the
# genome-wide negbinom DE analysis in figures/4e_pax5_target_enrichment
# (pax5_target_DE_RNA.csv), filters for ProB cells, and selects the top 50
# genes by adjusted p-value.
#
# The heatmap displays 6 columns — one per biological replicate — where each
# column represents the average normalised expression (SCT data slot) of all
# cells from that mouse within the ProB population:
#   Columns 1-3: WT replicates  (WT1, WT2, WT3)
#   Columns 4-6: MUT replicates (Mut1, Mut2, Mut3)
#
# Expression values are z-scored per row so colour reflects relative
# expression change rather than absolute magnitude.
#
# INPUTS
#   figures/4e_pax5_target_enrichment/pax5_target_DE_RNA.csv  — pre-computed DE
#   data/anselm_labeled.rds                                    — Seurat object
#
# OUTPUTS
#   figures/fig5a_proB_deg/proB_top50_DEG_heatmap.pdf
#   figures/fig5a_proB_deg/proB_top50_DEG_results.csv

# ── Packages ──────────────────────────────────────────────────────────────────

source("setup.R")

# Assumes working directory is the repository root
out_dir <- "figures/fig5a_proB_deg"
de_csv  <- "figures/4e_pax5_target_enrichment/pax5_target_DE_RNA.csv"

# ── Load DE results and select top 50 ProB genes ──────────────────────────────

if (!file.exists(de_csv)) {
  stop("DE results not found at: ", de_csv,
       "\nPlease run figures/4e_pax5_target_enrichment/pax5_target_enrichment.R first.")
}

de_all <- read.csv(de_csv)

top50 <- de_all %>%
  filter(cell_type == "ProB") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  head(50)

if (nrow(top50) == 0) stop("No ProB DE results found in ", de_csv)

message("Top 50 ProB DEGs loaded (", sum(top50$avg_log2FC > 0), " up, ",
        sum(top50$avg_log2FC < 0), " down in MUT)")

write.csv(top50, file.path(out_dir, "proB_top50_DEG_results.csv"), row.names = FALSE)

# ── Load Seurat object and subset to ProB ─────────────────────────────────────

seu <- readRDS("data/anselm_labeled.rds")
Idents(seu) <- "FinalLab3"
proB <- subset(seu, idents = "ProB")
meta <- proB@meta.data

# ── Build replicate-average expression matrix ─────────────────────────────────

sct_data <- GetAssayData(proB, assay = "SCT", layer = "data")

# Only keep top50 genes present in the assay
genes_present <- intersect(top50$gene, rownames(sct_data))
if (length(genes_present) < nrow(top50)) {
  message("Note: ", nrow(top50) - length(genes_present),
          " genes not found in SCT assay and will be excluded.")
}
sct_top <- sct_data[genes_present, , drop = FALSE]

# Replicate order: WT first, then MUT
rep_order <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")

avg_mat <- sapply(rep_order, function(rep) {
  cells <- rownames(meta[meta$Mouse == rep, ])
  cells <- intersect(cells, colnames(sct_top))
  if (length(cells) == 0) return(rep(NA_real_, nrow(sct_top)))
  if (length(cells) == 1) return(as.numeric(sct_top[, cells]))
  rowMeans(as.matrix(sct_top[, cells]))
})
rownames(avg_mat) <- genes_present

# Z-score per row
avg_scaled <- t(scale(t(avg_mat)))
avg_scaled[is.nan(avg_scaled)] <- 0

# ── Annotations ───────────────────────────────────────────────────────────────

col_anno <- data.frame(
  Genotype = c("WT", "WT", "WT", "MUT", "MUT", "MUT"),
  row.names = rep_order
)

row_anno <- data.frame(
  Direction = ifelse(top50$avg_log2FC[match(genes_present, top50$gene)] > 0,
                     "Up in MUT", "Down in MUT"),
  row.names = genes_present
)

anno_colours <- list(
  Genotype  = c(WT = "#2166AC", MUT = "#D6604D"),
  Direction = c("Up in MUT" = "#D6604D", "Down in MUT" = "#2166AC")
)

# ── Plot ──────────────────────────────────────────────────────────────────────

col_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pdf(file.path(out_dir, "proB_top50_DEG_heatmap.pdf"), width = 6, height = 12)
pheatmap(
  avg_scaled,
  color             = col_palette,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  annotation_col    = col_anno,
  annotation_row    = row_anno,
  annotation_colors = anno_colours,
  fontsize_row      = 7,
  fontsize_col      = 10,
  border_color      = NA,
  main              = "ProB top 50 DEGs — MUT vs WT\n(z-scored mean SCT expression per replicate)"
)
dev.off()

message("Done. Outputs saved to ", out_dir)
