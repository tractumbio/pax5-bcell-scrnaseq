# Heatmap of differentially expressed genes in ProB cells (MUT vs WT)
#
# This script reads pre-computed differential expression results from the
# genome-wide negbinom DE analysis in figures/4e_pax5_target_enrichment
# (pax5_target_DE_RNA.csv) and filters for ProB cells.
#
# ALL significant DEGs are shown in the heatmap. A row annotation marks
# the top 50 genes by adjusted p-value. Each column represents the average
# normalised expression (SCT data slot) across cells from one biological
# replicate (mouse), z-scored per row:
#   Columns 1-3: WT replicates  (WT1, WT2, WT3)
#   Columns 4-6: MUT replicates (Mut1, Mut2, Mut3)
#
# INPUTS
#   figures/4e_pax5_target_enrichment/pax5_target_DE_RNA.csv
#   data/anselm_labeled.rds
#
# OUTPUTS
#   figures/fig5a_proB_deg/proB_DEG_heatmap.pdf
#   figures/fig5a_proB_deg/proB_top50_DEG_results.csv

source("setup.R")

out_dir <- "figures/fig5a_proB_deg"
de_csv  <- "figures/4e_pax5_target_enrichment/pax5_target_DE_RNA.csv"

# ── Load DE results ────────────────────────────────────────────────────────────

if (!file.exists(de_csv)) {
  stop("DE results not found at: ", de_csv,
       "\nPlease run figures/4e_pax5_target_enrichment/pax5_target_enrichment.R first.")
}

de_proB <- read.csv(de_csv) %>%
  filter(cell_type == "ProB") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

if (nrow(de_proB) == 0) stop("No ProB DE results found in ", de_csv)

top50_genes <- de_proB$gene[1:min(50, nrow(de_proB))]

write.csv(de_proB %>% filter(gene %in% top50_genes),
          file.path(out_dir, "proB_top50_DEG_results.csv"), row.names = FALSE)

message(nrow(de_proB), " ProB DEGs loaded — top 50 will be annotated.")

# ── Load Seurat object and subset to ProB ─────────────────────────────────────

seu  <- readRDS("data/anselm_labeled.rds")
Idents(seu) <- "FinalLab3"
proB <- subset(seu, idents = "ProB")
meta <- proB@meta.data

# ── Build replicate-average expression matrix ─────────────────────────────────

sct_data <- GetAssayData(proB, assay = "SCT", layer = "data")

genes_use <- intersect(de_proB$gene, rownames(sct_data))
sct_mat   <- sct_data[genes_use, , drop = FALSE]

rep_order <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")

avg_mat <- sapply(rep_order, function(rep) {
  cells <- rownames(meta[meta$Mouse == rep, ])
  cells <- intersect(cells, colnames(sct_mat))
  if (length(cells) == 0) return(rep(NA_real_, nrow(sct_mat)))
  if (length(cells) == 1) return(as.numeric(sct_mat[, cells]))
  rowMeans(as.matrix(sct_mat[, cells]))
})
rownames(avg_mat) <- genes_use

# Z-score per row
avg_scaled <- t(scale(t(avg_mat)))
avg_scaled[is.nan(avg_scaled)] <- 0

# ── ComplexHeatmap annotations ────────────────────────────────────────────────

# Column annotation — genotype
col_anno <- HeatmapAnnotation(
  Genotype = c("WT", "WT", "WT", "MUT", "MUT", "MUT"),
  col      = list(Genotype = c(WT = "#2166AC", MUT = "#D6604D")),
  annotation_name_side = "left"
)

# Row annotation — direction + top50 highlight
is_top50    <- genes_use %in% top50_genes
direction   <- ifelse(de_proB$avg_log2FC[match(genes_use, de_proB$gene)] > 0,
                      "Up in MUT", "Down in MUT")

row_anno <- rowAnnotation(
  Top50 = anno_simple(
    as.integer(is_top50),
    col    = c("0" = "grey90", "1" = "#F1A340"),
    height = unit(3, "mm")
  ),
  Direction = direction,
  col = list(
    Direction = c("Up in MUT" = "#D6604D", "Down in MUT" = "#2166AC")
  ),
  annotation_name_side = "top",
  show_legend = c(Top50 = TRUE, Direction = TRUE)
)

# Colour scale
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#D6604D")
)

# Gene labels: only label top 50
gene_labels <- ifelse(is_top50, genes_use, "")

# Dynamic height: ~2px per gene, min 8 inches
plot_h <- max(8, length(genes_use) * 0.08)

# ── Draw heatmap ──────────────────────────────────────────────────────────────

ht <- Heatmap(
  avg_scaled,
  name               = "z-score",
  col                = col_fun,
  cluster_columns    = FALSE,
  cluster_rows       = TRUE,
  show_row_names     = TRUE,
  row_labels         = gene_labels,
  row_names_gp       = gpar(fontsize = 6),
  show_column_names  = TRUE,
  column_names_gp    = gpar(fontsize = 10),
  top_annotation     = col_anno,
  right_annotation   = row_anno,
  column_title       = "ProB DEGs — MUT vs WT\n(z-scored mean SCT per replicate)",
  column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
  border             = FALSE,
  use_raster         = TRUE,
  raster_quality     = 5
)
ht
pdf(file.path(out_dir, "proB_DEG_heatmap.pdf"), width = 7, height = plot_h)
draw(ht)
dev.off()

message("Done. Outputs saved to ", out_dir)
