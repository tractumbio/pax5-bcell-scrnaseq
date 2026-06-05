# Heatmap of top 60 differentially expressed genes in ProB cells (MUT vs WT)
#
# Reads pre-computed DE results from output/fig4/supporting_data/
# pax5_target_DE_RNA.csv and produces three heatmaps:
#
#   1. proB_heatmap_top30up_top30down.pdf
#      Top 30 by p-value with positive FC + top 30 by p-value with negative FC.
#      Ribosomal protein genes (Rps/Rpl) excluded.
#
#   2. proB_heatmap_top60_pval.pdf
#      Top 60 by p-value regardless of direction.
#      Ribosomal protein genes (Rps/Rpl) excluded.
#
#   3. proB_heatmap_top60_pval_withribo.pdf
#      Top 60 by p-value regardless of direction.
#      Ribosomal protein genes retained.
#
# Each column shows the average normalised SCT expression across cells from
# one biological replicate (mouse), z-scored per row.
#   Columns 1-3: WT replicates  (WT1, WT2, WT3)
#   Columns 4-6: MUT replicates (Mut1, Mut2, Mut3)
#
# INPUTS
#   output/fig4/supporting_data/pax5_target_DE_RNA.csv
#   data/bcells_annotated.rds
#
# OUTPUTS
#   output/fig5/plots/ (PDFs) and output/fig5/supporting_data/ (CSVs): proB_heatmap_top30up_top30down.pdf
#   output/fig5/plots/ (PDFs) and output/fig5/supporting_data/ (CSVs): proB_heatmap_top60_pval.pdf
#   output/fig5/plots/ (PDFs) and output/fig5/supporting_data/ (CSVs): proB_heatmap_top60_pval_withribo.pdf
#   output/fig5/plots/ (PDFs) and output/fig5/supporting_data/ (CSVs): proB_top60_DEG_results.csv

source("setup.R")
library(pheatmap)
library(RColorBrewer)

plots_dir <- "output/fig5/plots"
sdata_dir <- "output/fig5/supporting_data"
de_csv  <- "output/fig4/supporting_data/pax5_target_DE_RNA.csv"

# ── Load DE results ────────────────────────────────────────────────────────────

if (!file.exists(de_csv)) {
  stop("DE results not found at: ", de_csv,
       "\nPlease run output/fig4/scripts/fig4e_pax5_target_enrichment.R first.")
}

de_proB_all <- read.csv(de_csv) %>%
  filter(cell_type == "ProB") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

if (nrow(de_proB_all) == 0) stop("No ProB DE results found in ", de_csv)

# Ribosomal-filtered version
ribo_genes  <- de_proB_all$gene[grepl("^Rps|^Rpl", de_proB_all$gene)]
message("Found ", length(ribo_genes), " ribosomal protein genes (Rps/Rpl).")
de_proB_noribo <- de_proB_all %>% filter(!grepl("^Rps|^Rpl", gene))

# Gene selections
genes_30up_30down <- c(
  de_proB_noribo %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj) %>% head(30) %>% pull(gene),
  de_proB_noribo %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj) %>% head(30) %>% pull(gene)
)
genes_top60_noribo <- de_proB_noribo %>% arrange(p_val_adj) %>% head(60) %>% pull(gene)
genes_top60_withribo <- de_proB_all  %>% arrange(p_val_adj) %>% head(60) %>% pull(gene)

write.csv(de_proB_noribo %>% filter(gene %in% genes_30up_30down),
          file.path(sdata_dir, "proB_top60_DEG_results.csv"), row.names = FALSE)

# ── Load Seurat object and subset to ProB ─────────────────────────────────────

seu  <- readRDS("data/bcells_annotated.rds")
Idents(seu) <- "FinalLab3"
proB <- subset(seu, idents = "ProB")
meta <- proB@meta.data

sct_data  <- GetAssayData(proB, assay = "SCT", layer = "data")
rep_order <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")

col_anno_df <- data.frame(
  Genotype  = c("WT", "WT", "WT", "MUT", "MUT", "MUT"),
  row.names = rep_order
)
anno_colours <- list(
  Genotype  = c(WT = "#2166AC", MUT = "#D6604D"),
  Direction = c("Up in MUT" = "#D6604D", "Down in MUT" = "#2166AC")
)

# ── Helper functions ──────────────────────────────────────────────────────────

make_avg_matrix <- function(genes, de_ref) {
  genes  <- intersect(genes, rownames(sct_data))
  mat    <- sct_data[genes, , drop = FALSE]
  avg    <- sapply(rep_order, function(rep) {
    cells <- intersect(rownames(meta[meta$Mouse == rep, ]), colnames(mat))
    if (length(cells) == 0) return(rep(NA_real_, nrow(mat)))
    if (length(cells) == 1) return(as.numeric(mat[, cells]))
    rowMeans(as.matrix(mat[, cells]))
  })
  rownames(avg) <- genes
  scaled <- t(scale(t(avg)))
  scaled[is.nan(scaled)] <- 0
  direction <- ifelse(de_ref$avg_log2FC[match(genes, de_ref$gene)] > 0,
                      "Up in MUT", "Down in MUT")
  list(mat = scaled, direction = direction, genes = genes)
}

plot_heatmap <- function(res, title, filename) {
  row_anno_df <- data.frame(Direction = res$direction, row.names = res$genes)
  cell_h_pts  <- 10
  plot_h      <- (nrow(res$mat) * cell_h_pts / 72) + 2.5
  options(bitmapType = "quartz")
  pdf(file.path(plots_dir, filename), width = 6, height = plot_h)
  pt = pheatmap(
    res$mat,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    annotation_col    = col_anno_df,
    annotation_row    = row_anno_df,
    annotation_colors = anno_colours,
    cellheight        = cell_h_pts,
    fontsize_row      = 8,
    fontsize_col      = 10,
    border_color      = NA,
    main              = title
  )
  print(pt)
  dev.off()
  message("Saved: ", filename)
}

# ── Plot 1: top 30 up + top 30 down (no ribo) ─────────────────────────────────

message("Plotting heatmap 1: top 30 up + top 30 down (Rps/Rpl excluded) ...")
res1 <- make_avg_matrix(genes_30up_30down, de_proB_noribo)
plot_heatmap(res1,
  title    = "ProB top 30 up + top 30 down DEGs — MUT vs WT\n(Rps/Rpl excluded, z-scored SCT)",
  filename = "proB_heatmap_top30up_top30down.pdf")

# ── Plot 2: top 60 by p-value (no ribo) ──────────────────────────────────────

message("Plotting heatmap 2: top 60 by p-value (Rps/Rpl excluded) ...")
res2 <- make_avg_matrix(genes_top60_noribo, de_proB_noribo)
plot_heatmap(res2,
  title    = "ProB top 60 DEGs by p-value — MUT vs WT\n(Rps/Rpl excluded, z-scored SCT)",
  filename = "proB_heatmap_top60_pval.pdf")

# ── Plot 3: top 60 by p-value (with ribo) ────────────────────────────────────

message("Plotting heatmap 3: top 60 by p-value (Rps/Rpl retained) ...")
res3 <- make_avg_matrix(genes_top60_withribo, de_proB_all)
plot_heatmap(res3,
  title    = "ProB top 60 DEGs by p-value — MUT vs WT\n(z-scored SCT)",
  filename = "proB_heatmap_top60_pval_withribo.pdf")

message("\nDone. Fig5a outputs saved under output/fig5/")






