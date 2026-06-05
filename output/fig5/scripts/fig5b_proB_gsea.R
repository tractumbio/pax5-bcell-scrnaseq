# Gene set enrichment analysis across B cell developmental stages (MUT vs WT)
#
# Uses pre-computed genome-wide negbinom DE results from
# output/fig4/supporting_data/pax5_target_DE_RNA.csv. Genes are ranked
# by average log2FC (MUT vs WT) per cell type. fgsea is run against four
# MSigDB v2026.1 mouse gene set collections (gene symbols):
#
#   - Hallmark       (mh.all)         — 50 curated hallmark gene sets
#   - Reactome       (m2.cp.reactome) — curated pathway gene sets
#   - GO Biological  (m5.go.bp)       — GO biological process terms
#   - TF targets     (m3.gtrd)        — transcription factor binding site targets
#
# Each collection is tested twice:
#   - All genes (full ranked list)
#   - Ribosomal protein genes removed (Rps/Rpl excluded from ranked list)
#
# Results are visualised as lollipop plots (NES on x, pathway on y, colour by
# direction, dot size by -log10 FDR) — one plot per cell type per collection.
#
# INPUTS
#   output/fig4/supporting_data/pax5_target_DE_RNA.csv
#   data/msigdb_v2026.1.Mm_files_to_download_locally.zip
#
# OUTPUTS
#   output/fig5/supporting_data/  — fgsea CSVs per cell type per collection
#   output/fig5/plots/            — lollipop PDFs per cell type per collection,
#                                   plus {Stage}_top60_heatmap.pdf per B-cell stage

source("setup.R")

sdata_dir <- "output/fig5/supporting_data"
plots_dir <- "output/fig5/plots"
results_dir <- sdata_dir
dir.create(results_dir, showWarnings = FALSE)
dir.create(plots_dir,   showWarnings = FALSE)

de_csv   <- "output/fig4/supporting_data/pax5_target_DE_RNA.csv"
zip_file <- "data/msigdb_v2026.1.Mm_files_to_download_locally.zip"
gmt_dir  <- "data/msigdb_v2026.1.Mm_files_to_download_locally/msigdb_v2026.1.Mm_GMTs"

cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")

# ── Extract GMT files if needed ───────────────────────────────────────────────

gmt_files <- list(
  hallmark = file.path(gmt_dir, "mh.all.v2026.1.Mm.symbols.gmt"),
  reactome = file.path(gmt_dir, "m2.cp.reactome.v2026.1.Mm.symbols.gmt"),
  gobp     = file.path(gmt_dir, "m5.go.bp.v2026.1.Mm.symbols.gmt"),
  tf       = file.path(gmt_dir, "m3.gtrd.v2026.1.Mm.symbols.gmt")
)

needed <- names(gmt_files)[!file.exists(unlist(gmt_files))]
if (length(needed) > 0) {
  message("Extracting missing GMT files ...")
  zip_paths <- c(
    hallmark = "msigdb_v2026.1.Mm_files_to_download_locally/msigdb_v2026.1.Mm_GMTs/mh.all.v2026.1.Mm.symbols.gmt",
    reactome = "msigdb_v2026.1.Mm_files_to_download_locally/msigdb_v2026.1.Mm_GMTs/m2.cp.reactome.v2026.1.Mm.symbols.gmt",
    gobp     = "msigdb_v2026.1.Mm_files_to_download_locally/msigdb_v2026.1.Mm_GMTs/m5.go.bp.v2026.1.Mm.symbols.gmt",
    tf       = "msigdb_v2026.1.Mm_files_to_download_locally/msigdb_v2026.1.Mm_GMTs/m3.gtrd.v2026.1.Mm.symbols.gmt"
  )
  unzip(zip_file, files = zip_paths[needed], exdir = "data")
}
message("GMT files ready.")

# ── Load gene sets ────────────────────────────────────────────────────────────

message("Loading gene sets ...")
pathways <- lapply(gmt_files, gmtPathways)
message("  Hallmark: ", length(pathways$hallmark), " sets")
message("  Reactome: ", length(pathways$reactome), " sets")
message("  GO BP:    ", length(pathways$gobp),     " sets")
message("  TF:       ", length(pathways$tf),       " sets")

# ── Load DE results ───────────────────────────────────────────────────────────

if (!file.exists(de_csv)) {
  stop("DE results not found at: ", de_csv,
       "\nPlease run output/fig4/scripts/fig4e_pax5_target_enrichment.R first.")
}

de_all <- read.csv(de_csv)

# ── Helper: build ranked list ─────────────────────────────────────────────────

make_ranked <- function(de, remove_ribo = FALSE) {
  if (remove_ribo) de <- de %>% filter(!grepl("^Rps|^Rpl", gene))
  ranked <- setNames(de$avg_log2FC, de$gene)
  sort(ranked, decreasing = TRUE)
}

# ── Helper: run fgsea ─────────────────────────────────────────────────────────

run_fgsea_safe <- function(pw, ranked) {
  tryCatch(
    fgsea(pathways = pw, stats = ranked,
          minSize = 10, maxSize = 500, nPermSimple = 10000) %>%
      mutate(leading_edge_genes = sapply(leadingEdge, paste, collapse = " "),
             leading_edge_n     = lengths(leadingEdge)) %>%
      select(pathway, NES, pval, padj, size, ES, leading_edge_n, leading_edge_genes) %>%
      arrange(padj),
    error = function(e) { message("  fgsea error: ", e$message); NULL }
  )
}

# ── Helper: lollipop plot ─────────────────────────────────────────────────────

# Strip common prefixes from pathway names for cleaner labels
clean_name <- function(x) {
  x <- gsub("^(HALLMARK_|REACTOME_|GOBP_|GTRD_|MH_)", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_title(x)
}

plot_lollipop <- function(res, title, filename, n_top = 15) {
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))

  sig <- res %>% filter(padj < 0.05)
  if (nrow(sig) == 0) {
    message("  No significant pathways — skipping plot: ", filename)
    return(invisible(NULL))
  }

  top_up   <- sig %>% filter(NES > 0) %>% arrange(padj) %>% head(n_top)
  top_down <- sig %>% filter(NES < 0) %>% arrange(padj) %>% head(n_top)
  plot_df  <- bind_rows(top_up, top_down) %>%
    mutate(
      label     = clean_name(pathway),
      label     = factor(label, levels = rev(label)),
      neg_log10 = pmin(-log10(padj), 10),
      direction = ifelse(NES > 0, "Up in MUT", "Down in MUT")
    )

  x_lim <- max(abs(plot_df$NES), na.rm = TRUE) * 1.1

  p <- ggplot(plot_df, aes(x = NES, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, y = label, yend = label),
                 colour = "grey50", linewidth = 0.6) +
    geom_point(aes(colour = direction, size = neg_log10)) +
    scale_colour_manual(values = c("Up in MUT" = "#D6604D", "Down in MUT" = "#2166AC"),
                        name = NULL) +
    scale_size_continuous(range = c(2, 7), name = expression(-log[10](FDR))) +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    labs(title = title, x = "NES (MUT vs WT)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y    = element_text(size = 8),
      legend.position = "right",
      plot.title     = element_text(face = "bold", size = 11)
    )

  plot_h <- max(4, nrow(plot_df) * 0.28 + 1.5)
  options(bitmapType = "quartz")
  ggsave(file.path(plots_dir, filename), p, width = 11, height = plot_h)
  message("  Saved: ", filename)
}

# ── Main loop: all cell types × all collections × with/without ribo ───────────

collection_names <- c(hallmark = "Hallmark", reactome = "Reactome",
                       gobp = "GO BP", tf = "TF (GTRD)")

for (ct in cell_order) {

  de_ct <- de_all %>% filter(cell_type == ct)

  if (nrow(de_ct) < 10) {
    message("Skipping ", ct, ": insufficient DE genes.")
    next
  }

  message("\n── ", ct, " ──────────────────────────────────────────────")

  ranked_all   <- make_ranked(de_ct, remove_ribo = FALSE)
  ranked_noribo <- make_ranked(de_ct, remove_ribo = TRUE)

  n_ribo <- length(ranked_all) - length(ranked_noribo)
  message("  Genes: ", length(ranked_all), " total, ", n_ribo, " Rps/Rpl removed in noribo version.")

  ct_slug <- gsub(" ", "_", ct)

  for (col in names(pathways)) {

    col_name <- collection_names[col]
    message("  Running ", col_name, " ...")

    # With all genes
    res_all <- run_fgsea_safe(pathways[[col]], ranked_all)
    if (!is.null(res_all)) {
      write.csv(res_all,
                file.path(results_dir, paste0(ct_slug, "_", col, "_allgenes.csv")),
                row.names = FALSE)
      plot_lollipop(res_all,
                    title    = paste0(ct, " — ", col_name, " (all genes)"),
                    filename = paste0(ct_slug, "_", col, "_allgenes.pdf"))
    }

    # Without ribo genes
    res_noribo <- run_fgsea_safe(pathways[[col]], ranked_noribo)
    if (!is.null(res_noribo)) {
      write.csv(res_noribo,
                file.path(results_dir, paste0(ct_slug, "_", col, "_noribo.csv")),
                row.names = FALSE)
      plot_lollipop(res_noribo,
                    title    = paste0(ct, " — ", col_name, " (Rps/Rpl excluded)"),
                    filename = paste0(ct_slug, "_", col, "_noribo.pdf"))
    }
  }
}

message("\nGSEA done. Results in ", results_dir, "\nPlots in ", plots_dir)

# ── Per-stage top-60 DEG heatmaps ─────────────────────────────────────────────
# For each B-cell stage: top 30 up + top 30 down DEGs (by adjusted p-value within
# each direction, Rps/Rpl excluded). Columns = biological replicates (WT1-3,
# Mut1-3); values = mean SCT normalised expression per mouse, z-scored per gene.

message("\nGenerating per-stage top-60 DEG heatmaps ...")

library(pheatmap)
library(RColorBrewer)

seu <- readRDS("data/bcells_annotated.rds")
Idents(seu) <- "FinalLab3"

rep_order   <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")
sct_data    <- GetAssayData(seu, assay = "SCT", layer = "data")

col_anno_df <- data.frame(
  Genotype  = c("WT", "WT", "WT", "MUT", "MUT", "MUT"),
  row.names = rep_order
)
anno_colours <- list(
  Genotype  = c(WT = "#2166AC", MUT = "#D6604D"),
  Direction = c("Up in MUT" = "#D6604D", "Down in MUT" = "#2166AC")
)

for (ct in cell_order) {

  de_ct <- de_all %>%
    filter(cell_type == ct, !grepl("^Rps|^Rpl", gene)) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))

  if (nrow(de_ct) < 10) {
    message("  Skipping ", ct, ": insufficient DE genes.")
    next
  }

  top_up   <- de_ct %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj) %>% head(30) %>% pull(gene)
  top_down <- de_ct %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj) %>% head(30) %>% pull(gene)
  top_genes <- c(top_up, top_down)

  ct_cells  <- WhichCells(seu, idents = ct)
  meta_ct   <- seu@meta.data[ct_cells, ]
  genes_use <- intersect(top_genes, rownames(sct_data))

  avg_mat <- sapply(rep_order, function(rep) {
    cells <- intersect(rownames(meta_ct[meta_ct$Mouse == rep, ]), ct_cells)
    if (length(cells) == 0) return(rep(NA_real_, length(genes_use)))
    if (length(cells) == 1) return(as.numeric(sct_data[genes_use, cells]))
    rowMeans(as.matrix(sct_data[genes_use, cells]))
  })
  rownames(avg_mat) <- genes_use

  avg_scaled <- t(scale(t(avg_mat)))
  avg_scaled[is.nan(avg_scaled)] <- 0

  direction   <- ifelse(de_ct$avg_log2FC[match(genes_use, de_ct$gene)] > 0,
                        "Up in MUT", "Down in MUT")
  row_anno_df <- data.frame(Direction = direction, row.names = genes_use)

  cell_h_pts <- 10
  plot_h     <- (nrow(avg_scaled) * cell_h_pts / 72) + 2.5
  ct_slug    <- gsub(" ", "_", ct)

  options(bitmapType = "quartz")
  pdf(file.path(plots_dir, paste0(ct_slug, "_top60_heatmap.pdf")),
      width = 6, height = plot_h)
  pheatmap(
    avg_scaled,
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
    main              = paste0(ct, " top 60 DEGs (MUT vs WT, Rps/Rpl excluded)")
  )
  dev.off()
  message("  Saved: ", ct_slug, "_top60_heatmap.pdf")
}

message("\nDone. All fig5b outputs saved under output/fig5/")
