# Differential expression of Pax5 across B cell developmental stages (MUT vs WT)
#
# This script loads a labelled Seurat object (bcells_annotated.rds) and tests whether
# Pax5 expression differs between MUT and WT cells within each B cell type defined by
# the FinalLab3 metadata column. Cell types analysed: PreproB, ProB, PreB, Immature B,
# FoB, MzB. Two complementary DE methods are run in parallel:
#
#   1. SCT assay (Wilcoxon): PrepSCTFindMarkers() re-corrects the SCT normalisation
#      across all compared cells, then FindMarkers() is run on the normalised SCT data
#      slot. This is the Seurat-recommended workflow for SCT-normalised data.
#
#   2. RNA assay (negative binomial): FindMarkers() is run directly on raw UMI counts
#      from the RNA assay using the negbinom test, which models count overdispersion
#      explicitly and is appropriate for raw integer count data.
#
# For both methods, fold changes are expressed as MUT relative to WT (positive FC =
# higher in MUT). Results are exported as CSVs and visualised as lollipop plots where
# the x-axis shows log2FC, dot size encodes the mean % of cells expressing Pax5 (0-100),
# and adjusted p-values are labelled on each point. A per-cell CSV combining all metadata
# with Pax5 expression values from both assays is also exported.

source("setup.R")

plots_dir <- "output/fig4/plots"
sdata_dir <- "output/fig4/supporting_data"

# Load object
seu <- readRDS("data/bcells_annotated.rds")

# Set identities to FinalLab3
Idents(seu) <- "FinalLab3"

# Re-correct SCT model across all cells before DE (Seurat recommended workflow)
seu <- PrepSCTFindMarkers(seu)

cell_types <- levels(Idents(seu))

# --- Export per-cell metadata + Pax5 expression ---
pax5_sct <- FetchData(seu, vars = "Pax5", assay = "SCT", layer = "data")
pax5_rna <- FetchData(seu, vars = "Pax5", assay = "RNA", layer = "counts")
cell_meta <- seu@meta.data

per_cell <- cell_meta %>%
  mutate(
    cell_barcode      = rownames(cell_meta),
    Pax5_SCT_norm     = pax5_sct[rownames(cell_meta), "Pax5"],
    Pax5_RNA_counts   = pax5_rna[rownames(cell_meta), "Pax5"]
  )

write.csv(per_cell, file.path(sdata_dir, "pax5_per_cell_metadata_expression.csv"), row.names = FALSE)

# --- DE helper ---
run_de <- function(seu, cell_types, assay, slot, test) {
  results_list <- list()

  for (ct in cell_types) {
    ct_cells  <- WhichCells(seu, idents = ct)
    meta_ct   <- seu@meta.data[ct_cells, ]

    wt_cells  <- rownames(meta_ct[meta_ct$Genotype == "WT",  ])
    mut_cells <- rownames(meta_ct[meta_ct$Genotype == "MUT", ])

    if (length(wt_cells) < 3 || length(mut_cells) < 3) {
      message("Skipping ", ct, ": insufficient cells (WT=", length(wt_cells), ", MUT=", length(mut_cells), ")")
      next
    }

    tryCatch({
      res <- FindMarkers(
        seu,
        ident.1         = mut_cells,
        ident.2         = wt_cells,
        features        = "Pax5",
        test.use        = test,
        assay           = assay,
        slot            = slot,
        logfc.threshold = 0,
        min.pct         = 0
      )
      res$cell_type <- ct
      res$gene      <- rownames(res)
      results_list[[ct]] <- res
    }, error = function(e) {
      message("Error in ", ct, " [", assay, "]: ", e$message)
    })
  }

  results <- bind_rows(results_list)
  if (nrow(results) == 0) return(results)
  results <- results[, c("cell_type", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]
  results[order(results$p_val_adj), ]
}

# --- Method 1: SCT assay (PrepSCTFindMarkers + Wilcoxon on normalised data) ---
message("\n=== SCT assay (Wilcoxon) ===")
res_sct <- run_de(seu, cell_types, assay = "SCT", slot = "data",   test = "wilcox")
print(res_sct)
write.csv(res_sct, file.path(sdata_dir, "pax5_DE_MUT_vs_WT_SCT.csv"), row.names = FALSE)

# --- Method 2: RNA assay (negbinom on raw counts) ---
message("\n=== RNA assay (negbinom) ===")
res_rna <- run_de(seu, cell_types, assay = "RNA", slot = "counts", test = "negbinom")
print(res_rna)
write.csv(res_rna, file.path(sdata_dir, "pax5_DE_MUT_vs_WT_RNA.csv"), row.names = FALSE)

# --- Lollipop plots ---

cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")

format_pval <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
         ifelse(p < 0.05, formatC(p, format = "f", digits = 3),
                formatC(p, format = "f", digits = 3)))
}

make_lollipop <- function(res, title) {
  df <- res %>%
    filter(cell_type %in% cell_order) %>%
    mutate(
      cell_type = factor(cell_type, levels = cell_order),
      pct_expr  = (pct.1 + pct.2) / 2 * 100,
      p_label   = paste0("p=", format_pval(p_val_adj))
    )

  x_range <- max(abs(df$avg_log2FC), na.rm = TRUE)
  x_limit <- x_range * 1.6

  ggplot(df, aes(x = avg_log2FC, y = cell_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = avg_log2FC, y = cell_type, yend = cell_type),
                 colour = "grey40", linewidth = 0.7) +
    geom_point(aes(size = pct_expr, colour = avg_log2FC)) +
    geom_text(aes(x = avg_log2FC + sign(avg_log2FC) * x_range * 0.12,
                  label = p_label),
              size = 3, hjust = ifelse(df$avg_log2FC >= 0, 0, 1)) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey80", high = "#D6604D",
                           midpoint = 0, name = "log2FC\n(MUT/WT)") +
    scale_size_continuous(range = c(2, 9), limits = c(0, 100), name = "% cells\nexpressing") +
    scale_y_discrete(limits = rev(cell_order)) +
    xlim(-x_limit, x_limit) +
    labs(title = title,
         x     = "log2 fold change (MUT vs WT)",
         y     = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y     = element_text(size = 11),
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 13)
    )
}

p_sct <- make_lollipop(res_sct, "Pax5 MUT vs WT — SCT (Wilcoxon)")
p_rna <- make_lollipop(res_rna, "Pax5 MUT vs WT — RNA (negbinom)")

ggsave(file.path(plots_dir, "pax5_lollipop_SCT.pdf"), p_sct, width = 7, height = 4)
ggsave(file.path(plots_dir, "pax5_lollipop_RNA.pdf"), p_rna, width = 7, height = 4)

message("\nDone. Fig4d outputs saved under output/fig4/")
