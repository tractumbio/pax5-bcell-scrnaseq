# Pax5 target gene set enrichment across B cell developmental stages (MUT vs WT)
#
# This script tests whether genes activated or repressed by Pax5 are coordinately
# dysregulated in MUT compared to WT cells across six B cell developmental stages
# (PreproB, ProB, PreB, Immature B, FoB, MzB) defined by the FinalLab3 metadata column.
#
# The analysis proceeds in two stages per cell type:
#
# STAGE 1 — Differential expression (MUT vs WT)
#   FindMarkers() is run using the negative binomial test on raw UMI counts from the
#   RNA assay (assay = "RNA", slot = "counts"). This tests all genes expressed in
#   each cell type. Fold changes are expressed as MUT relative to WT (positive FC =
#   higher in MUT). Cell types with fewer than 3 cells in either genotype are skipped.
#
# STAGE 2 — Gene set enrichment (fgsea)
#   All genes tested in Stage 1 are ranked by average log2 fold change (MUT vs WT).
#   fgsea is then run using this full ranked list against two gene sets from
#   pax5_targets.rds: genes activated by Pax5 and genes repressed by Pax5. A negative
#   NES for the activated set indicates that Pax5 targets are down-regulated in MUT
#   (consistent with loss of Pax5 function); a positive NES for the repressed set
#   indicates de-repression.
#
# OUTPUTS
#   pax5_target_DE_RNA.csv          — full DE results for all genes, all cell types
#   pax5_fgsea_results.csv          — fgsea NES, p-value, FDR per gene set per cell type
#   pax5_fgsea_enrichment_plot.pdf  — dot/bar summary plot of NES across cell types

source("setup.R")

# Assumes working directory is the repository root (set by main_script.R or manually)

out_dir    <- "figures/4e_pax5_target_enrichment"
cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")

# ── Load data ──────────────────────────────────────────────────────────────────

seu          <- readRDS("data/anselm_labeled.rds")
pax5_targets <- readRDS("data/pax5_targets.rds")

Idents(seu) <- "FinalLab3"
cell_types  <- intersect(cell_order, levels(Idents(seu)))

# ── STAGE 1: DE for all genes per cell type (negbinom, RNA counts) ─────────────

de_cache <- file.path(out_dir, "pax5_target_DE_cache.rds")

if (file.exists(de_cache)) {
  warning("DE cache found at '", de_cache, "' — skipping differential expression and using previously computed results.")
  de_list <- readRDS(de_cache)
} else {
  de_list <- list()

  for (ct in cell_types) {
    ct_cells  <- WhichCells(seu, idents = ct)
    meta_ct   <- seu@meta.data[ct_cells, ]

    wt_cells  <- rownames(meta_ct[meta_ct$Genotype == "WT",  ])
    mut_cells <- rownames(meta_ct[meta_ct$Genotype == "MUT", ])

    if (length(wt_cells) < 3 || length(mut_cells) < 3) {
      message("Skipping ", ct, ": insufficient cells (WT=", length(wt_cells), ", MUT=", length(mut_cells), ")")
      next
    }

    message("Running DE for ", ct, " ...")

    tryCatch({
      res <- FindMarkers(
        seu,
        ident.1         = mut_cells,
        ident.2         = wt_cells,
        test.use        = "negbinom",
        assay           = "RNA",
        slot            = "counts",
        logfc.threshold = 0,
        min.pct         = 0,
        verbose         = FALSE
      )
      res$gene      <- rownames(res)
      res$cell_type <- ct
      de_list[[ct]] <- res
    }, error = function(e) {
      message("Error in ", ct, ": ", e$message)
    })
  }

  saveRDS(de_list, de_cache)
  message("DE results saved to cache: ", de_cache)
}

de_results <- bind_rows(de_list)
de_results <- de_results[, c("cell_type", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]

write.csv(de_results, file.path(out_dir, "pax5_target_DE_RNA.csv"), row.names = FALSE)

# ── STAGE 2: fgsea using avg_log2FC as ranking metric ─────────────────────────

fgsea_list <- list()

for (ct in names(de_list)) {
  res <- de_list[[ct]]

  # Full ranked list of all genes tested, ranked by avg_log2FC (MUT vs WT)
  ranked <- setNames(res$avg_log2FC, res$gene)
  ranked <- sort(ranked, decreasing = TRUE)

  tryCatch({
    fg <- fgsea(
      pathways   = pax5_targets,
      stats      = ranked,
      minSize    = 5,
      maxSize    = 500,
      nPermSimple = 10000
    )
    fg$cell_type <- ct
    fgsea_list[[ct]] <- fg
  }, error = function(e) {
    message("fgsea error in ", ct, ": ", e$message)
  })
}

fgsea_results <- bind_rows(fgsea_list)
fgsea_results <- fgsea_results %>%
  select(cell_type, pathway, NES, pval, padj, size, ES) %>%
  arrange(cell_type, padj)

write.csv(fgsea_results, file.path(out_dir, "pax5_fgsea_results.csv"), row.names = FALSE)

# ── Plot: NES per gene set per cell type ───────────────────────────────────────

plot_df <- fgsea_results %>%
  mutate(
    cell_type  = factor(cell_type, levels = cell_order),
    sig_label  = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

p <- ggplot(plot_df, aes(x = cell_type, y = NES, fill = pathway)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sig_label,
                y = NES + sign(NES) * 0.05),
            position = position_dodge(width = 0.7),
            size = 4, vjust = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c("#D6604D", "#2166AC"),
                    name   = "Pax5 target set") +
  scale_x_discrete(limits = cell_order) +
  labs(
    title = "Pax5 target gene set enrichment (MUT vs WT)",
    x     = NULL,
    y     = "Normalised Enrichment Score (NES)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13)
  )

ggsave(file.path(out_dir, "pax5_fgsea_enrichment_plot.pdf"), p, width = 8, height = 5)

message("\nDone. All outputs saved to ", out_dir)
