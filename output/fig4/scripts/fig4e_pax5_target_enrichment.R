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

# Assumes working directory is the repository root (set by main_script.R)
plots_dir <- "output/fig4/plots"
sdata_dir <- "output/fig4/supporting_data"
cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")

# ── Load data ──────────────────────────────────────────────────────────────────

seu          <- readRDS("data/bcells_annotated.rds")
pax5_targets <- readRDS("data/pax5_targets.rds")

Idents(seu) <- "FinalLab3"
cell_types  <- intersect(cell_order, levels(Idents(seu)))

# ── STAGE 1: DE for all genes per cell type (negbinom, RNA counts) ─────────────

de_cache <- file.path(sdata_dir, "pax5_target_DE_cache.rds")

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

write.csv(de_results, file.path(sdata_dir, "pax5_target_DE_RNA.csv"), row.names = FALSE)

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
  mutate(leading_edge_n = lengths(leadingEdge)) %>%
  select(cell_type, pathway, NES, pval, padj, size, ES, leading_edge_n) %>%
  arrange(cell_type, padj)

# CSV 1: enrichment results with leading edge genes as a semicolon-separated string
fgsea_export <- bind_rows(fgsea_list) %>%
  mutate(
    leading_edge_n     = lengths(leadingEdge),
    leading_edge_genes = sapply(leadingEdge, paste, collapse = " ")
  ) %>%
  select(cell_type, pathway, NES, pval, padj, size, ES, leading_edge_n, leading_edge_genes) %>%
  arrange(cell_type, padj)

write.csv(fgsea_export, file.path(sdata_dir, "pax5_fgsea_results_with_leading_edge.csv"), row.names = FALSE)

# CSV 2: leading edge genes × DE results
de_full <- bind_rows(de_list)

leading_edge_de <- fgsea_export %>%
  select(cell_type, pathway, leading_edge_genes) %>%
  tidyr::separate_rows(leading_edge_genes, sep = ";") %>%
  rename(gene = leading_edge_genes) %>%
  left_join(de_full, by = c("cell_type", "gene")) %>%
  select(cell_type, pathway, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2) %>%
  arrange(cell_type, pathway, p_val_adj)

write.csv(leading_edge_de, file.path(sdata_dir, "pax5_leading_edge_DE_results.csv"), row.names = FALSE)

# CSV 3: DE results subset to leading edge genes (one row per unique gene × cell type)
leading_edge_genes_per_ct <- fgsea_export %>%
  select(cell_type, leading_edge_genes) %>%
  tidyr::separate_rows(leading_edge_genes, sep = " ") %>%
  rename(gene = leading_edge_genes) %>%
  distinct(cell_type, gene)

leading_edge_de_subset <- de_full %>%
  inner_join(leading_edge_genes_per_ct, by = c("cell_type", "gene")) %>%
  select(cell_type, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2) %>%
  arrange(cell_type, p_val_adj)

write.csv(leading_edge_de_subset, file.path(sdata_dir, "pax5_leading_edge_DE_subset.csv"), row.names = FALSE)

# Rebuild clean fgsea_results for downstream plotting
fgsea_results <- fgsea_export

# ── Count genes per set expressed in >10 cells (RNA counts assay) ─────────────

counts_mat   <- GetAssayData(seu, assay = "RNA", layer = "counts")
cells_per_gene <- rowSums(counts_mat > 0)
expressed_genes <- names(cells_per_gene[cells_per_gene > 10])

geneset_sizes <- sapply(pax5_targets, function(genes) {
  sum(genes %in% expressed_genes)
})

# Build facet labels: "Activated (n=142)"
facet_labels <- setNames(
  paste0(names(geneset_sizes), " (n=", geneset_sizes, ")"),
  names(geneset_sizes)
)

# ── Plot: lollipop, NES on x, cell type on y, faceted by pathway ──────────────

x_limit <- max(abs(fgsea_results$NES), na.rm = TRUE) * 1.5

plot_df <- fgsea_results %>%
  mutate(
    cell_type     = factor(cell_type, levels = rev(cell_order)),
    pathway_label = facet_labels[pathway],
    p_label       = paste0(
      "p=", ifelse(padj < 0.001,
                   formatC(padj, format = "e", digits = 1),
                   formatC(padj, format = "f", digits = 3))
    )
  )

p <- ggplot(plot_df, aes(x = NES, y = cell_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES, y = cell_type, yend = cell_type),
               colour = "grey40", linewidth = 0.7) +
  geom_point(aes(colour = NES), size = 3) +
  geom_text(aes(x = NES + sign(NES) * x_limit * 0.04,
                label = p_label),
            size = 2.8,
            hjust = ifelse(plot_df$NES >= 0, 0, 1)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey70", high = "#D6604D",
                         midpoint = 0, name = "NES") +
  scale_x_continuous(limits = c(-x_limit, x_limit)) +
  facet_wrap(~ pathway_label, ncol = 2) +
  labs(
    title = "Pax5 target gene set enrichment (MUT vs WT)",
    x     = "Normalised Enrichment Score (NES)",
    y     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    panel.spacing   = unit(1.2, "lines")
  )

ggsave(file.path(plots_dir, "pax5_fgsea_enrichment_plot.pdf"), p, width = 10, height = 5)

# ── Per-cell-type leading edge barplots ────────────────────────────────────────

# Join leading edge genes with their pathway and DE stats
le_plot_df <- fgsea_export %>%
  select(cell_type, pathway, leading_edge_genes) %>%
  tidyr::separate_rows(leading_edge_genes, sep = " ") %>%
  rename(gene = leading_edge_genes) %>%
  left_join(
    bind_rows(de_list) %>% select(cell_type, gene, avg_log2FC, p_val_adj),
    by = c("cell_type", "gene")
  ) %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(
    neg_log10_p   = -log10(p_val_adj + .Machine$double.eps),
    pathway_label = facet_labels[pathway]
  )

for (ct in cell_order) {
  df_ct <- le_plot_df %>%
    filter(cell_type == ct) %>%
    group_by(pathway_label) %>%
    mutate(gene = reorder(gene, avg_log2FC)) %>%
    ungroup()

  if (nrow(df_ct) == 0) next

  n_genes   <- length(unique(df_ct$gene))
  plot_h    <- max(4, n_genes * 0.22 + 2)

  p_le <- ggplot(df_ct, aes(x = avg_log2FC, y = gene, fill = neg_log10_p)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_gradient(low = "grey85", high = "#B2182B",
                        name = expression(-log[10](p[adj]))) +
    facet_wrap(~ pathway_label, ncol = 2, scales = "free_y") +
    labs(
      title = paste0("Leading edge genes — ", ct, " (MUT vs WT)"),
      x     = "log2 fold change (MUT vs WT)",
      y     = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      axis.text.y      = element_text(size = 7),
      legend.position  = "right",
      plot.title       = element_text(face = "bold", size = 12),
      panel.spacing    = unit(1, "lines")
    )

  fname <- file.path(plots_dir, paste0("pax5_leading_edge_", gsub(" ", "_", ct), ".pdf"))
  ggsave(fname, p_le, width = 10, height = plot_h)
}

message("\nDone. Fig4e outputs saved under output/fig4/")
