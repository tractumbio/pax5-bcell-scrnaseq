# UMAP visualisation of B-cell label annotation — Figure 4a-c
#
# This script loads the fully annotated query Seurat object produced by
# data_processing/label_transfer_and_annotation.R (anselm_labeled.rds) and generates
# four UMAP plots that together constitute Figure 4a-c:
#
#   Fig 4 (clusters)  — cells coloured by unsupervised Seurat cluster identity,
#                        showing the raw clustering output before any label transfer.
#
#   Fig 4a (FinalLab) — cells coloured by fine-grained transferred labels (up to 13
#                        B-cell subtypes directly inherited from the Lee et al. reference
#                        atlas via majority-vote assignment per cluster).
#
#   Fig 4b (FinalLab2) — same as 4a but with related sub-types collapsed: all Pro B
#                         variants merged into "ProB", all Pre B / BCR variants into
#                         "PreB". Gives a cleaner view of major developmental stages.
#
#   Fig 4c (FinalLab3) — extends FinalLab2 by splitting "Mature B" clusters into
#                         functional sub-types (FoB, MzB, GerminalCenter, Bmem, DZ, LZ)
#                         based on which marker gene module score is highest in each
#                         cluster. This is the final label used in all downstream analyses.
#
# In addition to the PDF plots, a CSV is exported containing, for every cell:
#   - UMAP coordinates (UMAP_1, UMAP_2)
#   - All metadata columns from the Seurat object
#   - FinalLab, FinalLab2, FinalLab3 label columns
# This table can be used for statistical analyses or figure reproduction without
# needing to reload the full Seurat object.

source("setup.R")

# Assumes working directory is the repository root (set by main_script.R or manually)

out_dir <- "figures/4a-c_umap_plots"

b.query <- readRDS("data/anselm_labeled.rds")

# Set consistent cell type ordering for all legends and labels
cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")

b.query$FinalLab2 <- factor(b.query$FinalLab2,
                             levels = c(cell_order, setdiff(unique(b.query$FinalLab2), cell_order)))
b.query$FinalLab3 <- factor(b.query$FinalLab3,
                             levels = c(cell_order, setdiff(unique(b.query$FinalLab3), cell_order)))

# ── Plots ──────────────────────────────────────────────────────────────────────

p_clust <- DimPlot(b.query, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Seurat clusters")

p_fine <- DimPlot(b.query, group.by = "FinalLab", label = TRUE, repel = TRUE) +
  ggtitle("Fine-grained B-cell labels (FinalLab)")

p_merged <- DimPlot(b.query, group.by = "FinalLab2", label = TRUE, repel = TRUE) +
  ggtitle("Merged B-cell labels (FinalLab2)")

p_mature <- DimPlot(b.query, group.by = "FinalLab3", label = TRUE, repel = TRUE) +
  ggtitle("B-cell labels with Mature B module sub-classification (FinalLab3)")

ggsave(file.path(out_dir, "Fig4_clusters.pdf"),          p_clust,  width = 8, height = 6)
ggsave(file.path(out_dir, "Fig4a_finelabels.pdf"),       p_fine,   width = 8, height = 6)
ggsave(file.path(out_dir, "Fig4b_merged_labels.pdf"),    p_merged, width = 8, height = 6)
ggsave(file.path(out_dir, "Fig4c_mature_sublabels.pdf"), p_mature, width = 8, height = 6)

# ── Per-cell metadata + UMAP coordinates CSV ──────────────────────────────────

umap_coords <- as.data.frame(Embeddings(b.query, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

per_cell <- cbind(
  cell            = rownames(b.query@meta.data),
  umap_coords,
  b.query@meta.data
)

write.csv(per_cell, file.path(out_dir, "fig4_per_cell_metadata_umap.csv"), row.names = FALSE)

message("Done. Plots and metadata CSV saved to ", out_dir)
