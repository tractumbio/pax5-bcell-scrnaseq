# UMAP visualisation of B-cell label annotation — Figure 4a-c
#
# This script loads the fully annotated query Seurat object produced by
# data_processing/label_transfer_and_annotation.R (bcells_annotated.rds) and generates
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

plots_dir <- "output/fig4/plots"
sdata_dir <- "output/fig4/supporting_data"

b.query <- readRDS("data/bcells_annotated.rds")

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

ggsave(file.path(plots_dir, "Fig4_clusters.pdf"),          p_clust,  width = 8, height = 6)
ggsave(file.path(plots_dir, "Fig4a_finelabels.pdf"),       p_fine,   width = 8, height = 6)
ggsave(file.path(plots_dir, "Fig4b_merged_labels.pdf"),    p_merged, width = 8, height = 6)
ggsave(file.path(plots_dir, "Fig4c_mature_sublabels.pdf"), p_mature, width = 8, height = 6)

# ── Per-cell metadata + UMAP coordinates CSV ──────────────────────────────────

umap_coords <- as.data.frame(Embeddings(b.query, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

per_cell <- cbind(
  cell            = rownames(b.query@meta.data),
  umap_coords,
  b.query@meta.data
)

write.csv(per_cell, file.path(sdata_dir, "fig4_per_cell_metadata_umap.csv"), row.names = FALSE)

# ── Fig 4d: canonical marker dot plot ─────────────────────────────────────────
# Literature-backed identity markers (5 per stage), grouped by developmental
# stage. Cell types ordered PreproB -> MzB on the y-axis.

marker_list <- list(
  # Cd7/Runx2/Flt3 are PreproB-specific in this dataset; Ly6d dropped as it
  # actually peaks in Immature B / MzB here rather than PreproB.
  PreproB    = c("Cd7", "Runx2", "Flt3", "Irf8", "Tcf4"),
  ProB       = c("Rag1", "Dntt", "Vpreb1", "Igll1", "Ebf1"),
  PreB       = c("Ccnd2", "Vpreb3", "Cd2", "Cd24a", "Myb"),
  `Immature B` = c("Ighm", "Ms4a1", "Cd79a", "Cd79b", "Cd74"),
  FoB        = c("Fcer2a", "Ccr7", "Sell", "Klf2", "Gpr183"),
  MzB        = c("Cr2", "Cd1d1", "Cd9", "S1pr3", "Dtx1")
)

# Keep only genes present in the object (warn about any missing)
all_markers <- unlist(marker_list, use.names = FALSE)
present     <- all_markers[all_markers %in% rownames(b.query)]
missing     <- setdiff(all_markers, present)
if (length(missing) > 0)
  message("Marker genes not found in object (skipped): ", paste(missing, collapse = ", "))

Idents(b.query)    <- "FinalLab3"
b.query$FinalLab3  <- factor(b.query$FinalLab3,
                             levels = c(cell_order, setdiff(unique(b.query$FinalLab3), cell_order)))

p_dot <- DotPlot(
  b.query,
  features = marker_list,
  group.by = "FinalLab3",
  assay    = "SCT",
  cols     = c("lightgrey", "#D6604D")
) +
  scale_y_discrete(limits = rev(cell_order)) +
  labs(title = "Canonical B-cell stage markers", x = NULL, y = NULL) +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text       = element_text(size = 9, face = "bold"),
    plot.title       = element_text(face = "bold")
  )

ggsave(file.path(plots_dir, "Fig4d_marker_dotplot.pdf"), p_dot, width = 14, height = 5)

message("Done. Fig4a-d outputs saved under output/fig4/")
