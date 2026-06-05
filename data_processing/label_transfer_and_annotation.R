# B-cell label transfer and hierarchical annotation of the Pax5 query dataset
#
# This script takes two Seurat objects — a published reference atlas of mouse B-cell
# subtypes (Lee et al., Nature Communications) and an experimental query dataset
# (anselm_demux_filtered.rds) — and assigns B-cell identity labels to each cell in
# the query through the following steps:
#
# STEP 1 — Preprocessing
#   Both objects are normalised, scored for cell cycle phase (S and G2M scores), and
#   re-normalised with SCTransform while regressing out cell cycle variation. The query
#   is additionally clustered (Leiden algorithm, resolution 0.5) and embedded in UMAP
#   space using the top 40 PCs.
#
# STEP 2 — SCT integration
#   Reference and query are integrated into a shared embedding using Seurat's SCT-aware
#   anchor-based integration (FindIntegrationAnchors + IntegrateData). The combined
#   object is re-clustered and re-embedded in the integrated PCA space.
#
# STEP 3 — Label transfer
#   Cell type labels from the integrated reference are transferred to query cells using
#   FindTransferAnchors + TransferData. A majority-vote rule then assigns each Seurat
#   cluster its most frequently predicted label, producing FinalLab (13 fine-grained
#   subtypes). Sub-types are further collapsed: all "Pro B" variants → ProB, all
#   "Pre B / BCR" variants → PreB, yielding FinalLab2.
#
# STEP 4 — Mature B sub-classification
#   Marker gene module scores (FoB, MzB, GerminalCenter, Bmem, DZ, LZ) are added for
#   each cell. Clusters labelled "Mature B" in FinalLab are re-assigned to whichever
#   mature sub-type has the highest mean module score in that cluster, producing
#   FinalLab3 — the final label column used downstream.
#
# OUTPUT
#   data/anselm_labeled.rds   — fully annotated Seurat object with FinalLab/2/3
#   plots/Fig4_umap_data.csv  — UMAP coordinates + all label columns per cell

source("setup.R")

# Assumes working directory is the repository root (set by main_script.R or manually)

# ── STEP 1: Preprocessing ─────────────────────────────────────────────────────

reference <- readRDS("data/lee_etat_natcom_data.rds")
query     <- readRDS("data/anselm_demux_filtered.rds")

# Remove low-quality reference cluster
Idents(reference) <- "ident"
reference <- subset(reference, idents = "High Mitochondrial B", invert = TRUE)

# Mouse gene symbols use title case
s.genes   <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)

reference <- reference %>%
  NormalizeData() %>%
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
  SCTransform(vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)

query <- query %>%
  NormalizeData() %>%
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
  SCTransform(vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:40, verbose = FALSE) %>%
  FindClusters(algorithm = 4, resolution = 0.5, verbose = FALSE) %>%
  RunUMAP(dims = 1:40, verbose = FALSE)

# ── STEP 2: SCT integration ────────────────────────────────────────────────────

intList  <- list(ref = reference, q = query)
features <- SelectIntegrationFeatures(object.list = intList, nfeatures = 3000)
intList  <- PrepSCTIntegration(object.list = intList, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(
  object.list          = intList,
  normalization.method = "SCT",
  anchor.features      = features
)

combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
DefaultAssay(combined) <- "integrated"

combined <- combined %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:40, verbose = FALSE) %>%
  FindClusters(algorithm = 4, resolution = 0.5, verbose = FALSE) %>%
  RunUMAP(dims = 1:40, verbose = FALSE)

# ── STEP 3: Label transfer and majority-vote assignment ────────────────────────

b.query <- intList$q

b.anchors <- FindTransferAnchors(
  reference           = combined,
  query               = b.query,
  dims                = 1:40,
  reference.reduction = "pca"
)

predictions <- TransferData(anchorset = b.anchors, refdata = combined$ident, dims = 1:40)
b.query     <- AddMetaData(b.query, metadata = predictions)

# Majority vote: assign each cluster its most frequent predicted label
metaB  <- b.query@meta.data
newIDs <- table(metaB$predicted.id, metaB$seurat_clusters) %>%
  as.data.frame() %>%
  group_by(Var2) %>%
  top_n(n = 1, wt = Freq) %>%
  filter(Freq > 0)

newIdents        <- setNames(newIDs$Var1, newIDs$Var2)
b.query          <- RenameIdents(b.query, newIdents)
b.query$FinalLab <- Idents(b.query)

# Collapse sub-types → FinalLab2
b.query@meta.data <- b.query@meta.data %>%
  mutate(FinalLab2 = case_when(
    grepl("Pro B", FinalLab)     ~ "ProB",
    grepl("Pre B|BCR", FinalLab) ~ "PreB",
    TRUE                         ~ as.character(FinalLab)
  ))

saveRDS(b.query, "data/anselm_labeled.rds")

# ── STEP 4: Mature B sub-classification by module score ───────────────────────

marker_modules <- readRDS("data/marker_modules.rds")
safe_names     <- gsub("-", "_", names(marker_modules))

for (i in seq_along(marker_modules)) {
  b.query <- AddModuleScore(b.query, features = list(marker_modules[[i]]), name = safe_names[i])
}

mature_clusters   <- names(newIdents)[newIdents == "Mature B"]
mature_modules    <- c("FoB", "MzB", "GerminalCenter", "Bmem", "DZ.genes", "LZ.genes")
mature_score_cols <- paste0(gsub("-", "_", mature_modules), "1")

cluster_scores <- b.query@meta.data %>%
  filter(seurat_clusters %in% mature_clusters) %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(mature_score_cols), mean), .groups = "drop")

best_label <- mature_modules[apply(cluster_scores[, mature_score_cols], 1, which.max)]
mature_map <- setNames(best_label, as.character(cluster_scores$seurat_clusters))

b.query@meta.data <- b.query@meta.data %>%
  mutate(FinalLab3 = case_when(
    as.character(seurat_clusters) %in% names(mature_map) ~ mature_map[as.character(seurat_clusters)],
    TRUE ~ as.character(FinalLab2)
  ))

# Export UMAP coordinates + all label columns
umap_coords <- as.data.frame(Embeddings(b.query, "umap"))
umap_coords$cell            <- rownames(umap_coords)
umap_coords$seurat_clusters <- b.query$seurat_clusters
umap_coords$FinalLab        <- as.character(b.query$FinalLab)
umap_coords$FinalLab2       <- b.query$FinalLab2
umap_coords$FinalLab3       <- b.query$FinalLab3

write.csv(umap_coords, "plots/Fig4_umap_data.csv", row.names = FALSE)

saveRDS(b.query, "data/anselm_labeled.rds")

message("Done. Annotated object saved to data/anselm_labeled.rds")
