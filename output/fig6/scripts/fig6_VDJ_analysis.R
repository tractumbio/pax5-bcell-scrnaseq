# BCR VDJ analysis — Figure 6
#
# Analyses B-cell receptor (BCR) VDJ repertoire data from Cell Ranger
# filtered contig annotations alongside the annotated single-cell object.
#
# A cell is counted as "Productive" if it has ANY productive full-length contig
# (not necessarily a paired productive heavy + light chain).
#
# ANALYSIS 1 — UMAP: all cells and productive BCR cells
#   Two-panel UMAP coloured by B-cell label (FinalLab3).
#
# ANALYSIS 2 — Productive/unproductive ratio by B-cell type
#   Mean productive fraction per B-cell type, averaged across the 3 replicates
#   per genotype with SEM error bars. Mice with cells but no productive BCR
#   contribute an explicit 0 (not dropped), so means/SEMs are unbiased.
#
# ANALYSIS 3 — Dominant clonotype fraction per genotype
#   Mean dominant clonotype fraction per genotype (WT / MUT) with SEM error
#   bars and individual mouse points. Empty/None clonotype IDs are excluded.
#
# ANALYSIS 4 — BCR gene usage heatmaps split by chain
#   Per-gene V/D/J/C proportions per mouse, with explicit zeros for genes
#   absent in a mouse. A GLM (Gamma family) + emmeans tests WT vs MUT for
#   genes present in both genotypes; p-values are FDR-corrected (BH) across
#   all genes. Genotype-exclusive genes (present in only one genotype) cannot
#   be modelled and are recovered separately. Heatmaps (raw + z-scored) per
#   chain (IGH, IGK, IGL) show FDR-significant + genotype-exclusive genes,
#   annotated by segment and significance flag.
#
# INPUTS
#   data/vdj_contigs.rds   (combined Cell Ranger filtered contigs, all 6 mice;
#                           falls back to data/{sample}/filtered_contig_annotations.csv)
#   data/bcells_annotated.rds
#
# OUTPUTS
#   output/fig6/plots/fig6a_umap.pdf
#   output/fig6/plots/fig6b_productive_by_celltype.pdf
#   output/fig6/plots/fig6c_clonotype_expansion.pdf
#   output/fig6/plots/fig6d_gene_usage_{chain}_raw.pdf
#   output/fig6/plots/fig6d_gene_usage_{chain}_scaled.pdf
#   output/fig6/supporting_data/fig6_gene_usage_stats.csv

source("setup.R")
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(emmeans)
library(scales)
library(tibble)
library(purrr)
library(patchwork)

plots_dir <- "output/fig6/plots"
sdata_dir <- "output/fig6/supporting_data"
samples    <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")
cell_order <- c("PreproB", "ProB", "PreB", "Immature B", "FoB", "MzB")
options(bitmapType = "quartz")

# ── Load VDJ contigs ──────────────────────────────────────────────────────────

message("Loading VDJ contigs ...")
# Combined VDJ contig table (all 6 mice) with a Mouse column. Falls back to the
# per-sample Cell Ranger CSVs if the combined RDS is not present.
vdj_rds <- "data/vdj_contigs.rds"
if (file.exists(vdj_rds)) {
  contigs_df <- readRDS(vdj_rds)
} else {
  contigs_df <- lapply(samples, function(s) {
    path <- file.path("data", s, "filtered_contig_annotations.csv")
    if (!file.exists(path)) stop("Missing: ", path)
    read.csv(path, stringsAsFactors = FALSE) %>% mutate(Mouse = s)
  }) %>% bind_rows()
}

contigs_df <- contigs_df %>%
  filter(is_cell == TRUE, productive == TRUE, full_length == TRUE)

message("  ", nrow(contigs_df), " productive full-length contigs, ",
        length(unique(contigs_df$barcode)), " unique cells.")

# ── Load Seurat object ────────────────────────────────────────────────────────

message("Loading Seurat object ...")
seu  <- readRDS("data/bcells_annotated.rds")
umap <- as.data.frame(Embeddings(seu, "umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
umap <- umap[rownames(seu@meta.data), ]   # enforce row alignment

sc_meta <- seu@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  cbind(umap) %>%
  mutate(
    productive = ifelse(barcode %in% contigs_df$barcode, "Productive", "Unproductive"),
    Mouse      = factor(Mouse, levels = samples),
    Genotype   = factor(Genotype, levels = c("WT", "MUT")),
    FinalLab3  = factor(FinalLab3, levels = cell_order)
  )

label_pos <- sc_meta %>%
  group_by(FinalLab3) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")

# ── ANALYSIS 1: UMAP ─────────────────────────────────────────────────────────

message("Plotting UMAP ...")

umap_theme <- theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

p_all <- ggplot(sc_meta, aes(UMAP_1, UMAP_2, colour = FinalLab3)) +
  geom_point(alpha = 0.3, size = 0.4) +
  geom_text(data = label_pos, aes(label = FinalLab3), size = 3,
            colour = "black", fontface = "bold") +
  umap_theme + ggtitle("All cells")

p_prod <- ggplot(sc_meta %>% filter(productive == "Productive"),
                 aes(UMAP_1, UMAP_2, colour = FinalLab3)) +
  geom_point(alpha = 0.4, size = 0.4) +
  geom_text(data = label_pos, aes(label = FinalLab3), size = 3,
            colour = "black", fontface = "bold") +
  umap_theme + ggtitle("Productive BCR cells")

pdf(file.path(plots_dir, "fig6a_umap.pdf"), width = 10, height = 5)
print(p_all + p_prod)
dev.off()
message("  Saved: fig6a_umap.pdf")

# ── ANALYSIS 2: Productive fraction by cell type × genotype ──────────────────

message("Plotting productive fraction by cell type ...")

# Count productive and total per mouse x cell type. Zeros are kept: a mouse
# with cells but no productive BCR contributes frac = 0 (NOT dropped), so the
# mean and SEM are computed over all mice that have cells of that type.
prod_per_mouse <- sc_meta %>%
  filter(FinalLab3 %in% cell_order) %>%
  group_by(Mouse, Genotype, FinalLab3) %>%
  summarise(n_prod  = sum(productive == "Productive"),
            n_total = n(),
            .groups = "drop") %>%
  mutate(frac = n_prod / n_total)

write.csv(prod_per_mouse,
          file.path(sdata_dir, "fig6b_productive_per_mouse.csv"), row.names = FALSE)

prod_summ <- prod_per_mouse %>%
  group_by(Genotype, FinalLab3) %>%
  summarise(mean_frac = mean(frac),
            sem       = sd(frac) / sqrt(n()),
            n_mice    = n(),
            .groups   = "drop") %>%
  mutate(FinalLab3 = factor(FinalLab3, levels = cell_order))

p_ct <- ggplot(prod_summ, aes(x = FinalLab3, y = mean_frac, fill = Genotype)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_frac - sem, ymax = mean_frac + sem),
                position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = c(WT = "#2166AC", MUT = "#D6604D")) +
  scale_x_discrete(limits = cell_order) +
  scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
  labs(x = NULL, y = "% productive cells", fill = "Genotype",
       title = "BCR productive fraction by B-cell type (mean ± SEM)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold"))

pdf(file.path(plots_dir, "fig6b_productive_by_celltype.pdf"), width = 9, height = 5)
print(p_ct)
dev.off()
message("  Saved: fig6b_productive_by_celltype.pdf")

# ── ANALYSIS 3: Dominant clonotype fraction per genotype ─────────────────────

message("Plotting clonotype expansion ...")

comm_barcodes  <- intersect(contigs_df$barcode, colnames(seu))
contigs_filter <- contigs_df %>%
  filter(barcode %in% comm_barcodes) %>%
  left_join(seu@meta.data %>%
              mutate(barcode = rownames(.)) %>%
              select(barcode, Genotype),   # Mouse already in contigs_df
            by = "barcode")

clones <- contigs_filter %>%
  filter(!is.na(raw_clonotype_id), raw_clonotype_id != "", raw_clonotype_id != "None") %>%
  select(barcode, raw_clonotype_id, Mouse, Genotype) %>%
  distinct() %>%
  group_by(Mouse, Genotype, raw_clonotype_id) %>%
  summarise(nClones = n(), .groups = "drop") %>%
  group_by(Mouse, Genotype) %>%
  mutate(clone_frac = nClones / sum(nClones)) %>%
  ungroup()

clones_top <- clones %>%
  group_by(Mouse, Genotype) %>%
  slice_max(nClones, n = 1, with_ties = FALSE)

clones_summ <- clones_top %>%
  group_by(Genotype) %>%
  summarise(mean_frac = mean(clone_frac),
            sem       = sd(clone_frac) / sqrt(n()),
            .groups   = "drop")

p_clone <- ggplot(clones_summ, aes(x = Genotype, y = mean_frac, fill = Genotype)) +
  geom_col(width = 0.6, colour = "black") +
  geom_errorbar(aes(ymin = mean_frac - sem, ymax = mean_frac + sem), width = 0.2) +
  geom_jitter(data = clones_top, aes(y = clone_frac, colour = Mouse),
              width = 0.1, size = 3) +
  scale_fill_manual(values  = c(WT = "#2166AC", MUT = "#D6604D")) +
  scale_colour_manual(values = setNames(viridis(6), samples)) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Dominant clonotype fraction",
       fill = "Genotype", colour = "Mouse",
       title = "Dominant clonotype fraction (mean ± SEM)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

pdf(file.path(plots_dir, "fig6c_clonotype_expansion.pdf"), width = 4, height = 4)
print(p_clone)
dev.off()
message("  Saved: fig6c_clonotype_expansion.pdf")

# ── ANALYSIS 4: Gene usage heatmaps split by chain ───────────────────────────

message("Computing per-gene BCR segment proportions ...")

contigs_split <- contigs_filter %>%
  group_by(Mouse) %>%
  group_split() %>%
  setNames(sapply(., function(x) x$Mouse[1]))

# Per-mouse, per-chain proportion of each gene segment (observed genes only)
prop_obs <- imap(contigs_split, function(x, nm) {
  x %>%
    select(v_gene, d_gene, j_gene, c_gene, chain) %>%
    mutate(across(c(v_gene, d_gene, j_gene, c_gene),
                  ~ifelse(. == "", "NoGene", .))) %>%
    group_by(chain) %>%
    group_split() %>%
    setNames(sapply(., function(d) d$chain[1])) %>%
    imap(function(d, ch) {
      lapply(c("v_gene","d_gene","j_gene","c_gene"), function(seg) {
        tbl <- table(d[[seg]])
        as.data.frame(tbl / sum(tbl)) %>%
          rename(gene = Var1, Freq = Freq) %>%
          mutate(chain = ch, gene_type = seg)
      }) %>% bind_rows()
    }) %>%
    bind_rows() %>%
    mutate(Mouse = nm)
}) %>% bind_rows() %>%
  filter(gene != "NoGene") %>%
  mutate(gene = as.character(gene))

# Fill explicit zeros: every (chain, gene_type, gene) must have all 6 mice.
# A gene absent in a mouse is a true 0, not missing data. This makes the GLM
# means and the heatmap consistent and lets us detect genotype-exclusive genes.
prop_full <- prop_obs %>%
  group_by(chain, gene_type, gene) %>%
  complete(Mouse = samples, fill = list(Freq = 0)) %>%
  ungroup() %>%
  mutate(Genotype  = ifelse(startsWith(Mouse, "W"), "WT", "MUT"),
         gene_type = factor(gene_type, c("v_gene","d_gene","j_gene","c_gene")))

# ── Genotype-exclusive genes (present in exactly one genotype) ─────────────────
# These cannot be tested by the Gamma model (one group is all-zero) but are
# biologically important, so we recover them separately for the heatmaps.
gene_geno <- prop_full %>%
  group_by(chain, gene_type, gene, Genotype) %>%
  summarise(total = sum(Freq), .groups = "drop") %>%
  pivot_wider(names_from = Genotype, values_from = total, values_fill = 0)

exclusive_genes <- gene_geno %>%
  filter((WT == 0) != (MUT == 0)) %>%
  mutate(exclusive = ifelse(WT == 0, "MUT-only", "WT-only")) %>%
  select(chain, gene_type, gene, exclusive)

message("  ", nrow(exclusive_genes), " genotype-exclusive genes recovered.")

# ── GLM + emmeans on testable genes (present in BOTH genotypes; Freq > 0) ──────
message("  Running GLM statistics ...")
stats <- prop_full %>%
  filter(Freq > 0) %>%                    # Gamma requires strictly positive y
  group_by(gene_type) %>%
  group_split() %>%
  map(function(x) {
    tryCatch({
      glm(Freq ~ Genotype * gene, data = x, family = Gamma(link = "log")) %>%
        emmeans(~ Genotype | gene) %>%
        pairs() %>% summary() %>% as.data.frame() %>%
        filter(!is.na(p.value)) %>%
        mutate(gene_type = as.character(x$gene_type[1]))
    }, error = function(e) NULL)
  }) %>%
  bind_rows()

# Multiple-testing correction across ALL gene contrasts (Benjamini–Hochberg)
stats <- stats %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  left_join(prop_full %>% select(gene, chain, gene_type) %>%
              mutate(gene_type = as.character(gene_type)) %>% distinct(),
            by = c("gene", "gene_type")) %>%
  arrange(padj)

write.csv(stats, file.path(sdata_dir, "fig6_gene_usage_stats.csv"), row.names = FALSE)
write.csv(exclusive_genes,
          file.path(sdata_dir, "fig6_gene_usage_exclusive.csv"), row.names = FALSE)

sig_genes <- stats %>% filter(padj < 0.05)
message("  ", sum(stats$p.value < 0.05, na.rm = TRUE), " raw p<0.05; ",
        nrow(sig_genes), " significant after BH FDR<0.05.")

# Genes shown in heatmaps: FDR-significant + genotype-exclusive
heatmap_genes <- bind_rows(
  sig_genes        %>% select(chain, gene) %>% mutate(flag = "FDR sig"),
  exclusive_genes  %>% select(chain, gene, flag = exclusive)
) %>% distinct(chain, gene, .keep_all = TRUE)

chains_present <- intersect(c("IGH", "IGK", "IGL"), unique(heatmap_genes$chain))

col_ha <- HeatmapAnnotation(
  Genotype = ifelse(startsWith(samples, "W"), "WT", "MUT"),
  col = list(Genotype = c(WT = "#2166AC", MUT = "#D6604D"))
)

for (ch in chains_present) {

  genes_ch <- heatmap_genes %>% filter(chain == ch)
  if (nrow(genes_ch) == 0) next

  prop_mat <- prop_full %>%
    filter(chain == ch, gene %in% genes_ch$gene) %>%
    select(gene, Freq, Mouse) %>%
    pivot_wider(names_from = Mouse, values_from = Freq, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  col_order <- samples[samples %in% colnames(prop_mat)]
  prop_mat  <- prop_mat[, col_order, drop = FALSE]

  # Row annotations: segment + significance flag
  ann <- prop_full %>% filter(chain == ch) %>%
    select(gene, gene_type) %>% distinct() %>%
    right_join(genes_ch %>% select(gene, flag), by = "gene") %>%
    distinct(gene, .keep_all = TRUE) %>%
    column_to_rownames("gene")
  ann <- ann[rownames(prop_mat), , drop = FALSE]

  row_ha <- rowAnnotation(
    Segment = as.character(ann$gene_type),
    Flag    = ann$flag,
    col = list(
      Segment = setNames(viridis(length(unique(na.omit(ann$gene_type))), option = "A"),
                         unique(na.omit(as.character(ann$gene_type)))),
      Flag    = c("FDR sig" = "grey20", "WT-only" = "#2166AC", "MUT-only" = "#D6604D")
    )
  )

  ht_h <- max(5, nrow(prop_mat) * 0.15 + 2)

  ht_raw <- Heatmap(prop_mat,
    name               = "Proportion",
    col                = colorRamp2(c(0, max(prop_mat)), c("white", "#1A6B3C")),
    top_annotation     = col_ha,
    left_annotation    = row_ha,
    cluster_columns    = FALSE,
    row_names_gp       = gpar(fontsize = 7),
    column_title       = paste0(ch, " — gene usage (raw proportions)"),
    heatmap_legend_param = list(title = "Proportion")
  )

  prop_scaled <- t(scale(t(prop_mat)))
  prop_scaled[is.nan(prop_scaled)] <- 0

  ht_scaled <- Heatmap(prop_scaled,
    name               = "z-score",
    col                = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D")),
    top_annotation     = col_ha,
    left_annotation    = row_ha,
    cluster_columns    = FALSE,
    row_names_gp       = gpar(fontsize = 7),
    column_title       = paste0(ch, " — gene usage (z-scored)"),
    heatmap_legend_param = list(title = "z-score")
  )

  pdf(file.path(plots_dir, paste0("fig6d_gene_usage_", ch, "_raw.pdf")),
      width = 7, height = ht_h)
  draw(ht_raw)
  dev.off()

  pdf(file.path(plots_dir, paste0("fig6d_gene_usage_", ch, "_scaled.pdf")),
      width = 7, height = ht_h)
  draw(ht_scaled)
  dev.off()

  message("  Saved heatmaps for chain: ", ch, " (", nrow(prop_mat), " genes)")
}

message("\nDone. Fig6 outputs saved under output/fig6/")
