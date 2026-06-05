# Quantify the effect of the fig6 correctness fixes (old vs new logic).
# Run from repo root: Rscript output/fig6/scripts/fig6_qc_test_fix_effects.R

suppressMessages({
  library(Seurat); library(dplyr); library(tidyr); library(purrr)
  library(emmeans); library(tibble)
})

samples    <- c("WT1","WT2","WT3","Mut1","Mut2","Mut3")
cell_order <- c("PreproB","ProB","PreB","Immature B","FoB","MzB")

contigs_df <- lapply(samples, function(s)
  read.csv(file.path("data", s, "filtered_contig_annotations.csv"),
           stringsAsFactors = FALSE) %>% mutate(Mouse = s)) %>%
  bind_rows() %>%
  filter(is_cell == TRUE, productive == TRUE, full_length == TRUE)

seu <- readRDS("data/bcells_annotated.rds")

sc_meta <- seu@meta.data %>%
  mutate(barcode = rownames(.),
         productive = ifelse(barcode %in% contigs_df$barcode, "Productive", "Unproductive"),
         Genotype = factor(Genotype, levels = c("WT","MUT"))) %>%
  filter(FinalLab3 %in% cell_order)

cat("\n=================  FIX 1: productive fraction by cell type  =================\n")

# OLD: drop all-unproductive mouse x celltype, average survivors
old <- sc_meta %>%
  group_by(Mouse, Genotype, FinalLab3, productive) %>% summarise(n = n(), .groups="drop") %>%
  group_by(Mouse, Genotype, FinalLab3) %>% mutate(frac = n/sum(n)) %>%
  filter(productive == "Productive") %>% ungroup() %>%
  group_by(Genotype, FinalLab3) %>%
  summarise(old_mean = mean(frac), old_sem = sd(frac)/sqrt(n()), old_n = n(), .groups="drop")

# NEW: keep zeros
new <- sc_meta %>%
  group_by(Mouse, Genotype, FinalLab3) %>%
  summarise(n_prod = sum(productive=="Productive"), n_total = n(), .groups="drop") %>%
  mutate(frac = n_prod/n_total) %>%
  group_by(Genotype, FinalLab3) %>%
  summarise(new_mean = mean(frac), new_sem = sd(frac)/sqrt(n()), new_n = n(), .groups="drop")

cmp <- full_join(old, new, by = c("Genotype","FinalLab3")) %>%
  mutate(FinalLab3 = factor(FinalLab3, cell_order)) %>%
  arrange(Genotype, FinalLab3) %>%
  mutate(across(c(old_mean,old_sem,new_mean,new_sem), ~round(.*100,1)),
         delta_mean_pp = round(new_mean - old_mean, 1))
print(as.data.frame(cmp), row.names = FALSE)
cat("\n(values in %, delta_mean_pp = new - old in percentage points; ",
    "old_n/new_n = mice averaged)\n")

cat("\n=================  FIX 3: multiple-testing correction  =================\n")

contigs_filter <- contigs_df %>%
  filter(barcode %in% colnames(seu)) %>%
  left_join(seu@meta.data %>% mutate(barcode = rownames(.)) %>% select(barcode, Genotype),
            by = "barcode")
contigs_split <- contigs_filter %>% group_by(Mouse) %>% group_split() %>%
  setNames(sapply(., \(x) x$Mouse[1]))

prop_obs <- imap(contigs_split, function(x, nm) {
  x %>% select(v_gene,d_gene,j_gene,c_gene,chain) %>%
    mutate(across(c(v_gene,d_gene,j_gene,c_gene), ~ifelse(.=="","NoGene",.))) %>%
    group_by(chain) %>% group_split() %>% setNames(sapply(., \(d) d$chain[1])) %>%
    imap(function(d, ch) lapply(c("v_gene","d_gene","j_gene","c_gene"), function(seg){
      tbl <- table(d[[seg]]); as.data.frame(tbl/sum(tbl)) %>%
        rename(gene=Var1, Freq=Freq) %>% mutate(chain=ch, gene_type=seg)
    }) %>% bind_rows()) %>% bind_rows() %>% mutate(Mouse = nm)
}) %>% bind_rows() %>% filter(gene != "NoGene") %>% mutate(gene = as.character(gene))

prop_full <- prop_obs %>% group_by(chain, gene_type, gene) %>%
  complete(Mouse = samples, fill = list(Freq = 0)) %>% ungroup() %>%
  mutate(Genotype = ifelse(startsWith(Mouse,"W"),"WT","MUT"))

stats <- prop_full %>% filter(Freq > 0) %>% group_by(gene_type) %>% group_split() %>%
  map(function(x) tryCatch(
    glm(Freq ~ Genotype*gene, data=x, family=Gamma(link="log")) %>%
      emmeans(~Genotype|gene) %>% pairs() %>% summary() %>% as.data.frame() %>%
      filter(!is.na(p.value)) %>% mutate(gene_type = as.character(x$gene_type[1])),
    error=function(e) NULL)) %>% bind_rows() %>%
  mutate(padj = p.adjust(p.value, "BH"))

cat("Total gene contrasts tested:        ", nrow(stats), "\n")
cat("Significant at raw p < 0.05:        ", sum(stats$p.value < 0.05), "\n")
cat("Significant at BH FDR < 0.05:       ", sum(stats$padj < 0.05), "\n")

cat("\n=================  FIX 2: genotype-exclusive genes recovered  =============\n")
gene_geno <- prop_full %>% group_by(chain, gene_type, gene, Genotype) %>%
  summarise(total = sum(Freq), .groups="drop") %>%
  pivot_wider(names_from = Genotype, values_from = total, values_fill = 0)
excl <- gene_geno %>% filter((WT==0) != (MUT==0))
cat("Genotype-exclusive genes (dropped by old code, now recovered): ", nrow(excl), "\n")
cat("  WT-only: ", sum(excl$MUT==0), " | MUT-only: ", sum(excl$WT==0), "\n")

cat("\n=================  FIX 5: empty clonotype IDs  =============================\n")
n_empty <- contigs_filter %>%
  filter(is.na(raw_clonotype_id) | raw_clonotype_id %in% c("","None")) %>% nrow()
cat("Contig rows with empty/None clonotype_id (now excluded): ", n_empty, "\n\n")
