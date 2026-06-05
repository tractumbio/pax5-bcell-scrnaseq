# Data

The input data for this analysis are **not stored in git** (they are large
binary files). They are archived on Zenodo and downloaded automatically by
`main_script.R` on the first run.

> **Download:** https://doi.org/10.5281/zenodo.20553303

To set up manually instead, download the archive from the DOI above and place
the files in this `data/` directory so it matches the manifest below.

## File manifest

| File | Description | Used by |
|------|-------------|---------|
| `bcells_raw_filtered.rds` | Raw, quality-filtered query Seurat object (WT + MUT B cells) | data processing |
| `reference_atlas_lee.rds` | Published mouse B-cell reference atlas (Lee et al., *Nat. Commun.*) | data processing |
| `bcells_annotated.rds` | Fully annotated query object (FinalLab/2/3). Produced by the data-processing script; also provided pre-computed on Zenodo so figures can be run without re-running Stage 1 | all figure scripts |
| `marker_modules.rds` | Named list of marker-gene modules for B-cell subtypes | data processing |
| `pax5_targets.rds` | Named list: genes activated and genes repressed by Pax5 | Figure 4e |
| `vdj_contigs.rds` | Combined 10x Cell Ranger filtered BCR contig annotations for all 6 mice (WT1–3, Mut1–3), with a `Mouse` column | Figure 6 |
| `msigdb_genesets.rds` | The four MSigDB v2026.1 mouse gene-set collections used here — `hallmark`, `reactome`, `gobp`, `tf` (gene-symbol). A named list of named gene-set lists. Derived from MSigDB and subject to the [MSigDB terms of use](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) | Figure 5b |

## Notes

- `bcells_annotated.rds` uses standardised stage labels in `FinalLab3`:
  `PreproB, ProB, PreB, Immature B, FoB, MzB`.
- The per-sample Cell Ranger VDJ folders (`WT1/ … Mut3/`) are optional — the
  combined `vdj_contigs.rds` supersedes them; Figure 6 falls back to the
  per-sample CSVs only if the combined file is absent.
