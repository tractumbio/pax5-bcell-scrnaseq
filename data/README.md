# Data

The `.rds` files required to run this analysis are not tracked in git due to their size.
Place the following files in this directory before running any scripts:

| File | Description |
|------|-------------|
| `anselm_demux_filtered.rds` | Raw demultiplexed query Seurat object (input to data processing) |
| `lee_etat_natcom_data.rds` | Published B-cell reference atlas (Lee et al., *Nature Communications*) |
| `anselm_labeled.rds` | Fully annotated query object — produced by `data_processing/label_transfer_and_annotation.R`, required by all figure scripts |
| `marker_modules.rds` | Named list of marker gene sets used to sub-classify Mature B clusters |
| `pax5_targets.rds` | Named list with two elements: genes activated and genes repressed by Pax5 |
| `nat_comm_markers.xlsx` | Marker gene table from the Lee et al. reference (supplementary) |
