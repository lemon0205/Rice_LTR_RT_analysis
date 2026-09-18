# TE-PAV population analyses

This module generates the TE presence/absence genotype matrix, calculates dTE differentiation, identifies high-FST regions, evaluates gene proximity, and prepares GO-slim analysis inputs.

## Inputs

- Per-accession TE annotations and PAV/synteny files.
- MSU gene BED annotation.
- dTE FST results and GO-slim annotation file.

## Core order

1. `generate_te_matrix_synteny.py`
2. `dte_fst_calc.R`
3. `build_fst1_regions.py`
4. `build_te_fate_group_distance.py`
5. `te_nearby_gene_for_go.py` and `export_te_four_groups_go_lists.py`
6. `build_four_TE_class_2kb_gene_lists.R` and `run_four_TE_class_GOSlim_enrichment.R`

Set `DTE_FST_RESULTS`, `FST1_REGION_OUTDIR`, `PROJECT_DATA_ROOT`, and `GO_2KB_DIR` as appropriate for the later steps. The module produces tabular analysis outputs; figure rendering is not included.
