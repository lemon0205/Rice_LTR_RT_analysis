# WGBS analyses

This module calculates TE-level CHH methylation contrasts and derives centromeric/pericentromeric small-RNA summary tables.

## Inputs

- CHH coverage files and TE BED annotations.
- DeepTools matrices for MSU 21–22-nt and 24-nt small-RNA signals.
- Centromere coordinate file.

## Scripts

- `compute_TE_CHH_loss_from_cov.R`: primary TE-level CHH coverage summary.
- `recompute_sRNA_CHH_loss_from_cov.R`: optional coverage-based recalculation pathway.
- `extract_MSU_centromeric_sRNA_profiles.py`: quantitative profile and context-table extraction.

Set `POLIV_ANALYSIS_DIR`, `DEEPTOOLS_MATRIX_DIR`, `CENTROMERE_BED`, and the centromeric output variables in the local path configuration before use.
