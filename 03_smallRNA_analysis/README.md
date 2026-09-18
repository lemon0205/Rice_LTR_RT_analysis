# Small-RNA analyses

This module quantifies TE-associated 21–22-nt and 24-nt siRNA signals and tests their relationships with methylation, TE length, and genomic context.

## Inputs

- 21–22-nt and 24-nt BigWigs or deepTools matrices.
- `15_TE_signals_tissue_averaged.tsv` from the processed data archive.

## Script order

1. `compute_TE_sRNA_log2ratio_from_bw.R` when TE-level BigWig summaries are required.
2. `extract_TE_coverage_metric_columns.py` and `extract_24nt_boundary_enrichment.py` for matrix-derived summaries.
3. `analyze_deeptools_matrix_spatial_peaks.py` for spatial metrics and matrix metadata.
4. `analyze_21_22nt_integrated_evidence.R` and `analyze_21_22nt_length_interaction.R` for the final quantitative tests.

Set `TE_SIGNAL_TABLE`, `TE_SIGNAL_OUTDIR`, and `DEEPTOOLS_MATRIX_DIR` before execution.
