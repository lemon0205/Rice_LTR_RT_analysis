# TE annotation and landscape analyses

This module derives intact-LTR counts and TE-level presence, insertion-age, retention, and distance-to-gene summary tables used in the early manuscript results.

## Inputs

- Per-accession TE GFF/GTF annotations.
- MSU gene BED annotation.
- MSU-anchored TE genotype matrix.

## Script order

1. `extract_intact_ltr_counts.py`
2. `build_ltr_time_presence_distance.py`
3. `build_te_distance_frequency.py`
4. `build_te_retention_by_distance.py`
5. `analyze_ltr_insertion_age_group_diff.R`

The Python scripts expose their required inputs through command-line arguments. The final R script consumes the generated insertion-age tables. No figure-generation code is included.
