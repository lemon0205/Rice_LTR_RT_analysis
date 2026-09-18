# Solo-LTR family analyses

This module contains the table construction, model fitting, family-level evidence integration, age-composition tests, and robustness analyses supporting the solo-LTR/family-dynamics results.

## Inputs

- Family-by-accession abundance table.
- Intact-LTR and solo-LTR element tables.
- Context and epigenetic summary tables from the processed data archive.

## Execution order

Run numbered scripts in ascending order where dependencies are present. The core sequence is `01_build_analysis_tables.R`, `02_fit_context_adjusted_turnover_model.R`, `03_build_family_evidence_table.R`, `04_test_expansion_context_association.R`, `06_decompose_GJ_XI_context_turnover.R`, and `08_integrate_family_trajectories.R`; higher-numbered scripts provide age, robustness, and epigenetic sensitivity analyses.

All model formulas, thresholds, and randomisation settings are retained in the source scripts. This module produces analysis tables only.
