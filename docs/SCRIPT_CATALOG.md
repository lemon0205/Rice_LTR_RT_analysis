# Analysis script catalog

This catalog summarizes the analysis scripts distributed with this repository. Exact filenames, thresholds, sample sets, and model formulas are defined in the source scripts.

| Module | Scripts | Principal inputs | Principal analysis outputs |
|---|---|---|---|
| TE annotation and landscape | `extract_intact_ltr_counts.py`; `build_ltr_time_presence_distance.py`; `build_te_distance_frequency.py`; `build_te_retention_by_distance.py`; `analyze_ltr_insertion_age_group_diff.R` | TE GFF/GTF annotation, MSU gene BED, and MSU-anchored TE genotype matrix | Intact-LTR count tables, TE presence/identity/distance tables, frequency/retention summaries, and group-difference statistics |
| WGBS | `compute_TE_CHH_loss_from_cov.R`; `recompute_sRNA_CHH_loss_from_cov.R`; `extract_MSU_centromeric_sRNA_profiles.py` | CHH coverage files, TE BEDs, deepTools matrices, and centromere coordinates | TE-level CHH loss table; centromeric context annotation, profile summary, and TE count tables |
| Small RNA | `compute_TE_sRNA_log2ratio_from_bw.R`; `extract_24nt_boundary_enrichment.py`; `extract_TE_coverage_metric_columns.py`; `analyze_deeptools_matrix_spatial_peaks.py`; `analyze_21_22nt_integrated_evidence.R`; `analyze_21_22nt_length_interaction.R` | 21–22-nt/24-nt BigWigs or deepTools matrices and `15_TE_signals_tissue_averaged.tsv` | TE-level siRNA ratios, boundary/spatial metrics, matrix metadata, integrated-evidence tables, and length-interaction statistics |
| Pol IV perturbation | `build_mutant_minus_WT_bigwigs.bash` | GSE131319 WT and `osnrpd1ab` CHH and small-RNA BigWig replicates | WT mean, mutant mean, and WT-minus-mutant CHH/small-RNA BigWigs |
| Solo-LTR family analysis | `01_...R`, `02_...R`, `03_...R`, `04_...R`, `06_...R`, `08_...R`, `15_...R`, `22_...R`, `23_...R`, `38_...R`, `39_...R`, `47_...R`, `50_...R` | Family-by-accession table, intact/solo LTR element tables, and quantitative epigenetic tables | Eligibility tables, context-adjusted turnover model results, evidence and robustness tables, age-composition tests, and CHH/dynamic-class statistics |
| TE-PAV population analysis | `generate_te_matrix_synteny.py`; `dte_fst_calc.R`; `build_fst1_regions.py`; `build_te_fate_group_distance.py`; `te_nearby_gene_for_go.py`; `export_te_four_groups_go_lists.py`; `build_four_TE_class_2kb_gene_lists.R`; `run_four_TE_class_GOSlim_enrichment.R` | Per-accession TE/PAV annotations, MSU gene BED, dTE FST table, and GO-slim assignment | TE genotype matrix, dTE FST table, FST=1 regions, gene-distance summaries, nearby-gene and GO input/result tables |

Figure-generation, rendering, and layout scripts are intentionally excluded. This repository distributes analysis code and tabular/statistical result generation only.
