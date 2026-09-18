# Items requiring author confirmation before public release

No original files were moved or deleted. The items below are deliberately not copied into this analysis-only release until their status is confirmed.

## Multiple or mixed-purpose versions

| Topic | Candidate source files | Current treatment | Why confirmation is needed |
|---|---|---|---|
| dTE FST | `Figure7/pan-TE/output/dte_fst_calc.R`; `Figure7/pan-TE/output/dte_fst_analysis.R` | Retained `dte_fst_calc.R` only | `dte_fst_analysis.R` contains both calculation and PDF rendering; the retained script is the calculation-only version. Confirm it produced the submitted values. |
| Family-dynamics classification | `Figure6/.../17_define_context_aware_family_dynamics.R`; `Figure6/.../42_define_raw_observed_family_dynamics.R` | Neither copied | Both mix analysis with figure rendering and describe different primary frameworks. Confirm which submitted analysis underlies Figure 6 before splitting only its non-plot calculation portion. |
| CHH/sRNA loss | `compute_TE_CHH_loss_from_cov.R`; `recompute_sRNA_CHH_loss_from_cov.R`; `compute_TE_sRNA_log2ratio_from_bw.R` | All retained | The `recompute` script contains an optional large-coverage recalculation switch. Confirm whether the submitted table came from the direct coverage calculation or the recomputation path. |
| Figure S9 | `build_FigureS9_final_logic.R`; `build_FigureS9_replanned.R`; `build_FigureS9_TE_PAV_dTE_robustness.R`; `build_FigureS9_TE_SNP_convergence.R` | Not copied | These are combined analysis-and-plot scripts with explicit alternative/replanned versions. Only the non-plot GO list/enrichment helpers are included. |

## Software-version record

No lockfile, environment file, or command log recording exact software versions was found in the scanned script set. Do **not** invent an `environment.yml` or `requirements.txt`. Before submission, run and save `R --version`, `python --version`, `deepTools --version`, `samtools --version`, `bedtools --version`, and the R package versions from the original analysis environment; then add the verified values to the root README.
