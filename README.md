# Rice LTR-RT, RdDM, and TE-PAV analysis code

This repository is an analysis-only public release supporting the main results of a rice LTR retrotransposon (LTR-RT), RNA-directed DNA methylation (RdDM), and transposable-element presence/absence variation (TE-PAV) study. It intentionally excludes figure-generation, layout, and rendering code, as well as raw sequencing data, large processed matrices, caches, and local paths.

The release was prepared as a copy from the working project. It does not change the original project or modify the scientific logic or analysis parameters in the retained scripts. The only code edits replace hard-coded local paths with environment variables or positional command-line arguments, and remove a mixed plotting section from the Pol IV track-construction shell script.

## Repository layout

```text
01_TE_annotation_and_landscape/scripts/  Intact-LTR extraction and TE age/distance summary tables
02_WGBS_analysis/scripts/                CHH/centromeric quantitative summaries
03_smallRNA_analysis/scripts/            21–22-nt and 24-nt siRNA quantitative analyses
04_PolIV_perturbation/scripts/           osnrpd1ab minus WT signal-track construction
05_soloLTR_family_analysis/scripts/      Solo-LTR/family turnover, robustness, and statistics
06_TE_PAV_population_analysis/scripts/   TE-PAV matrix, dTE FST, regions, and GO inputs
config/                                  Local-path template (not data)
docs/                                    Provenance and release-review notes
```

## Data access and placement

Raw sequencing reads are **not** included and must never be committed to this repository.

| Data source | Use in this repository | Expected local placement |
|---|---|---|
| GSA `PRJCA071954` | Study-generated sequencing data used in methylation/small-RNA related analyses | Download under a local `data/raw/PRJCA071954/` directory; generate alignment and signal files outside version control. |
| GEO `GSE131319` | Public Pol IV perturbation data used for WT versus `osnrpd1ab` CHH and small-RNA comparisons | Download under a local `data/raw/GSE131319/`; derived BigWig files are supplied to scripts through `POLIV_ANALYSIS_DIR` or shell arguments. |
| Zenodo (processed archive) | Large derived TE signal tables, deepTools matrices, coverage files, TE-PAV matrices, and result tables | Download the manuscript-associated archive after its DOI is finalized; place files in a local `data/processed/` directory and point the environment variables below to their actual locations. |

Required filename patterns are preserved in scripts. Key examples are `*_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz`, `15_TE_signals_tissue_averaged.tsv`, `TE_CHH_loss_from_cov.tsv`, `te_genotype_matrix_strict_noNA.csv`, and `dte_fst_results.csv`.

## Configuration and execution

Copy `config/paths.env.example` to a local untracked file, fill in paths to downloaded/derived data, and source it before running the scripts. All placeholders are intentionally generic; no personal path, account name, credential, token, or SSH information is present in this release.

```bash
source config/paths.env.local
python 03_smallRNA_analysis/scripts/analyze_deeptools_matrix_spatial_peaks.py
Rscript 03_smallRNA_analysis/scripts/analyze_21_22nt_integrated_evidence.R
```

The retained scripts preserve their original thresholds and model definitions. Run modules only after their stated input tables have been generated or obtained from Zenodo:

1. `01_TE_annotation_and_landscape`: extract intact LTR counts and construct TE presence/gene-distance tables.
2. `02_WGBS_analysis`: calculate TE-level CHH loss and extract centromeric small-RNA matrix summaries.
3. `03_smallRNA_analysis`: calculate TE-level siRNA ratios, extract deepTools metrics, and fit integrated 21–22-nt/24-nt analyses.
4. `04_PolIV_perturbation`: run `build_mutant_minus_WT_bigwigs.bash <chh_bw_dir> <srna_bw_dir> <outdir> [threads]` to construct the six retained signal tracks.
5. `05_soloLTR_family_analysis`: run scripts in numeric order where applicable; their original working-directory assumptions are preserved.
6. `06_TE_PAV_population_analysis`: generate the TE matrix, calculate dTE FST, merge FST=1 regions, then export gene lists and GO inputs.

## Figure-result linkage

No drawing scripts are included. The following map identifies the retained analysis scripts that generate data/statistics used by the corresponding figures; final rendering occurred outside this release.

| Result | Retained analysis modules |
|---|---|
| Figure 1 | `01_TE_annotation_and_landscape` intact-LTR extraction and distance/presence summaries |
| Figure 2 | `01_TE_annotation_and_landscape` LTR time, TE retention, and gene-distance summaries |
| Figure 3 | `02_WGBS_analysis` and `03_smallRNA_analysis` deepTools/siRNA/CHH quantitative tables |
| Figure 4 | `03_smallRNA_analysis` and `04_PolIV_perturbation` quantitative signal analyses |
| Figure 5 | Pol IV perturbation quantitative summaries in modules 02–04 |
| Figure 6 | `05_soloLTR_family_analysis` family turnover, context, age, and CHH statistics |
| Figure 7 | `06_TE_PAV_population_analysis` TE-PAV, dTE FST, hotspot, distance, and GO-input analyses |
| Supplementary figures | Modules 02, 03, 05, and 06; see `docs/SCRIPT_PROVENANCE.tsv` and `docs/REVIEW_NEEDED.md` |

## Software dependencies

Dependencies are reported only when they are imported or invoked by retained code: R (with `data.table`, `GenomicRanges`, `rtracklayer`, `glmmTMB`, and `lme4` where called); Python (with `numpy` and `pandas` where called); and deepTools (`bigwigAverage`, `bigwigCompare`) for Pol IV tracks. Exact versions were not recorded in the scanned files and are therefore intentionally not guessed. See `docs/REVIEW_NEEDED.md` for the required version-capture step before submission.

## Scope and provenance

`docs/SCRIPT_CATALOG.md` records input and output types for every retained module. `docs/SCRIPT_PROVENANCE.tsv` maps each copied script to its original project path. `docs/REVIEW_NEEDED.md` lists duplicated or mixed analysis/plot versions that require author confirmation. Temporary files, failed/replanned versions, rendered figures, cached bytecode, intermediate data, and personal configuration are excluded.
