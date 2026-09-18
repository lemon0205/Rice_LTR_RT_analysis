# Computational analyses for rice LTR-RT, RdDM, and TE-PAV

This repository contains the computational analysis scripts supporting a study of long terminal repeat retrotransposons (LTR-RTs), RNA-directed DNA methylation (RdDM), and transposable-element presence/absence variation (TE-PAV) in rice. The repository is prepared for public release alongside the manuscript and is limited to the code required to generate the quantitative analysis tables and statistical results underlying the principal and supplementary results.

Figure rendering, graphic layout, raw sequencing reads, alignment files, BigWig tracks, deepTools matrices, and other large intermediate files are intentionally not distributed here.

## Repository contents

```text
01_TE_annotation_and_landscape/   LTR-RT annotation, insertion-age, retention, and gene-distance analyses
02_WGBS_analysis/                 TE-level CHH methylation and centromeric signal summaries
03_smallRNA_analysis/             21–22-nt and 24-nt siRNA quantitative analyses
04_PolIV_perturbation/            Pol IV mutant versus wild-type signal-track preparation
05_soloLTR_family_analysis/       Solo-LTR, family-dynamics, and epigenetic statistical analyses
06_TE_PAV_population_analysis/    TE-PAV matrix, dTE FST, hotspot, and GO analyses
config/                           Example local-path configuration
docs/                             Script catalog and repository documentation
```

Each module contains a module-specific README describing its analytical scope, required inputs, expected outputs, and script order.

## Data availability

Raw sequencing data are not included in this repository.

| Dataset | Accession or archive | Role in the analyses |
|---|---|---|
| Study-generated sequencing data | GSA `PRJCA071954` | Methylation and small-RNA analyses across rice accessions |
| Public Pol IV perturbation data | GEO `GSE131319` | Wild-type and `osnrpd1ab` methylation/small-RNA comparisons |
| Large processed analysis inputs and source tables | Zenodo deposition associated with the manuscript | DeepTools matrices, coverage files, TE-PAV matrices, and derived analysis tables |

After the Zenodo record is available, download its processed files to a local data directory and define the paths in `config/paths.env.example`. Large files remain excluded by `.gitignore`.

## Software requirements

The analyses use R, Python, and deepTools. R packages invoked by the retained scripts include `data.table`, `GenomicRanges`, `rtracklayer`, `glmmTMB`, and `lme4`; Python scripts use `numpy` and `pandas`; and the Pol IV track-preparation step uses `bigwigAverage` and `bigwigCompare` from deepTools.

Exact software versions should be recorded from the original computational environment before the archival release is finalized. No dependency lockfile is provided because an unverified environment specification would compromise reproducibility.

## Running the analyses

1. Obtain raw or processed inputs from the repositories listed above.
2. Copy `config/paths.env.example` to an untracked local configuration file and update each path.
3. Run the required module(s) according to their module README. Scripts retain the thresholds, sample definitions, and model specifications used in the analysis.

For example:

```bash
source config/paths.env.local
python 03_smallRNA_analysis/scripts/analyze_deeptools_matrix_spatial_peaks.py
Rscript 03_smallRNA_analysis/scripts/analyze_21_22nt_integrated_evidence.R
```

## Relationship to manuscript results

| Manuscript result | Supporting analysis module(s) |
|---|---|
| Figure 1 | `01_TE_annotation_and_landscape` |
| Figure 2 | `01_TE_annotation_and_landscape` |
| Figure 3 | `02_WGBS_analysis`, `03_smallRNA_analysis` |
| Figure 4 | `03_smallRNA_analysis`, `04_PolIV_perturbation` |
| Figure 5 | `02_WGBS_analysis`, `03_smallRNA_analysis`, `04_PolIV_perturbation` |
| Figure 6 | `05_soloLTR_family_analysis` |
| Figure 7 | `06_TE_PAV_population_analysis` |
| Supplementary figures | `02_WGBS_analysis`, `03_smallRNA_analysis`, `05_soloLTR_family_analysis`, and `06_TE_PAV_population_analysis` |

The script-level input and output catalog is available in `docs/SCRIPT_CATALOG.md`.

## Code availability

The complete analysis code used to generate the quantitative results reported in the manuscript is available from this repository. Raw data are available through GSA `PRJCA071954` and GEO `GSE131319`; large processed inputs and source data are available from the associated Zenodo record.

## License

This repository is distributed under the terms of the included [LICENSE](LICENSE).
