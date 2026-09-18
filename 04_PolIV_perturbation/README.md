# Pol IV perturbation analysis

This module prepares wild-type and `osnrpd1ab` difference tracks used as inputs to downstream quantitative analyses.

Run:

```bash
bash scripts/build_mutant_minus_WT_bigwigs.bash <chh_bigwig_dir> <small_rna_bigwig_dir> <output_dir> [threads]
```

The script retains the original replicate structure, a 50-bp bin size, and the WT-minus-mutant subtraction. It produces mean and contrast BigWig tracks only; plotting is intentionally outside the repository scope.
