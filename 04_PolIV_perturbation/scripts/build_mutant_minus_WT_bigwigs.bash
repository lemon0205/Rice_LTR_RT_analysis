#!/usr/bin/env bash
# Build the mutant-minus-WT signal tracks used for quantitative TE summaries.
# This release contains signal-track construction only.
# Usage: bash build_mutant_minus_WT_bigwigs.bash <chh_bw_dir> <srna_bw_dir> <outdir> [threads]
set -euo pipefail

chh_bw_dir=$1
srna_bw_dir=$2
outdir=$3
threads=${4:-8}
mkdir -p "$outdir"

bigwigAverage -b "$chh_bw_dir/CHH_WT_base_BSseq180209_rep1_1.clean.bw" "$chh_bw_dir/CHH_WT_base_BSseq180209_rep2_1.clean.bw" -o "$outdir/CHH_WT_mean.bw" --binSize 50 --numberOfProcessors "$threads"
bigwigAverage -b "$chh_bw_dir/CHH_osnrpd1ab_base_BSseq180209_rep1_1.clean.bw" "$chh_bw_dir/CHH_osnrpd1ab_base_BSseq180209_rep2_1.clean.bw" -o "$outdir/CHH_osnrpd1ab_mean.bw" --binSize 50 --numberOfProcessors "$threads"
bigwigCompare -b1 "$outdir/CHH_WT_mean.bw" -b2 "$outdir/CHH_osnrpd1ab_mean.bw" --binSize 50 --operation subtract -o "$outdir/CHH_WT_minus_osnrpd1ab.bw" --numberOfProcessors "$threads"

bigwigAverage -b "$srna_bw_dir/WT_base_sRNA180209_rep1_log2cpm.bw" "$srna_bw_dir/WT_base_sRNA180209_rep2_log2cpm.bw" "$srna_bw_dir/WT_base_sRNA180209_rep3_log2cpm.bw" -o "$outdir/sRNA_WT_mean.bw" --binSize 50 --numberOfProcessors "$threads"
bigwigAverage -b "$srna_bw_dir/osnrpd1ab_base_sRNA180209_rep1_log2cpm.bw" "$srna_bw_dir/osnrpd1ab_base_sRNA180209_rep2_log2cpm.bw" "$srna_bw_dir/osnrpd1ab_base_sRNA180209_rep3_log2cpm.bw" -o "$outdir/sRNA_osnrpd1ab_mean.bw" --binSize 50 --numberOfProcessors "$threads"
bigwigCompare -b1 "$outdir/sRNA_WT_mean.bw" -b2 "$outdir/sRNA_osnrpd1ab_mean.bw" --binSize 50 --operation subtract -o "$outdir/sRNA_WT_minus_osnrpd1ab.bw" --numberOfProcessors "$threads"
