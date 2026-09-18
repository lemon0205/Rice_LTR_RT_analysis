suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

base_dir <- Sys.getenv("POLIV_ANALYSIS_DIR")
if (!nzchar(base_dir)) stop("Set POLIV_ANALYSIS_DIR to the directory containing GSE131319/ and diff/.")
bw_dir <- file.path(base_dir, "GSE131319/bw/24nt_siRNA_bw")
out_dir <- file.path(base_dir, "diff/Scatter_plot")
chh_file <- file.path(out_dir, "TE_CHH_loss_from_cov.tsv")

read_te_from_chh <- function(path) {
  te <- fread(path, select = c("chr", "start", "end", "TE_id", "family", "strand", "TE_length"))
  te[!strand %in% c("+", "-"), strand := "*"]
  te[, te_index := .I]
  te
}

summarize_log2cpm_bw_for_te <- function(bw_file, sample_name, te, te_gr) {
  message("Reading ", sample_name, ": ", basename(bw_file))
  bw <- rtracklayer::import(BigWigFile(bw_file), which = te_gr)

  out <- data.table(
    te_index = te$te_index,
    signal_width = 0,
    weighted_signal_sum = 0
  )

  if (length(bw) > 0) {
    hits <- findOverlaps(bw, te_gr, ignore.strand = TRUE)
    if (length(hits) > 0) {
      overlaps <- pintersect(bw[queryHits(hits)], te_gr[subjectHits(hits)])
      hit_dt <- data.table(
        te_index = subjectHits(hits),
        score = mcols(bw)$score[queryHits(hits)],
        overlap_width = width(overlaps)
      )
      agg <- hit_dt[, .(
        signal_width = sum(overlap_width),
        weighted_signal_sum = sum(score * overlap_width)
      ), by = te_index]
      out[agg, `:=`(
        signal_width = i.signal_width,
        weighted_signal_sum = i.weighted_signal_sum
      ), on = "te_index"]
    }
  }

  # Missing bigWig bases are treated as 0, so divide by full TE length.
  out[, mean_log2CPM := weighted_signal_sum / te$TE_length]
  setnames(
    out,
    old = c("signal_width", "weighted_signal_sum", "mean_log2CPM"),
    new = paste0(sample_name, c("_sRNA_signal_width", "_sRNA_weighted_log2CPM_sum", "_sRNA_mean_log2CPM"))
  )
  out
}

te <- read_te_from_chh(chh_file)
te_gr <- GRanges(
  seqnames = te$chr,
  ranges = IRanges(start = te$start + 1L, end = te$end),
  strand = te$strand
)
mcols(te_gr)$te_index <- te$te_index

bw_files <- list(
  WT_rep1 = file.path(bw_dir, "WT_base_sRNA180209_rep1_log2cpm.bw"),
  WT_rep2 = file.path(bw_dir, "WT_base_sRNA180209_rep2_log2cpm.bw"),
  nrpd1ab_rep1 = file.path(bw_dir, "osnrpd1ab_base_sRNA180209_rep1_log2cpm.bw"),
  nrpd1ab_rep2 = file.path(bw_dir, "osnrpd1ab_base_sRNA180209_rep2_log2cpm.bw")
)

sample_tables <- lapply(names(bw_files), function(sample_name) {
  summarize_log2cpm_bw_for_te(bw_files[[sample_name]], sample_name, te, te_gr)
})

srna <- te[, .(chr, start, end, TE_id, family, strand, TE_length, te_index)]
for (sample_table in sample_tables) {
  srna <- merge(srna, sample_table, by = "te_index", all.x = TRUE, sort = FALSE)
}

srna[, WT_mean_sRNA_log2CPM := rowMeans(.SD, na.rm = TRUE),
     .SDcols = c("WT_rep1_sRNA_mean_log2CPM", "WT_rep2_sRNA_mean_log2CPM")]
srna[, nrpd1ab_mean_sRNA_log2CPM := rowMeans(.SD, na.rm = TRUE),
     .SDcols = c("nrpd1ab_rep1_sRNA_mean_log2CPM", "nrpd1ab_rep2_sRNA_mean_log2CPM")]
srna[, sRNA_log2_ratio_WT_over_nrpd1ab := WT_mean_sRNA_log2CPM - nrpd1ab_mean_sRNA_log2CPM]
srna[, has_sRNA_signal := rowSums(.SD > 0, na.rm = TRUE) > 0,
     .SDcols = c(
       "WT_rep1_sRNA_signal_width", "WT_rep2_sRNA_signal_width",
       "nrpd1ab_rep1_sRNA_signal_width", "nrpd1ab_rep2_sRNA_signal_width"
     )]

setcolorder(srna, c(
  "chr", "start", "end", "TE_id", "family", "strand", "TE_length",
  "WT_rep1_sRNA_mean_log2CPM", "WT_rep2_sRNA_mean_log2CPM",
  "nrpd1ab_rep1_sRNA_mean_log2CPM", "nrpd1ab_rep2_sRNA_mean_log2CPM",
  "WT_mean_sRNA_log2CPM", "nrpd1ab_mean_sRNA_log2CPM",
  "sRNA_log2_ratio_WT_over_nrpd1ab", "has_sRNA_signal"
))
srna[, te_index := NULL]

srna_file <- file.path(out_dir, "TE_sRNA_log2ratio_from_bw.tsv")
fwrite(srna, srna_file, sep = "\t", quote = FALSE, na = "NA")

chh <- fread(chh_file)
combined <- merge(
  chh,
  srna[, .(
    chr, start, end, TE_id,
    WT_rep1_sRNA_mean_log2CPM, WT_rep2_sRNA_mean_log2CPM,
    nrpd1ab_rep1_sRNA_mean_log2CPM, nrpd1ab_rep2_sRNA_mean_log2CPM,
    WT_mean_sRNA_log2CPM, nrpd1ab_mean_sRNA_log2CPM,
    sRNA_log2_ratio_WT_over_nrpd1ab, has_sRNA_signal
  )],
  by = c("chr", "start", "end", "TE_id"),
  all.x = TRUE,
  sort = FALSE
)

combined_file <- file.path(out_dir, "TE_CHH_loss_sRNA_log2ratio_combined.tsv")
fwrite(combined, combined_file, sep = "\t", quote = FALSE, na = "NA")

summary_tbl <- combined[, .(
  n_TE = .N,
  n_CHH_covered = sum(has_CHH_coverage == TRUE, na.rm = TRUE),
  n_with_sRNA_signal = sum(has_sRNA_signal == TRUE, na.rm = TRUE),
  median_CHH_loss = median(CHH_loss_WT_minus_nrpd1ab, na.rm = TRUE),
  median_sRNA_log2_ratio = median(sRNA_log2_ratio_WT_over_nrpd1ab, na.rm = TRUE),
  spearman_CHH_sRNA = suppressWarnings(cor(
    CHH_loss_WT_minus_nrpd1ab,
    sRNA_log2_ratio_WT_over_nrpd1ab,
    method = "spearman",
    use = "complete.obs"
  ))
), by = family]

summary_file <- file.path(out_dir, "TE_CHH_loss_sRNA_log2ratio_summary.tsv")
fwrite(summary_tbl, summary_file, sep = "\t", quote = FALSE, na = "NA")

message("Wrote: ", srna_file)
message("Wrote: ", combined_file)
message("Wrote: ", summary_file)
