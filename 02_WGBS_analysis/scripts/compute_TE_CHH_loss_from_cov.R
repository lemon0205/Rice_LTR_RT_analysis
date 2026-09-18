suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

base_dir <- Sys.getenv("POLIV_ANALYSIS_DIR")
if (!nzchar(base_dir)) stop("Set POLIV_ANALYSIS_DIR to the directory containing GSE131319/ and diff/.")
cov_dir <- file.path(base_dir, "GSE131319/bw/CHH_cov")
bed_dir <- file.path(base_dir, "GSE131319/bw/bed")
out_dir <- file.path(base_dir, "diff/Scatter_plot")

min_chh_coverage <- 1
chunk_size <- 2000000

read_te_bed <- function(path, family) {
  bed <- fread(path, header = FALSE)
  setnames(bed, c("chr", "start0", "end", "TE_id", "score", "strand"))
  bed[!strand %in% c("+", "-"), strand := "*"]
  bed[, `:=`(
    start1 = start0 + 1L,
    family = family,
    TE_length = end - start0
  )]
  bed
}

te <- rbindlist(list(
  read_te_bed(file.path(bed_dir, "MSU_TE_Gypsy.bed"), "Gypsy"),
  read_te_bed(file.path(bed_dir, "MSU_TE_Copia.bed"), "Copia")
), use.names = TRUE)
te[, te_index := .I]

te_gr <- GRanges(
  seqnames = te$chr,
  ranges = IRanges(start = te$start1, end = te$end),
  strand = te$strand
)
mcols(te_gr)$te_index <- te$te_index

summarize_cov_for_te <- function(cov_file, sample_name) {
  message("Reading ", sample_name, ": ", basename(cov_file))

  n_te <- nrow(te)
  meth_sum <- numeric(n_te)
  total_sum <- numeric(n_te)
  chunk_i <- 0L

  con <- gzfile(cov_file, "rt")
  on.exit(close(con), add = TRUE)

  repeat {
    chunk <- tryCatch(
      read.table(
        con,
        sep = "\t",
        header = FALSE,
        nrows = chunk_size,
        quote = "",
        comment.char = "",
        colClasses = c("character", "integer", "integer", "numeric", "numeric", "numeric")
      ),
      error = function(e) {
        if (grepl("no lines available", conditionMessage(e), fixed = TRUE)) {
          return(NULL)
        }
        stop(e)
      }
    )

    if (is.null(chunk) || nrow(chunk) == 0) {
      break
    }

    chunk_i <- chunk_i + 1L
    if (chunk_i %% 10L == 0L) {
      message("  processed chunks: ", chunk_i)
    }

    setDT(chunk)
    setnames(chunk, c("chr", "start", "end", "meth_pct", "meth", "unmeth"))
    chunk[, total := meth + unmeth]
    chunk <- chunk[total >= min_chh_coverage]
    if (nrow(chunk) == 0) {
      next
    }

    for (chr_i in unique(chunk$chr)) {
      chr_chunk <- chunk[chr == chr_i]
      te_chr <- te_gr[as.character(seqnames(te_gr)) == chr_i]
      if (nrow(chr_chunk) == 0 || length(te_chr) == 0) {
        next
      }

      cyt_gr <- GRanges(
        seqnames = chr_i,
        ranges = IRanges(start = chr_chunk$start, end = chr_chunk$end)
      )
      hits <- findOverlaps(cyt_gr, te_chr, ignore.strand = TRUE)
      if (length(hits) == 0) {
        next
      }

      te_idx <- mcols(te_chr)$te_index[subjectHits(hits)]
      hit_dt <- data.table(
        te_index = te_idx,
        meth = chr_chunk$meth[queryHits(hits)],
        total = chr_chunk$total[queryHits(hits)]
      )
      agg <- hit_dt[, .(meth = sum(meth), total = sum(total)), by = te_index]
      meth_sum[agg$te_index] <- meth_sum[agg$te_index] + agg$meth
      total_sum[agg$te_index] <- total_sum[agg$te_index] + agg$total
    }
  }

  data.table(
    te_index = seq_len(n_te),
    methylated_reads = meth_sum,
    total_CHH_reads = total_sum,
    CHH_methylation_pct = fifelse(total_sum > 0, 100 * meth_sum / total_sum, NA_real_)
  )
}

cov_files <- list(
  WT_rep1 = file.path(cov_dir, "CHH_WT_base_BSseq180209_rep1_1.clean.bedgraph.gz.bismark.cov.gz"),
  WT_rep2 = file.path(cov_dir, "CHH_WT_base_BSseq180209_rep2_1.clean.bedgraph.gz.bismark.cov.gz"),
  nrpd1ab_rep1 = file.path(cov_dir, "CHH_osnrpd1ab_base_BSseq180209_rep1_1.clean.bedgraph.gz.bismark.cov.gz"),
  nrpd1ab_rep2 = file.path(cov_dir, "CHH_osnrpd1ab_base_BSseq180209_rep2_1.clean.bedgraph.gz.bismark.cov.gz")
)

sample_tables <- lapply(names(cov_files), function(sample_name) {
  x <- summarize_cov_for_te(cov_files[[sample_name]], sample_name)
  setnames(
    x,
    old = c("methylated_reads", "total_CHH_reads", "CHH_methylation_pct"),
    new = paste0(sample_name, c("_methylated_reads", "_total_CHH_reads", "_CHH_pct"))
  )
  x
})

out <- te[, .(chr, start = start0, end, TE_id, family, strand, TE_length, te_index)]
for (sample_table in sample_tables) {
  out <- merge(out, sample_table, by = "te_index", all.x = TRUE, sort = FALSE)
}

out[, WT_mean_CHH_pct := rowMeans(.SD, na.rm = TRUE), .SDcols = c("WT_rep1_CHH_pct", "WT_rep2_CHH_pct")]
out[, nrpd1ab_mean_CHH_pct := rowMeans(.SD, na.rm = TRUE), .SDcols = c("nrpd1ab_rep1_CHH_pct", "nrpd1ab_rep2_CHH_pct")]
out[, CHH_loss_WT_minus_nrpd1ab := WT_mean_CHH_pct - nrpd1ab_mean_CHH_pct]
out[, has_CHH_coverage := is.finite(WT_mean_CHH_pct) & is.finite(nrpd1ab_mean_CHH_pct)]

setcolorder(out, c(
  "chr", "start", "end", "TE_id", "family", "strand", "TE_length",
  "WT_rep1_CHH_pct", "WT_rep2_CHH_pct",
  "nrpd1ab_rep1_CHH_pct", "nrpd1ab_rep2_CHH_pct",
  "WT_mean_CHH_pct", "nrpd1ab_mean_CHH_pct", "CHH_loss_WT_minus_nrpd1ab",
  "has_CHH_coverage"
))
out[, te_index := NULL]

output_file <- file.path(out_dir, "TE_CHH_loss_from_cov.tsv")
fwrite(out, output_file, sep = "\t", quote = FALSE, na = "NA")

summary_file <- file.path(out_dir, "TE_CHH_loss_from_cov_summary.tsv")
summary_tbl <- out[has_CHH_coverage == TRUE, .(
  n_TE = .N,
  median_WT_CHH_pct = median(WT_mean_CHH_pct, na.rm = TRUE),
  median_nrpd1ab_CHH_pct = median(nrpd1ab_mean_CHH_pct, na.rm = TRUE),
  median_CHH_loss = median(CHH_loss_WT_minus_nrpd1ab, na.rm = TRUE),
  mean_CHH_loss = mean(CHH_loss_WT_minus_nrpd1ab, na.rm = TRUE)
), by = family]
fwrite(summary_tbl, summary_file, sep = "\t", quote = FALSE)

message("Wrote: ", output_file)
message("Wrote: ", summary_file)
