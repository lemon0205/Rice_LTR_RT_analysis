suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
})

base_dir <- Sys.getenv("POLIV_ANALYSIS_DIR")
if (!nzchar(base_dir)) stop("Set POLIV_ANALYSIS_DIR to the directory containing GSE131319/ and diff/.")
scatter_dir <- file.path(base_dir, "diff/Scatter_plot")
bw_dir <- file.path(base_dir, "GSE131319/bw")
cov_dir <- file.path(bw_dir, "CHH_cov")
bed_dir <- file.path(bw_dir, "bed")

# Set to TRUE only when you want to recompute CHH from the large cov.gz files.
RUN_COV_RECALC <- FALSE
min_chh_coverage <- 1

read_te_bed <- function(path, family) {
  bed <- fread(path, header = FALSE)
  setnames(bed, c("chr", "start0", "end", "TE_id", "score", "strand"))
  bed[, `:=`(
    start = start0 + 1L,
    family = family,
    te_index = .I
  )]
  gr <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = bed$start, end = bed$end),
    strand = bed$strand
  )
  mcols(gr)$TE_id <- bed$TE_id
  mcols(gr)$family <- bed$family
  mcols(gr)$te_index <- bed$te_index
  list(table = bed, gr = gr)
}

summarize_cov_for_te <- function(cov_file, te_gr, min_cov = 1, chunk_size = 2000000) {
  n_te <- length(te_gr)
  meth_sum <- numeric(n_te)
  total_sum <- numeric(n_te)

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
      error = function(e) NULL
    )

    if (is.null(chunk) || nrow(chunk) == 0) {
      break
    }

    setDT(chunk)
    setnames(chunk, c("chr", "start", "end", "meth_pct", "meth", "unmeth"))
    chunk[, total := meth + unmeth]
    chunk <- chunk[total >= min_cov]

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

  ifelse(total_sum > 0, 100 * meth_sum / total_sum, NA_real_)
}

summarize_log2cpm_bw_for_te <- function(bw_file, te_gr) {
  bw <- rtracklayer::import(BigWigFile(bw_file), which = te_gr)
  if (length(bw) == 0) {
    return(rep(NA_real_, length(te_gr)))
  }

  hits <- findOverlaps(bw, te_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(rep(NA_real_, length(te_gr)))
  }

  overlaps <- pintersect(bw[queryHits(hits)], te_gr[subjectHits(hits)])
  hit_dt <- data.table(
    te_index = subjectHits(hits),
    score = mcols(bw)$score[queryHits(hits)],
    width = width(overlaps)
  )
  agg <- hit_dt[, .(weighted_signal = sum(score * width) / sum(width)), by = te_index]
  out <- rep(0, length(te_gr))
  out[agg$te_index] <- agg$weighted_signal
  out
}

diagnose_existing_loss_tables <- function() {
  for (family in c("Gypsy", "Copia")) {
    tab_path <- file.path(scatter_dir, paste0(family, "_sRNA_CHH_loss.tab"))
    dat <- fread(tab_path)
    setnames(dat, gsub("^#|'|\"", "", names(dat)))
    message("\n", family, " existing table")
    print(summary(dat[, .(sRNA_loss, CHH_loss)]))
  }
}

diagnose_existing_loss_tables()

message("\nCurrent interpretation:")
message("- CHH_loss was generated from CHH bigWig means as WT - osnrpd1ab.")
message("- sRNA_loss was generated from files named *_log2cpm.bw, so it is a log2CPM-scale WT - osnrpd1ab difference.")
message("- For CHH, cov.gz read counts allow a more defensible weighted methylation estimate: sum(methylated reads) / sum(total CHH reads) per TE.")

if (RUN_COV_RECALC) {
  te_sets <- list(
    Gypsy = read_te_bed(file.path(bed_dir, "MSU_TE_Gypsy.bed"), "Gypsy"),
    Copia = read_te_bed(file.path(bed_dir, "MSU_TE_Copia.bed"), "Copia")
  )

  cov_files <- list(
    WT_rep1 = file.path(cov_dir, "CHH_WT_base_BSseq180209_rep1_1.clean.bedgraph.gz.bismark.cov.gz"),
    WT_rep2 = file.path(cov_dir, "CHH_WT_base_BSseq180209_rep2_1.clean.bedgraph.gz.bismark.cov.gz"),
    nrpd1ab_rep1 = file.path(cov_dir, "CHH_osnrpd1ab_base_BSseq180209_rep1_1.clean.bedgraph.gz.bismark.cov.gz"),
    nrpd1ab_rep2 = file.path(cov_dir, "CHH_osnrpd1ab_base_BSseq180209_rep2_1.clean.bedgraph.gz.bismark.cov.gz")
  )

  sRNA_files <- list(
    WT_rep1 = file.path(bw_dir, "24nt_siRNA_bw/WT_base_sRNA180209_rep1_log2cpm.bw"),
    WT_rep2 = file.path(bw_dir, "24nt_siRNA_bw/WT_base_sRNA180209_rep2_log2cpm.bw"),
    nrpd1ab_rep1 = file.path(bw_dir, "24nt_siRNA_bw/osnrpd1ab_base_sRNA180209_rep1_log2cpm.bw"),
    nrpd1ab_rep2 = file.path(bw_dir, "24nt_siRNA_bw/osnrpd1ab_base_sRNA180209_rep2_log2cpm.bw")
  )

  for (family in names(te_sets)) {
    te_gr <- te_sets[[family]]$gr
    te_tbl <- te_sets[[family]]$table

    chh <- lapply(cov_files, summarize_cov_for_te, te_gr = te_gr, min_cov = min_chh_coverage)
    srna <- lapply(sRNA_files, summarize_log2cpm_bw_for_te, te_gr = te_gr)

    out <- data.table(
      chr = as.character(seqnames(te_gr)),
      start = start(te_gr) - 1L,
      end = end(te_gr),
      TE_id = te_tbl$TE_id,
      sRNA_log2CPM_diff = rowMeans(cbind(srna$WT_rep1, srna$WT_rep2), na.rm = TRUE) -
        rowMeans(cbind(srna$nrpd1ab_rep1, srna$nrpd1ab_rep2), na.rm = TRUE),
      CHH_methylation_loss = rowMeans(cbind(chh$WT_rep1, chh$WT_rep2), na.rm = TRUE) -
        rowMeans(cbind(chh$nrpd1ab_rep1, chh$nrpd1ab_rep2), na.rm = TRUE)
    )

    fwrite(
      out,
      file.path(scatter_dir, paste0(family, "_sRNA_log2CPM_CHH_cov_loss.tsv")),
      sep = "\t"
    )
  }
}
