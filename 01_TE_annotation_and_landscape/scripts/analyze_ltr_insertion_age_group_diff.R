#!/usr/bin/env Rscript

# Summarize the youngest LTR insertion-age bin and test XL/GL group differences.
# Run from this directory, or pass --indir/--outdir explicitly.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

indir <- get_arg("--indir", ".")
outdir <- get_arg("--outdir", indir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hist_file <- file.path(indir, "ltr_insertion_age_histogram.tsv")
ages_file <- file.path(indir, "ltr_insertion_ages.tsv")

hist_df <- read.delim(hist_file, stringsAsFactors = FALSE, check.names = FALSE)
ages_df <- read.delim(ages_file, stringsAsFactors = FALSE, check.names = FALSE)

# Same grouping as Figure2/XJ_GJ_TE/plots/far_distance_boxplot.R
indica_varieties <- c("IR64", "WSSM", "G46", "Tumba", "TM", "N22", "Basmati1")
japonica_varieties <- c("LJ", "MSU", "NamRoo", "Lemont")

add_group <- function(x) {
  ifelse(x %in% indica_varieties, "indica",
         ifelse(x %in% japonica_varieties, "japonica", NA_character_))
}

wide_family_counts <- function(df, value_col) {
  out <- reshape(
    df[, c("Accession", "Group", "Family", value_col)],
    idvar = c("Accession", "Group"),
    timevar = "Family",
    direction = "wide"
  )
  names(out) <- sub(paste0("^", value_col, "\\."), "", names(out))
  for (fam in c("Copia", "Gypsy")) {
    if (!fam %in% names(out)) out[[fam]] <- 0
    out[[fam]][is.na(out[[fam]])] <- 0
  }
  out$Total <- out$Copia + out$Gypsy
  out
}

bin0 <- hist_df[hist_df$bin_start_Ma == 0, c("Accession", "Family", "count")]
bin0$Group <- add_group(bin0$Accession)
bin0 <- bin0[!is.na(bin0$Group), c("Accession", "Group", "Family", "count")]

exact0 <- ages_df[ages_df$age_Ma == 0, c("Accession", "Family", "TE_ID")]
exact0 <- aggregate(TE_ID ~ Accession + Family, exact0, length)
names(exact0)[3] <- "count"
exact0$Group <- add_group(exact0$Accession)
exact0 <- exact0[!is.na(exact0$Group), c("Accession", "Group", "Family", "count")]

bin0_wide <- wide_family_counts(bin0, "count")
names(bin0_wide)[names(bin0_wide) %in% c("Copia", "Gypsy", "Total")] <-
  paste0(c("Copia", "Gypsy", "Total"), "_bin0_0_0.05Ma")

exact0_wide <- wide_family_counts(exact0, "count")
names(exact0_wide)[names(exact0_wide) %in% c("Copia", "Gypsy", "Total")] <-
  paste0(c("Copia", "Gypsy", "Total"), "_exact0")

rank_df <- merge(bin0_wide, exact0_wide, by = c("Accession", "Group"), all = TRUE)
count_cols <- setdiff(names(rank_df), c("Accession", "Group"))
rank_df[count_cols][is.na(rank_df[count_cols])] <- 0
rank_df <- rank_df[order(-rank_df$Total_bin0_0_0.05Ma), ]

write.table(
  rank_df,
  file.path(outdir, "ltr_insertion_age_bin0_rank.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

make_test_rows <- function(df, metric) {
  rows <- list()
  for (fam in c("Copia", "Gypsy", "Total")) {
    d <- if (fam == "Total") {
      aggregate(count ~ Accession + Group, df, sum)
    } else {
      df[df$Family == fam, c("Accession", "Group", "count")]
    }
    indica <- d$count[d$Group == "indica"]
    japonica <- d$count[d$Group == "japonica"]
    wt <- wilcox.test(count ~ Group, data = d, exact = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      metric = metric,
      family = fam,
      n_indica = length(indica),
      n_japonica = length(japonica),
      mean_indica = mean(indica),
      mean_japonica = mean(japonica),
      median_indica = median(indica),
      median_japonica = median(japonica),
      W = unname(wt$statistic),
      p_value = wt$p.value
    )
  }
  do.call(rbind, rows)
}

test_df <- rbind(
  make_test_rows(bin0, "bin_0_0.05Ma"),
  make_test_rows(exact0, "exact_age_0")
)

write.table(
  test_df,
  file.path(outdir, "ltr_insertion_age_group_wilcox.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Wrote: ", file.path(outdir, "ltr_insertion_age_bin0_rank.tsv"))
message("Wrote: ", file.path(outdir, "ltr_insertion_age_group_wilcox.tsv"))
