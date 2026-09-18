options(stringsAsFactors = FALSE)

input_file <- "output/te_genotype_matrix_strict_noNA.csv"
output_dir <- "output/dte_fst"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

japonica <- c("MSU", "Lemont", "LJ", "NamRoo")
indica <- c("Basmati1", "G46", "IR64", "N22", "TM", "Tumba", "WSSM")
outgroup <- "CG14"

meta_cols <- c("TE_ID", "Family", "Chromosome", "Start", "End", "Length")
sample_order <- c(japonica, indica, outgroup)

chr_rank <- function(x) {
  x <- sub("^Chr", "", x)
  suppressWarnings(num <- as.numeric(x))
  ifelse(is.na(num), Inf, num)
}

jaccard_distance <- function(mat) {
  n <- nrow(mat)
  out <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      a <- mat[i, ]
      b <- mat[j, ]
      inter <- sum(a == 1 & b == 1)
      union <- sum(a == 1 | b == 1)
      d <- if (union == 0) 0 else 1 - inter / union
      out[i, j] <- d
      out[j, i] <- d
    }
  }
  out
}

calc_biallelic_fst <- function(p1, p2, n1, n2) {
  p_bar <- (n1 * p1 + n2 * p2) / (n1 + n2)
  hs <- (n1 * (2 * p1 * (1 - p1)) + n2 * (2 * p2 * (1 - p2))) / (n1 + n2)
  ht <- 2 * p_bar * (1 - p_bar)
  if (ht <= 0) {
    return(0)
  }
  fst <- (ht - hs) / ht
  max(0, fst)
}

calc_group_bias <- function(p_jap, p_ind) {
  if (p_jap > p_ind) {
    "japonica"
  } else if (p_ind > p_jap) {
    "indica"
  } else {
    "shared"
  }
}

te <- read.csv(input_file, check.names = FALSE)

missing_cols <- setdiff(c(meta_cols, sample_order), colnames(te))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

sample_mat <- te[, sample_order]
sample_mat[] <- lapply(sample_mat, as.numeric)

jap_mat <- sample_mat[, japonica, drop = FALSE]
ind_mat <- sample_mat[, indica, drop = FALSE]
all_main <- c(japonica, indica)

te$japonica_freq <- rowMeans(jap_mat)
te$indica_freq <- rowMeans(ind_mat)
te$delta_freq <- abs(te$japonica_freq - te$indica_freq)
te$group_bias <- mapply(calc_group_bias, te$japonica_freq, te$indica_freq)

n_jap <- length(japonica)
n_ind <- length(indica)
te$fst <- mapply(calc_biallelic_fst, te$japonica_freq, te$indica_freq,
                 MoreArgs = list(n1 = n_jap, n2 = n_ind))

te$fisher_p <- apply(sample_mat[, all_main, drop = FALSE], 1, function(x) {
  jap_present <- sum(x[japonica] == 1)
  jap_absent <- sum(x[japonica] == 0)
  ind_present <- sum(x[indica] == 1)
  ind_absent <- sum(x[indica] == 0)
  contingency <- matrix(c(jap_present, jap_absent, ind_present, ind_absent), nrow = 2, byrow = TRUE)
  fisher.test(contingency)$p.value
})
te$FDR <- p.adjust(te$fisher_p, method = "BH")

main_polymorphic <- apply(sample_mat[, all_main, drop = FALSE], 1, function(x) length(unique(x)) > 1)
te$main_group_polymorphic <- main_polymorphic
te$is_high_fst_025 <- te$fst >= 0.25
te$is_top_5pct_fst <- te$fst >= as.numeric(quantile(te$fst, 0.95, na.rm = TRUE))
te$is_top_1pct_fst <- te$fst >= as.numeric(quantile(te$fst, 0.99, na.rm = TRUE))

te$ChromRank <- chr_rank(te$Chromosome)
te <- te[order(te$ChromRank, te$Start, te$End), ]
te$genome_order <- seq_len(nrow(te))
te <- te[, c(
  "TE_ID", "Family", "Chromosome", "Start", "End", "Length",
  sample_order,
  "japonica_freq", "indica_freq", "delta_freq", "fst",
  "fisher_p", "FDR", "group_bias", "main_group_polymorphic",
  "is_high_fst_025", "is_top_5pct_fst", "is_top_1pct_fst", "genome_order"
)]

write.csv(te, file.path(output_dir, "dte_fst_results.csv"), row.names = FALSE)

candidate_idx <- with(te, main_group_polymorphic & fst >= 0.25 & delta_freq >= 0.5)
candidates <- te[candidate_idx, ]
candidates <- candidates[order(-candidates$fst, -candidates$delta_freq, candidates$genome_order), ]
write.csv(candidates, file.path(output_dir, "high_fst_te_candidates.csv"), row.names = FALSE)

family_total <- table(te$Family)
family_high <- table(candidates$Family)
family_names <- sort(unique(te$Family))
family_summary <- data.frame(
  Family = family_names,
  total_sites = as.integer(family_total[family_names]),
  high_fst_sites = as.integer(ifelse(is.na(family_high[family_names]), 0, family_high[family_names]))
)
family_summary$high_fst_sites[is.na(family_summary$high_fst_sites)] <- 0
family_summary$high_fst_fraction <- ifelse(family_summary$total_sites > 0,
                                           family_summary$high_fst_sites / family_summary$total_sites,
                                           0)
family_summary <- family_summary[order(-family_summary$high_fst_sites, -family_summary$high_fst_fraction), ]
write.csv(family_summary, file.path(output_dir, "family_level_differentiation_summary.csv"), row.names = FALSE)

sample_binary <- t(as.matrix(sample_mat[, sample_order, drop = FALSE]))
sample_binary <- sample_binary[sample_order, , drop = FALSE]
sample_dist <- jaccard_distance(sample_binary)
sample_sim <- 1 - sample_dist
sample_distance_df <- as.data.frame(as.table(sample_dist))
colnames(sample_distance_df) <- c("Sample1", "Sample2", "JaccardDistance")
sample_distance_df$JaccardSimilarity <- 1 - sample_distance_df$JaccardDistance
write.csv(sample_distance_df, file.path(output_dir, "sample_distance_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(sample_sim), file.path(output_dir, "sample_similarity_matrix.csv"))

summary_lines <- c(
  sprintf("Total TE sites: %d", nrow(te)),
  sprintf("Main-group polymorphic TE sites: %d", sum(te$main_group_polymorphic)),
  sprintf("High-FST candidates (fst >= 0.25 and delta_freq >= 0.5): %d", nrow(candidates)),
  sprintf("Top 5%% FST cutoff: %.4f", as.numeric(quantile(te$fst, 0.95, na.rm = TRUE))),
  sprintf("Top 1%% FST cutoff: %.4f", as.numeric(quantile(te$fst, 0.99, na.rm = TRUE)))
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
