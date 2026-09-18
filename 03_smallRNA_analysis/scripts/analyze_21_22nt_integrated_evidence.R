#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

infile <- Sys.getenv("TE_SIGNAL_TABLE")
outdir <- Sys.getenv("TE_SIGNAL_OUTDIR")
if (!nzchar(infile) || !nzchar(outdir)) stop("Set TE_SIGNAL_TABLE and TE_SIGNAL_OUTDIR.")

dt <- fread(infile, showProgress = TRUE)
dt[, category := paste(family, gene_proximity, sep = "_")]
category_order <- c(
  "Copia_Overlap", "Copia_Near", "Copia_Far",
  "Gypsy_Overlap", "Gypsy_Near", "Gypsy_Far"
)

rank_z <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (sum(keep) < 3L) return(out)
  r <- rank(x[keep], ties.method = "average")
  s <- sd(r)
  if (!is.finite(s) || s == 0) return(out)
  out[keep] <- (r - mean(r)) / s
  out
}

signflip_test <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) {
    return(list(
      n = n, mean = NA_real_, median = NA_real_, lower = NA_real_,
      upper = NA_real_, positive = NA_integer_, p = NA_real_
    ))
  }
  obs <- mean(x)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  perm <- as.vector(signs %*% x / n)
  p <- mean(abs(perm) >= abs(obs) - 1e-12)
  se <- sd(x) / sqrt(n)
  crit <- qt(0.975, df = n - 1)
  list(
    n = n,
    mean = obs,
    median = median(x),
    lower = obs - crit * se,
    upper = obs + crit * se,
    positive = sum(x > 0),
    p = p
  )
}

message("1/4: Detection-state distribution...")
dt[, state := fifelse(
  mean_24nt_sRNA_log2CPM > 0 & mean_21_22nt_sRNA_log2CPM > 0,
  "24+ / 21-22+",
  fifelse(
    mean_24nt_sRNA_log2CPM > 0 & mean_21_22nt_sRNA_log2CPM <= 0,
    "24+ / 21-22-",
    fifelse(
      mean_24nt_sRNA_log2CPM <= 0 & mean_21_22nt_sRNA_log2CPM > 0,
      "24- / 21-22+",
      "24- / 21-22-"
    )
  )
)]

state_by_accession <- dt[, .N, by = .(accession, subspecies, family, gene_proximity, category, state)]
state_by_accession[, proportion := N / sum(N), by = .(accession, category)]
fwrite(
  state_by_accession,
  file.path(outdir, "19_small_RNA_detection_states_by_accession.tsv"),
  sep = "\t", quote = FALSE
)

state_summary <- state_by_accession[, .(
  n_accessions = .N,
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  min_proportion = min(proportion),
  max_proportion = max(proportion)
), by = .(family, gene_proximity, category, state)]
fwrite(
  state_summary,
  file.path(outdir, "20_small_RNA_detection_state_summary.tsv"),
  sep = "\t", quote = FALSE
)

message("2/4: Four-state methylation contrasts...")
meth_cols <- c(
  CHH = "mean_CHH_pct",
  CHG = "mean_CHG_pct",
  CpG = "mean_CpG_pct"
)
contrast_defs <- list(
  combined_targeting = c("24+ / 21-22+", "24+ / 21-22-"),
  candidate_noncanonical = c("24- / 21-22+", "24- / 21-22-")
)

state_effects <- rbindlist(lapply(names(meth_cols), function(context) {
  response <- meth_cols[[context]]
  rbindlist(lapply(names(contrast_defs), function(contrast) {
    lev <- contrast_defs[[contrast]]
    dt[, {
      a <- get(response)[state == lev[1] & is.finite(get(response))]
      b <- get(response)[state == lev[2] & is.finite(get(response))]
      list(
        context = context,
        contrast = contrast,
        numerator_state = lev[1],
        denominator_state = lev[2],
        n_numerator = length(a),
        n_denominator = length(b),
        median_numerator = if (length(a)) median(a) else NA_real_,
        median_denominator = if (length(b)) median(b) else NA_real_,
        median_difference_pct_points = if (length(a) >= 30L && length(b) >= 30L) {
          median(a) - median(b)
        } else {
          NA_real_
        }
      )
    }, by = .(accession, subspecies, family, gene_proximity, category)]
  }))
}))
fwrite(
  state_effects,
  file.path(outdir, "21_four_state_methylation_effects_by_accession.tsv"),
  sep = "\t", quote = FALSE
)

state_tests <- state_effects[, {
  z <- signflip_test(median_difference_pct_points)
  list(
    n_accessions = z$n,
    mean_difference_pct_points = z$mean,
    median_difference_pct_points = z$median,
    CI95_lower = z$lower,
    CI95_upper = z$upper,
    positive_accessions = z$positive,
    exact_signflip_p = z$p
  )
}, by = .(family, gene_proximity, category, context, contrast)]
state_tests[, FDR := p.adjust(exact_signflip_p, method = "BH"), by = .(contrast, context)]
fwrite(
  state_tests,
  file.path(outdir, "22_four_state_methylation_tests.tsv"),
  sep = "\t", quote = FALSE
)

message("3/4: Continuous 24-nt x 21-22-nt models...")
fit_interaction <- function(d, response) {
  y <- rank_z(d[[response]])
  x24 <- rank_z(d$mean_24nt_sRNA_log2CPM)
  x21 <- rank_z(d$mean_21_22nt_sRNA_log2CPM)
  xlen <- rank_z(d$log10_TE_length)
  keep <- is.finite(y) & is.finite(x24) & is.finite(x21) & is.finite(xlen)
  n <- sum(keep)
  if (n < 100L) {
    return(list(n = n, beta_24 = NA_real_, beta_21_22 = NA_real_, beta_interaction = NA_real_))
  }
  X <- cbind(1, x24[keep], x21[keep], x24[keep] * x21[keep], xlen[keep])
  fit <- lm.fit(X, y[keep])
  list(
    n = n,
    beta_24 = unname(fit$coefficients[2]),
    beta_21_22 = unname(fit$coefficients[3]),
    beta_interaction = unname(fit$coefficients[4])
  )
}

model_effects <- rbindlist(lapply(names(meth_cols), function(context) {
  response <- meth_cols[[context]]
  dt[, c(list(context = context), fit_interaction(.SD, response)),
    by = .(accession, subspecies, family, gene_proximity, category)
  ]
}))
fwrite(
  model_effects,
  file.path(outdir, "23_continuous_interaction_models_by_accession.tsv"),
  sep = "\t", quote = FALSE
)

model_long <- melt(
  model_effects,
  id.vars = c(
    "accession", "subspecies", "family", "gene_proximity",
    "category", "context", "n"
  ),
  measure.vars = c("beta_24", "beta_21_22", "beta_interaction"),
  variable.name = "term",
  value.name = "standardized_rank_beta"
)
model_tests <- model_long[, {
  z <- signflip_test(standardized_rank_beta)
  list(
    n_accessions = z$n,
    mean_beta = z$mean,
    median_beta = z$median,
    CI95_lower = z$lower,
    CI95_upper = z$upper,
    positive_accessions = z$positive,
    exact_signflip_p = z$p
  )
}, by = .(family, gene_proximity, category, context, term)]
model_tests[, FDR := p.adjust(exact_signflip_p, method = "BH"), by = .(context, term)]
fwrite(
  model_tests,
  file.path(outdir, "24_continuous_interaction_model_tests.tsv"),
  sep = "\t", quote = FALSE
)

context_wide <- dcast(
  model_long[term %chin% c("beta_21_22", "beta_interaction")],
  accession + subspecies + family + gene_proximity + category + term ~ context,
  value.var = "standardized_rank_beta"
)
context_tests <- rbindlist(lapply(c("CHG", "CpG"), function(comparator) {
  context_wide[, {
    z <- signflip_test(CHH - get(comparator))
    list(
      comparator_context = comparator,
      n_accessions = z$n,
      mean_beta_difference = z$mean,
      median_beta_difference = z$median,
      CI95_lower = z$lower,
      CI95_upper = z$upper,
      positive_accessions = z$positive,
      exact_signflip_p = z$p
    )
  }, by = .(family, gene_proximity, category, term)]
}))
context_tests[, FDR := p.adjust(exact_signflip_p, method = "BH"), by = term]
fwrite(
  context_tests,
  file.path(outdir, "24b_CHH_vs_maintenance_context_model_tests.tsv"),
  sep = "\t", quote = FALSE
)

message("4/4: Residual 21-22-nt enrichment by TE category...")
residual_dt <- dt[, {
  y <- rank_z(mean_21_22nt_sRNA_log2CPM)
  x24 <- rank_z(mean_24nt_sRNA_log2CPM)
  xchh <- rank_z(mean_CHH_pct)
  xchg <- rank_z(mean_CHG_pct)
  xcpg <- rank_z(mean_CpG_pct)
  xlen <- rank_z(log10_TE_length)
  keep <- is.finite(y) & is.finite(x24) & is.finite(xchh) &
    is.finite(xchg) & is.finite(xcpg) & is.finite(xlen)
  residual <- rep(NA_real_, .N)
  if (sum(keep) >= 100L) {
    fit <- lm.fit(
      cbind(1, x24[keep], xchh[keep], xchg[keep], xcpg[keep], xlen[keep]),
      y[keep]
    )
    residual[keep] <- fit$residuals
  }
  .(
    family, gene_proximity, category,
    residual_21_22 = residual
  )
}, by = .(accession, subspecies)]

residual_category <- residual_dt[is.finite(residual_21_22), .(
  n_TEs = .N,
  mean_residual = mean(residual_21_22),
  median_residual = median(residual_21_22)
), by = .(accession, subspecies, family, gene_proximity, category)]
fwrite(
  residual_category,
  file.path(outdir, "25_residual_21_22_enrichment_by_accession.tsv"),
  sep = "\t", quote = FALSE
)

wide <- dcast(
  residual_category,
  accession + subspecies ~ category,
  value.var = "median_residual"
)
comparators <- setdiff(category_order, "Copia_Overlap")
residual_contrasts <- rbindlist(lapply(comparators, function(comp) {
  diffs <- wide[["Copia_Overlap"]] - wide[[comp]]
  z <- signflip_test(diffs)
  data.table(
    reference = "Copia_Overlap",
    comparator = comp,
    n_accessions = z$n,
    mean_paired_difference = z$mean,
    median_paired_difference = z$median,
    CI95_lower = z$lower,
    CI95_upper = z$upper,
    positive_accessions = z$positive,
    exact_signflip_p = z$p
  )
}))
residual_contrasts[, FDR := p.adjust(exact_signflip_p, method = "BH")]
fwrite(
  residual_contrasts,
  file.path(outdir, "26_Copia_overlap_residual_enrichment_tests.tsv"),
  sep = "\t", quote = FALSE
)

message("Integrated evidence analysis complete.")
