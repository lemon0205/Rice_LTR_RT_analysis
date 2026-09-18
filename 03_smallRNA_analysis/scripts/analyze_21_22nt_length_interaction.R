#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

infile <- Sys.getenv("TE_SIGNAL_TABLE")
outdir <- Sys.getenv("TE_SIGNAL_OUTDIR", dirname(infile))
if (!nzchar(infile)) stop("Set TE_SIGNAL_TABLE; optionally set TE_SIGNAL_OUTDIR.")
dt <- fread(infile, showProgress = FALSE)
dt[, length_class := fifelse(TE_length <= 500, "Short (<=500 bp)",
                       fifelse(TE_length <= 1500, "Medium (501-1500 bp)", "Long (>1500 bp)"))]
dt[, category := paste(family, length_class, sep = "_")]

rank_z <- function(x) {
  out <- rep(NA_real_, length(x)); keep <- is.finite(x)
  if (sum(keep) < 3L) return(out)
  r <- rank(x[keep], ties.method = "average"); s <- sd(r)
  if (!is.finite(s) || s == 0) return(out)
  out[keep] <- (r - mean(r)) / s; out
}
signflip_test <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 3L) return(list(n=n, mean=NA_real_, median=NA_real_, lower=NA_real_, upper=NA_real_, positive=NA_integer_, p=NA_real_))
  obs <- mean(x); signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  perm <- as.vector(signs %*% x / n); se <- sd(x) / sqrt(n); crit <- qt(0.975, df=n-1)
  list(n=n, mean=obs, median=median(x), lower=obs-crit*se, upper=obs+crit*se,
       positive=sum(x > 0), p=mean(abs(perm) >= abs(obs)-1e-12))
}
fit_one <- function(d, response) {
  y <- rank_z(d[[response]]); x24 <- rank_z(d$mean_24nt_sRNA_log2CPM)
  x21 <- rank_z(d$mean_21_22nt_sRNA_log2CPM); xlen <- rank_z(d$log10_TE_length)
  keep <- is.finite(y) & is.finite(x24) & is.finite(x21) & is.finite(xlen)
  if (sum(keep) < 100L) return(c(n=sum(keep), beta_24=NA, beta_21_22=NA, beta_interaction=NA))
  X <- cbind(1, x24[keep], x21[keep], x24[keep]*x21[keep], xlen[keep])
  fit <- lm.fit(X, y[keep])
  c(n=sum(keep), beta_24=fit$coefficients[2], beta_21_22=fit$coefficients[3], beta_interaction=fit$coefficients[4])
}

meth_cols <- c(CHH="mean_CHH_pct", CHG="mean_CHG_pct", CpG="mean_CpG_pct")
effects <- rbindlist(lapply(names(meth_cols), function(ctx) {
  rbindlist(lapply(split(dt, list(dt$accession, dt$family, dt$length_class), drop=TRUE), function(x) {
    b <- fit_one(x, meth_cols[[ctx]])
    data.table(accession=x$accession[1], subspecies=x$subspecies[1], family=x$family[1],
      length_class=x$length_class[1], category=x$category[1], context=ctx, n=unname(b[1]),
      beta_24=unname(b[2]), beta_21_22=unname(b[3]), beta_interaction=unname(b[4]))
  }))
}))
fwrite(effects, file.path(outdir, "25_length_interaction_models_by_accession.tsv"), sep="\t", quote=FALSE)
tests <- effects[, { z <- signflip_test(beta_interaction); .(n_accessions=z$n, mean_beta_interaction=z$mean,
  median_beta_interaction=z$median, CI95_lower=z$lower, CI95_upper=z$upper, positive_accessions=z$positive,
  exact_signflip_p=z$p) }, by=.(family, length_class, category, context)]
tests[, FDR := p.adjust(exact_signflip_p, method="BH"), by=.(context)]
fwrite(tests, file.path(outdir, "26_length_interaction_model_tests.tsv"), sep="\t", quote=FALSE)
chh <- tests[context == "CHH"]
fwrite(chh, file.path(outdir, "26_length_interaction_CHH_summary.tsv"), sep="\t", quote=FALSE)
print(chh[, .(family, length_class, n_accessions, mean_beta_interaction, CI95_lower, CI95_upper, exact_signflip_p, FDR)])
