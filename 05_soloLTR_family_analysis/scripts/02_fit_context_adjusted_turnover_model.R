#!/usr/bin/env Rscript

# Estimate family-specific XI-versus-GJ solo/intact shifts before and after
# adjustment for gene proximity. These are structural-turnover associations,
# not estimates of deletion rates.

suppressPackageStartupMessages({library(data.table); library(glmmTMB)})

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
counts_path <- file.path(intermediate_dir, "05_family_accession_context_model_counts.tsv.gz")
out_path <- file.path(intermediate_dir, "07_family_GJ_XI_turnover_raw_vs_context_adjusted.tsv")
qc_path <- file.path(intermediate_dir, "07_turnover_model_QC.tsv")

d <- as.data.table(read.delim(gzfile(counts_path), check.names = FALSE))
d <- d[group %in% c("GJ", "XI") & N_intact_context > 0]
availability <- d[, .(N_intact_context = sum(N_intact_context)), by = .(family_id, group)]
availability <- dcast(availability, family_id ~ group, value.var = "N_intact_context", fill = 0)
model_keep <- availability[GJ > 10 & XI > 10, family_id]
d <- d[family_id %in% model_keep]
d[, `:=`(
  group = factor(group, levels = c("GJ", "XI")),
  family_id = factor(family_id),
  accession = factor(accession),
  gene_proximity = factor(gene_proximity, levels = c("far_1kb", "near_1kb", "overlap")),
  log_intact_exposure = log(N_intact_context)
)]

# A no-intercept parameterisation produces one solo-per-intact log rate for
# every family/group combination, making the XI-GJ contrast explicit.
fit_one <- function(include_context) {
  rhs <- if (include_context) {
    "0 + family_id:group + gene_proximity + (1|accession) + offset(log_intact_exposure)"
  } else {
    "0 + family_id:group + (1|accession) + offset(log_intact_exposure)"
  }
  glmmTMB(as.formula(paste("N_solo_context ~", rhs)), family = nbinom2, data = d)
}

adjusted_fit <- fit_one(TRUE)

contrast_by_family <- function(fit, include_context) {
  beta <- fixef(fit)$cond
  V <- vcov(fit)$cond
  fixed_formula <- if (include_context) {
    ~ 0 + family_id:group + gene_proximity
  } else {
    ~ 0 + family_id:group
  }
  fam <- levels(d$family_id)
  ans <- rbindlist(lapply(fam, function(f) {
    nd <- data.frame(
      family_id = factor(c(f, f), levels = levels(d$family_id)),
      group = factor(c("GJ", "XI"), levels = levels(d$group)),
      gene_proximity = factor(c("far_1kb", "far_1kb"), levels = levels(d$gene_proximity))
    )
    X <- model.matrix(fixed_formula, nd)
    # Align columns defensively in case a factor level is absent from a fit.
    Xfull <- matrix(0, nrow = 2, ncol = length(beta), dimnames = list(NULL, names(beta)))
    common <- intersect(colnames(X), names(beta))
    Xfull[, common] <- X[, common, drop = FALSE]
    cvec <- Xfull[2, ] - Xfull[1, ]
    nz <- which(cvec != 0)
    estimate <- sum(cvec[nz] * beta[nz])
    se <- sqrt(drop(t(cvec[nz]) %*% V[nz, nz, drop = FALSE] %*% cvec[nz]))
    z <- estimate / se
    data.table(family_id = f, log_IRR_XI_vs_GJ = estimate, SE = se, z = z,
               p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
               IRR_XI_vs_GJ = exp(estimate))
  }))
  ans[, q_BH := p.adjust(p_value, method = "BH")]
  ans
}

# Unadjusted values are descriptive rate ratios.  Inference is deliberately
# reserved for the converged context-adjusted model below.
raw_rates <- d[, .(N_solo = sum(N_solo_context), N_intact = sum(N_intact_context)), by = .(family_id, group)]
raw_rates <- dcast(raw_rates, family_id ~ group, value.var = c("N_solo", "N_intact"), fill = 0)
raw <- raw_rates[, .(
  family_id,
  raw_log_IRR_XI_vs_GJ = log((N_solo_XI / N_intact_XI) / (N_solo_GJ / N_intact_GJ)),
  raw_IRR_XI_vs_GJ = (N_solo_XI / N_intact_XI) / (N_solo_GJ / N_intact_GJ)
)]
adjusted <- contrast_by_family(adjusted_fit, TRUE)
setnames(adjusted, setdiff(names(adjusted), "family_id"), paste0("adjusted_", setdiff(names(adjusted), "family_id")))

out <- merge(raw, adjusted, by = "family_id")
out[, `:=`(
  attenuation_abs_logIRR = abs(raw_log_IRR_XI_vs_GJ) - abs(adjusted_log_IRR_XI_vs_GJ),
  attenuation_fraction = fifelse(abs(raw_log_IRR_XI_vs_GJ) > 0,
                                1 - abs(adjusted_log_IRR_XI_vs_GJ) / abs(raw_log_IRR_XI_vs_GJ), NA_real_),
  raw_direction = fifelse(raw_log_IRR_XI_vs_GJ > 0, "XI higher", "GJ higher"),
  adjusted_direction = fifelse(adjusted_log_IRR_XI_vs_GJ > 0, "XI higher", "GJ higher")
)]
meta <- fread(file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"))
out <- merge(meta[eligible_Nintact_gt10_both_groups == TRUE], out, by = "family_id", all.y = TRUE)
setorder(out, adjusted_q_BH, adjusted_p_value)
fwrite(out, out_path, sep = "\t")

qc <- data.table(
  model = "context_adjusted",
  n_rows = nrow(d),
  n_families = uniqueN(d$family_id),
  convergence_code = adjusted_fit$fit$convergence,
  pdHess = adjusted_fit$sdr$pdHess,
  AIC = AIC(adjusted_fit)
)
fwrite(qc, qc_path, sep = "\t")
fwrite(availability, file.path(intermediate_dir, "07_context_model_family_availability.tsv"), sep = "\t")
saveRDS(list(context_adjusted = adjusted_fit), file.path(intermediate_dir, "07_turnover_models.rds"))
