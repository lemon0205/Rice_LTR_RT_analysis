#!/usr/bin/env Rscript

# Leave-one-accession-out robustness check for occupancy and context-adjusted
# solo:intact dynamics. Each fold repeats the family-wise occupancy test and
# the same negative-binomial context model used for the primary classification.
suppressPackageStartupMessages({library(data.table); library(glmmTMB)})

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
eligible <- fread(file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"))[eligible_Nintact_gt10_both_groups == TRUE, family_id]
family <- fread(file.path(project_root, "data", "solo-LTR_family", "family_by_accession.tsv"))[group %in% c("GJ", "XI") & family_id %in% eligible]
counts <- as.data.table(read.delim(gzfile(file.path(intermediate_dir, "05_family_accession_context_model_counts.tsv.gz")), check.names = FALSE))
counts <- counts[group %in% c("GJ", "XI") & N_intact_context > 0 & family_id %in% eligible]
baseline <- fread(file.path(intermediate_dir, "20_context_aware_family_dynamics.tsv"))[, .(family_id, baseline_class = dynamics_class)]

fit_structural <- function(x) {
  avail <- x[, .(N_intact_context = sum(N_intact_context)), by = .(family_id, group)]
  avail <- dcast(avail, family_id ~ group, value.var = "N_intact_context", fill = 0)
  keep <- avail[GJ > 10 & XI > 10, family_id]
  z <- copy(x[family_id %in% keep])
  z[, `:=`(group = factor(group, levels = c("GJ", "XI")), family_id = factor(family_id),
            accession = factor(accession), gene_proximity = factor(gene_proximity, levels = c("far_1kb", "near_1kb", "overlap")),
            log_intact_exposure = log(N_intact_context))]
  fit <- glmmTMB(N_solo_context ~ 0 + family_id:group + gene_proximity + (1 | accession) + offset(log_intact_exposure),
                 family = nbinom2, data = z)
  beta <- fixef(fit)$cond; V <- vcov(fit)$cond; fam <- levels(z$family_id)
  ans <- rbindlist(lapply(fam, function(f) {
    nd <- data.frame(family_id = factor(c(f, f), levels = fam), group = factor(c("GJ", "XI"), levels = c("GJ", "XI")),
                     gene_proximity = factor(c("far_1kb", "far_1kb"), levels = c("far_1kb", "near_1kb", "overlap")))
    X <- model.matrix(~ 0 + family_id:group + gene_proximity, nd)
    Xfull <- matrix(0, nrow = 2, ncol = length(beta), dimnames = list(NULL, names(beta)))
    common <- intersect(colnames(X), names(beta)); Xfull[, common] <- X[, common, drop = FALSE]
    cv <- Xfull[2, ] - Xfull[1, ]; nz <- which(cv != 0)
    est <- sum(cv[nz] * beta[nz]); se <- sqrt(drop(t(cv[nz]) %*% V[nz, nz, drop = FALSE] %*% cv[nz]))
    data.table(family_id = f, adjusted_log_IRR_XI_vs_GJ = est, adjusted_p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
  }))
  ans[, adjusted_q := p.adjust(adjusted_p, method = "BH")]
  ans
}

fold_results <- rbindlist(lapply(sort(unique(family$accession)), function(drop_acc) {
  f <- family[accession != drop_acc]
  occ <- f[, .(mean_bp = mean(family_size_bp)), by = .(family_id, superfamily, group)]
  occ <- dcast(occ, family_id + superfamily ~ group, value.var = "mean_bp")
  occ[, `:=`(delta_bp_XI_minus_GJ = XI - GJ, delta_log2_size_XI_minus_GJ = log2(XI / GJ))]
  occ[, occupancy_p := vapply(family_id, function(id) {
    q <- f[family_id == id]
    wilcox.test(family_size_bp ~ group, data = q, exact = TRUE)$p.value
  }, numeric(1))]
  occ[, occupancy_q := p.adjust(occupancy_p, method = "BH")]
  st <- fit_structural(counts[accession != drop_acc])
  out <- merge(occ, st, by = "family_id", all.x = TRUE)
  out[, `:=`(occupancy_significant = occupancy_q < 0.05,
             structural_significant = adjusted_q < 0.05)]
  out[, dynamics_class := fcase(
    occupancy_significant & !structural_significant, "Amplification-like",
    !occupancy_significant & structural_significant, "Balanced structural divergence",
    occupancy_significant & structural_significant & sign(delta_log2_size_XI_minus_GJ) * sign(adjusted_log_IRR_XI_vs_GJ) < 0, "Removal-like",
    occupancy_significant & structural_significant, "Mixed / high-turnover",
    default = "Drifting")]
  out[, excluded_accession := drop_acc]
  out
}), fill = TRUE)

out <- merge(fold_results, baseline, by = "family_id", all.x = TRUE)
out[, `:=`(direction_retained = sign(delta_bp_XI_minus_GJ) == sign(delta_bp_XI_minus_GJ[excluded_accession == excluded_accession][1]),
           class_retained = dynamics_class == baseline_class)]
# Direction stability must be evaluated relative to the full-data estimate.
base_direction <- fread(file.path(intermediate_dir, "20_context_aware_family_dynamics.tsv"))[, .(family_id, baseline_delta_bp = delta_bp_XI_minus_GJ)]
out <- merge(out, base_direction, by = "family_id", all.x = TRUE)
out[, direction_retained := sign(delta_bp_XI_minus_GJ) == sign(baseline_delta_bp)]
summary <- out[, .(
  n_folds = .N,
  direction_stable_folds = sum(direction_retained, na.rm = TRUE),
  direction_stability = mean(direction_retained, na.rm = TRUE),
  class_stable_folds = sum(class_retained, na.rm = TRUE),
  class_stability = mean(class_retained, na.rm = TRUE),
  min_abs_delta_mb = min(abs(delta_bp_XI_minus_GJ), na.rm = TRUE),
  max_abs_delta_mb = max(abs(delta_bp_XI_minus_GJ), na.rm = TRUE)
), by = .(family_id, baseline_class)]

fwrite(out, file.path(intermediate_dir, "26_leave_one_accession_out_family_dynamics.tsv"), sep = "\t")
fwrite(summary, file.path(intermediate_dir, "26_leave_one_accession_out_stability_summary.tsv"), sep = "\t")
fwrite(summary, file.path(supp_dir, "Supplementary_Table_F6_leave_one_accession_out_stability.tsv"), sep = "\t")
message("Completed ", uniqueN(out$excluded_accession), " leave-one-accession-out folds.")
