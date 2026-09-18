#!/usr/bin/env Rscript

# Reassess recent-expansion evidence without using young-intact counts as an
# age proxy.  The young component is tested as young:older intact composition
# (quasi-binomial log odds) within each family, and LTR identity is tested as
# an accession-level GJ-XI difference by an exact label permutation test.
suppressPackageStartupMessages(library(data.table))

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
data_root <- file.path(project_root, "data", "solo-LTR_family")
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

exact_label_test <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]; group <- group[ok]
  n_gj <- sum(group == "GJ")
  if (length(x) < 4L || n_gj < 2L || sum(group == "XI") < 2L) return(c(effect = NA_real_, p_value = NA_real_))
  observed <- mean(x[group == "XI"]) - mean(x[group == "GJ"])
  combos <- combn(seq_along(x), n_gj)
  null <- apply(combos, 2L, function(gj) mean(x[-gj]) - mean(x[gj]))
  c(effect = observed, p_value = mean(abs(null) >= abs(observed) - 1e-12))
}

fit_young_odds <- function(x) {
  x <- x[N_intact > 0L]
  if (uniqueN(x$group) != 2L || min(x[, .N, by = group]$N) < 2L) return(c(effect = NA_real_, p_value = NA_real_))
  x[, N_old_intact := pmax(N_intact - N_young_intact, 0L)]
  fit <- tryCatch(glm(cbind(N_young_intact, N_old_intact) ~ group, family = quasibinomial(link = "logit"), data = x), error = function(e) NULL)
  if (is.null(fit) || !"groupXI" %in% names(coef(fit))) return(c(effect = NA_real_, p_value = NA_real_))
  co <- summary(fit)$coefficients
  if (!is.finite(co["groupXI", "Estimate"]) || !is.finite(co["groupXI", "Pr(>|t|)"])) return(c(effect = NA_real_, p_value = NA_real_))
  c(effect = co["groupXI", "Estimate"], p_value = co["groupXI", "Pr(>|t|)"])
}

family <- fread(file.path(data_root, "family_by_accession.tsv"))[group %in% c("GJ", "XI")]
eligibility <- fread(file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"))[eligible_Nintact_gt10_both_groups == TRUE]
family <- family[family_id %in% eligibility$family_id]

age <- family[, {
  young <- fit_young_odds(.SD)
  identity <- exact_label_test(median_ltr_identity, group)
  .(young_odds_logOR_XI_vs_GJ = young["effect"], p_young_odds = young["p_value"],
    delta_mean_young_fraction_XI_minus_GJ = mean(young_intact_fraction[group == "XI"], na.rm = TRUE) - mean(young_intact_fraction[group == "GJ"], na.rm = TRUE),
    delta_median_identity_XI_minus_GJ = identity["effect"], p_identity_exact = identity["p_value"])
}, by = .(family_id, superfamily)]
age[, `:=`(q_young_odds_BH = p.adjust(p_young_odds, method = "BH"),
           q_identity_BH = p.adjust(p_identity_exact, method = "BH")), by = superfamily]

dyn <- fread(file.path(project_root, "family_dynamics_nintact10", "family_dynamics_GJ_XI_Nintact10.tsv"),
             select = c("family_id", "delta_log2_size_XI_minus_GJ", "q_size_BH"))
out <- merge(eligibility[, .(family_id, superfamily, N_intact_GJ, N_intact_XI)], dyn, by = "family_id", all.x = TRUE)
out <- merge(out, age, by = c("family_id", "superfamily"), all.x = TRUE)
out[, `:=`(
  abundance_shift = is.finite(q_size_BH) & q_size_BH < 0.05,
  young_composition_support = is.finite(q_young_odds_BH) & q_young_odds_BH < 0.05 & sign(young_odds_logOR_XI_vs_GJ) == sign(delta_log2_size_XI_minus_GJ),
  identity_support = is.finite(q_identity_BH) & q_identity_BH < 0.05 & sign(delta_median_identity_XI_minus_GJ) == sign(delta_log2_size_XI_minus_GJ)
)]
out[, recent_expansion_supported := abundance_shift & young_composition_support & identity_support]
out[, abs_abundance_effect := abs(delta_log2_size_XI_minus_GJ)]
setorder(out, superfamily, -recent_expansion_supported, -abs_abundance_effect)
out[, abs_abundance_effect := NULL]

fwrite(out, file.path(intermediate_dir, "19_recent_expansion_age_composition_reassessment.tsv"), sep = "\t")
fwrite(out, file.path(supp_dir, "Supplementary_Table_F6_recent_expansion_age_composition_reassessment.tsv"), sep = "\t")
summary <- out[, .N, by = .(superfamily, abundance_shift, young_composition_support, identity_support, recent_expansion_supported)]
fwrite(summary, file.path(intermediate_dir, "19_recent_expansion_age_composition_summary.tsv"), sep = "\t")
message("Eligible families: ", nrow(out), "; abundance shifts: ", sum(out$abundance_shift),
        "; young-composition-supported: ", sum(out$abundance_shift & out$young_composition_support),
        "; recent-expansion-supported: ", sum(out$recent_expansion_supported))
