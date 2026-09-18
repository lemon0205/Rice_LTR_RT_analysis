#!/usr/bin/env Rscript

# Independent family-level solo-LTR abundance test.  This is deliberately
# separate from solo:intact: a changed ratio can arise from either component,
# whereas this table identifies families with a population contrast in solo
# copy number itself.
suppressPackageStartupMessages(library(data.table))

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
data_dir <- file.path(project_root, "data", "solo-LTR_family")
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")

exact_group_test <- function(values, groups) {
  obs <- mean(values[groups == "XI"]) - mean(values[groups == "GJ"])
  gj_sets <- combn(seq_along(values), sum(groups == "GJ"), simplify = FALSE)
  nul <- vapply(gj_sets, function(gj) {
    xi <- setdiff(seq_along(values), gj)
    mean(values[xi]) - mean(values[gj])
  }, numeric(1))
  c(effect = obs, p = mean(abs(nul) >= abs(obs) - 1e-12))
}

family <- fread(file.path(data_dir, "family_by_accession.tsv"))[group %in% c("GJ", "XI")]
solo <- fread(file.path(data_dir, "soloLTR_TSD_supported.best_Os_family.keep_conflict_drop_nohit.tsv"))[
  group %in% c("GJ", "XI"), .(N_solo = .N), by = .(accession, group, family_id, superfamily)]
meta <- unique(family[, .(family_id, superfamily)])
acc <- unique(family[, .(accession, group)])
d <- CJ(family_id = meta$family_id, accession = acc$accession, unique = TRUE)
d <- merge(d, meta, by = "family_id", all.x = TRUE)
d <- merge(d, acc, by = "accession", all.x = TRUE)
d <- merge(d, solo[, .(accession, family_id, N_solo)], by = c("accession", "family_id"), all.x = TRUE)
d[is.na(N_solo), N_solo := 0]

raw <- fread(file.path(intermediate_dir, "20_raw_observed_family_dynamics.tsv"))[
  , .(family_id, eligible_Nintact_gt10_both_groups, dynamics_class,
      occupancy_direction, delta_bp_XI_minus_GJ, raw_log_SI_XI_minus_GJ)]
res <- d[, {
  t <- exact_group_test(log2(N_solo + 0.5), group)
  .(mean_Nsolo_GJ = mean(N_solo[group == "GJ"]),
    mean_Nsolo_XI = mean(N_solo[group == "XI"]),
    log2_solo_abundance_XI_GJ = t[["effect"]], p_solo_abundance_exact = t[["p"]])
}, by = .(family_id, superfamily)]
res <- merge(res, raw, by = "family_id")
res[, q_solo_abundance_BH := NA_real_]
for (sf in unique(res$superfamily)) {
  ii <- res$superfamily == sf & res$eligible_Nintact_gt10_both_groups
  res[ii, q_solo_abundance_BH := p.adjust(p_solo_abundance_exact, "BH")]
}
res[, solo_abundance_significant := q_solo_abundance_BH < 0.05]
res[, solo_abundance_direction := fifelse(log2_solo_abundance_XI_GJ > 0, "XI higher", "GJ higher")]
res[, abs_log2_solo_abundance_contrast := abs(log2_solo_abundance_XI_GJ)]
setorder(res, superfamily, -solo_abundance_significant, -abs_log2_solo_abundance_contrast, family_id)

fwrite(res, file.path(intermediate_dir, "36_family_solo_LTR_abundance_differences.tsv"), sep = "\t")
fwrite(res, file.path(supp_dir, "Supplementary_Table_F7_family_solo_LTR_abundance_differences.tsv"), sep = "\t")
summary <- res[eligible_Nintact_gt10_both_groups == TRUE,
  .(n_families = .N, n_solo_abundance_significant = sum(solo_abundance_significant),
    XI_higher = sum(solo_abundance_significant & solo_abundance_direction == "XI higher"),
    GJ_higher = sum(solo_abundance_significant & solo_abundance_direction == "GJ higher")),
  by = .(superfamily, dynamics_class)]
fwrite(summary, file.path(intermediate_dir, "36_family_solo_LTR_abundance_differences_summary.tsv"), sep = "\t")
message("Eligible solo-LTR abundance differences: ", sum(res$eligible_Nintact_gt10_both_groups & res$solo_abundance_significant), "/", sum(res$eligible_Nintact_gt10_both_groups))
print(summary)
