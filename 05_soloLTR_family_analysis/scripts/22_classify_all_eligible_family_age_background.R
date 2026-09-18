#!/usr/bin/env Rscript

# Apply the same 0.002-divergence-bin age definition to every eligible family,
# then evaluate whether an age class carries more absolute GJ-XI occupancy
# divergence than expected from family count or baseline genomic occupancy.
suppressPackageStartupMessages({library(data.table)})

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
bin_width <- 0.002
set.seed(20260810)

eligible <- fread(file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"))[eligible_Nintact_gt10_both_groups == TRUE, family_id]
elements <- fread(file.path(project_root, "data", "solo-LTR_family", "intact_ltr_elements.tsv"))[group %in% c("GJ", "XI") & family_id %in% eligible & is.finite(ltr_identity)]
elements[, age_bin := floor(pmax(0, 1 - ltr_identity) / bin_width) * bin_width]
acc_total <- elements[, .(n_intact = .N), by = .(family_id, accession, group)]
bins <- elements[, .(n = .N), by = .(family_id, accession, group, age_bin)]
bins <- merge(bins, acc_total, by = c("family_id", "accession", "group"), all.x = TRUE)
bins[, frequency := n / n_intact]
max_bin <- elements[, .(max_bin = max(age_bin)), by = family_id]
grid <- merge(acc_total, max_bin, by = "family_id")[, .(age_bin = seq(0, max_bin, by = bin_width)), by = .(family_id, accession, group)]
land <- merge(grid, bins[, .(family_id, accession, group, age_bin, frequency)], by = c("family_id", "accession", "group", "age_bin"), all.x = TRUE)
land[is.na(frequency), frequency := 0]
land <- land[, .(mean_frequency = mean(frequency)), by = .(family_id, group, age_bin)]
pooled <- land[, .(family_frequency = mean(mean_frequency)), by = .(family_id, age_bin)]

age_classes <- pooled[, {
  first <- family_frequency[age_bin == min(age_bin)][1]
  peak <- max(family_frequency)
  .(first_bin_frequency = first, peak_frequency = peak,
    age_landscape_class = fifelse(abs(first - peak) < 1e-12, "Young-peaked", fifelse(first >= 0.05, "Moderate", "Old")))
}, by = family_id]

family_stats <- fread(file.path(project_root, "family_dynamics_nintact10", "family_dynamics_GJ_XI_Nintact10.tsv"),
                      select = c("family_id", "superfamily", "mean_family_size_bp_GJ", "mean_family_size_bp_XI", "delta_bp_XI_minus_GJ"))
d <- merge(family_stats[family_id %in% eligible], age_classes, by = "family_id")
d[, `:=`(
  baseline_occupancy_mb = (mean_family_size_bp_GJ + mean_family_size_bp_XI) / 2e6,
  absolute_occupancy_difference_mb = abs(delta_bp_XI_minus_GJ) / 1e6,
  age_landscape_class = factor(age_landscape_class, levels = c("Young-peaked", "Moderate", "Old"))
)]

summary <- d[, .(
  n_families = .N,
  baseline_occupancy_mb = sum(baseline_occupancy_mb),
  absolute_occupancy_difference_mb = sum(absolute_occupancy_difference_mb)
), by = age_landscape_class]
summary[, `:=`(
  baseline_occupancy_share = baseline_occupancy_mb / sum(baseline_occupancy_mb),
  absolute_difference_share = absolute_occupancy_difference_mb / sum(absolute_occupancy_difference_mb)
)]
summary[, enrichment_vs_baseline := absolute_difference_share / baseline_occupancy_share]

# Exchangeability null: preserve age-class membership and all family-level
# absolute differences, then permute the latter among families.
B <- 10000L
null <- matrix(NA_real_, nrow = B, ncol = nrow(summary))
for (b in seq_len(B)) {
  permuted <- sample(d$absolute_occupancy_difference_mb, replace = FALSE)
  null[b, ] <- tapply(permuted, d$age_landscape_class, sum)[as.character(summary$age_landscape_class)]
}
summary[, `:=`(
  null_mean_absolute_difference_mb = colMeans(null),
  permutation_p_enrichment = vapply(seq_len(.N), function(i) mean(null[, i] >= absolute_occupancy_difference_mb[i]), numeric(1)),
  permutation_p_depletion = vapply(seq_len(.N), function(i) mean(null[, i] <= absolute_occupancy_difference_mb[i]), numeric(1))
)]
summary[, permutation_q_enrichment_BH := p.adjust(permutation_p_enrichment, method = "BH")]
summary[, permutation_q_depletion_BH := p.adjust(permutation_p_depletion, method = "BH")]

fwrite(d, file.path(intermediate_dir, "25_all_eligible_family_age_classes.tsv"), sep = "\t")
fwrite(summary, file.path(intermediate_dir, "25_all_eligible_age_background_contribution.tsv"), sep = "\t")
fwrite(d, file.path(supp_dir, "Supplementary_Table_F6_all_eligible_family_age_classes.tsv"), sep = "\t")
fwrite(summary, file.path(supp_dir, "Supplementary_Table_F6_all_eligible_age_background_contribution.tsv"), sep = "\t")
message("All eligible age classes: ", paste(names(table(d$age_landscape_class)), as.integer(table(d$age_landscape_class)), collapse = "; "))
