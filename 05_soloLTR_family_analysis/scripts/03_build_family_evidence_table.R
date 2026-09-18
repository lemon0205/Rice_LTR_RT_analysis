#!/usr/bin/env Rscript

# Combine abundance, age, and context-adjusted structural-turnover evidence.
# This table exposes every decision flag; it deliberately does not force a
# biological class when a family has insufficient or conflicting evidence.

suppressPackageStartupMessages(library(data.table))

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

dyn <- fread(file.path(project_root, "family_dynamics_nintact10", "family_dynamics_GJ_XI_Nintact10.tsv"))
turnover <- fread(file.path(intermediate_dir, "07_family_GJ_XI_turnover_raw_vs_context_adjusted.tsv"))

d <- merge(turnover, dyn[, .(family_id, delta_log2_size_XI_minus_GJ, q_size_BH,
                              delta_N_young_intact_XI_minus_GJ, q_N_young_intact_BH,
                              mean_median_identity_GJ, mean_median_identity_XI)],
           by = "family_id", all.x = TRUE)
d[, `:=`(
  abundance_shift_supported = !is.na(q_size_BH) & q_size_BH < 0.05,
  abundance_direction = fifelse(delta_log2_size_XI_minus_GJ > 0, "XI higher", "GJ higher"),
  young_support_directional = !is.na(q_N_young_intact_BH) & q_N_young_intact_BH < 0.05 &
    sign(delta_N_young_intact_XI_minus_GJ) == sign(delta_log2_size_XI_minus_GJ),
  identity_directional = !is.na(mean_median_identity_XI) &
    sign(mean_median_identity_XI - mean_median_identity_GJ) == sign(delta_log2_size_XI_minus_GJ),
  adjusted_turnover_supported = !is.na(adjusted_q_BH) & adjusted_q_BH < 0.05 &
    sign(adjusted_log_IRR_XI_vs_GJ) != sign(delta_log2_size_XI_minus_GJ),
  context_attenuated = !is.na(attenuation_fraction) & attenuation_fraction >= 0.50
)]
d[, expansion_evidence := abundance_shift_supported & young_support_directional & identity_directional]
d[, provisional_interpretation := fcase(
  !abundance_shift_supported, "No supported abundance divergence",
  expansion_evidence & adjusted_turnover_supported, "Mixed: expansion and adjusted turnover evidence",
  expansion_evidence, "Expansion-driven candidate",
  adjusted_turnover_supported, "Removal-consistent candidate",
  context_attenuated, "Context-attenuated candidate; no adjusted turnover support",
  default = "Abundance-shift candidate; mechanism unresolved"
)]

keep_cols <- c("family_id", "superfamily", "N_intact_GJ", "N_intact_XI",
  "delta_log2_size_XI_minus_GJ", "q_size_BH", "delta_N_young_intact_XI_minus_GJ", "q_N_young_intact_BH",
  "mean_median_identity_GJ", "mean_median_identity_XI", "raw_IRR_XI_vs_GJ",
  "adjusted_IRR_XI_vs_GJ", "adjusted_q_BH", "attenuation_fraction",
  "abundance_shift_supported", "young_support_directional", "identity_directional",
  "adjusted_turnover_supported", "context_attenuated", "provisional_interpretation")
setcolorder(d, keep_cols)
setorder(d, superfamily, provisional_interpretation, family_id)
fwrite(d, file.path(intermediate_dir, "08_family_expansion_turnover_evidence.tsv"), sep = "\t")
fwrite(d, file.path(supp_dir, "Supplementary_Table_F6_family_expansion_turnover_evidence.tsv"), sep = "\t")

summary <- d[, .N, by = .(superfamily, provisional_interpretation)][order(superfamily, provisional_interpretation)]
fwrite(summary, file.path(intermediate_dir, "08_family_evidence_class_summary.tsv"), sep = "\t")
