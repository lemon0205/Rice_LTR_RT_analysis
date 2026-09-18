#!/usr/bin/env Rscript

# Integrate abundance/age evidence with the GJ-XI context decomposition.
# Thresholds are effect-size rules for transparent descriptive classification,
# not per-family hypothesis tests: 1.5-fold structural contrast and >=50%
# attenuation after common-composition standardisation.

suppressPackageStartupMessages(library(data.table))

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

dec <- fread(file.path(intermediate_dir, "10_family_GJ_XI_context_decomposition.tsv"))
dyn <- fread(file.path(project_root, "family_dynamics_nintact10", "family_dynamics_GJ_XI_Nintact10.tsv"))
dyn <- dyn[, .(family_id, delta_log2_size_XI_minus_GJ, q_size_BH,
               delta_N_young_intact_XI_minus_GJ, q_N_young_intact_BH,
               mean_median_identity_GJ, mean_median_identity_XI)]
d <- merge(dec, dyn, by = "family_id", all.x = TRUE)

structural_cutoff <- log(1.5)
d[, `:=`(
  abundance_shift = !is.na(q_size_BH) & q_size_BH < 0.05,
  age_shift = !is.na(q_N_young_intact_BH) & q_N_young_intact_BH < 0.05 &
    sign(delta_N_young_intact_XI_minus_GJ) == sign(delta_log2_size_XI_minus_GJ),
  identity_shift = !is.na(mean_median_identity_XI) &
    sign(mean_median_identity_XI - mean_median_identity_GJ) == sign(delta_log2_size_XI_minus_GJ),
  raw_structural_divergence = is.finite(raw_log_IRR_XI_vs_GJ) & abs(raw_log_IRR_XI_vs_GJ) >= structural_cutoff
)]
d[, `:=`(
  context_major = is.finite(context_associated_fraction) & raw_structural_divergence & context_associated_fraction >= 0.50,
  within_context_major = is.finite(standardized_log_IRR_XI_vs_GJ) & abs(standardized_log_IRR_XI_vs_GJ) >= structural_cutoff
)]
d[, expansion_evidence := abundance_shift & age_shift & identity_shift]
d[, provisional_trajectory := fcase(
  expansion_evidence & context_major, "Expansion plus context-associated turnover",
  expansion_evidence & within_context_major, "Expansion plus within-context turnover",
  expansion_evidence, "Expansion-driven",
  context_major, "Context-associated turnover",
  within_context_major, "Within-context turnover",
  abundance_shift, "Abundance divergence; mechanism unresolved",
  default = "No supported abundance divergence"
)]

keep <- c("family_id", "superfamily", "N_intact_GJ.x", "N_intact_XI.x",
          "delta_log2_size_XI_minus_GJ", "q_size_BH", "delta_N_young_intact_XI_minus_GJ", "q_N_young_intact_BH",
          "mean_median_identity_GJ", "mean_median_identity_XI", "raw_log_IRR_XI_vs_GJ",
          "standardized_log_IRR_XI_vs_GJ", "context_associated_log_component", "context_associated_fraction",
          "abundance_shift", "age_shift", "identity_shift", "expansion_evidence",
          "raw_structural_divergence", "context_major", "within_context_major", "provisional_trajectory")
d <- d[, ..keep]
setnames(d, c("N_intact_GJ.x", "N_intact_XI.x"), c("N_intact_GJ", "N_intact_XI"))
setorder(d, superfamily, provisional_trajectory, family_id)
fwrite(d, file.path(intermediate_dir, "12_integrated_family_trajectories.tsv"), sep = "\t")
fwrite(d, file.path(supp_dir, "Supplementary_Table_F6_integrated_family_trajectories.tsv"), sep = "\t")
summary <- d[, .N, by = .(superfamily, provisional_trajectory)][order(superfamily, provisional_trajectory)]
fwrite(summary, file.path(intermediate_dir, "12_integrated_trajectory_summary.tsv"), sep = "\t")
