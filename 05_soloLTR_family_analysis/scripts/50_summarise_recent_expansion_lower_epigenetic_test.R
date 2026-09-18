#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
inter <- file.path(root, "02_intermediate_data")
files <- c(
  intact_all = file.path(inter, "48_amplification_age_stratified_5prime_epigenetic_family_effects.tsv"),
  intact_young = file.path(inter, "48_amplification_age_stratified_5prime_epigenetic_family_effects_young_intact_only.tsv"),
  solo = file.path(inter, "49_recent_amplification_solo_CHH_sRNA_family_effects.tsv")
)
out <- rbindlist(lapply(names(files), function(nm) {
  d <- fread(files[[nm]])[age_landscape_class == "Young-peaked"]
  d[, {
    x <- higher_minus_other[is.finite(higher_minus_other)]
    list(analysis = nm, n_families = length(x),
         median_higher_minus_other = median(x), n_positive = sum(x > 0),
         p_lower = wilcox.test(x, alternative = "less", exact = TRUE)$p.value)
  }, by = metric]
}))
out[, q_BH_CHH_and_sRNA := p.adjust(p_lower, method = "BH"), by = analysis]
fwrite(out, file.path(root, "05_supplementary_tables", "Supplementary_Table_S7_recent_expansion_lower_epigenetic_test.tsv"), sep = "\t")
print(out)
