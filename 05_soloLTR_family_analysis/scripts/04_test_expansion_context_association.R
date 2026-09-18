#!/usr/bin/env Rscript

# Test whether expansion-driven candidates occupy distinct present-day genomic
# contexts. Formal inference is limited to expansion evidence (8 families) vs
# all remaining eligible families; removal (n=1) and mixed (n=2) are exported
# descriptively only because their sample sizes do not support a stable model.

suppressPackageStartupMessages({library(data.table); library(glmmTMB); library(lme4)})

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

evidence <- fread(file.path(intermediate_dir, "08_family_expansion_turnover_evidence.tsv"))
element <- as.data.table(read.delim(gzfile(file.path(intermediate_dir, "03_eligible_element_context_table.tsv.gz")), check.names = FALSE))
element <- merge(element, evidence[, .(family_id, provisional_interpretation, expansion_evidence)], by = "family_id", all.x = TRUE)
element[, expansion_group := factor(ifelse(expansion_evidence, "Expansion-driven", "All other eligible families"),
                                    levels = c("All other eligible families", "Expansion-driven"))]

# Descriptive table retained for every provisional class, including rare classes.
desc <- element[, .(
  n_elements = uniqueN(element_id), n_families = uniqueN(family_id),
  median_gene_distance_bp = as.numeric(median(nearest_gene_distance_bp, na.rm = TRUE)),
  median_TE_length_bp = as.numeric(median(TE_length_bp, na.rm = TRUE)),
  mean_R_24nt = mean(R_24nt_signal, na.rm = TRUE), mean_Y_24nt = mean(Y_24nt_signal, na.rm = TRUE),
  R_methylated = sum(R_methylated_reads, na.rm = TRUE), R_coverage = sum(R_total_coverage, na.rm = TRUE),
  Y_methylated = sum(Y_methylated_reads, na.rm = TRUE), Y_coverage = sum(Y_total_coverage, na.rm = TRUE)
), by = .(provisional_interpretation, superfamily, group, structure, gene_proximity)]
desc[, `:=`(R_weighted_CHH = fifelse(R_coverage > 0, R_methylated / R_coverage, NA_real_),
            Y_weighted_CHH = fifelse(Y_coverage > 0, Y_methylated / Y_coverage, NA_real_))]
fwrite(desc, file.path(intermediate_dir, "09_all_classes_context_epigenetic_descriptives.tsv"), sep = "\t")

# One observation per family/accession/structure: avoids element-level
# pseudoreplication in the formal expansion-versus-other models.
fa <- element[, .(
  n_elements = uniqueN(element_id),
  n_far = uniqueN(element_id[gene_proximity == "far_1kb"]),
  median_log_distance = as.numeric(median(log1p(nearest_gene_distance_bp), na.rm = TRUE)),
  median_log_length = as.numeric(median(log(TE_length_bp), na.rm = TRUE)),
  mean_R_24nt = mean(R_24nt_signal, na.rm = TRUE), mean_Y_24nt = mean(Y_24nt_signal, na.rm = TRUE),
  R_methylated = sum(R_methylated_reads, na.rm = TRUE), R_coverage = sum(R_total_coverage, na.rm = TRUE),
  Y_methylated = sum(Y_methylated_reads, na.rm = TRUE), Y_coverage = sum(Y_total_coverage, na.rm = TRUE)
), by = .(family_id, accession, group, superfamily, structure, expansion_group)]
fa[, `:=`(family_id = factor(family_id), accession = factor(accession), group = factor(group),
          superfamily = factor(superfamily), structure = factor(structure))]
fwrite(fa, file.path(intermediate_dir, "09_family_accession_expansion_context_model_data.tsv.gz"), sep = "\t")

extract_lmer <- function(fit, outcome) {
  z <- as.data.frame(summary(fit)$coefficients)
  term <- "expansion_groupExpansion-driven"
  data.table(outcome = outcome, estimate = z[term, "Estimate"], SE = z[term, "Std. Error"],
             statistic = z[term, "t value"], p_value = 2 * pnorm(abs(z[term, "t value"]), lower.tail = FALSE))
}
extract_glmm <- function(fit, outcome) {
  z <- as.data.frame(summary(fit)$coefficients$cond)
  term <- "expansion_groupExpansion-driven"
  data.table(outcome = outcome, estimate = z[term, "Estimate"], SE = z[term, "Std. Error"],
             statistic = z[term, "z value"], p_value = z[term, "Pr(>|z|)"])
}

models <- list(); results <- list()
models$far <- glmmTMB(cbind(n_far, n_elements - n_far) ~ expansion_group + group + superfamily + structure +
                        (1 | family_id) + (1 | accession), family = betabinomial(link = "logit"), data = fa)
results$far <- extract_glmm(models$far, "far_fraction")
models$distance <- lmer(median_log_distance ~ expansion_group + group + superfamily + structure +
                          (1 | family_id), data = fa, REML = FALSE)
results$distance <- extract_lmer(models$distance, "log1p_gene_distance")
models$length <- lmer(median_log_length ~ expansion_group + group + superfamily + structure +
                        (1 | family_id), data = fa, REML = FALSE)
results$length <- extract_lmer(models$length, "log_TE_length")
for (condition in c("R", "Y")) {
  signal <- paste0("mean_", condition, "_24nt")
  meth <- paste0(condition, "_methylated")
  cov <- paste0(condition, "_coverage")
  x <- copy(fa)
  x[, log_signal := log1p(get(signal))]
  random_term <- if (condition == "R") "(1 | accession)" else "(1 | family_id) + (1 | accession)"
  models[[paste0("s24_", condition)]] <- lmer(as.formula(paste("log_signal ~ expansion_group + group + superfamily + structure +", random_term)), data = x, REML = FALSE)
  results[[paste0("s24_", condition)]] <- extract_lmer(models[[paste0("s24_", condition)]], paste0("24nt_", condition))
  y <- x[get(cov) > 0]
  y[, `:=`(meth = get(meth), unmeth = get(cov) - get(meth))]
  models[[paste0("chh_", condition)]] <- glmmTMB(cbind(meth, unmeth) ~ expansion_group + group + superfamily + structure +
                                                    (1 | family_id) + (1 | accession), family = betabinomial(link = "logit"), data = y)
  results[[paste0("chh_", condition)]] <- extract_glmm(models[[paste0("chh_", condition)]], paste0("CHH_", condition))
}

out <- rbindlist(results)
out[, `:=`(q_BH = p.adjust(p_value, method = "BH"), effect_ratio = exp(estimate))]
fwrite(out, file.path(intermediate_dir, "09_expansion_vs_other_context_epigenetic_effects.tsv"), sep = "\t")
fwrite(out, file.path(supp_dir, "Supplementary_Table_F6_expansion_context_models.tsv"), sep = "\t")
saveRDS(models, file.path(intermediate_dir, "09_expansion_context_models.rds"))
