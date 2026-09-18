#!/usr/bin/env Rscript

# Primary Figure 6 analysis: decompose each family's GJ-XI solo/intact
# difference into a context-associated component and a within-context component.
# Counts are modelled as solo-LTR events per intact parent element (NB with an
# intact-count offset); this is a rate model, not a deletion-rate model.

suppressPackageStartupMessages({library(data.table); library(glmmTMB)})

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
supp_dir <- file.path(analysis_root, "05_supplementary_tables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

in_path <- file.path(intermediate_dir, "04_family_accession_context_raw_counts.tsv.gz")
d <- as.data.table(read.delim(gzfile(in_path), check.names = FALSE))
d <- d[group %in% c("GJ", "XI") & N_intact_context > 0]
d[, `:=`(group = factor(group, levels = c("GJ", "XI")), family_id = factor(family_id),
          accession = factor(accession),
          gene_proximity = factor(gene_proximity, levels = c("far_1kb", "near_1kb", "overlap")),
          log_intact_exposure = log(N_intact_context))]

# Family-specific GJ-XI slopes are partially pooled, avoiding unstable separate
# fixed effects for all 86 families. A shared group-by-proximity term supplies
# the within-context standardisation surface; accession is a random intercept.
fit <- glmmTMB(
  N_solo_context ~ group * gene_proximity + (1 + group || family_id) + (1 | accession) +
    offset(log_intact_exposure),
  family = nbinom2, data = d
)

beta <- fixef(fit)$cond
V <- vcov(fit)$cond
families <- levels(d$family_id)
contexts <- levels(d$gene_proximity)

make_X <- function(family, group, context) {
  nd <- data.frame(
    family_id = factor(family, levels = families),
    group = factor(group, levels = levels(d$group)),
    gene_proximity = factor(context, levels = contexts)
  )
  X <- model.matrix(~ group * gene_proximity, nd)
  full <- matrix(0, nrow = 1, ncol = length(beta), dimnames = list(NULL, names(beta)))
  common <- intersect(colnames(X), names(beta))
  full[, common] <- X[, common, drop = FALSE]
  full
}

family_re <- ranef(fit)$cond$family_id
family_re_intercept <- function(f) {
  if ("(Intercept)" %in% colnames(family_re)) family_re[f, "(Intercept)"] else 0
}
family_re_groupXI <- function(f) {
  candidate <- grep("groupXI", colnames(family_re), value = TRUE)
  if (length(candidate) > 0) family_re[f, candidate[1]] else 0
}

# Common weights are family-specific but group-invariant: pooled GJ+XI intact
# composition. Therefore the standardised contrast answers the counterfactual
# question "what if both groups had the same proximity composition?".
composition <- d[, .(N_intact_context = sum(N_intact_context)), by = .(family_id, group, gene_proximity)]
composition <- composition[, .(N_intact_context = sum(N_intact_context)), by = .(family_id, gene_proximity)]
composition[, weight := N_intact_context / sum(N_intact_context), by = family_id]

raw <- d[, .(N_solo = sum(N_solo_context), N_intact = sum(N_intact_context)), by = .(family_id, group)]
raw <- dcast(raw, family_id ~ group, value.var = c("N_solo", "N_intact"), fill = 0)
raw[, `:=`(raw_rate_GJ = N_solo_GJ / N_intact_GJ,
            raw_rate_XI = N_solo_XI / N_intact_XI)]
raw[, raw_log_IRR_XI_vs_GJ := fifelse(raw_rate_GJ > 0 & raw_rate_XI > 0,
                                      log(raw_rate_XI / raw_rate_GJ), NA_real_)]

decompose_one <- function(f) {
  w <- composition[family_id == f][match(contexts, gene_proximity), weight]
  w[is.na(w)] <- 0
  if (sum(w) == 0) return(data.table(family_id = f))
  w <- w / sum(w)
  Xg <- do.call(rbind, lapply(contexts, function(c) make_X(f, "GJ", c)))
  Xx <- do.call(rbind, lapply(contexts, function(c) make_X(f, "XI", c)))
  eta_g <- drop(Xg %*% beta) + family_re_intercept(f)
  eta_x <- drop(Xx %*% beta) + family_re_intercept(f) + family_re_groupXI(f)
  rate_g <- sum(w * exp(eta_g)); rate_x <- sum(w * exp(eta_x))
  standardized <- log(rate_x / rate_g)
  data.table(family_id = f, standardized_log_IRR_XI_vs_GJ = standardized,
             standardized_IRR_XI_vs_GJ = exp(standardized))
}

standardized <- rbindlist(lapply(families, decompose_one), fill = TRUE)
out <- merge(raw, standardized, by = "family_id", all.x = TRUE)
out[, `:=`(
  context_associated_log_component = raw_log_IRR_XI_vs_GJ - standardized_log_IRR_XI_vs_GJ,
  context_associated_fraction = fifelse(!is.na(raw_log_IRR_XI_vs_GJ) & abs(raw_log_IRR_XI_vs_GJ) > 0,
                                       (raw_log_IRR_XI_vs_GJ - standardized_log_IRR_XI_vs_GJ) / raw_log_IRR_XI_vs_GJ,
                                       NA_real_)
)]

meta <- fread(file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"))
out <- merge(meta[eligible_Nintact_gt10_both_groups == TRUE], out, by = "family_id", all.y = TRUE)
out[, sort_abs_context_component := abs(context_associated_log_component)]
setorder(out, -sort_abs_context_component)
out[, sort_abs_context_component := NULL]
fwrite(out, file.path(intermediate_dir, "10_family_GJ_XI_context_decomposition.tsv"), sep = "\t")
fwrite(out, file.path(supp_dir, "Supplementary_Table_F6_family_context_decomposition.tsv"), sep = "\t")

weights_out <- dcast(composition, family_id ~ gene_proximity, value.var = "weight", fill = 0)
fwrite(weights_out, file.path(intermediate_dir, "10_family_pooled_context_weights.tsv"), sep = "\t")
qc <- data.table(metric = c("count_rows", "families_in_model", "convergence_code", "pdHess", "AIC", "families_standardized"),
                 value = c(nrow(d), uniqueN(d$family_id), fit$fit$convergence, fit$sdr$pdHess, AIC(fit), nrow(standardized)))
fwrite(qc, file.path(intermediate_dir, "10_context_decomposition_model_QC.tsv"), sep = "\t")
saveRDS(fit, file.path(intermediate_dir, "10_context_decomposition_model.rds"))
