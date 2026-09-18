#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(glmmTMB)})
root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
panel <- Sys.getenv("SOLO_LTR_CHH_INPUT_DIR")
if (!nzchar(panel)) stop("Set SOLO_LTR_CHH_INPUT_DIR to the directory containing S7D_dynamic_family_solo_intact_CHH_within_accession_data.tsv.")
d <- fread(file.path(panel, "S7D_dynamic_family_solo_intact_CHH_within_accession_data.tsv"))
d[, `:=`(family_id = factor(family_id), accession = factor(accession),
          gene_proximity = factor(gene_proximity), structure = factor(structure, levels = c("Intact LTR end", "Solo-LTR")),
          superfamily = factor(superfamily), dynamics_class = factor(dynamics_class, levels = c("Amplification-like", "Removal-like")))]
m <- glmmTMB(cbind(meth, unmeth) ~ structure * dynamics_class + superfamily + gene_proximity + (1 | family_id) + (1 | accession),
             family = betabinomial(link = "logit"), data = d)
co <- summary(m)$coefficients$cond; V <- vcov(m)$cond
effect <- function(cl) {
  v <- setNames(rep(0, nrow(co)), rownames(co)); v["structureSolo-LTR"] <- 1
  if (cl == "Removal-like") v["structureSolo-LTR:dynamics_classRemoval-like"] <- 1
  est <- sum(v * co[, "Estimate"]); se <- sqrt(as.numeric(t(v) %*% V %*% v))
  data.table(dynamics_class = cl, estimate_solo_minus_intact = est, SE = se,
             lower = est - 1.96 * se, upper = est + 1.96 * se,
             p_value = 2 * pnorm(-abs(est / se)), odds_ratio = exp(est))
}
out <- rbindlist(lapply(c("Amplification-like", "Removal-like"), effect))
out[, q_BH := p.adjust(p_value, method = "BH")]
interaction <- co["structureSolo-LTR:dynamics_classRemoval-like", ]
extra <- data.table(term = "solo_LTR_by_removal_interaction", estimate = interaction["Estimate"], SE = interaction["Std. Error"], p_value = interaction["Pr(>|z|)"])
fwrite(out, file.path(panel, "S7D_solo_intact_CHH_by_dynamic_class.tsv"), sep = "\t")
fwrite(extra, file.path(panel, "S7D_solo_intact_CHH_dynamic_class_interaction.tsv"), sep = "\t")
print(out); print(extra)
