#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(glmmTMB)})
root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
d <- as.data.table(read.delim(gzfile(file.path(root, "02_intermediate_data", "30_new_solo_intact_context_epigenetic_model_data.tsv.gz")), check.names = FALSE))
dyn <- fread(file.path(root, "02_intermediate_data", "20_context_aware_family_dynamics.tsv"))[dynamics_class %in% c("Amplification-like", "Removal-like"), family_id]
d <- d[family_id %in% dyn]
d[, `:=`(family_id = factor(family_id), accession = factor(accession), group = factor(group),
          gene_proximity = factor(gene_proximity), structure = factor(structure, levels = c("Intact LTR end", "Solo-LTR")),
          superfamily = factor(superfamily, levels = c("Copia", "Gypsy")))]
m <- glmmTMB(cbind(meth, unmeth) ~ structure * superfamily + group + gene_proximity + (1 | family_id) + (1 | accession),
             family = betabinomial(link = "logit"), data = d)
co <- summary(m)$coefficients$cond; V <- vcov(m)$cond
out <- rbindlist(lapply(c("Copia", "Gypsy"), function(sf) {
  v <- setNames(rep(0, nrow(co)), rownames(co)); v["structureSolo-LTR"] <- 1
  if (sf == "Gypsy") v["structureSolo-LTR:superfamilyGypsy"] <- 1
  est <- sum(v * co[, "Estimate"]); se <- sqrt(as.numeric(t(v) %*% V %*% v))
  data.table(superfamily = sf, estimate_solo_minus_intact = est, SE = se, p_value = 2 * pnorm(-abs(est / se)))
}))
out[, q_BH := p.adjust(p_value, "BH")]
print(out)
fwrite(out, file.path(root, "02_intermediate_data", "32_dynamic_solo_intact_CHH_old_data_group_adjusted_verification.tsv"), sep = "\t")
