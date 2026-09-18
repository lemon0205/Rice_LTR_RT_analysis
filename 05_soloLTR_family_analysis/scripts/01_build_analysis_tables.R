#!/usr/bin/env Rscript

# Build the fixed N_intact > 10 eligibility set and the analysis-ready
# family × accession × gene-proximity tables for Figure 6.
# Source files are read only; all outputs are written under this analysis folder.

suppressPackageStartupMessages(library(data.table))

analysis_root <- normalizePath(file.path(getwd(), "Figure6_context_adjusted_analysis"), winslash = "/", mustWork = TRUE)
project_root <- dirname(analysis_root)
data_root <- file.path(project_root, "data")
input_dir <- file.path(analysis_root, "01_input_data")
intermediate_dir <- file.path(analysis_root, "02_intermediate_data")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

family_path <- file.path(data_root, "solo-LTR_family", "family_by_accession.tsv")
element_path <- file.path(data_root, "LTR_CHH_siRNA", "analysis_tables", "LTR_epigenetic_element_master.tsv.gz")

manifest <- data.table(
  input_name = c("family_by_accession", "LTR_epigenetic_element_master"),
  path = c(normalizePath(family_path, winslash = "/"), normalizePath(element_path, winslash = "/")),
  purpose = c("Eligibility, abundance, young-intact, and LTR-identity metrics",
              "Element context, 24-nt siRNA, CHH, and solo/intact counts")
)
fwrite(manifest, file.path(input_dir, "input_file_manifest.tsv"), sep = "\t")

family <- fread(family_path)[group %in% c("GJ", "XI")]
setorder(family, family_id, accession)

# The user-specified criterion applies to total intact elements within each group.
eligible <- family[, .(
  N_intact_GJ = sum(N_intact[group == "GJ"], na.rm = TRUE),
  N_intact_XI = sum(N_intact[group == "XI"], na.rm = TRUE),
  n_accessions_GJ = uniqueN(accession[group == "GJ"]),
  n_accessions_XI = uniqueN(accession[group == "XI"])
), by = .(family_id, superfamily)]
eligible[, eligible_Nintact_gt10_both_groups := N_intact_GJ > 10 & N_intact_XI > 10]
setorder(eligible, superfamily, family_id)
fwrite(eligible, file.path(intermediate_dir, "01_family_eligibility_Nintact_gt10.tsv"), sep = "\t")

keep <- eligible[eligible_Nintact_gt10_both_groups == TRUE, family_id]
family_keep <- family[family_id %in% keep]
fwrite(family_keep, file.path(intermediate_dir, "02_family_accession_metrics_eligible.tsv.gz"), sep = "\t")

# Collapse paired intact LTR ends to one parent element. The context is assigned
# from the nearest end (minimum distance); this avoids double-counting intact TEs.
element <- as.data.table(read.delim(gzfile(element_path), check.names = FALSE))[group %in% c("GJ", "XI") & family_id %in% keep]
prox_from_distance <- function(x) fifelse(x == 0, "overlap", fifelse(x < 1000, "near_1kb", "far_1kb"))

intact_parent <- element[structure == "intact_end", .(
  superfamily = first(superfamily),
  nearest_gene_distance_bp = min(nearest_gene_distance_bp, na.rm = TRUE),
  TE_length_bp = median(length_bp, na.rm = TRUE),
  young_status = first(young_status),
  R_24nt_signal = mean(R_24nt_signal, na.rm = TRUE),
  Y_24nt_signal = mean(Y_24nt_signal, na.rm = TRUE),
  R_methylated_reads = sum(R_methylated_reads, na.rm = TRUE),
  R_total_coverage = sum(R_total_coverage, na.rm = TRUE),
  Y_methylated_reads = sum(Y_methylated_reads, na.rm = TRUE),
  Y_total_coverage = sum(Y_total_coverage, na.rm = TRUE)
), by = .(accession, group, family_id, parent_intact_id)]
intact_parent[, `:=`(structure = "intact", gene_proximity = prox_from_distance(nearest_gene_distance_bp))]

solo <- element[structure == "solo", .(
  superfamily = first(superfamily),
  nearest_gene_distance_bp = first(nearest_gene_distance_bp),
  TE_length_bp = first(length_bp),
  R_24nt_signal = first(R_24nt_signal),
  Y_24nt_signal = first(Y_24nt_signal),
  R_methylated_reads = first(R_methylated_reads),
  R_total_coverage = first(R_total_coverage),
  Y_methylated_reads = first(Y_methylated_reads),
  Y_total_coverage = first(Y_total_coverage)
), by = .(accession, group, family_id, solo_id = region_id)]
solo[, `:=`(structure = "solo", young_status = NA_character_, gene_proximity = prox_from_distance(nearest_gene_distance_bp))]

element_level <- rbindlist(list(
  intact_parent[, .(accession, group, family_id, superfamily, structure, element_id = parent_intact_id,
                    gene_proximity, nearest_gene_distance_bp, TE_length_bp, young_status,
                    R_24nt_signal, Y_24nt_signal, R_methylated_reads, R_total_coverage,
                    Y_methylated_reads, Y_total_coverage)],
  solo[, .(accession, group, family_id, superfamily, structure, element_id = solo_id,
           gene_proximity, nearest_gene_distance_bp, TE_length_bp, young_status,
           R_24nt_signal, Y_24nt_signal, R_methylated_reads, R_total_coverage,
           Y_methylated_reads, Y_total_coverage)]
), use.names = TRUE)
fwrite(element_level, file.path(intermediate_dir, "03_eligible_element_context_table.tsv.gz"), sep = "\t")

# Count table for the context-adjusted structural-turnover model.  Each
# denominator is a count of unique intact parent elements in the same context.
counts <- element_level[, .(N_elements = uniqueN(element_id)), by = .(accession, group, family_id, superfamily, structure, gene_proximity)]
counts_wide <- dcast(counts, accession + group + family_id + superfamily + gene_proximity ~ structure,
                     value.var = "N_elements", fill = 0)
for (nm in c("intact", "solo")) if (!nm %in% names(counts_wide)) counts_wide[, (nm) := 0L]
setnames(counts_wide, c("intact", "solo"), c("N_intact_context", "N_solo_context"))

# Exposure rows without intact elements cannot contribute to an offset model.
counts_model <- counts_wide[N_intact_context > 0]
fwrite(counts_wide, file.path(intermediate_dir, "04_family_accession_context_raw_counts.tsv.gz"), sep = "\t")
fwrite(counts_model, file.path(intermediate_dir, "05_family_accession_context_model_counts.tsv.gz"), sep = "\t")

# Context-stratified summaries for 24-nt siRNA and CHH; these are separate from
# the count model because they represent present-day signals, not removal rates.
signals <- element_level[, .(
  n_elements = uniqueN(element_id),
  median_gene_distance_bp = as.numeric(median(nearest_gene_distance_bp, na.rm = TRUE)),
  median_TE_length_bp = as.numeric(median(TE_length_bp, na.rm = TRUE)),
  mean_R_24nt_signal = mean(R_24nt_signal, na.rm = TRUE),
  mean_Y_24nt_signal = mean(Y_24nt_signal, na.rm = TRUE),
  R_methylated_reads = sum(R_methylated_reads, na.rm = TRUE),
  R_total_coverage = sum(R_total_coverage, na.rm = TRUE),
  Y_methylated_reads = sum(Y_methylated_reads, na.rm = TRUE),
  Y_total_coverage = sum(Y_total_coverage, na.rm = TRUE)
), by = .(accession, group, family_id, superfamily, structure, gene_proximity)]
signals[, `:=`(R_weighted_CHH = fifelse(R_total_coverage > 0, R_methylated_reads / R_total_coverage, NA_real_),
               Y_weighted_CHH = fifelse(Y_total_coverage > 0, Y_methylated_reads / Y_total_coverage, NA_real_))]
fwrite(signals, file.path(intermediate_dir, "06_family_accession_context_epigenetic_summary.tsv.gz"), sep = "\t")

qc <- data.table(
  metric = c("families_total_GJ_XI", "families_eligible_Nintact_gt10", "eligible_family_accession_rows",
             "element_rows_eligible", "intact_parent_elements", "solo_elements", "count_model_rows"),
  value = c(uniqueN(family$family_id), length(keep), nrow(family_keep), nrow(element_level),
            nrow(intact_parent), nrow(solo), nrow(counts_model))
)
fwrite(qc, file.path(intermediate_dir, "00_build_table_QC.tsv"), sep = "\t")
