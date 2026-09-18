options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("PROJECT_DATA_ROOT")
if (!nzchar(root)) stop("Set PROJECT_DATA_ROOT to the root containing Supplementary figure/ and Figure7/.")
s8_dir <- file.path(root, "Supplementary figure", "S8")
out_dir <- file.path(s8_dir, "GO_2kb_gene_lists")
list_dir <- file.path(out_dir, "gene_lists")
table_dir <- file.path(out_dir, "tables")
dir.create(list_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

te_file <- file.path(s8_dir, "data", "te_four_groups.csv")
gene_bed <- file.path(root, "Figure7", "pan-TE", "genome", "MSU_genes.bed")

class_levels <- c("Shared-TE", "Low-Fst TE", "Mid-Fst TE", "High-Fst TE")
safe_stem <- c(
  "Shared-TE" = "Shared_TE",
  "Low-Fst TE" = "Low_Fst_TE",
  "Mid-Fst TE" = "Mid_Fst_TE",
  "High-Fst TE" = "High_Fst_TE"
)

te <- fread(te_file)
te <- te[Group != "msu_specific_te" & !is.na(fst)]
te[, TE_Class := fifelse(
  Group == "all_shared_te", "Shared-TE",
  fifelse(fst == 1, "High-Fst TE",
    fifelse(fst >= 0.25, "Mid-Fst TE", "Low-Fst TE")
  )
)]
te[, `:=`(
  Chromosome = ifelse(grepl("^Chr", Chromosome), Chromosome, paste0("Chr", Chromosome)),
  TE_Start = as.integer(Start),
  TE_End = as.integer(End),
  start = pmax(0L, as.integer(Start) - 2000L),
  end = as.integer(End) + 2000L
)]
te_int <- te[, .(Chromosome, start, end, TE_ID, TE_Class, Family, TE_Start, TE_End)]

genes <- fread(gene_bed, header = FALSE, sep = "\t")
setnames(genes, c("Chromosome", "Gene_Start", "Gene_End", "Attributes", "Score", "Strand"))
genes[, Gene_ID := sub("^.*ID=([^;]+).*$", "\\1", Attributes)]
genes[, Gene_Name := sub("^.*Name=([^;]+).*$", "\\1", Attributes)]
genes[, Gene_Note := sub("^.*Note=([^;]+).*$", "\\1", Attributes)]
genes[!grepl("Note=", Attributes), Gene_Note := ""]
genes[, Gene_Note := URLdecode(Gene_Note)]
genes[, `:=`(start = as.integer(Gene_Start), end = as.integer(Gene_End))]
gene_int <- genes[, .(Chromosome, start, end, Gene_ID, Gene_Name, Gene_Note,
                      Gene_Start = as.integer(Gene_Start), Gene_End = as.integer(Gene_End), Strand)]

setkey(te_int, Chromosome, start, end)
setkey(gene_int, Chromosome, start, end)
pairs <- foverlaps(gene_int, te_int, type = "any", nomatch = 0L)
pairs[, Distance_bp := fifelse(
  Gene_End < TE_Start, TE_Start - Gene_End,
  fifelse(Gene_Start > TE_End, Gene_Start - TE_End, 0L)
)]
pairs <- pairs[Distance_bp <= 2000]
pairs[, Relation := fifelse(
  Distance_bp == 0, "overlap",
  fifelse(Gene_End < TE_Start, "left_within_2kb", "right_within_2kb")
)]
pairs <- unique(pairs[, .(
  TE_ID, TE_Class, Family, Chromosome,
  TE_Start, TE_End, Gene_ID, Gene_Name, Gene_Start, Gene_End, Strand,
  Distance_bp, Relation, Gene_Note
)])
setorder(pairs, TE_Class, Chromosome, TE_Start, Gene_ID)

fwrite(pairs, file.path(table_dir, "all_four_TE_classes_gene_pairs_within_2kb.csv"), bom = TRUE)

summary_rows <- list()
for (cls in class_levels) {
  cls_pairs <- pairs[TE_Class == cls]
  gene_table <- unique(cls_pairs[, .(Gene_ID, Gene_Name, Chromosome, Gene_Start,
                                     Gene_End, Strand, Gene_Note)])
  setorder(gene_table, Chromosome, Gene_Start, Gene_ID)
  stem <- safe_stem[[cls]]
  fwrite(gene_table, file.path(table_dir, paste0(stem, "_genes_within_2kb.csv")), bom = TRUE)
  writeLines(gene_table$Gene_ID, file.path(list_dir, paste0(stem, "_gene_ids_2kb.txt")))
  summary_rows[[cls]] <- data.table(
    TE_Class = cls,
    TE_loci = uniqueN(te_int[TE_Class == cls, TE_ID]),
    TE_gene_pairs = nrow(cls_pairs),
    Unique_genes_2kb = nrow(gene_table)
  )
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, "four_TE_class_2kb_gene_summary.csv"), bom = TRUE)

all_genes <- unique(genes[, .(Gene_ID, Gene_Name, Chromosome, Gene_Start, Gene_End, Strand, Gene_Note)])
setorder(all_genes, Chromosome, Gene_Start, Gene_ID)
fwrite(all_genes, file.path(table_dir, "background_all_MSU_genes.csv"), bom = TRUE)
writeLines(all_genes$Gene_ID, file.path(list_dir, "background_all_MSU_gene_ids.txt"))

te_nearby_background <- unique(pairs[, .(Gene_ID, Gene_Name, Chromosome, Gene_Start,
                                         Gene_End, Strand, Gene_Note)])
setorder(te_nearby_background, Chromosome, Gene_Start, Gene_ID)
fwrite(te_nearby_background,
       file.path(table_dir, "background_all_four_class_TE_nearby_genes.csv"), bom = TRUE)
writeLines(te_nearby_background$Gene_ID,
           file.path(list_dir, "background_all_four_class_TE_nearby_gene_ids.txt"))

writeLines(c(
  "Four TE-class gene lists for GO analysis",
  "Definition: every MSU gene overlapping or within 2 kb of an analyzed MSU-reference TE locus.",
  "TE classes match Figure 7: Shared-TE, Low-Fst TE, Mid-Fst TE, High-Fst TE.",
  "Recommended enrichment background: background_all_four_class_TE_nearby_gene_ids.txt.",
  "Alternative broad background: background_all_MSU_gene_ids.txt.",
  "The TE-nearby background is preferred because it controls the genomic-context bias introduced by selecting genes near TEs."
), file.path(out_dir, "README_GO_gene_lists.txt"))

print(summary_dt)
message("Saved four-class 2-kb gene lists to: ", out_dir)
