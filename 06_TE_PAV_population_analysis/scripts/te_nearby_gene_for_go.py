import csv
import os
from collections import defaultdict


HIGH_FST_FILE = "output/dte_fst/high_fst_te_candidates.csv"
SPECIFIC_TE_FILE = "output/table2_specific_TEs.csv"
GENE_BED_FILE = "genome/MSU_genes.bed"
OUTPUT_DIR = "GO_result"
MAX_DISTANCE = 2000


def ensure_output_dir(path):
    os.makedirs(path, exist_ok=True)


def parse_gene_id(attr_field):
    for item in attr_field.split(";"):
        if item.startswith("ID="):
            return item.split("=", 1)[1]
    return attr_field


def load_genes(gene_bed_file):
    genes_by_chrom = defaultdict(list)
    with open(gene_bed_file, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 4:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            gene_id = parse_gene_id(row[3])
            genes_by_chrom[chrom].append({
                "Gene_ID": gene_id,
                "Gene_Chrom": chrom,
                "Gene_Start": start,
                "Gene_End": end,
            })
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda x: (x["Gene_Start"], x["Gene_End"], x["Gene_ID"]))
    return genes_by_chrom


def relation_and_distance(te_start, te_end, gene_start, gene_end):
    if gene_end >= te_start and gene_start <= te_end:
        return "overlap", 0
    if gene_end < te_start:
        return "upstream_within_2kb", te_start - gene_end
    return "downstream_within_2kb", gene_start - te_end


def find_best_gene(te_chrom, te_start, te_end, genes_by_chrom):
    best = None
    for gene in genes_by_chrom.get(te_chrom, []):
        relation, distance = relation_and_distance(te_start, te_end, gene["Gene_Start"], gene["Gene_End"])
        if distance > MAX_DISTANCE:
            if gene["Gene_Start"] > te_end and (gene["Gene_Start"] - te_end) > MAX_DISTANCE:
                break
            continue
        candidate = {
            **gene,
            "Distance_bp": distance,
            "Relation": relation,
        }
        if best is None:
            best = candidate
            continue
        if candidate["Distance_bp"] < best["Distance_bp"]:
            best = candidate
            continue
        if candidate["Distance_bp"] == best["Distance_bp"] and candidate["Gene_Start"] < best["Gene_Start"]:
            best = candidate
    return best


def annotate_high_fst(high_fst_file, genes_by_chrom):
    records = []
    with open(high_fst_file, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            te_chrom = row["Chromosome"]
            te_start = int(row["Start"])
            te_end = int(row["End"])
            best_gene = find_best_gene(te_chrom, te_start, te_end, genes_by_chrom)
            if not best_gene:
                continue
            records.append({
                "TE_ID": row["TE_ID"],
                "Family": row["Family"],
                "Chromosome": te_chrom,
                "Start": te_start,
                "End": te_end,
                "fst": row["fst"],
                "delta_freq": row["delta_freq"],
                "group_bias": row["group_bias"],
                "Gene_ID": best_gene["Gene_ID"],
                "Gene_Chrom": best_gene["Gene_Chrom"],
                "Gene_Start": best_gene["Gene_Start"],
                "Gene_End": best_gene["Gene_End"],
                "Distance_bp": best_gene["Distance_bp"],
                "Relation": best_gene["Relation"],
            })
    return records


def annotate_specific(specific_te_file, genes_by_chrom):
    records = []
    with open(specific_te_file, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            te_chrom = row["Chrom"]
            te_start = int(row["Start"])
            te_end = int(row["End"])
            best_gene = find_best_gene(te_chrom, te_start, te_end, genes_by_chrom)
            if not best_gene:
                continue
            records.append({
                "TE_ID": row["TE_ID"],
                "Family": row["Family"],
                "Chrom": te_chrom,
                "Start": te_start,
                "End": te_end,
                "Source_Variety": row["Source_Variety"],
                "Type": row["Type"],
                "Gene_ID": best_gene["Gene_ID"],
                "Gene_Chrom": best_gene["Gene_Chrom"],
                "Gene_Start": best_gene["Gene_Start"],
                "Gene_End": best_gene["Gene_End"],
                "Distance_bp": best_gene["Distance_bp"],
                "Relation": best_gene["Relation"],
            })
    return records


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_gene_list(path, rows):
    seen = set()
    ordered = []
    for row in rows:
        gene_id = row["Gene_ID"]
        if gene_id not in seen:
            seen.add(gene_id)
            ordered.append(gene_id)
    with open(path, "w", newline="") as handle:
        for gene_id in ordered:
            handle.write(f"{gene_id}\n")


def main():
    ensure_output_dir(OUTPUT_DIR)
    genes_by_chrom = load_genes(GENE_BED_FILE)

    high_fst_rows = annotate_high_fst(HIGH_FST_FILE, genes_by_chrom)
    specific_rows = annotate_specific(SPECIFIC_TE_FILE, genes_by_chrom)

    high_fst_fields = [
        "TE_ID", "Family", "Chromosome", "Start", "End",
        "fst", "delta_freq", "group_bias",
        "Gene_ID", "Gene_Chrom", "Gene_Start", "Gene_End",
        "Distance_bp", "Relation"
    ]
    specific_fields = [
        "TE_ID", "Family", "Chrom", "Start", "End",
        "Source_Variety", "Type",
        "Gene_ID", "Gene_Chrom", "Gene_Start", "Gene_End",
        "Distance_bp", "Relation"
    ]

    write_csv(os.path.join(OUTPUT_DIR, "high_fst_te_nearby_genes.csv"), high_fst_fields, high_fst_rows)
    write_csv(os.path.join(OUTPUT_DIR, "specific_te_nearby_genes.csv"), specific_fields, specific_rows)
    write_gene_list(os.path.join(OUTPUT_DIR, "high_fst_te_gene_list_unique.txt"), high_fst_rows)
    write_gene_list(os.path.join(OUTPUT_DIR, "specific_te_gene_list_unique.txt"), specific_rows)

    print(f"High-FST TE with nearby genes: {len(high_fst_rows)}")
    print(f"Specific TE with nearby genes: {len(specific_rows)}")


if __name__ == "__main__":
    main()
