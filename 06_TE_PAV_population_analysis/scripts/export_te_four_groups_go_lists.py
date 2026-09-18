import csv
import os
from collections import defaultdict


INPUT_FILE = "output/te_four_groups/te_four_groups_distance_summary.csv"
OUTPUT_DIR = "GO_result"
GROUP_ORDER = ["fst_1_te", "msu_specific_te", "all_shared_te", "other_te"]


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def main():
    ensure_dir(OUTPUT_DIR)

    grouped_rows = defaultdict(list)
    grouped_gene_ids = defaultdict(list)
    seen_ids = defaultdict(set)

    with open(INPUT_FILE, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            group = row["Group"]
            if group not in GROUP_ORDER:
                continue
            grouped_rows[group].append(row)
            gene_id = row["Nearest_Gene_ID"]
            if gene_id and gene_id not in seen_ids[group]:
                seen_ids[group].add(gene_id)
                grouped_gene_ids[group].append(gene_id)

    for group in GROUP_ORDER:
        csv_path = os.path.join(OUTPUT_DIR, f"{group}_nearby_genes.csv")
        txt_path = os.path.join(OUTPUT_DIR, f"{group}_nearby_gene_ids.txt")

        with open(csv_path, "w", newline="") as handle:
            fieldnames = [
                "TE_ID",
                "Group",
                "Family",
                "Nearest_Gene_ID",
                "Distance_bp",
                "Distance_Bin",
                "Relation",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(grouped_rows[group])

        with open(txt_path, "w", newline="") as handle:
            for gene_id in grouped_gene_ids[group]:
                handle.write(f"{gene_id}\n")

        print(f"{group}: {len(grouped_rows[group])} TE-gene pairs, {len(grouped_gene_ids[group])} unique genes")


if __name__ == "__main__":
    main()
