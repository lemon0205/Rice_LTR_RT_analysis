#!/usr/bin/env python3
"""Build distance-to-gene summaries for four TE fate groups."""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path


GROUP_ORDER = ["all_shared_te", "msu_specific_te", "fst_1_te", "other_te"]
DIST_ORDER = ["Overlap", "Near_1kb", "Far_1kb"]


def load_gene_index(path: Path):
    by_chr = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end < start:
                start, end = end, start
            by_chr[chrom].append((start, end))
    indexed = {}
    for chrom, vals in by_chr.items():
        vals.sort()
        indexed[chrom] = {"intervals": vals, "starts": [v[0] for v in vals]}
    return indexed


def nearest_distance(chrom: str, start: int, end: int, genes) -> int | None:
    if chrom not in genes:
        return None
    intervals = genes[chrom]["intervals"]
    starts = genes[chrom]["starts"]
    idx = bisect.bisect_right(starts, end)
    best = None
    for j in range(max(0, idx - 3), min(len(intervals), idx + 4)):
        gs, ge = intervals[j]
        if end >= gs and start <= ge:
            return 0
        dist = gs - end if end < gs else start - ge
        if dist >= 0 and (best is None or dist < best):
            best = dist
    return best


def dist_class(distance: int | None) -> str:
    if distance is None:
        return "No_gene"
    if distance == 0:
        return "Overlap"
    if distance <= 1000:
        return "Near_1kb"
    return "Far_1kb"


def write_summary(rows: list[dict[str, object]], outdir: Path):
    count = Counter()
    totals = Counter()
    count_family = Counter()
    totals_family = Counter()
    for row in rows:
        g = row["Group"]
        d = row["Distance_class"]
        f = row["Family"]
        count[(g, d)] += 1
        totals[g] += 1
        count_family[(g, f, d)] += 1
        totals_family[(g, f)] += 1

    with (outdir / "fate_group_distance_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["Group", "Distance_class", "count", "proportion"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for g in GROUP_ORDER:
            for d in DIST_ORDER:
                c = count[(g, d)]
                writer.writerow({"Group": g, "Distance_class": d, "count": c, "proportion": c / totals[g] if totals[g] else 0})

    with (outdir / "fate_group_distance_by_family_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["Group", "Family", "Distance_class", "count", "proportion"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for g in GROUP_ORDER:
            for f in ["Copia", "Gypsy"]:
                for d in DIST_ORDER:
                    c = count_family[(g, f, d)]
                    writer.writerow({"Group": g, "Family": f, "Distance_class": d, "count": c, "proportion": c / totals_family[(g, f)] if totals_family[(g, f)] else 0})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--te-four-groups", required=True, type=Path)
    parser.add_argument("--genes-bed", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = load_gene_index(args.genes_bed)

    rows = []
    with args.te_four_groups.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            group = row["Group"]
            family = row["Family"]
            if group not in GROUP_ORDER or family not in ("Copia", "Gypsy"):
                continue
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
            if end < start:
                start, end = end, start
            distance = nearest_distance(chrom, start, end, genes)
            out = dict(row)
            out["Distance_to_gene_bp"] = "" if distance is None else distance
            out["Distance_class"] = dist_class(distance)
            rows.append(out)

    fields = list(rows[0].keys())
    with (args.outdir / "te_fate_group_distance_master.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    write_summary(rows, args.outdir)
    print("rows", len(rows))
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()

