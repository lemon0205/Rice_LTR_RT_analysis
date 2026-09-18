#!/usr/bin/env python3
"""Build TE distance-to-gene versus presence-frequency table."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import defaultdict
from pathlib import Path

SAMPLES = ["MSU", "Basmati1", "CG14", "G46", "IR64", "Lemont", "LJ", "N22", "NamRoo", "TM", "Tumba", "WSSM"]


def load_gene_index(path: Path):
    by_chr = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            if end < start:
                start, end = end, start
            by_chr[chrom].append((start, end))
    idx = {}
    for chrom, vals in by_chr.items():
        vals.sort()
        idx[chrom] = {"intervals": vals, "starts": [v[0] for v in vals]}
    return idx


def nearest_distance(chrom: str, start: int, end: int, genes) -> int | None:
    if chrom not in genes:
        return None
    intervals = genes[chrom]["intervals"]
    starts = genes[chrom]["starts"]
    i = bisect.bisect_right(starts, end)
    best = None
    for j in range(max(0, i - 3), min(len(intervals), i + 4)):
        gs, ge = intervals[j]
        if end >= gs and start <= ge:
            return 0
        d = gs - end if end < gs else start - ge
        if d >= 0 and (best is None or d < best):
            best = d
    return best


def distance_bin(distance: int | None) -> str:
    if distance is None:
        return "No_gene"
    if distance == 0:
        return "Overlap"
    if distance <= 100:
        return "1-100bp"
    if distance <= 250:
        return "100-250bp"
    if distance <= 500:
        return "250-500bp"
    if distance <= 1000:
        return "500bp-1kb"
    if distance <= 2000:
        return "1-2kb"
    if distance <= 5000:
        return "2-5kb"
    if distance <= 10000:
        return "5-10kb"
    if distance <= 20000:
        return "10-20kb"
    if distance <= 50000:
        return "20-50kb"
    return ">50kb"


def parse_young_intervals(copia_young: Path, gypsy_young: Path):
    by_key = defaultdict(list)
    for fam, path in [("Copia", copia_young), ("Gypsy", gypsy_young)]:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "repeat_region":
                    continue
                chrom, start, end = p[0], int(p[3]), int(p[4])
                if end < start:
                    start, end = end, start
                by_key[(fam, chrom)].append((start, end))
    idx = {}
    for key, vals in by_key.items():
        vals.sort()
        idx[key] = {"intervals": vals, "starts": [v[0] for v in vals]}
    return idx


def is_young(fam: str, chrom: str, start: int, end: int, young_idx) -> bool:
    key = (fam, chrom)
    if key not in young_idx:
        return False
    intervals = young_idx[key]["intervals"]
    starts = young_idx[key]["starts"]
    i = bisect.bisect_right(starts, end)
    for j in range(max(0, i - 5), min(len(intervals), i + 5)):
        ys, ye = intervals[j]
        if end >= ys and start <= ye:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--genes-bed", required=True, type=Path)
    ap.add_argument("--copia-young", required=True, type=Path)
    ap.add_argument("--gypsy-young", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    genes = load_gene_index(args.genes_bed)
    young_idx = parse_young_intervals(args.copia_young, args.gypsy_young)
    rows = []
    with args.matrix.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            fam = row["Family"]
            if fam not in ("Copia", "Gypsy"):
                continue
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
            if end < start:
                start, end = end, start
            pc = sum(int(row[s]) for s in SAMPLES)
            dist = nearest_distance(chrom, start, end, genes)
            y = is_young(fam, chrom, start, end, young_idx)
            rows.append({
                "TE_ID": row["TE_ID"],
                "Family": fam,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Length": row["Length"],
                "presence_count": pc,
                "presence_frequency": pc / len(SAMPLES),
                "distance_to_gene_bp": "" if dist is None else dist,
                "log10_distance_plus1": "" if dist is None else math.log10(dist + 1),
                "distance_bin": distance_bin(dist),
                "Age_group": "Young" if y else "All",
                "low_frequency": 1 if pc <= 3 else 0,
            })

    fields = list(rows[0].keys())
    with (args.outdir / "te_distance_frequency_master.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print("rows", len(rows))
    print("young", sum(1 for r in rows if r["Age_group"] == "Young"))


if __name__ == "__main__":
    main()

