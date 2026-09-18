#!/usr/bin/env python3
"""Build TE retention table from MSU-anchored genotype matrix and gene distance."""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path

SAMPLES = ["MSU", "Basmati1", "CG14", "G46", "IR64", "Lemont", "LJ", "N22", "NamRoo", "TM", "Tumba", "WSSM"]


def parse_attrs(text: str) -> dict[str, str]:
    attrs = {}
    for part in text.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def load_gene_index(path: Path):
    by_chr = defaultdict(list)
    with path.open(encoding="utf-8") as f:
        for line in f:
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
        idx[chrom] = {"intervals": vals, "starts": [x[0] for x in vals]}
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


def distance_class(d: int | None) -> str:
    if d is None:
        return "No_gene"
    if d == 0:
        return "Overlap"
    if d <= 1000:
        return "Near_1kb"
    return "Far_1kb"


def freq_class(c: int) -> str:
    if c == 1:
        return "private"
    if c <= 3:
        return "low_frequency"
    if c <= 8:
        return "intermediate"
    if c <= 11:
        return "shared"
    return "core"


def load_young_intervals(copia_young: Path, gypsy_young: Path):
    by_key = defaultdict(list)
    for family, path in [("Copia", copia_young), ("Gypsy", gypsy_young)]:
        with path.open(encoding="utf-8") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "repeat_region":
                    continue
                chrom, start, end = p[0], int(p[3]), int(p[4])
                if end < start:
                    start, end = end, start
                by_key[(family, chrom)].append((start, end))
    idx = {}
    for key, vals in by_key.items():
        vals.sort()
        idx[key] = {"intervals": vals, "starts": [x[0] for x in vals]}
    return idx


def overlaps_young(family: str, chrom: str, start: int, end: int, young_idx) -> bool:
    key = (family, chrom)
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


def summarize(rows, outdir: Path):
    # Presence summary
    groups = defaultdict(list)
    for r in rows:
        for panel in ["All", r["Age_group"]]:
            if panel == "NonYoung":
                continue
            groups[(panel, r["Family"], r["Distance_class"])].append(float(r["presence_frequency"]))
    with (outdir / "presence_frequency_summary.tsv").open("w", newline="", encoding="utf-8") as f:
        fields = ["Panel", "Family", "Distance_class", "n", "mean", "median"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for key, vals in sorted(groups.items()):
            vals2 = sorted(vals)
            med = vals2[len(vals2)//2] if len(vals2) % 2 else (vals2[len(vals2)//2-1] + vals2[len(vals2)//2]) / 2
            w.writerow({"Panel": key[0], "Family": key[1], "Distance_class": key[2], "n": len(vals), "mean": sum(vals)/len(vals), "median": med})

    # Frequency class summary
    counts = defaultdict(Counter)
    totals = Counter()
    for r in rows:
        for panel in ["All", r["Age_group"]]:
            if panel == "NonYoung":
                continue
            key = (panel, r["Family"], r["Distance_class"])
            counts[key][r["frequency_class"]] += 1
            totals[key] += 1
    with (outdir / "frequency_class_summary.tsv").open("w", newline="", encoding="utf-8") as f:
        fields = ["Panel", "Family", "Distance_class", "frequency_class", "count", "proportion"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for key in sorted(totals):
            for cls in ["private", "low_frequency", "intermediate", "shared", "core"]:
                c = counts[key][cls]
                w.writerow({"Panel": key[0], "Family": key[1], "Distance_class": key[2], "frequency_class": cls, "count": c, "proportion": c / totals[key]})


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
    young = load_young_intervals(args.copia_young, args.gypsy_young)

    rows = []
    with args.matrix.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            fam = row["Family"]
            if fam not in ("Copia", "Gypsy"):
                continue
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
            if end < start:
                start, end = end, start
            pc = sum(int(row[s]) for s in SAMPLES)
            dist = nearest_distance(chrom, start, end, genes)
            age_group = "Young" if overlaps_young(fam, chrom, start, end, young) else "NonYoung"
            rows.append({
                "TE_ID": row["TE_ID"],
                "Family": fam,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Length": row["Length"],
                "presence_count": pc,
                "presence_frequency": pc / len(SAMPLES),
                "frequency_class": freq_class(pc),
                "distance_to_gene_bp": "" if dist is None else dist,
                "Distance_class": distance_class(dist),
                "Age_group": age_group,
            })

    fields = list(rows[0].keys())
    with (args.outdir / "te_retention_master_table.tsv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summarize(rows, args.outdir)

    print("rows", len(rows))
    print("young", sum(1 for r in rows if r["Age_group"] == "Young"))
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()

