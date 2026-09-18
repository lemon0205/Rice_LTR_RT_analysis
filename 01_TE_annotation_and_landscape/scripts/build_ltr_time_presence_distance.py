#!/usr/bin/env python3
"""
Build MSU-anchored TE presence, gene-distance, and identity summaries.

Inputs:
  - te_genotype_matrix_strict_noNA.csv
  - MSU_genes.bed
  - MSU_TE.gff3

Outputs:
  - te_presence_distance_identity.tsv
  - summary_by_family_position.tsv
  - frequency_class_by_family_position.tsv
  - identity_presence_long.tsv
  - README_LTR_time_analysis.md
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median


SAMPLES = [
    "MSU",
    "Basmati1",
    "CG14",
    "G46",
    "IR64",
    "Lemont",
    "LJ",
    "N22",
    "NamRoo",
    "TM",
    "Tumba",
    "WSSM",
]


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs = {}
    for item in attr_text.strip().split(";"):
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        attrs[k] = v
    return attrs


def safe_float(value: str | None) -> float | None:
    if value is None or value == "" or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_te_gff(gff_path: Path) -> dict[str, dict[str, object]]:
    te = {}
    with gff_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_text = parts
            attrs = parse_attrs(attr_text)
            te_id = attrs.get("ID")
            if not te_id:
                continue
            te[te_id] = {
                "gff_chrom": chrom,
                "gff_start": int(start),
                "gff_end": int(end),
                "feature": feature,
                "classification": attrs.get("classification", ""),
                "identity": safe_float(attrs.get("identity")),
                "ltr_identity": safe_float(attrs.get("ltr_identity")),
                "method": attrs.get("method", ""),
                "name": attrs.get("Name", ""),
            }
    return te


def load_genes_bed(bed_path: Path) -> dict[str, dict[str, list[int]]]:
    intervals = defaultdict(list)
    with bed_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if end < start:
                start, end = end, start
            intervals[chrom].append((start, end))

    indexed = {}
    for chrom, vals in intervals.items():
        vals.sort()
        starts = [x[0] for x in vals]
        ends = [x[1] for x in vals]
        indexed[chrom] = {"intervals": vals, "starts": starts, "ends": ends}
    return indexed


def nearest_gene_distance(chrom: str, start: int, end: int, genes) -> int | None:
    if chrom not in genes:
        return None
    data = genes[chrom]
    intervals = data["intervals"]
    starts = data["starts"]
    idx = bisect.bisect_right(starts, end)

    best = None
    for j in (idx - 2, idx - 1, idx, idx + 1):
        if j < 0 or j >= len(intervals):
            continue
        gs, ge = intervals[j]
        if end >= gs and start <= ge:
            return 0
        if end < gs:
            dist = gs - end
        else:
            dist = start - ge
        if dist >= 0 and (best is None or dist < best):
            best = dist
    return best


def position_class(distance: int | None) -> str:
    if distance is None:
        return "no_gene_on_chrom"
    if distance == 0:
        return "overlap"
    if distance <= 1000:
        return "near_1kb"
    return "far_1kb"


def frequency_class(presence_count: int) -> str:
    if presence_count == 1:
        return "private"
    if presence_count <= 3:
        return "low_frequency"
    if presence_count <= 8:
        return "intermediate"
    if presence_count <= 11:
        return "shared"
    return "core"


def num_summary(values: list[float]) -> dict[str, object]:
    vals = [v for v in values if v is not None and not math.isnan(v)]
    if not vals:
        return {
            "n": 0,
            "mean": "",
            "median": "",
            "min": "",
            "max": "",
        }
    return {
        "n": len(vals),
        "mean": f"{mean(vals):.6g}",
        "median": f"{median(vals):.6g}",
        "min": f"{min(vals):.6g}",
        "max": f"{max(vals):.6g}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--genes-bed", required=True, type=Path)
    parser.add_argument("--te-gff", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    te_gff = load_te_gff(args.te_gff)
    genes = load_genes_bed(args.genes_bed)

    detailed_path = args.outdir / "te_presence_distance_identity.tsv"
    summary_path = args.outdir / "summary_by_family_position.tsv"
    freq_path = args.outdir / "frequency_class_by_family_position.tsv"
    long_path = args.outdir / "identity_presence_long.tsv"
    readme_path = args.outdir / "README_LTR_time_analysis.md"

    detailed_rows = []
    with args.matrix.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            te_id = row["TE_ID"]
            chrom = row["Chromosome"]
            start = int(row["Start"])
            end = int(row["End"])
            if end < start:
                start, end = end, start
            presence_count = sum(int(row[s]) for s in SAMPLES)
            presence_freq = presence_count / len(SAMPLES)
            dist = nearest_gene_distance(chrom, start, end, genes)
            gff_info = te_gff.get(te_id, {})
            identity = gff_info.get("identity")
            ltr_identity = gff_info.get("ltr_identity")
            detailed_rows.append(
                {
                    "TE_ID": te_id,
                    "Family": row["Family"],
                    "Chromosome": chrom,
                    "Start": start,
                    "End": end,
                    "Length": row["Length"],
                    "presence_count": presence_count,
                    "presence_frequency": f"{presence_freq:.6g}",
                    "frequency_class": frequency_class(presence_count),
                    "distance_to_nearest_gene_bp": "" if dist is None else dist,
                    "position_class": position_class(dist),
                    "identity": "" if identity is None else f"{identity:.6g}",
                    "ltr_identity": "" if ltr_identity is None else f"{ltr_identity:.6g}",
                    "classification": gff_info.get("classification", ""),
                    "method": gff_info.get("method", ""),
                }
            )

    fieldnames = list(detailed_rows[0].keys()) if detailed_rows else []
    with detailed_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(detailed_rows)

    grouped = defaultdict(list)
    freq_counts = defaultdict(Counter)
    for row in detailed_rows:
        key = (row["Family"], row["position_class"])
        grouped[key].append(row)
        freq_counts[key][row["frequency_class"]] += 1

    summary_fields = [
        "Family",
        "position_class",
        "n",
        "presence_count_mean",
        "presence_count_median",
        "presence_frequency_mean",
        "identity_n",
        "identity_mean",
        "identity_median",
        "ltr_identity_n",
        "ltr_identity_mean",
        "ltr_identity_median",
        "distance_median_bp",
    ]
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for (fam, pos), rows in sorted(grouped.items()):
            pc = [float(r["presence_count"]) for r in rows]
            pf = [float(r["presence_frequency"]) for r in rows]
            identity_vals = [float(r["identity"]) for r in rows if r["identity"] != ""]
            ltr_vals = [float(r["ltr_identity"]) for r in rows if r["ltr_identity"] != ""]
            dist_vals = [float(r["distance_to_nearest_gene_bp"]) for r in rows if r["distance_to_nearest_gene_bp"] != ""]
            writer.writerow(
                {
                    "Family": fam,
                    "position_class": pos,
                    "n": len(rows),
                    "presence_count_mean": f"{mean(pc):.6g}",
                    "presence_count_median": f"{median(pc):.6g}",
                    "presence_frequency_mean": f"{mean(pf):.6g}",
                    "identity_n": len(identity_vals),
                    "identity_mean": "" if not identity_vals else f"{mean(identity_vals):.6g}",
                    "identity_median": "" if not identity_vals else f"{median(identity_vals):.6g}",
                    "ltr_identity_n": len(ltr_vals),
                    "ltr_identity_mean": "" if not ltr_vals else f"{mean(ltr_vals):.6g}",
                    "ltr_identity_median": "" if not ltr_vals else f"{median(ltr_vals):.6g}",
                    "distance_median_bp": "" if not dist_vals else f"{median(dist_vals):.6g}",
                }
            )

    freq_fields = ["Family", "position_class", "frequency_class", "count", "proportion"]
    with freq_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=freq_fields, delimiter="\t")
        writer.writeheader()
        for key, counts in sorted(freq_counts.items()):
            fam, pos = key
            total = sum(counts.values())
            for cls in ["private", "low_frequency", "intermediate", "shared", "core"]:
                writer.writerow(
                    {
                        "Family": fam,
                        "position_class": pos,
                        "frequency_class": cls,
                        "count": counts.get(cls, 0),
                        "proportion": f"{counts.get(cls, 0) / total:.6g}" if total else "",
                    }
                )

    with long_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["TE_ID", "Family", "position_class", "metric", "value"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in detailed_rows:
            for metric in ["presence_frequency", "identity", "ltr_identity"]:
                if row[metric] != "":
                    writer.writerow(
                        {
                            "TE_ID": row["TE_ID"],
                            "Family": row["Family"],
                            "position_class": row["position_class"],
                            "metric": metric,
                            "value": row[metric],
                        }
                    )

    readme_path.write_text(
        "\n".join(
            [
                "# LTR time / retention analysis",
                "",
                "This folder summarizes MSU-anchored TE presence frequency, distance to the nearest MSU gene, and TE identity values.",
                "",
                "## Inputs",
                f"- Genotype matrix: `{args.matrix}`",
                f"- MSU genes BED: `{args.genes_bed}`",
                f"- MSU TE GFF3: `{args.te_gff}`",
                "",
                "## Position classes",
                "- `overlap`: TE overlaps an MSU gene interval.",
                "- `near_1kb`: TE does not overlap a gene but is within 1 kb of the nearest gene.",
                "- `far_1kb`: TE is more than 1 kb away from the nearest gene.",
                "",
                "## Frequency classes",
                "- `private`: present in 1 of 12 accessions.",
                "- `low_frequency`: present in 2-3 accessions.",
                "- `intermediate`: present in 4-8 accessions.",
                "- `shared`: present in 9-11 accessions.",
                "- `core`: present in all 12 accessions.",
                "",
                "## Interpretation notes",
                "- `identity` is the EDTA library-to-target sequence identity and can be used as a relative age/conservation proxy.",
                "- `ltr_identity` is preferable for intact LTR insertion age when available, but it may be sparse.",
                "- The key test is whether gene-proximal TEs have lower presence frequency and/or higher identity, supporting lower long-term retention near genes.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"wrote {detailed_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {freq_path}")
    print(f"wrote {long_path}")
    print(f"wrote {readme_path}")


if __name__ == "__main__":
    main()

