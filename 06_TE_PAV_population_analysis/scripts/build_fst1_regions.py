from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter
from pathlib import Path


MERGE_DISTANCE = 2000
INPUT_CSV = Path(os.environ["DTE_FST_RESULTS"])
OUT_DIR = Path(os.environ["FST1_REGION_OUTDIR"])


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def classify_region(te_count: int) -> str:
    if te_count == 1:
        return "single_te"
    if te_count == 2:
        return "multi_te_2"
    return "hotspot_ge3"


def finalize_region(region_index: int, members: list[dict[str, str]]) -> tuple[dict[str, object], list[dict[str, object]]]:
    starts = [int(row["Start"]) for row in members]
    ends = [int(row["End"]) for row in members]
    families = [row["Family"] for row in members]
    family_counts = Counter(families)
    dominant_family = sorted(family_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]
    biases = {row["group_bias"] for row in members}
    group_bias = sorted(biases)[0] if len(biases) == 1 else "mixed"
    fst_values = [safe_float(row["fst"]) for row in members]
    delta_values = [safe_float(row["delta_freq"]) for row in members]
    region_id = f"FST1R_{region_index:06d}"

    region_row = {
        "Region_ID": region_id,
        "Chromosome": members[0]["Chromosome"],
        "Start": min(starts),
        "End": max(ends),
        "Region_Length": max(ends) - min(starts),
        "TE_Count": len(members),
        "Region_Class": classify_region(len(members)),
        "TE_IDs": ";".join(row["TE_ID"] for row in members),
        "Families": ";".join(sorted(f"{fam}:{cnt}" for fam, cnt in family_counts.items())),
        "Dominant_Family": dominant_family,
        "Group_Bias": group_bias,
        "Max_FST": max(fst_values),
        "Mean_FST": sum(fst_values) / len(fst_values),
        "Mean_Delta_Freq": sum(delta_values) / len(delta_values),
    }

    member_rows: list[dict[str, object]] = []
    for order, row in enumerate(members, start=1):
        member_rows.append(
            {
                "Region_ID": region_id,
                "Member_Order": order,
                "TE_ID": row["TE_ID"],
                "Family": row["Family"],
                "Chromosome": row["Chromosome"],
                "Start": int(row["Start"]),
                "End": int(row["End"]),
                "Length": int(row["Length"]),
                "FST": safe_float(row["fst"]),
                "Delta_Freq": safe_float(row["delta_freq"]),
                "Group_Bias": row["group_bias"],
            }
        )
    return region_row, member_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fst1_rows = [row for row in reader if safe_float(row["fst"]) == 1.0]

    fst1_rows.sort(key=lambda row: (row["Chromosome"], int(row["Start"]), int(row["End"]), row["TE_ID"]))

    regions: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    current: list[dict[str, str]] = []
    previous_chr = None
    previous_end = None
    region_index = 1

    for row in fst1_rows:
        chrom = row["Chromosome"]
        start = int(row["Start"])
        end = int(row["End"])

        if not current:
            current = [row]
        elif chrom == previous_chr and start - int(previous_end) <= MERGE_DISTANCE:
            current.append(row)
        else:
            region_row, member_rows = finalize_region(region_index, current)
            regions.append(region_row)
            members.extend(member_rows)
            region_index += 1
            current = [row]

        previous_chr = chrom
        previous_end = end

    if current:
        region_row, member_rows = finalize_region(region_index, current)
        regions.append(region_row)
        members.extend(member_rows)

    regions_csv = OUT_DIR / "fst1_high_fst_regions.csv"
    with regions_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Region_ID",
                "Chromosome",
                "Start",
                "End",
                "Region_Length",
                "TE_Count",
                "Region_Class",
                "TE_IDs",
                "Families",
                "Dominant_Family",
                "Group_Bias",
                "Max_FST",
                "Mean_FST",
                "Mean_Delta_Freq",
            ],
        )
        writer.writeheader()
        writer.writerows(regions)

    members_csv = OUT_DIR / "fst1_high_fst_region_te_members.csv"
    with members_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Region_ID",
                "Member_Order",
                "TE_ID",
                "Family",
                "Chromosome",
                "Start",
                "End",
                "Length",
                "FST",
                "Delta_Freq",
                "Group_Bias",
            ],
        )
        writer.writeheader()
        writer.writerows(members)

    bed_path = OUT_DIR / "fst1_high_fst_regions.bed"
    with bed_path.open("w", encoding="utf-8", newline="") as handle:
        for row in regions:
            handle.write(
                f"{row['Chromosome']}\t{row['Start']}\t{row['End']}\t{row['Region_ID']}|{row['Region_Class']}|{row['TE_Count']}\n"
            )

    class_counts = Counter(row["Region_Class"] for row in regions)
    hotspot_regions = [row for row in regions if row["Region_Class"] == "hotspot_ge3"]
    hotspot_regions_sorted = sorted(
        hotspot_regions,
        key=lambda row: (-int(row["TE_Count"]), -int(row["Region_Length"]), row["Region_ID"]),
    )
    region_lengths = sorted(int(row["Region_Length"]) for row in regions)
    te_counts = sorted(int(row["TE_Count"]) for row in regions)

    def median(values: list[int]) -> float:
        if not values:
            return 0.0
        n = len(values)
        mid = n // 2
        if n % 2 == 1:
            return float(values[mid])
        return (values[mid - 1] + values[mid]) / 2.0

    summary_lines = [
        f"Merge distance: {MERGE_DISTANCE} bp",
        f"fst=1 TE count: {len(fst1_rows)}",
        f"Total fst1 regions: {len(regions)}",
        f"single_te regions: {class_counts.get('single_te', 0)}",
        f"TE_Count = 2 regions: {class_counts.get('multi_te_2', 0)}",
        f"TE_Count >= 3 hotspots: {class_counts.get('hotspot_ge3', 0)}",
        f"hotspot_ge3 proportion among fst=1 TEs: {sum(int(row['TE_Count']) for row in hotspot_regions) / len(fst1_rows):.4f}",
        f"Median region length: {median(region_lengths):.1f}",
        f"Median TE per region: {median(te_counts):.1f}",
        "",
        "Top 10 hotspot_ge3 regions:",
    ]

    for row in hotspot_regions_sorted[:10]:
        summary_lines.append(
            f"{row['Region_ID']}\t{row['Chromosome']}:{row['Start']}-{row['End']}\tTE_Count={row['TE_Count']}\tLength={row['Region_Length']}\tFamily={row['Dominant_Family']}\tBias={row['Group_Bias']}"
        )

    summary_path = OUT_DIR / "fst1_region_summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"fst=1 TE count: {len(fst1_rows)}")
    print(f"Total fst1 regions: {len(regions)}")
    print(f"single_te regions: {class_counts.get('single_te', 0)}")
    print(f"TE_Count = 2 regions: {class_counts.get('multi_te_2', 0)}")
    print(f"TE_Count >= 3 hotspots: {class_counts.get('hotspot_ge3', 0)}")


if __name__ == "__main__":
    main()
