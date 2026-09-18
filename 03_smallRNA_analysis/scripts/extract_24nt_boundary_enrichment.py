#!/usr/bin/env python3
"""Quantify 24-nt siRNA boundary enrichment from deepTools matrices.

For each accession/family/TE-length/gene-proximity group, the script averages
root and young-panicle matrices and reports a boundary enrichment index:

    log2((mean of the outer 20% at each TE-body edge + 1e-4) /
         (mean of the inner 60% of the TE body + 1e-4))

The metric is calculated from the original TE rows, not from plotted profiles.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path


SCRIPT = Path(__file__).resolve()
FIGURE_ROOT = SCRIPT.parents[3]
DATA_ROOT = FIGURE_ROOT.parent / "data" / "TE-CHH-siRNA"
MATRIX_DIR = DATA_ROOT / "matrices"
SIGNAL_TABLE = (
    FIGURE_ROOT
    / "Figure4"
    / "noncanonical_RdDM_21_22nt"
    / "data"
    / "15_TE_signals_tissue_averaged.tsv"
)
OUT_DIR = SCRIPT.parent.parent / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE = OUT_DIR / "24nt_boundary_enrichment_by_accession.tsv"

PSEUDOCOUNT = 1e-4


def length_class(length: float) -> str:
    if length <= 500:
        return "Short (<=500 bp)"
    if length <= 1500:
        return "Medium (501-1500 bp)"
    return "Long (>1500 bp)"


def load_metadata() -> dict[tuple[str, str, str, int, int], tuple[float, str]]:
    metadata: dict[tuple[str, str, str, int, int], tuple[float, str]] = {}
    with SIGNAL_TABLE.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (
                row["accession"],
                row["family"],
                row["chr"],
                int(row["start"]),
                int(row["end"]),
            )
            metadata[key] = (float(row["TE_length"]), row["gene_proximity"])
    return metadata


def family_for_row(row_index: int, header: dict) -> str:
    boundaries = header["group_boundaries"]
    labels = header["group_labels"]
    for idx, label in enumerate(labels):
        if boundaries[idx] <= row_index < boundaries[idx + 1]:
            if "Copia" in label:
                return "Copia"
            if "Gypsy" in label:
                return "Gypsy"
    raise ValueError(f"Unable to assign family for row {row_index}")


def process_matrix(
    matrix_path: Path,
    metadata: dict[tuple[str, str, str, int, int], tuple[float, str]],
    accum: dict[tuple[str, str, str, str], list[float]],
) -> tuple[int, int]:
    accession_match = re.match(r"(.+?)_T\d+_24nt_", matrix_path.name)
    if not accession_match:
        raise ValueError(f"Unexpected matrix name: {matrix_path.name}")
    accession = accession_match.group(1)
    matched = 0
    missing = 0

    with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
        header_line = next(handle)
        if not header_line.startswith("@"):
            raise ValueError(f"Missing deepTools JSON header: {matrix_path}")
        header = json.loads(header_line[1:])
        upstream_bins = int(header["upstream"][0] / header["bin size"][0])
        body_bins = int(header["body"][0] / header["bin size"][0])
        sample_boundaries = header["sample_boundaries"]
        edge_bins = max(1, int(round(body_bins * 0.20)))

        for row_index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            family = family_for_row(row_index, header)
            key = (accession, family, chrom, start, end)
            meta = metadata.get(key)
            if meta is None:
                missing += 1
                continue
            te_length, proximity = meta
            values = [float(x) if x not in {"nan", "NA", ""} else math.nan for x in fields[6:]]

            replicate_edges: list[float] = []
            replicate_inners: list[float] = []
            for sample_start, sample_end in zip(sample_boundaries[:-1], sample_boundaries[1:]):
                block = values[sample_start:sample_end]
                body = block[upstream_bins : upstream_bins + body_bins]
                edges = body[:edge_bins] + body[-edge_bins:]
                inner = body[edge_bins:-edge_bins]
                edge_valid = [x for x in edges if math.isfinite(x)]
                inner_valid = [x for x in inner if math.isfinite(x)]
                if edge_valid and inner_valid:
                    replicate_edges.append(sum(edge_valid) / len(edge_valid))
                    replicate_inners.append(sum(inner_valid) / len(inner_valid))

            if not replicate_edges or not replicate_inners:
                continue
            edge_mean = sum(replicate_edges) / len(replicate_edges)
            inner_mean = sum(replicate_inners) / len(replicate_inners)
            group = (accession, family, length_class(te_length), proximity)
            accum[group][0] += edge_mean
            accum[group][1] += inner_mean
            accum[group][2] += 1
            matched += 1

    return matched, missing


def main() -> None:
    metadata = load_metadata()
    accum: dict[tuple[str, str, str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0])
    matrix_files = sorted(MATRIX_DIR.glob("*_24nt_[RY]_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz"))
    if not matrix_files:
        raise FileNotFoundError(f"No 24-nt matrices found in {MATRIX_DIR}")

    matched_total = 0
    missing_total = 0
    for matrix_path in matrix_files:
        matched, missing = process_matrix(matrix_path, metadata, accum)
        matched_total += matched
        missing_total += missing

    with OUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "accession",
                "family",
                "length_class",
                "gene_proximity",
                "n_TE_tissue_rows",
                "mean_boundary_signal",
                "mean_inner_body_signal",
                "boundary_enrichment_log2",
            ]
        )
        for group in sorted(accum):
            edge_sum, inner_sum, count = accum[group]
            edge_mean = edge_sum / count
            inner_mean = inner_sum / count
            enrichment = math.log2((edge_mean + PSEUDOCOUNT) / (inner_mean + PSEUDOCOUNT))
            writer.writerow(
                [
                    *group,
                    int(count),
                    f"{edge_mean:.10g}",
                    f"{inner_mean:.10g}",
                    f"{enrichment:.10g}",
                ]
            )

    print(f"Wrote {OUT_FILE}")
    print(f"Matrices: {len(matrix_files)}; matched rows: {matched_total}; unmatched rows: {missing_total}")


if __name__ == "__main__":
    main()
