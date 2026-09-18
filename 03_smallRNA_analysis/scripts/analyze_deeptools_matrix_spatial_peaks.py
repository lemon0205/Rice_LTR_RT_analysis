#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


MATRIX_DIR = Path(os.environ["DEEPTOOLS_MATRIX_DIR"])
ANNOTATION_FILE = Path(os.environ["TE_SIGNAL_TABLE"])
OUT_DIR = Path(os.environ["TE_SIGNAL_OUTDIR"])

FILE_RE = re.compile(
    r"^(?P<accession>.+?)_(?P<tcode>T\d+)_(?P<signal>21-22nt|24nt)_"
    r"(?P<tissue>[RY])_TE_scale2kb_up2kb_down2kb_bin50\.matrix\.gz$"
)


def read_matrix(path: Path):
    with gzip.open(path, "rt") as handle:
        header_line = handle.readline().rstrip("\n")
        if not header_line.startswith("@"):
            raise ValueError(f"Missing deepTools JSON header: {path}")
        header = json.loads(header_line[1:])
        sample_boundaries = header["sample_boundaries"]
        group_boundaries = header["group_boundaries"]
        group_labels = header["group_labels"]
        n_bins = sample_boundaries[1] - sample_boundaries[0]
        n_samples = len(sample_boundaries) - 1

        metadata = []
        values = []
        row_index = 0
        group_index = 0
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            while (
                group_index + 1 < len(group_boundaries) - 1
                and row_index >= group_boundaries[group_index + 1]
            ):
                group_index += 1
            label = group_labels[group_index]
            family = "Copia" if "Copia" in label else "Gypsy"
            numeric = np.asarray(fields[6:], dtype=np.float32)
            if numeric.size != n_samples * n_bins:
                raise ValueError(
                    f"Unexpected value count in {path.name}, row {row_index}: "
                    f"{numeric.size} != {n_samples * n_bins}"
                )
            numeric = numeric.reshape(n_samples, n_bins)
            row_profile = np.nanmean(numeric, axis=0).astype(np.float32)
            metadata.append(
                (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    fields[3],
                    fields[5],
                    family,
                )
            )
            values.append(row_profile)
            row_index += 1

    meta = pd.DataFrame(
        metadata, columns=["chr", "start", "end", "TE_id", "strand", "family"]
    )
    matrix = np.vstack(values).astype(np.float32)
    return header, meta, matrix


def assert_same_rows(reference: pd.DataFrame, other: pd.DataFrame, label: str):
    columns = ["chr", "start", "end", "TE_id", "strand", "family"]
    if len(reference) != len(other) or not reference[columns].equals(other[columns]):
        raise ValueError(f"TE row order mismatch for {label}")


def safe_mean(matrix: np.ndarray, start: int, stop: int):
    return np.nanmean(matrix[:, start:stop], axis=1)


def matrix_metrics(matrix: np.ndarray, prefix: str):
    # 0:40 = 2-kb upstream; 40:80 = scaled TE body; 80:120 = 2-kb downstream.
    upstream = safe_mean(matrix, 0, 40)
    body = safe_mean(matrix, 40, 80)
    downstream = safe_mean(matrix, 80, 120)
    five_boundary = safe_mean(matrix, 35, 45)
    body_interior = safe_mean(matrix, 45, 75)
    three_boundary = safe_mean(matrix, 75, 85)
    boundary = np.nanmean(
        np.column_stack([five_boundary, three_boundary]), axis=1
    )
    flank = np.nanmean(np.column_stack([upstream, downstream]), axis=1)
    window = np.nanmean(matrix, axis=1)
    maximum = np.nanmax(matrix, axis=1)
    max_bin = np.nanargmax(matrix, axis=1)
    detected_fraction = np.mean(matrix > 0, axis=1)
    total = np.nansum(matrix, axis=1)
    sorted_matrix = np.sort(np.nan_to_num(matrix, nan=0.0), axis=1)
    top10_share = np.divide(
        np.sum(sorted_matrix[:, -10:], axis=1),
        total,
        out=np.zeros_like(total),
        where=total > 0,
    )
    eps = np.float32(1e-6)
    boundary_vs_interior = np.log2((boundary + eps) / (body_interior + eps))
    body_vs_flank = np.log2((body + eps) / (flank + eps))

    peak_region = np.full(len(max_bin), "body_interior", dtype=object)
    peak_region[max_bin < 35] = "upstream_flank"
    peak_region[(max_bin >= 35) & (max_bin < 45)] = "5prime_boundary"
    peak_region[(max_bin >= 45) & (max_bin < 75)] = "body_interior"
    peak_region[(max_bin >= 75) & (max_bin < 85)] = "3prime_boundary"
    peak_region[max_bin >= 85] = "downstream_flank"

    return pd.DataFrame(
        {
            f"{prefix}_window_mean": window,
            f"{prefix}_upstream_mean": upstream,
            f"{prefix}_5prime_boundary_mean": five_boundary,
            f"{prefix}_body_mean": body,
            f"{prefix}_body_interior_mean": body_interior,
            f"{prefix}_3prime_boundary_mean": three_boundary,
            f"{prefix}_downstream_mean": downstream,
            f"{prefix}_boundary_mean": boundary,
            f"{prefix}_flank_mean": flank,
            f"{prefix}_max": maximum,
            f"{prefix}_max_bin": max_bin + 1,
            f"{prefix}_peak_region": peak_region,
            f"{prefix}_detected_bin_fraction": detected_fraction,
            f"{prefix}_top10_bin_signal_share": top10_share,
            f"{prefix}_boundary_vs_interior_log2ratio": boundary_vs_interior,
            f"{prefix}_body_vs_flank_log2ratio": body_vs_flank,
        }
    )


def load_annotations():
    usecols = [
        "accession",
        "subspecies",
        "chr",
        "start",
        "end",
        "TE_id",
        "family",
        "TE_length",
        "gene_proximity",
        "distance_to_gene",
        "mean_CHH_pct",
        "mean_CHG_pct",
        "mean_CpG_pct",
        "mean_24nt_sRNA_log2CPM",
        "mean_21_22nt_sRNA_log2CPM",
    ]
    annotation = pd.read_csv(
        ANNOTATION_FILE, sep="\t", usecols=usecols, low_memory=False
    )
    annotation["length_group"] = pd.cut(
        annotation["TE_length"],
        bins=[-np.inf, 500, 1500, np.inf],
        labels=["Short (<=500 bp)", "Medium (501-1500 bp)", "Long (>1500 bp)"],
    ).astype(str)
    return annotation


def aggregate_profiles(
    metadata: pd.DataFrame,
    matrix21: np.ndarray,
    matrix24: np.ndarray,
):
    rows = []
    accession = metadata["accession"].iloc[0]
    group_cols = ["family", "gene_proximity", "length_group"]
    for keys, index in metadata.groupby(group_cols, observed=True).groups.items():
        index = np.asarray(list(index), dtype=int)
        if index.size == 0:
            continue
        for signal, matrix in (("21-22nt", matrix21), ("24nt", matrix24)):
            sub = matrix[index, :]
            for statistic, profile in (
                ("mean", np.nanmean(sub, axis=0)),
                ("median", np.nanmedian(sub, axis=0)),
            ):
                for bin_index, value in enumerate(profile, start=1):
                    if bin_index <= 40:
                        segment = "upstream"
                        relative_position = -2 + (bin_index - 0.5) * 0.05
                    elif bin_index <= 80:
                        segment = "TE_body"
                        relative_position = (bin_index - 40 - 0.5) / 40
                    else:
                        segment = "downstream"
                        relative_position = (bin_index - 80 - 0.5) * 0.05
                    rows.append(
                        {
                            "accession": accession,
                            "family": keys[0],
                            "gene_proximity": keys[1],
                            "length_group": keys[2],
                            "signal": signal,
                            "statistic": statistic,
                            "bin": bin_index,
                            "segment": segment,
                            "relative_position": relative_position,
                            "profile_value": float(value),
                            "n_TEs": int(index.size),
                        }
                    )
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    annotation = load_annotations()

    files = {}
    for path in MATRIX_DIR.glob("*.matrix.gz"):
        match = FILE_RE.match(path.name)
        if not match:
            continue
        info = match.groupdict()
        key = (info["accession"], info["tcode"], info["signal"], info["tissue"])
        files[key] = path

    accessions = sorted({(key[0], key[1]) for key in files})
    all_metrics = []
    all_profiles = []
    header_records = []

    for accession, tcode in accessions:
        print(f"Processing {accession} ({tcode})...", flush=True)
        loaded = {}
        for signal in ("21-22nt", "24nt"):
            for tissue in ("R", "Y"):
                path = files[(accession, tcode, signal, tissue)]
                header, meta, matrix = read_matrix(path)
                loaded[(signal, tissue)] = (header, meta, matrix)
                header_records.append(
                    {
                        "accession": accession,
                        "t_code": tcode,
                        "signal": signal,
                        "tissue": tissue,
                        "file": str(path),
                        "upstream_bp": header["upstream"][0],
                        "body_scaled_bp": header["body"][0],
                        "downstream_bp": header["downstream"][0],
                        "bin_size": header["bin size"][0],
                        "bin_average_type": header["bin avg type"],
                        "missing_data_as_zero": header["missing data as zero"],
                        "sort_regions": header["sort regions"],
                        "sort_using": header["sort using"],
                        "n_rows": len(meta),
                        "n_samples": len(header["sample_labels"]),
                        "sample_labels": ",".join(header["sample_labels"]),
                    }
                )

        reference_meta = loaded[("21-22nt", "R")][1]
        for key, (_, meta, _) in loaded.items():
            assert_same_rows(reference_meta, meta, f"{accession} {key}")

        matrix21 = np.nanmean(
            np.stack(
                [
                    loaded[("21-22nt", "R")][2],
                    loaded[("21-22nt", "Y")][2],
                ]
            ),
            axis=0,
        )
        matrix24 = np.nanmean(
            np.stack(
                [loaded[("24nt", "R")][2], loaded[("24nt", "Y")][2]]
            ),
            axis=0,
        )

        metrics = reference_meta.copy()
        metrics.insert(0, "t_code", tcode)
        metrics.insert(0, "accession", accession)
        metrics["TE_length_matrix"] = metrics["end"] - metrics["start"]
        metrics = pd.concat(
            [
                metrics.reset_index(drop=True),
                matrix_metrics(matrix21, "sRNA21_22"),
                matrix_metrics(matrix24, "sRNA24"),
            ],
            axis=1,
        )

        annotation_accession = annotation[annotation["accession"] == accession]
        merged = metrics.merge(
            annotation_accession,
            on=["accession", "chr", "start", "end", "TE_id", "family"],
            how="left",
            validate="one_to_one",
        )
        match_rate = merged["gene_proximity"].notna().mean()
        if match_rate < 0.98:
            raise ValueError(
                f"Annotation match rate for {accession} is only {match_rate:.3%}"
            )

        all_metrics.append(merged)
        all_profiles.append(aggregate_profiles(merged, matrix21, matrix24))

    metrics_all = pd.concat(all_metrics, ignore_index=True)
    profiles_all = pd.concat(all_profiles, ignore_index=True)
    headers = pd.DataFrame(header_records)

    metrics_all.to_csv(
        OUT_DIR / "33_matrix_TE_spatial_metrics.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    profiles_all.to_csv(
        OUT_DIR / "34_matrix_spatial_profiles.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    headers.to_csv(
        OUT_DIR / "36_deeptools_matrix_metadata.tsv", sep="\t", index=False
    )

    print("Matrix profile extraction complete.", flush=True)


if __name__ == "__main__":
    main()
