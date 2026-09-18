#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ["CENTROMERIC_SRNA_ANALYSIS_DIR"])
MATRIX_DIR = Path(os.environ["DEEPTOOLS_MATRIX_DIR"])
CENTROMERE_FILE = Path(os.environ["CENTROMERE_BED"])
OUT_DIR = Path(os.environ.get("CENTROMERIC_SRNA_OUTDIR", ROOT / "data"))

MATRIX_FILES = {
    ("24nt", "R"): MATRIX_DIR / "MSU_T6_24nt_R_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz",
    ("24nt", "Y"): MATRIX_DIR / "MSU_T6_24nt_Y_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz",
    ("21-22nt", "R"): MATRIX_DIR / "MSU_T6_21-22nt_R_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz",
    ("21-22nt", "Y"): MATRIX_DIR / "MSU_T6_21-22nt_Y_TE_scale2kb_up2kb_down2kb_bin50.matrix.gz",
}

LENGTH_LEVELS = [
    "Short (<=500 bp)",
    "Medium (501-1500 bp)",
    "Long (>1500 bp)",
]


def read_matrix(path: Path):
    with gzip.open(path, "rt") as handle:
        first = handle.readline().rstrip("\n")
        if not first.startswith("@"):
            raise ValueError(f"Missing deepTools JSON header in {path}")
        header = json.loads(first[1:])
        sample_boundaries = header["sample_boundaries"]
        group_boundaries = header["group_boundaries"]
        group_labels = header["group_labels"]
        n_bins = sample_boundaries[1] - sample_boundaries[0]
        n_samples = len(sample_boundaries) - 1

        metadata = []
        profiles = []
        row_index = 0
        group_index = 0
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            while (
                group_index + 1 < len(group_boundaries) - 1
                and row_index >= group_boundaries[group_index + 1]
            ):
                group_index += 1
            family = "Copia" if "Copia" in group_labels[group_index] else "Gypsy"
            values = np.asarray(fields[6:], dtype=np.float32)
            if values.size != n_samples * n_bins:
                raise ValueError(
                    f"Unexpected number of values in {path.name}, row {row_index}"
                )
            values = values.reshape(n_samples, n_bins)
            profiles.append(np.nanmean(values, axis=0).astype(np.float32))
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
            row_index += 1

    meta = pd.DataFrame(
        metadata,
        columns=["chr", "start", "end", "TE_id", "strand", "family"],
    )
    return header, meta, np.vstack(profiles).astype(np.float32)


def assert_same_rows(reference: pd.DataFrame, other: pd.DataFrame, label: str):
    columns = ["chr", "start", "end", "TE_id", "strand", "family"]
    if len(reference) != len(other) or not reference[columns].equals(other[columns]):
        raise ValueError(f"TE row order mismatch: {label}")


def load_centromeres():
    cent = pd.read_csv(
        CENTROMERE_FILE,
        sep="\t",
        header=None,
        names=["chr", "centromere_start", "centromere_end"],
    )
    return cent.set_index("chr")


def annotate_chromatin_context(meta: pd.DataFrame, centromeres: pd.DataFrame):
    out = meta.copy()
    out["TE_length"] = out["end"] - out["start"]
    out["midpoint"] = (out["start"] + out["end"]) / 2
    out["distance_to_centromere"] = np.nan

    for chrom, index in out.groupby("chr").groups.items():
        if chrom not in centromeres.index:
            continue
        row = centromeres.loc[chrom]
        cent_start = float(row["centromere_start"])
        cent_end = float(row["centromere_end"])
        midpoint = out.loc[index, "midpoint"].to_numpy()
        distance = np.where(
            (midpoint >= cent_start) & (midpoint <= cent_end),
            0.0,
            np.minimum(np.abs(midpoint - cent_start), np.abs(midpoint - cent_end)),
        )
        out.loc[index, "distance_to_centromere"] = distance

    out["chromatin_context"] = np.select(
        [
            out["distance_to_centromere"] <= 2_000_000,
            out["distance_to_centromere"] <= 5_000_000,
            out["distance_to_centromere"] > 5_000_000,
        ],
        ["Centromeric/pericentromeric", "Transition", "Distal arms"],
        default="Unclassified",
    )
    out["region_group"] = np.where(
        out["chromatin_context"] == "Centromeric/pericentromeric",
        "Centromeric/pericentromeric (<=2 Mb)",
        np.where(
            out["chromatin_context"] == "Unclassified",
            "Unclassified",
            "Other chromosome regions (>2 Mb)",
        ),
    )
    out["length_group"] = pd.cut(
        out["TE_length"],
        bins=[-np.inf, 500, 1500, np.inf],
        labels=LENGTH_LEVELS,
    ).astype(str)
    return out


def bin_annotation(bin_index: int):
    if bin_index < 40:
        return "Upstream", -2.0 + (bin_index + 0.5) * 0.05
    if bin_index < 80:
        return "TE body", (bin_index - 40 + 0.5) / 40
    return "Downstream", 1.0 + (bin_index - 80 + 0.5) * 0.05


def summarize_profiles(meta: pd.DataFrame, matrices: dict[str, np.ndarray]):
    rows = []
    selected = meta[meta["region_group"] != "Unclassified"].copy()
    for (region_group, family, length_group), index in selected.groupby(
        ["region_group", "family", "length_group"], observed=True
    ).groups.items():
        index = np.asarray(list(index), dtype=int)
        for signal, matrix in matrices.items():
            values = matrix[index, :]
            n_te = values.shape[0]
            mean = np.nanmean(values, axis=0)
            valid_n = np.sum(np.isfinite(values), axis=0)
            sd = np.nanstd(values, axis=0, ddof=1)
            se = np.divide(
                sd,
                np.sqrt(valid_n),
                out=np.full_like(sd, np.nan, dtype=float),
                where=valid_n > 1,
            )
            for i in range(values.shape[1]):
                segment, relative_position = bin_annotation(i)
                rows.append(
                    {
                        "accession": "MSU",
                        "signal": signal,
                        "region_group": region_group,
                        "family": family,
                        "length_group": length_group,
                        "bin": i + 1,
                        "segment": segment,
                        "relative_position": relative_position,
                        "n_TEs": n_te,
                        "mean_profile": float(mean[i]),
                        "SE_across_TEs": float(se[i]),
                        "CI95_lower": float(max(0.0, mean[i] - 1.96 * se[i])),
                        "CI95_upper": float(mean[i] + 1.96 * se[i]),
                    }
                )
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in MATRIX_FILES.values():
        if not path.exists():
            raise FileNotFoundError(path)
    if not CENTROMERE_FILE.exists():
        raise FileNotFoundError(CENTROMERE_FILE)

    loaded = {}
    for key, path in MATRIX_FILES.items():
        print(f"Reading {path.name}", flush=True)
        loaded[key] = read_matrix(path)

    reference_meta = loaded[("24nt", "R")][1]
    for key, (_, meta, _) in loaded.items():
        assert_same_rows(reference_meta, meta, str(key))

    # Each tissue contains two replicates. Average replicates within tissue,
    # then average root and young panicle so the two tissues have equal weight.
    matrices = {
        signal: np.nanmean(
            np.stack(
                [loaded[(signal, "R")][2], loaded[(signal, "Y")][2]]
            ),
            axis=0,
        )
        for signal in ("24nt", "21-22nt")
    }

    meta = annotate_chromatin_context(reference_meta, load_centromeres())
    profiles = summarize_profiles(meta, matrices)

    counts = (
        meta[meta["region_group"] != "Unclassified"]
        .groupby(["region_group", "family", "length_group"], observed=True)
        .size()
        .reset_index(name="n_TEs")
    )

    profiles.to_csv(
        OUT_DIR / "MSU_centromere2Mb_sRNA_mean_profiles.tsv",
        sep="\t",
        index=False,
    )
    counts.to_csv(
        OUT_DIR / "MSU_centromere_distance_TE_counts.tsv",
        sep="\t",
        index=False,
    )
    meta.to_csv(
        OUT_DIR / "MSU_TE_centromere_context_annotation.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print("MSU centromeric/pericentromeric profile extraction complete.")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
