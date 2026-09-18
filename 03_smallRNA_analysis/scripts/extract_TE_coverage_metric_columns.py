#!/usr/bin/env python3

import csv
import gzip
import sys


COMMON_COLUMNS = ["accession", "family", "gene_proximity", "TE_length"]


def main(input_path: str, output_path: str, signal_column: str) -> None:
    columns = COMMON_COLUMNS + [signal_column]
    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        missing = [column for column in columns if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        with open(output_path, "w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(
                target,
                fieldnames=columns,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                writer.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: extract_TE_coverage_metric_columns.py INPUT.tsv.gz OUTPUT.tsv SIGNAL_COLUMN"
        )
    main(sys.argv[1], sys.argv[2], sys.argv[3])
