#!/usr/bin/env python3
from pathlib import Path
import csv
import re
from collections import Counter, defaultdict


ROOT = Path(__file__).resolve().parents[2]
TE_GFF = ROOT / "Figure6" / "pan-TE" / "TE-GFF"
OUTDIR = Path(__file__).resolve().parent


def parse_attrs(text):
    attrs = {}
    for item in text.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def parse_gtf(path):
    elements = defaultdict(lambda: {
        "features": Counter(),
        "repeat": None,
        "ltrrt": None,
        "attrs": {},
        "ltr_ids": [],
        "tsd_ids": [],
    })

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attrs(attr_text)
            feature_id = attrs.get("ID", "")
            parent = attrs.get("Parent")

            if feature == "repeat_region":
                parent = feature_id
            if not parent:
                continue

            record = elements[parent]
            record["features"][feature] += 1
            if feature == "repeat_region":
                record["repeat"] = fields
                record["attrs"] = attrs
            elif feature.endswith("LTR_retrotransposon") and attrs.get("method") == "structural":
                record["ltrrt"] = fields
            elif feature == "long_terminal_repeat":
                record["ltr_ids"].append(feature_id)
            elif feature == "target_site_duplication":
                record["tsd_ids"].append(feature_id)

    return elements


def young_ids_for(family):
    ids = defaultdict(set)
    for path in sorted((TE_GFF / family).glob(f"*_{family}_young.gtf")):
        accession = path.name.split(f"_TE_{family}_young.gtf")[0]
        for parent, rec in parse_gtf(path).items():
            if is_intact_basic(rec):
                ids[accession].add(parent)
    return ids


def is_intact_basic(rec):
    features = rec["features"]
    ltrrt = sum(v for k, v in features.items() if k.endswith("LTR_retrotransposon"))
    return (
        features["repeat_region"] == 1
        and features["long_terminal_repeat"] == 2
        and ltrrt == 1
        and rec["ltrrt"] is not None
    )


def is_intact_strict_tsd(rec):
    return is_intact_basic(rec) and rec["features"]["target_site_duplication"] == 2


def accession_from_name(path, family):
    suffix = f"_TE_{family}.gtf"
    return path.name[:-len(suffix)]


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return ""


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    young_by_family = {
        "Copia": young_ids_for("Copia"),
        "Gypsy": young_ids_for("Gypsy"),
    }

    rows = []
    for family in ("Copia", "Gypsy"):
        for path in sorted((TE_GFF / family).glob(f"*_TE_{family}.gtf")):
            if path.name.endswith("_young.gtf"):
                continue
            accession = accession_from_name(path, family)
            elements = parse_gtf(path)
            for repeat_id, rec in sorted(elements.items()):
                if not is_intact_basic(rec):
                    continue
                repeat = rec["repeat"]
                ltrrt = rec["ltrrt"]
                attrs = rec["attrs"]
                seqid = ltrrt[0]
                start = int(ltrrt[3])
                end = int(ltrrt[4])
                ltrrt_attrs = parse_attrs(ltrrt[8])
                ltrrt_id = ltrrt_attrs.get("ID", "")
                rows.append({
                    "Accession": accession,
                    "Family": family,
                    "Repeat_ID": repeat_id,
                    "LTRRT_ID": ltrrt_id,
                    "Chr": seqid,
                    "Start": start,
                    "End": end,
                    "Strand": ltrrt[6],
                    "Length_bp": end - start + 1,
                    "LTR_identity": to_float(attrs.get("ltr_identity")),
                    "Motif": attrs.get("motif", ""),
                    "TSD": attrs.get("tsd", ""),
                    "LTR_count": rec["features"]["long_terminal_repeat"],
                    "TSD_count": rec["features"]["target_site_duplication"],
                    "Intact_basic": "TRUE",
                    "Intact_strict_TSD": "TRUE" if is_intact_strict_tsd(rec) else "FALSE",
                    "Young_intact": "TRUE" if repeat_id in young_by_family[family][accession] else "FALSE",
                })

    element_path = OUTDIR / "intact_ltr_elements.tsv"
    with element_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouped = defaultdict(lambda: {
        "Intact_count": 0,
        "Strict_TSD_count": 0,
        "Young_intact_count": 0,
        "Young_strict_TSD_count": 0,
        "Intact_length_bp": 0,
        "Young_intact_length_bp": 0,
    })
    for row in rows:
        key = (row["Accession"], row["Family"])
        grouped[key]["Intact_count"] += 1
        grouped[key]["Strict_TSD_count"] += row["Intact_strict_TSD"] == "TRUE"
        grouped[key]["Young_intact_count"] += row["Young_intact"] == "TRUE"
        grouped[key]["Young_strict_TSD_count"] += (
            row["Young_intact"] == "TRUE" and row["Intact_strict_TSD"] == "TRUE"
        )
        grouped[key]["Intact_length_bp"] += row["Length_bp"]
        if row["Young_intact"] == "TRUE":
            grouped[key]["Young_intact_length_bp"] += row["Length_bp"]

    summary_rows = []
    accession_order = [
        "Basmati1", "CG14", "G46", "IR64", "Lemont", "LJ",
        "MSU", "N22", "NamRoo", "TM", "Tumba", "WSSM",
    ]
    for accession in accession_order:
        for family in ("Copia", "Gypsy"):
            data = grouped[(accession, family)]
            intact = data["Intact_count"]
            young = data["Young_intact_count"]
            summary_rows.append({
                "Accession": accession,
                "Family": family,
                **data,
                "Young_fraction": young / intact if intact else "",
                "Strict_TSD_fraction": data["Strict_TSD_count"] / intact if intact else "",
            })

    summary_path = OUTDIR / "intact_ltr_count_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {element_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
