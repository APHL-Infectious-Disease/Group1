#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import re

COLUMN_MAP = {
    "Norovirus": ["norovirus"],
    "Rotavirus": ["rotavirus"],
    "Parvovirus": ["parvovirus"],
    "Enterovirus": ["enterovirus"],
    "Astrovirus": ["astrovirus"],
    "Poliovirus": ["poliovirus"],
    "Sapovirus": ["sapovirus"],
    "Adenovirus": ["adenovirus"],
    "Morbillivirus": ["morbillivirus", "measles virus"],
    "Unclassified": ["unclassified"],
}

TARGET_COLUMNS = [
    "Norovirus", "Rotavirus", "Parvovirus", "Enterovirus", "Astrovirus",
    "Poliovirus", "Sapovirus", "Adenovirus", "Morbillivirus", "Unclassified"
]

VIRUS_TARGET_COLUMNS = [
    "Norovirus", "Rotavirus", "Parvovirus", "Enterovirus", "Astrovirus",
    "Poliovirus", "Sapovirus", "Adenovirus", "Morbillivirus"
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra-meta", required=True)
    parser.add_argument("--kraken-reports", nargs="+", required=True)
    return parser.parse_args()


def infer_sra_id(filepath: str) -> str:
    """
    Robustly infer SRA accession from Kraken report filename.
    Works even if the filename has prefixes/suffixes around the accession.
    """
    name = Path(filepath).name

    # First try to find a real SRA accession anywhere in the filename
    m = re.search(r'([SED]RR\d+)', name, flags=re.IGNORECASE)
    if m:
        return m.group(1).upper()

    # Fallback patterns
    patterns = [
        r"(.+?)\.kraken2\.report\.txt$",
        r"(.+?)\.report\.txt$",
        r"(.+?)\.txt$",
    ]
    for pattern in patterns:
        m = re.match(pattern, name)
        if m:
            return m.group(1)

    return Path(filepath).stem


def normalize(text: str) -> str:
    return str(text).strip().lower()


def parse_kraken_report(report_file: str):
    sra_id = infer_sra_id(report_file)

    best_hits = {
        col: {"percent": 0.0, "taxid": "0", "organism": "", "display": "0"}
        for col in TARGET_COLUMNS
    }

    with open(report_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            try:
                percent = float(parts[0].strip())
            except ValueError:
                continue

            taxid = parts[4].strip()
            organism = parts[5].strip()
            org_norm = normalize(organism)

            for col, keywords in COLUMN_MAP.items():
                if any(k in org_norm for k in keywords):
                    if percent > best_hits[col]["percent"]:
                        best_hits[col] = {
                            "percent": percent,
                            "taxid": taxid,
                            "organism": organism,
                            "display": f"{percent};{taxid};{organism}"
                        }
                    break

    row = {"SRA accession": sra_id}

    detected_count = 0
    for col in TARGET_COLUMNS:
        row[col] = best_hits[col]["display"]
        row[f"{col}_pct"] = best_hits[col]["percent"]
        row[f"{col}_detected"] = 1 if best_hits[col]["percent"] > 0 else 0
        if best_hits[col]["percent"] > 0 and col != "Unclassified":
            detected_count += 1

    row["total_detected_targets"] = detected_count
    return row


def load_metadata(meta_file: str) -> pd.DataFrame:
    df = pd.read_csv(meta_file, dtype=str).fillna("")

    rename_map = {}

    if "Run" in df.columns:
        rename_map["Run"] = "SRA accession"
    elif "run_accession" in df.columns:
        rename_map["run_accession"] = "SRA accession"

    if "ScientificName" in df.columns:
        rename_map["ScientificName"] = "source"
    elif "scientific_name" in df.columns:
        rename_map["scientific_name"] = "source"

    if "Collection_Date" in df.columns:
        rename_map["Collection_Date"] = "collection date"
    elif "collection_date" in df.columns:
        rename_map["collection_date"] = "collection date"

    if "geo_loc_name_country_calc" in df.columns:
        rename_map["geo_loc_name_country_calc"] = "location"
    elif "geo_loc_name" in df.columns:
        rename_map["geo_loc_name"] = "location"
    elif "location" in df.columns:
        rename_map["location"] = "location"

    if "LibraryStrategy" in df.columns:
        rename_map["LibraryStrategy"] = "library_strategy"
    if "LibrarySource" in df.columns:
        rename_map["LibrarySource"] = "library_source"
    if "SampleName" in df.columns:
        rename_map["SampleName"] = "sample_name"
    if "BioProject" in df.columns:
        rename_map["BioProject"] = "bioproject"

    df = df.rename(columns=rename_map)

    if "SRA accession" not in df.columns:
        raise ValueError("Could not find a Run/run_accession/SRA accession column in metadata CSV")

    for col in [
        "source", "collection date", "location",
        "library_strategy", "library_source", "sample_name", "bioproject"
    ]:
        if col not in df.columns:
            df[col] = ""

    keep_cols = [
        "SRA accession", "source", "collection date", "location",
        "library_strategy", "library_source", "sample_name", "bioproject"
    ]

    # Normalize accession formatting
    df["SRA accession"] = df["SRA accession"].astype(str).str.strip().str.upper()

    return df[keep_cols].copy()


def main():
    args = parse_args()

    metadata_df = load_metadata(args.sra_meta)

    kraken_rows = [parse_kraken_report(fp) for fp in args.kraken_reports]
    kraken_df = pd.DataFrame(kraken_rows)

    if kraken_df.empty:
        kraken_df = pd.DataFrame(columns=["SRA accession"] + TARGET_COLUMNS)

    kraken_df["SRA accession"] = kraken_df["SRA accession"].astype(str).str.strip().str.upper()

    processed_accessions = set(kraken_df["SRA accession"].tolist())

    # LEFT merge so metadata rows are preserved even if something didn't process
    final_df = metadata_df.merge(kraken_df, how="left", on="SRA accession")

    final_df["processed"] = final_df["SRA accession"].isin(processed_accessions)

    for col in TARGET_COLUMNS:
        final_df[col] = final_df.get(col, "0").fillna("0")
        final_df[f"{col}_pct"] = pd.to_numeric(
            final_df.get(f"{col}_pct", 0), errors="coerce"
        ).fillna(0.0)
        final_df[f"{col}_detected"] = pd.to_numeric(
            final_df.get(f"{col}_detected", 0), errors="coerce"
        ).fillna(0).astype(int)

    final_df["total_detected_targets"] = pd.to_numeric(
        final_df.get("total_detected_targets", 0), errors="coerce"
    ).fillna(0).astype(int)

    final_df = final_df.sort_values("SRA accession")

    # Full table
    final_df.to_csv("metadata_postkraken.csv", index=False)

    # Hits-only table
    virus_detection_cols = [f"{col}_detected" for col in VIRUS_TARGET_COLUMNS]
    hits_only_df = final_df[
        final_df[virus_detection_cols].sum(axis=1) > 0
    ].copy()

    hits_only_df.to_csv("metadata_postkraken_hits_only.csv", index=False)


if __name__ == "__main__":
    main()