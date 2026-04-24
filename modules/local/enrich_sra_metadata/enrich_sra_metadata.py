#!/usr/bin/env python3

import argparse
import sys
import time
from typing import Dict, List

import pandas as pd
import requests
from lxml import etree

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runinfo", required=True, help="Input SRA RunInfo CSV")
    parser.add_argument("--out", required=True, help="Output enriched metadata CSV")
    parser.add_argument("--candidate-limit", type=int, default=1000, help="Number of raw SRA candidates to inspect before filtering")
    parser.add_argument("--max-runs", type=int, default=15, help="Number of qualified runs to keep")
    parser.add_argument("--max-sra-mb", type=float, default=500, help="Maximum SRA/FASTQ size in MB to keep. Use 0 to disable.")
    parser.add_argument("--usa-only", action="store_true", help="Keep only USA BioSample locations")
    return parser.parse_args()


def chunked(items: List[str], n: int):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def safe_series(df: pd.DataFrame, candidates: List[str], default: str = "") -> pd.Series:
    for col in candidates:
        if col in df.columns:
            return df[col].fillna("")
    return pd.Series([default] * len(df), index=df.index, dtype="string")


def is_usa_location(value: str) -> bool:
    v = str(value).strip().lower()
    return (
        v.startswith("usa")
        or v.startswith("u.s.a")
        or v.startswith("us:")
        or v.startswith("united states")
        or v.startswith("united states of america")
    )


def fetch_biosample_xml(biosample_ids: List[str]) -> bytes:
    ids = ",".join(biosample_ids)
    resp = requests.get(
        NCBI_EFETCH,
        params={
            "db": "biosample",
            "id": ids,
            "retmode": "xml",
        },
        timeout=120,
    )
    resp.raise_for_status()
    return resp.content


def parse_biosample_attributes(xml_bytes: bytes) -> Dict[str, Dict[str, str]]:
    root = etree.fromstring(xml_bytes)
    results: Dict[str, Dict[str, str]] = {}

    for bs in root.xpath(".//BioSample"):
        accession = bs.get("accession", "").strip()
        if not accession:
            continue

        record = {
            "biosample_collection_date": "",
            "biosample_location": "",
            "biosample_isolation_source": "",
            "biosample_lat_lon": "",
        }

        for attr in bs.xpath(".//Attributes/Attribute"):
            name = (attr.get("attribute_name") or "").strip().lower()
            value = "".join(attr.itertext()).strip()

            if name == "collection_date" and not record["biosample_collection_date"]:
                record["biosample_collection_date"] = value
            elif name == "geo_loc_name" and not record["biosample_location"]:
                record["biosample_location"] = value
            elif name == "isolation_source" and not record["biosample_isolation_source"]:
                record["biosample_isolation_source"] = value
            elif name == "lat_lon" and not record["biosample_lat_lon"]:
                record["biosample_lat_lon"] = value

        results[accession] = record

    return results


def enrich_runinfo(
    runinfo_path: str,
    candidate_limit: int,
    max_runs: int,
    max_sra_mb: float,
    usa_only: bool,
) -> pd.DataFrame:
    df = pd.read_csv(runinfo_path, dtype=str).fillna("")

    if "Run" not in df.columns:
        raise ValueError("RunInfo CSV must contain a 'Run' column")

    df["Run"] = df["Run"].astype(str).str.strip().str.upper()

    # Keep only SRR accessions, since fastq-dl was working best with SRR in this workflow.
    df = df[df["Run"].str.startswith("SRR", na=False)].copy()

    # Keep only records with downloadable paths when RunInfo provides this field.
    if "download_path" in df.columns:
        before = len(df)
        df = df[df["download_path"].astype(str).str.strip() != ""].copy()
        print(f"INFO: download_path filter retained {len(df)} of {before} SRR rows", file=sys.stderr)

    # Filter by SRA size before BioSample enrichment/download.
    # RunInfo commonly has a size_MB column.
    if max_sra_mb and max_sra_mb > 0:
        if "size_MB" in df.columns:
            before = len(df)
            df["size_MB_numeric"] = pd.to_numeric(df["size_MB"], errors="coerce")
            df = df[
                df["size_MB_numeric"].notna()
                & (df["size_MB_numeric"] > 0)
                & (df["size_MB_numeric"] <= max_sra_mb)
            ].copy()
            print(
                f"INFO: max_sra_mb filter <= {max_sra_mb} MB retained {len(df)} of {before} rows",
                file=sys.stderr,
            )
        else:
            print(
                "WARNING: RunInfo CSV does not contain size_MB column; max_sra_mb filter was skipped",
                file=sys.stderr,
            )

    if candidate_limit and candidate_limit > 0:
        df = df.head(candidate_limit).copy()

    if df.empty:
        raise ValueError(
            "No candidate SRR accessions remained after initial filters. "
            "Try increasing --max_sra_mb, increasing --sra_candidate_limit, or broadening --sra_query."
        )

    if "BioSample" not in df.columns:
        df["BioSample"] = ""

    biosamples = sorted({x.strip() for x in df["BioSample"].tolist() if x.strip()})
    biosample_map: Dict[str, Dict[str, str]] = {}

    for batch in chunked(biosamples, 200):
        try:
            xml_bytes = fetch_biosample_xml(batch)
            biosample_map.update(parse_biosample_attributes(xml_bytes))
        except Exception as e:
            print(f"WARNING: failed to fetch BioSample batch: {e}", file=sys.stderr)
        time.sleep(0.34)

    def lookup(bs_id: str, key: str) -> str:
        bs_id = str(bs_id or "").strip()
        if not bs_id:
            return ""
        return biosample_map.get(bs_id, {}).get(key, "")

    out = pd.DataFrame()
    out["SRA accession"] = safe_series(df, ["Run"])
    out["BioSample"] = safe_series(df, ["BioSample"])
    out["source"] = safe_series(df, ["source", "Source", "ScientificName", "scientific_name"])
    out["collection date"] = safe_series(df, ["collection_date", "Collection_Date"])
    out["location"] = safe_series(df, ["location", "geo_loc_name", "geo_loc_name_country_calc"])
    out["library_strategy"] = safe_series(df, ["LibraryStrategy", "library_strategy"])
    out["library_source"] = safe_series(df, ["LibrarySource", "library_source"])
    out["library_layout"] = safe_series(df, ["LibraryLayout", "library_layout"])
    out["sample_name"] = safe_series(df, ["SampleName", "sample_name"])
    out["bioproject"] = safe_series(df, ["BioProject", "bioproject"])
    out["download_path"] = safe_series(df, ["download_path"])
    out["size_MB"] = safe_series(df, ["size_MB", "size_mb"])

    out["biosample_collection_date"] = out["BioSample"].map(lambda x: lookup(x, "biosample_collection_date"))
    out["biosample_location"] = out["BioSample"].map(lambda x: lookup(x, "biosample_location"))
    out["biosample_isolation_source"] = out["BioSample"].map(lambda x: lookup(x, "biosample_isolation_source"))
    out["biosample_lat_lon"] = out["BioSample"].map(lambda x: lookup(x, "biosample_lat_lon"))

    out["collection date"] = out["collection date"].mask(
        out["collection date"].astype(str).str.strip() == "",
        out["biosample_collection_date"],
    )

    out["location"] = out["location"].mask(
        out["location"].astype(str).str.strip() == "",
        out["biosample_location"],
    )

    out["source"] = out["source"].mask(
        out["source"].astype(str).str.strip() == "",
        out["biosample_isolation_source"],
    )

    out["usa_sample"] = out["location"].map(is_usa_location)

    if usa_only:
        before = len(out)
        out = out[out["usa_sample"]].copy()
        print(f"INFO: USA BioSample filter retained {len(out)} of {before} candidate runs", file=sys.stderr)

    if max_runs and max_runs > 0:
        if len(out) < max_runs:
            print(
                f"WARNING: Requested {max_runs} qualified USA runs but only found {len(out)} "
                f"within candidate limit {candidate_limit} and max_sra_mb {max_sra_mb}. Proceeding with {len(out)}.",
                file=sys.stderr,
            )
        out = out.head(max_runs).copy()

    if out.empty:
        raise ValueError(
            "No SRA accessions remained after size + BioSample USA filtering. "
            "Try increasing --max_sra_mb, increasing --sra_candidate_limit, or broadening --sra_query."
        )

    out["SRA accession"] = out["SRA accession"].astype(str).str.strip().str.upper()

    return out


def main():
    args = parse_args()
    out_df = enrich_runinfo(
        runinfo_path=args.runinfo,
        candidate_limit=args.candidate_limit,
        max_runs=args.max_runs,
        max_sra_mb=args.max_sra_mb,
        usa_only=args.usa_only,
    )
    out_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
