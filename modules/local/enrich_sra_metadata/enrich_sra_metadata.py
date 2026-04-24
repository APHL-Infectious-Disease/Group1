#!/usr/bin/env python3

import argparse
import io
import sys
import time
from typing import Dict, List
import pandas as pd
import requests
from lxml import etree

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runinfo", required=True, help="Input SRA RunInfo CSV")
    parser.add_argument("--out", required=True, help="Output enriched metadata CSV")
    return parser.parse_args()


def chunked(items: List[str], n: int) -> List[List[str]]:
    for i in range(0, len(items), n):
        yield items[i:i+n]


def safe_series(df: pd.DataFrame, candidates: List[str], default: str = "") -> pd.Series:
    for col in candidates:
        if col in df.columns:
            return df[col].fillna("")
    return pd.Series([default] * len(df), index=df.index, dtype="string")


def fetch_biosample_xml(biosample_ids: List[str]) -> bytes:
    ids = ",".join(biosample_ids)
    resp = requests.get(
        NCBI_EFETCH,
        params={
            "db": "biosample",
            "id": ids,
            "retmode": "xml"
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

        attrs = bs.xpath(".//Attributes/Attribute")
        for attr in attrs:
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


def enrich_runinfo(runinfo_path: str) -> pd.DataFrame:
    df = pd.read_csv(runinfo_path, dtype=str).fillna("")

    if "Run" not in df.columns:
        raise ValueError("RunInfo CSV must contain a 'Run' column")

    biosample_col = "BioSample" if "BioSample" in df.columns else None
    if biosample_col is None:
        df["BioSample"] = ""
        biosample_col = "BioSample"

    biosamples = sorted({x.strip() for x in df[biosample_col].tolist() if x.strip()})
    biosample_map: Dict[str, Dict[str, str]] = {}

    # Batch BioSample requests to be polite / robust
    for batch in chunked(biosamples, 200):
        try:
            xml_bytes = fetch_biosample_xml(batch)
            biosample_map.update(parse_biosample_attributes(xml_bytes))
        except Exception as e:
            print(f"WARNING: failed to fetch BioSample batch: {e}", file=sys.stderr)
        time.sleep(0.34)

    def lookup(bs_id: str, key: str) -> str:
        bs_id = (bs_id or "").strip()
        if not bs_id:
            return ""
        return biosample_map.get(bs_id, {}).get(key, "")

    # Build a clean output table for downstream use
    out = pd.DataFrame()
    out["SRA accession"] = safe_series(df, ["Run"])
    out["BioSample"] = safe_series(df, ["BioSample"])
    out["source"] = safe_series(df, ["source", "Source", "ScientificName", "scientific_name"])
    out["collection date"] = safe_series(df, ["collection_date", "Collection_Date"])
    out["location"] = safe_series(df, ["location", "geo_loc_name", "geo_loc_name_country_calc"])
    out["library_strategy"] = safe_series(df, ["LibraryStrategy", "library_strategy"])
    out["library_source"] = safe_series(df, ["LibrarySource", "library_source"])
    out["sample_name"] = safe_series(df, ["SampleName", "sample_name"])
    out["bioproject"] = safe_series(df, ["BioProject", "bioproject"])
    out["biosample_collection_date"] = out["BioSample"].map(lambda x: lookup(x, "biosample_collection_date"))
    out["biosample_location"] = out["BioSample"].map(lambda x: lookup(x, "biosample_location"))
    out["biosample_isolation_source"] = out["BioSample"].map(lambda x: lookup(x, "biosample_isolation_source"))
    out["biosample_lat_lon"] = out["BioSample"].map(lambda x: lookup(x, "biosample_lat_lon"))

    # Prefer BioSample metadata when SRA RunInfo fields are blank
    out["collection date"] = out["collection date"].mask(
        out["collection date"].astype(str).str.strip() == "",
        out["biosample_collection_date"]
    )
    out["location"] = out["location"].mask(
        out["location"].astype(str).str.strip() == "",
        out["biosample_location"]
    )
    out["source"] = out["source"].mask(
        out["source"].astype(str).str.strip() == "",
        out["biosample_isolation_source"]
    )

    def is_usa_location(value: str) -> bool:
        v = str(value).strip().lower()
        return (
            v.startswith("usa") or
            v.startswith("united states") or
            v.startswith("united states of america")
        )

    out["usa_sample"] = out["location"].map(is_usa_location)

    return out

def main():
    args = parse_args()
    out_df = enrich_runinfo(args.runinfo)
    out_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
