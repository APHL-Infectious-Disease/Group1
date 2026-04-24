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
    parser.add_argument("--accessions", required=True, help="Text file of SRR/ERR/DRR accessions")
    parser.add_argument("--out", required=True, help="Output metadata CSV")
    return parser.parse_args()


def chunked(items: List[str], n: int):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def is_usa_location(value: str) -> bool:
    v = str(value).strip().lower()
    return (
        v.startswith("usa")
        or v.startswith("u.s.a")
        or v.startswith("us:")
        or v.startswith("united states")
        or v.startswith("united states of america")
    )


def fetch_runinfo(accessions: List[str]) -> pd.DataFrame:
    ids = ",".join(accessions)
    resp = requests.get(
        NCBI_EFETCH,
        params={
            "db": "sra",
            "id": ids,
            "rettype": "runinfo",
            "retmode": "text",
        },
        timeout=120,
    )
    resp.raise_for_status()

    text = resp.text.strip()
    if not text:
        return pd.DataFrame(
            columns=[
                "Run", "BioSample", "BioProject", "LibraryStrategy",
                "LibrarySource", "LibraryLayout", "SampleName", "ScientificName"
            ]
        )

    from io import StringIO
    return pd.read_csv(StringIO(text), dtype=str).fillna("")


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


def main():
    args = parse_args()

    with open(args.accessions, "r", encoding="utf-8") as handle:
        accession_list = [
            line.strip().upper()
            for line in handle
            if line.strip() and not line.strip().startswith("#")
        ]

    if not accession_list:
        raise ValueError("No accessions found in accessions file")

    runinfo_batches = []

    for batch in chunked(accession_list, 200):
        try:
            runinfo_batches.append(fetch_runinfo(batch))
        except Exception as e:
            print(f"WARNING: failed to fetch runinfo batch: {e}", file=sys.stderr)
        time.sleep(0.34)

    if runinfo_batches:
        runinfo_df = pd.concat(runinfo_batches, ignore_index=True).fillna("")
    else:
        runinfo_df = pd.DataFrame(
            columns=[
                "Run", "BioSample", "BioProject", "LibraryStrategy",
                "LibrarySource", "LibraryLayout", "SampleName", "ScientificName"
            ]
        )

    if "Run" in runinfo_df.columns:
        runinfo_df["Run"] = runinfo_df["Run"].astype(str).str.strip().str.upper()

    biosamples = sorted({
        x.strip()
        for x in runinfo_df.get("BioSample", pd.Series([], dtype=str)).tolist()
        if x.strip()
    })

    biosample_map: Dict[str, Dict[str, str]] = {}

    for batch in chunked(biosamples, 200):
        try:
            xml_bytes = fetch_biosample_xml(batch)
            biosample_map.update(parse_biosample_attributes(xml_bytes))
        except Exception as e:
            print(f"WARNING: failed to fetch BioSample batch: {e}", file=sys.stderr)
        time.sleep(0.34)

    requested_df = pd.DataFrame({"Run": accession_list})
    merged = requested_df.merge(runinfo_df, how="left", on="Run").fillna("")

    def lookup(bs_id: str, key: str) -> str:
        bs_id = str(bs_id or "").strip()
        if not bs_id:
            return ""
        return biosample_map.get(bs_id, {}).get(key, "")

    out = pd.DataFrame()
    out["SRA accession"] = merged["Run"]
    out["BioSample"] = merged.get("BioSample", "")
    out["source"] = merged.get("ScientificName", "")
    out["collection date"] = out["BioSample"].map(lambda x: lookup(x, "biosample_collection_date"))
    out["location"] = out["BioSample"].map(lambda x: lookup(x, "biosample_location"))
    out["library_strategy"] = merged.get("LibraryStrategy", "")
    out["library_source"] = merged.get("LibrarySource", "")
    out["library_layout"] = merged.get("LibraryLayout", "")
    out["sample_name"] = merged.get("SampleName", "")
    out["bioproject"] = merged.get("BioProject", "")
    out["biosample_isolation_source"] = out["BioSample"].map(lambda x: lookup(x, "biosample_isolation_source"))
    out["biosample_lat_lon"] = out["BioSample"].map(lambda x: lookup(x, "biosample_lat_lon"))

    out["source"] = out["source"].mask(
        out["source"].astype(str).str.strip() == "",
        out["biosample_isolation_source"],
    )

    out["usa_sample"] = out["location"].map(is_usa_location)

    out.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()

