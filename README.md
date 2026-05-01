# APHL-Infectious-Disease/group1

> **APHL-Infectious-Disease/group1**
>
> This repository contains a prototype pipeline developed for the 2026 APHL Hackathon by Group 1. All code, data products, and visualizations are for demonstration and development purposes only and should not be used for production, clinical, regulatory, or public health decision-making at this time.

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/APHL-Infectious-Disease/group1)
[![GitHub Actions CI Status](https://github.com/APHL-Infectious-Disease/group1/actions/workflows/nf-test.yml/badge.svg)](https://github.com/APHL-Infectious-Disease/group1/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/APHL-Infectious-Disease/group1/actions/workflows/linting.yml/badge.svg)](https://github.com/APHL-Infectious-Disease/group1/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

---

## Introduction

**APHL-Infectious-Disease/group1** is a Nextflow-based bioinformatics pipeline developed for the 2026 APHL Hackathon.

The pipeline enables re-analysis of publicly available sequencing data—particularly **wastewater metagenomic datasets**—to detect **measles and enteric viruses of public health interest**, including:

- Norovirus  
- Rotavirus  
- Enterovirus  
- Astrovirus  
- Adenovirus  
- Sapovirus  
- Poliovirus  
- Morbillivirus  

---

## Key Features

### Multiple Input Modes

| Mode | Description |
|------|------------|
| `--mode sra` | Automatically discover datasets from NCBI SRA |
| `--mode accessions` | Use a user-provided list of SRR/ERR accessions |
| `--mode samplesheet` | Provide local FASTQ files |

---

### Metadata Enrichment

- Pulls **RunInfo metadata from SRA**
- Enriches with **BioSample metadata via NCBI E-utilities**
- Extracts:
  - collection date
  - geographic location
  - isolation source
  - lat/lon (when available)
- Filters to **U.S.-based samples**

---

### Data Processing

- Parallel FASTQ download via `fastq-dl`
- Optional preprocessing (QC + trimming)
- Kraken2 classification against viral database

---

### Post-classification Analysis

Generates:

- `metadata_postkraken.csv`
- `metadata_postkraken_hits_only.csv`

---

### Optional Steps

| Feature | Parameter |
|--------|----------|
| Preprocessing | `--run_preprocessing` |
| MultiQC report | `--run_multiqc` |

---

## Workflow Overview

```text
SRA / Accessions / Samplesheet
              ↓
         FASTQ Download
              ↓
        (Optional QC)
              ↓
           Kraken2
              ↓
     Metadata Enrichment
              ↓
    PostKraken Matrix Build
              ↓
       Results Dashboard
```

## Subworkflows and Modules

### SRA_DISCOVERY
- Uses Entrez to search SRA  
- Retrieves RunInfo metadata  
- Outputs accession list + metadata  

### FETCH_READS (FASTQDL)
- Downloads FASTQ files from SRA  
- Supports parallel retrieval  

### MAKE_ACCESSIONS_METADATA
- Builds metadata for user-provided accessions  
- Pulls SRA + BioSample metadata  

### ENRICH_SRA_METADATA
- Enhances SRA metadata with BioSample fields  
- Improves completeness (dates, locations, etc.)  

### PREPROCESS_READS (Optional)
- Tools: `fastp`, `FastQC`, `seqkit`  
- Performs trimming + QC  

### PREPARE_KRAKEN_DB
- Accepts:
  - local DB (`--kraken2_db`)
  - remote DB (`--kraken2_db_url`)  
- Downloads and prepares database if needed  

### CLASSIFY_READS (Kraken2)
- Runs Kraken2 classification  
- Produces per-sample reports  

### BUILD_POSTKRAKEN_MATRIX
- Combines:
  - Kraken reports  
  - enriched metadata  
- Generates final matrices  

### POSTKRAKEN_MATRIX
- Parses Kraken outputs  
- Produces:
  - full matrix  
  - hits-only matrix  

### QC REPORT (Optional)
- Generates MultiQC report  

### RESULTS DASHBOARD
```bash
bash shiny.sh
```
- shiny.sh launches app.R, which combines and saves kraken files as a single summary file and generates an associated dashboard from the kraken summary file and SRA metadata file 
- This should automatically open the associated dashboard html. If not, got to Ports, select 8080 and click "open in browser" (the globe icon)

In a new codespace, it may be necissary to first run the following commands in the terminal:
```bash
conda install -c conda-forge r-base
conda install -c conda-forge r-tidyverse r-shiny r-leaflet r-thematic r-DT
```


---

## Usage

### SRA Mode

```bash
nextflow run main.nf \
  -profile conda \
  --mode sra \
  --kraken2_db assets/kraken2db_v2 \
  --max_runs 20 \
  --outdir results
```

### Accessions Mode

```bash
nextflow run main.nf \
  -profile conda \
  --mode accessions \
  --accessions docs/accessions.txt \
  --kraken2_db assets/kraken2db_v2 \
  --outdir results
```

### Samplesheet Mode

```bash
nextflow run main.nf \
  -profile conda \
  --mode samplesheet \
  --input path/to/samplesheet.csv \
  --kraken2_db assets/kraken2db_v2 \
  --outdir results
```

### Using a Remote Kraken2 Database

```bash
nextflow run main.nf \
  -profile conda \
  --mode samplesheet \
  --input path/to/samplesheet.csv \
  --kraken2_db_url https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz \
  --outdir results
```

## Optional Parameters

| Parameter | Description |
|----------|------------|
| `--run_preprocessing` | Enable QC + trimming |
| `--run_multiqc` | Generate MultiQC report |
| `--max_runs` | Limit SRA runs |
| `--subsample_size` | Subsample reads |

---

## Output

| File | Description |
|------|------------|
| `metadata_postkraken.csv` | Full detection matrix |
| `metadata_postkraken_hits_only.csv` | Filtered detections |
| `multiqc_report.html` | QC summary (optional) |
| `pipeline_info/` | Execution logs + versions |

## Known Limitations

- SRA metadata may be incomplete or inconsistent  
- BioSample metadata availability varies  
- Detection thresholds not yet standardized  
- Prototype only — not production ready  

---

## Future Directions

- Improve SRA query strategy  
- Add detection confidence scoring  
- Expand pathogen panel  
- Integrate dashboards  
- Optimize performance at scale  

## Credits

APHL-Infectious-Disease/group1 was written by Group 1.


This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
