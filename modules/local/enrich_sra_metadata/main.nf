process ENRICH_SRA_METADATA {
    tag "enrich_sra_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    path runinfo_csv

    output:
    path "sra_metadata_enriched.csv", emit: csv

    script:
    """
    python ${moduleDir}/enrich_sra_metadata.py \
        --runinfo ${runinfo_csv} \
        --out sra_metadata_enriched.csv
    """
}
