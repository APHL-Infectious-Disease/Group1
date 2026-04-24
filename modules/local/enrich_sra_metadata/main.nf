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
    def candidate_limit = params.sra_candidate_limit ?: 1000
    def max_runs = params.max_runs ?: 15
    def max_sra_mb = params.max_sra_mb ?: 500

    """
    python ${moduleDir}/enrich_sra_metadata.py \
        --runinfo ${runinfo_csv} \
        --out sra_metadata_enriched.csv \
        --candidate-limit ${candidate_limit} \
        --max-runs ${max_runs} \
        --max-sra-mb ${max_sra_mb} \
        --usa-only
    """
}

