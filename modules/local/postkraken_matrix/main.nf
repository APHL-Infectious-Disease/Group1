process POSTKRAKEN_MATRIX {
    tag "postkraken_matrix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.11' :
        'python:3.11' }"

    input:
    path sra_meta_tsv
    path kraken_reports

    output:
    path "metadata_postkraken.csv", emit: matrix
    path "metadata_postkraken_hits_only.csv", emit: hits_only
    path "versions.yml", emit: versions

    script:
    """
    python ${moduleDir}/postkraken_matrix.py \
        --sra-meta ${sra_meta_tsv} \
        --kraken-reports ${kraken_reports.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}