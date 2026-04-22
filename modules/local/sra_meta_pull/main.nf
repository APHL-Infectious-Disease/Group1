process SRA_META {
tag "$meta"
    label 'process_single'

    errorStrategy 'terminate'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(term)
    val database

    output:
    tuple val(meta), path("*.xml") , emit: xml
    tuple val(meta), path("sra_meta.tsv") , emit: sra_meta
    path("sra_meta_top3.tsv") , emit: sra_meta_top
    tuple val("${task.process}"), val('ENTREZDIRECT'), eval('esearch -version 2>&1'), emit: versions_esearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta}"
    def args = task.ext.args ?: '| efetch -format native'

    """
    esearch -db $database -query $term $args > ${prefix}.xml


    # Build header row
    {
        printf "experiment_id\tlibrary_strategy\tlibrary_selection\ttaxon_id\tscientific_name\tcollection_date\tlocation\tlat_lon"
        printf "\n"
    } > sra_meta.tsv

    xtract -input ${prefix}.xml -pattern EXPERIMENT_PACKAGE \\
        -first PRIMARY_ID \\
        -element LIBRARY_STRATEGY \\
        -element LIBRARY_SELECTION \\
        -element TAXON_ID \\
        -element SCIENTIFIC_NAME \\
        -block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE \\
        -if TAG -equals collection_date -element VALUE \\
        -block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE \\
        -if TAG -equals geo_loc_name -element VALUE \\
        -block SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE \\
        -if TAG -equals lat_lon -element VALUE | grep "USA" >> sra_meta.tsv
    
    cat sra_meta.tsv | grep -v Wisconsin | sort -R | tail -n20 | cut -f1 > sra_meta_top3.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch ${prefix}.xml

    """
}
