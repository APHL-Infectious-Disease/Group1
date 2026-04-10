process SRA_META {
    tag "Gathering SRA Runs"
    label 'process_single'

    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:24.0--he881be0_0':
        'biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    val run

    output:
    path("*.xml") , emit: xml
    path("*.tsv") , emit: tsv
    tuple val("${task.process}"), val('ENTREZDIRECT'), eval('esearch -version 2>&1'), emit: versions_esearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """

    bash /workspaces/Group1/modules/local/sra_meta_pull/sra_meta_pull.sh

    """

    stub:
    """
    touch sra.xml

    """
}
