/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML }   from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SRA_DISCOVERY }            from '../subworkflows/local/sra_discovery'
include { FETCH_READS }              from '../subworkflows/local/fetch_reads'
include { PREPROCESS_READS }         from '../subworkflows/local/preprocess_reads'
include { CLASSIFY_READS }           from '../subworkflows/local/classify_reads'
include { BUILD_POSTKRAKEN_MATRIX }  from '../subworkflows/local/build_postkraken_matrix'
include { REPORTING }                from '../subworkflows/local/reporting'
include { PREPARE_KRAKEN_DB }        from '../subworkflows/local/prepare_kraken_db'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAKE_ACCESSIONS_METADATA } from '../modules/local/make_accessions_metadata/main'
include { ENRICH_SRA_METADATA }      from '../modules/local/enrich_sra_metadata/main'

workflow GROUP1 {

    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUT SELECTION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'sra') {
        SRA_DISCOVERY()
        ch_input = SRA_DISCOVERY.out.reads

        FETCH_READS(
            ch_input
        )
        ch_reads = FETCH_READS.out.reads

        ENRICH_SRA_METADATA(
            SRA_DISCOVERY.out.tsv
        )

    } else if (params.mode == 'accessions') {

        if (!params.accessions) {
            error "When --mode accessions, you must provide --accessions <file>"
        }

        ch_accessions_file = Channel.fromPath(params.accessions, checkIfExists: true)

        ch_input = ch_accessions_file
            .splitText()
            .map { it.trim() }
            .filter { it && !it.startsWith('#') }
            .map { accession ->
                def meta = [
                    id         : accession,
                    sample     : accession,
                    single_end : false
                ]
                tuple(meta, accession)
            }

        FETCH_READS(
            ch_input
        )
        ch_reads = FETCH_READS.out.reads

        MAKE_ACCESSIONS_METADATA(
            ch_accessions_file
        )

    } else if (params.mode == 'samplesheet') {
        ch_reads = ch_samplesheet

    } else {
        error "Unsupported params.mode: ${params.mode}. Use 'sra', 'accessions', or 'samplesheet'."
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OPTIONAL PREPROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.run_preprocessing) {
        PREPROCESS_READS(
            ch_reads
        )

        ch_reads_for_kraken = PREPROCESS_READS.out.reads

    } else {
        ch_reads_for_kraken = ch_reads
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREPARE KRAKEN2 DB
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    PREPARE_KRAKEN_DB()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CLASSIFICATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    CLASSIFY_READS(
        ch_reads_for_kraken,
        PREPARE_KRAKEN_DB.out.db
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        POST-KRAKEN MATRIX
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'sra') {
        BUILD_POSTKRAKEN_MATRIX(
            ENRICH_SRA_METADATA.out.csv,
            CLASSIFY_READS.out.kraken2_report
        )
    } else if (params.mode == 'accessions') {
        BUILD_POSTKRAKEN_MATRIX(
            MAKE_ACCESSIONS_METADATA.out.csv,
            CLASSIFY_READS.out.kraken2_report
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SOFTWARE VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'group1_software_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        REPORTING / MULTIQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.run_multiqc) {

        ch_multiqc_files = Channel.empty()

        if (params.run_preprocessing) {
            ch_multiqc_files = ch_multiqc_files
                .mix(PREPROCESS_READS.out.fastqc_raw)
                .mix(PREPROCESS_READS.out.fastqc_trimmed)
                .mix(PREPROCESS_READS.out.fastp_json)
                .mix(PREPROCESS_READS.out.seqkit_stats)
        }

        ch_multiqc_files = ch_multiqc_files
            .mix(CLASSIFY_READS.out.kraken2_report)

        REPORTING(
            ch_multiqc_files,
            ch_collated_versions
        )
    }

    emit:
    reads               = ch_reads_for_kraken
    kraken2_report      = CLASSIFY_READS.out.kraken2_report
    metadata_postkraken = (params.mode == 'sra' || params.mode == 'accessions') ? BUILD_POSTKRAKEN_MATRIX.out.matrix : Channel.empty()
    metadata_postkraken_hits_only = (params.mode == 'sra' || params.mode == 'accessions') ? BUILD_POSTKRAKEN_MATRIX.out.hits_only : Channel.empty()
    multiqc_report      = params.run_multiqc ? REPORTING.out.report : Channel.empty()
    versions            = ch_versions
}



    