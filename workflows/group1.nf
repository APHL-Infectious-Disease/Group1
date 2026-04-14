/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_group1_pipeline'

// Retrieve metadata for SRA runs and get data
include { FASTQDL } from '../modules/nf-core/fastqdl/main'                                                              
                                                            
include {SRA_META} from '../modules/local/sra_meta_pull'
include {ENTREZDIRECT_ESEARCH} from '../modules/nf-core/entrezdirect/esearch'

                                                       
 include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'   
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GROUP1 {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    
    ch_versions = channel.empty()

    ENTREZDIRECT_ESEARCH(
        tuple('group1_esearch','"WGS[Strategy] AND USA AND Wastewater AND Metagenome"'
        ),
        "sra"
    )

    SRA_META(
        true
    )
    ch_versions = ch_versions.mix(SRA_META.out.versions_esearch)
    // ch_meta = SRA_META.out.tsv
    // ch_xml = SRA_META.out.xml
    
    
    ch_fastqdl = tuple([meta: 'SRX26273713', id: 'SRX26273713'],'SRX26273713') // example SRA run ID, replace with actual IDs as needed
    
    FASTQDL(
        ch_fastqdl
    )
    ch_fastq = FASTQDL.out.fastq
    ch_kraken_db = Channel.fromPath('/workspaces/Group1/assets/kraken2db/')
    KRAKEN2_KRAKEN2(
        ch_fastq,
        ch_kraken_db,
        false,
        true
    )
   

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'group1_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    // meta = ch_meta
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
