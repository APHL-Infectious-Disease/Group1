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

    SRA_META(
         tuple('group1_esearch','"WGS[Strategy] AND USA AND Wastewater AND Metagenome"'
        ),
        "sra"
    )
    // ch_versions = ch_versions.mix(SRA_META.out.versions_esearch)
    // ch_meta = SRA_META.out.tsv
    // ch_xml = SRA_META.out.xml
    
    // Read TSV and convert to a channel from sra_meta_top3.tsv, which contains the top 3 SRA run IDs
    
    tsvChannel = SRA_META.out.sra_meta_top
    .flatMap { file ->
        file.readLines().collect { line ->
            def cols = line.split('\t')
            cols        }
    }
     ch_fastqdl = tuple([meta: 'SRX26273713', id: 'SRX26273713'],'SRX26273713') // example SRA run ID
    
  
    FASTQDL(
        tsvChannel.map { cols -> tuple([meta: cols[0], id: cols[0]], cols[0]) }, // meta and id are the same in this case
    )
    ch_fastq = FASTQDL.out.fastq
    ch_kraken_db = Channel.fromPath('/workspaces/Group1/assets/kraken2db_v2/')
    .collect()

  
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
