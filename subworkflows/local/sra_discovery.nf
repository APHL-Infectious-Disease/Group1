/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ENTREZDIRECT_ESEARCH } from '../../modules/nf-core/entrezdirect/esearch'
include { ENRICH_SRA_METADATA }  from '../../modules/local/enrich_sra_metadata/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBWORKFLOW: DISCOVER SRA RUNS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SRA_DISCOVERY {

    main:

    ch_versions = Channel.empty()

    /*
    Query SRA broadly. The resulting RunInfo CSV may contain many candidate rows.
    USA filtering and max_runs limiting happen after BioSample enrichment.
    */
    ENTREZDIRECT_ESEARCH(
        tuple('group1_esearch', params.sra_query),
        "sra"
    )

    ch_versions = ch_versions.mix(ENTREZDIRECT_ESEARCH.out.versions_esearch)

    ch_runinfo_checked = ENTREZDIRECT_ESEARCH.out.runinfo
        .map { meta, csv ->
            if (!csv.exists() || csv.size() == 0) {
                error "ENTREZDIRECT_ESEARCH produced an empty runinfo CSV. Check params.sra_query: ${params.sra_query}"
            }
            csv
        }

    /*
    Enrich with BioSample metadata, keep USA samples, then apply max_runs.
    This makes --max_runs mean usable USA runs, not raw SRA search hits.
    */
    ENRICH_SRA_METADATA(
        ch_runinfo_checked
    )

    ch_accessions = ENRICH_SRA_METADATA.out.csv
        .splitCsv(header: true)
        .filter { row ->
            row['SRA accession'] && row['SRA accession'].toString().trim()
        }
        .map { row ->
            def accession = row['SRA accession'].toString().trim()

            def layout = ''
            if (row.containsKey('library_layout') && row['library_layout']) {
                layout = row['library_layout'].toString().toUpperCase()
            }

            def meta = [
                id         : accession,
                sample     : accession,
                single_end : (layout != 'PAIRED')
            ]

            tuple(meta, accession)
        }
        .ifEmpty {
            error "No USA SRA accessions were available after BioSample enrichment/filtering."
        }

    emit:
    reads     = ch_accessions
    tsv       = ENRICH_SRA_METADATA.out.csv
    xml       = Channel.empty()
    versions  = ch_versions
}
