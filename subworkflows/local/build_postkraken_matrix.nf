/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { POSTKRAKEN_MATRIX } from '../../modules/local/postkraken_matrix/main'

workflow BUILD_POSTKRAKEN_MATRIX {

    take:
    ch_sra_meta_tsv
    ch_kraken_reports

    main:
    ch_versions = Channel.empty()

    ch_meta_csv = ch_sra_meta_tsv.first()

    ch_report_files = ch_kraken_reports
        .map { meta, report -> report }
        .collect()

    POSTKRAKEN_MATRIX(
        ch_meta_csv,
        ch_report_files
    )

    ch_versions = ch_versions.mix(POSTKRAKEN_MATRIX.out.versions)

    emit:
    matrix    = POSTKRAKEN_MATRIX.out.matrix
    hits_only = POSTKRAKEN_MATRIX.out.hits_only
    versions  = ch_versions
}
