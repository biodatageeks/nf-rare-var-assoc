include { CALCULATE_F_OUTLIERS   } from '../../../modules/local/python/calc_f_outliers'
include { PLINK2_HET             } from '../../../modules/local/plink2/het'
include { PLINK2_INDEP_PAIRWISE  } from '../../../modules/local/plink2/indep_pairwise'
include { PLINK2_MAKEPGEN        } from '../../../modules/local/plink2/makepgen'

workflow F_COEFFICIENT_FILTERING {

    take:
    ch_pgen_pvar_psam
    ch_tracking_in

    main:

    ch_versions = Channel.empty()
    empty_tracking_file = file("${projectDir}/assets/empty_tracking.json", checkIfExists: true)
    trackingFirstOrEmpty = { ch -> ch.first().ifEmpty(empty_tracking_file) }
    // Tracking is auxiliary: keep a single record and fallback when upstream optional branches produce none.
    ch_tracking_single = trackingFirstOrEmpty(ch_tracking_in)
    
    PLINK2_INDEP_PAIRWISE (
        ch_pgen_pvar_psam.map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, []) }
            .combine(ch_tracking_single),
        Channel.value(params.plink2_indep_pairwise_window),
        Channel.value('indep_pairwise'),
        Channel.value(params.plink2_indep_pairwise_options)
    )
    ch_indep_pairwise_prune_in = PLINK2_INDEP_PAIRWISE.out.out_prune_in
    ch_indep_pairwise_prune_out = PLINK2_INDEP_PAIRWISE.out.out_prune_out
    ch_versions = ch_versions.mix(PLINK2_INDEP_PAIRWISE.out.versions.first())
    ch_tracking = PLINK2_INDEP_PAIRWISE.out.tracking_out.first()

    PLINK2_HET (
        ch_pgen_pvar_psam
            .join(ch_indep_pairwise_prune_in, by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, het_file -> tuple(meta, pgen_file, pvar_file, psam_file, het_file) }
            .combine(trackingFirstOrEmpty(PLINK2_INDEP_PAIRWISE.out.tracking_out)),
        Channel.value('het'),
        Channel.value('')
    )
    ch_het  = PLINK2_HET.out.out_het
    ch_versions = ch_versions.mix(PLINK2_HET.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_HET.out.tracking_out.first())

    CALCULATE_F_OUTLIERS (
        ch_het,
        Channel.value(params.inbreeding_outliers_range_stds),
        Channel.value('inbreeding_outliers')
    )
    ch_inbreeding_outliers  = CALCULATE_F_OUTLIERS.out.outliers
    ch_versions = ch_versions.mix(CALCULATE_F_OUTLIERS.out.versions.first())

    PLINK2_MAKEPGEN (
        ch_pgen_pvar_psam
            .join(ch_inbreeding_outliers, by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, outliers_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], outliers_file, []) }
            .combine(trackingFirstOrEmpty(PLINK2_HET.out.tracking_out)),
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('remove_inbreeding_outliers'),
        Channel.value('')
    )
    ch_pgen_pvar_psam_out  = PLINK2_MAKEPGEN.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN.out.tracking_out.first())


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    pgen_pvar_psam_out = ch_pgen_pvar_psam_out
    versions        = ch_versions
    tracking        = ch_tracking
}


def trackingFirstOrEmpty(ch) {
    ch.first().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}

def trackingLastOrEmpty(ch) {
    ch.last().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}
