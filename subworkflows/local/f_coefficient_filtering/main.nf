include { CALCULATE_F_OUTLIERS   } from '../../../modules/local/python/calc_f_outliers'
include { PLINK2_HET             } from '../../../modules/local/plink2/het'
include { PLINK2_INDEP_PAIRWISE  } from '../../../modules/local/plink2/indep_pairwise'
include { PLINK2_MAKEBED         } from '../../../modules/local/plink2/makebed'

workflow F_COEFFICIENT_FILTERING {

    take:
    ch_bed_bim_fam
    ch_tracking_in

    main:

    ch_versions = Channel.empty()
    
    PLINK2_INDEP_PAIRWISE (
        ch_bed_bim_fam.map { meta, bed_file, bim_file, fam_file -> tuple(meta, bed_file, bim_file, fam_file, []) },
        Channel.value(params.plink2_indep_pairwise_window),
        Channel.value('indep_pairwise'),
        Channel.value(params.plink2_indep_pairwise_options),
        ch_tracking_in.first()
    )
    ch_indep_pairwise_prune_in = PLINK2_INDEP_PAIRWISE.out.out_prune_in
    ch_indep_pairwise_prune_out = PLINK2_INDEP_PAIRWISE.out.out_prune_out
    ch_versions = ch_versions.mix(PLINK2_INDEP_PAIRWISE.out.versions.first())
    ch_tracking = PLINK2_INDEP_PAIRWISE.out.tracking_out.first()

    PLINK2_HET (
        ch_bed_bim_fam
            .join(ch_indep_pairwise_prune_in, by: 0)
            .map { meta, bed_file, bim_file, fam_file, het_file -> tuple(meta, bed_file, bim_file, fam_file, het_file) },
        Channel.value('het'),
        Channel.value(''),
        PLINK2_INDEP_PAIRWISE.out.tracking_out.first()
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

    PLINK2_MAKEBED (
        ch_bed_bim_fam
            .join(ch_inbreeding_outliers, by: 0)
            .map { meta, bed_file, bim_file, fam_file, outliers_file -> tuple(meta, bed_file, bim_file, fam_file, [], [], outliers_file, []) },
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value('remove_inbreeding_outliers'),
        Channel.value(''),
        PLINK2_HET.out.tracking_out.first()
    )
    ch_bed_bim_fam_out  = PLINK2_MAKEBED.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK2_MAKEBED.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEBED.out.tracking_out.first())


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    bed_bim_fam_out = ch_bed_bim_fam_out
    versions        = ch_versions
    tracking        = ch_tracking
}