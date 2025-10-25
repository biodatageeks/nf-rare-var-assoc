include { PLINK19_MAKESET           } from '../../../modules/local/plink19/makeset'
include { PLINK2_MAKEBED            } from '../../../modules/local/plink2/makebed'
include { PLINK2_INDEP_PAIRWISE     } from '../../../modules/local/plink2/indep_pairwise'
include { PLINK2_KING_CUTOFF        } from '../../../modules/local/plink2/king_cutoff'
include { PLINK2_PCA                } from '../../../modules/local/plink2/pca'
include { PLINK2_PROJECTION_SCORE   } from '../../../modules/local/plink2/projection_score'
include { DRAW_PC_PLOT              } from '../../../modules/local/python/draw_pc_plot'

workflow PCA {

    take:
    ch_pgen_pvar_psam
    ch_pheno
    ch_frq
    ch_tracking_in

    main:

    ch_hild = Channel.fromPath(params.hild_path, checkIfExists: true)
    ch_versions = Channel.empty()


    PLINK2_MAKEBED (
        ch_pgen_pvar_psam
            .map { meta, pgen_file, pvar_file, psam_file ->
                tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [])
            }
            .combine(ch_tracking_in.first()),
        Channel.value(''),
        Channel.value(''),
        Channel.value('makebed'),
        Channel.value('')
    )
    ch_bed_bim_fam = PLINK2_MAKEBED.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK2_MAKEBED.out.versions.first())
    ch_tracking = PLINK2_MAKEBED.out.tracking_out.first()

    PLINK19_MAKESET (
        ch_bed_bim_fam
            .combine(ch_hild)
            .combine(PLINK2_MAKEBED.out.tracking_out.first())
    )
    ch_hildset = PLINK19_MAKESET.out.out_set
    ch_versions = ch_versions.mix(PLINK19_MAKESET.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK19_MAKESET.out.tracking_out.first())

    PLINK2_INDEP_PAIRWISE (
        ch_pgen_pvar_psam
            .join(ch_hildset, by: 0)
            .combine(PLINK19_MAKESET.out.tracking_out.first()),
        Channel.value(params.plink2_indep_pairwise_window_pca),
        Channel.value('indep_pairwise'),
        Channel.value(params.plink2_indep_pairwise_options)
    )
    ch_indep_pairwise_prune_in = PLINK2_INDEP_PAIRWISE.out.out_prune_in
    ch_indep_pairwise_prune_out = PLINK2_INDEP_PAIRWISE.out.out_prune_out
    ch_versions = ch_versions.mix(PLINK2_INDEP_PAIRWISE.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_INDEP_PAIRWISE.out.tracking_out.first())

    PLINK2_KING_CUTOFF (
        ch_pgen_pvar_psam
            .join(ch_indep_pairwise_prune_in, by: 0)
            .combine(PLINK2_INDEP_PAIRWISE.out.tracking_out.first()),
        Channel.value(params.plink2_king_cutoff_threshold_pca),
        Channel.value('king_cutoff'),
        Channel.value(params.plink2_king_cutoff_options)
    )
    ch_king_cutoff_prune_in = PLINK2_KING_CUTOFF.out.out_prune_in
    ch_king_cutoff_prune_out = PLINK2_KING_CUTOFF.out.out_prune_out
    ch_versions = ch_versions.mix(PLINK2_KING_CUTOFF.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_KING_CUTOFF.out.tracking_out.first())

    PLINK2_PCA (
        ch_pgen_pvar_psam
            .join(ch_indep_pairwise_prune_in, by: 0)
            .join(ch_king_cutoff_prune_in, by: 0)
            .join(ch_frq, by: 0)
            .combine(PLINK2_KING_CUTOFF.out.tracking_out.first()),
        Channel.value(params.plink2_pca_settings),
        Channel.value('pca'),
        Channel.value(params.plink2_pca_options)
    )
    // ch_acount = PLINK2_PCA.out.acount
    ch_eigenval = PLINK2_PCA.out.eigenval
    ch_eigenvec = PLINK2_PCA.out.eigenvec
    ch_eigenvec_allele = PLINK2_PCA.out.eigenvec_allele
    ch_versions = ch_versions.mix(PLINK2_PCA.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_PCA.out.tracking_out.first())

    PLINK2_PROJECTION_SCORE (
        ch_pgen_pvar_psam
            .join(ch_frq, by: 0)
            .join(ch_eigenvec_allele, by: 0)
            .combine(PLINK2_PCA.out.tracking_out.first()),
        Channel.value(params.plink2_projection_score_settings),
        Channel.value('projection'),
        Channel.value(params.plink2_projection_score_options)
    )
    ch_sscore = PLINK2_PROJECTION_SCORE.out.sscore
    ch_versions = ch_versions.mix(PLINK2_PROJECTION_SCORE.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_PROJECTION_SCORE.out.tracking_out.first())

    if (params.skip_reporting == 'false') {
        DRAW_PC_PLOT (
            ch_sscore.join(ch_pheno, by: 0),
            Channel.value('pc_plot')
        )
        ch_plot_file = DRAW_PC_PLOT.out.plot_file
        ch_versions = ch_versions.mix(DRAW_PC_PLOT.out.versions.first())
    } else {
        ch_plot_file = Channel.empty()
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    sscore     = ch_sscore
    plot_file  = ch_plot_file
    versions   = ch_versions
    tracking   = ch_tracking
}