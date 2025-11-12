
include { BCFTOOLS_REPLACE_SAMPLE_NAMES      } from '../../../modules/local/bcftools/replace_sample_names'
include { BCFTOOLS_VCF2PSAM      } from '../../../modules/local/bcftools/vcf2psam'
include { BCFTOOLS_VCF2FRQ       } from '../../../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_TAG2TAG       } from '../../../modules/local/bcftools/tag2tag'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_1   } from '../../../modules/local/bcftools/view'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/local/bcftools/index'
include { BCFTOOLS_NORM          } from '../../../modules/local/bcftools/norm'
include { BCFTOOLS_ANNOTATE      } from '../../../modules/local/bcftools/annotate'
include { FIX_ZERO_PL            } from '../../../modules/local/python/fix_zero_PL'

workflow PREPARE {

    take:
    ch_input_vcf
    ch_vep_cachesubdir
    ch_all_samples

    main:
    ch_rename_chr = Channel.fromPath("${projectDir}/assets/rename_chr.txt", checkIfExists: true)
    ch_versions = Channel.empty()

    
    BCFTOOLS_REPLACE_SAMPLE_NAMES (
        ch_input_vcf,
        Channel.value(params.bcftools_replace_sample_names_sed_arg),
        Channel.value('replace_sample_names')
    )
    ch_vcf_with_sample_names_corrected = BCFTOOLS_REPLACE_SAMPLE_NAMES.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REPLACE_SAMPLE_NAMES.out.versions.first())

    BCFTOOLS_INDEX_1 (
        ch_vcf_with_sample_names_corrected
    )
    ch_vcf_with_sample_names_corrected_tbi = BCFTOOLS_INDEX_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())


    ch_view_1_input = ch_vcf_with_sample_names_corrected
        .join(ch_vcf_with_sample_names_corrected_tbi, by: 0)
        .join(ch_all_samples, by: 0)
        .map { meta, vcf_file, tbi_file, samples_file -> tuple(meta, vcf_file, tbi_file, [], [], samples_file, [], []) }

    BCFTOOLS_VIEW_1 (
        ch_view_1_input,
        Channel.value(params.bcftools_view_1_options),
        Channel.value("view1")
    )
    ch_all_samples_vcf  = BCFTOOLS_VIEW_1.out.vcf
    ch_all_samples_vcf_tbi  = BCFTOOLS_VIEW_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_1.out.versions.first())
    ch_tracking = BCFTOOLS_VIEW_1.out.tracking_out

    if (params.skip_preparation == false) {
        // NORM must be before ANNOTATE (where we assign variant ids) because we must first split multiallelic sites before assigning variant ids.
        // Otherwise we'll end up with duplicated ids with a comma within them and this will cause subsequent plink write-snplist steps to fail
        BCFTOOLS_NORM (
            ch_all_samples_vcf
                .join(ch_all_samples_vcf_tbi, by: 0)
                .combine(ch_vep_cachesubdir.map { t -> "${t}/${params.vep_fasta_path}" })
                .combine(BCFTOOLS_VIEW_1.out.tracking_out.first()),
            Channel.value("norm")
        )
        ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
        ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
        ch_tracking = ch_tracking.mix(BCFTOOLS_NORM.out.tracking_out.first())

        BCFTOOLS_ANNOTATE (
            ch_normalized_vcf
                .join(ch_normalized_vcf_tbi, by: 0)
                .combine(ch_rename_chr)
                .map { meta, vcf_file, tbi_file, rename_chr -> tuple(meta, vcf_file, tbi_file, [], [], [], rename_chr) },
            Channel.value("rename_chr")
        )
        ch_annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf
        ch_annotated_vcf_tbi = BCFTOOLS_ANNOTATE.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    } else {
        ch_annotated_vcf = ch_all_samples_vcf
        ch_annotated_vcf_tbi = ch_all_samples_vcf_tbi
    }

    if (params.use_dosage == 'true') {
        FIX_ZERO_PL (
            ch_annotated_vcf,
            Channel.value(params.fix_zero_pl_min_gq),
            Channel.value('calc_dosage')
        )
        ch_vcf_with_pl_corrected = FIX_ZERO_PL.out.vcf
        ch_versions = ch_versions.mix(FIX_ZERO_PL.out.versions.first())

        BCFTOOLS_INDEX_2 (
            ch_vcf_with_pl_corrected
        )
        ch_vcf_with_pl_corrected_tbi = BCFTOOLS_INDEX_2.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())
    } else {
        ch_vcf_with_pl_corrected = ch_annotated_vcf
        ch_vcf_with_pl_corrected_tbi = ch_annotated_vcf_tbi
    }


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    all_samples_vcf      = ch_vcf_with_pl_corrected
    all_samples_vcf_tbi  = ch_vcf_with_pl_corrected_tbi
    versions   = ch_versions
    tracking   = ch_tracking
}