
include { BCFTOOLS_REPLACE_SAMPLE_NAMES      } from '../../../modules/local/bcftools/replace_sample_names'
include { BCFTOOLS_VCF2PSAM      } from '../../../modules/local/bcftools/vcf2psam'
include { BCFTOOLS_VCF2FRQ       } from '../../../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_TAG2TAG       } from '../../../modules/local/bcftools/tag2tag'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_1   } from '../../../modules/local/bcftools/view'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/local/bcftools/index'
include { BCFTOOLS_NORM          } from '../../../modules/local/bcftools/norm'
include { BCFTOOLS_ANNOTATE      } from '../../../modules/local/bcftools/annotate'
//include { FIX_ZERO_PL            } from '../../../modules/local/python/fix_zero_PL'
include { FILTER_AND_ENHANCE_VCF } from '../../../modules/local/combo/filter_and_enhance_vcf'
include { VIEW_AND_FILTER2_POLARSBIO         } from '../../../modules/local/python/view_and_filter2_polarsbio'
include { FILTER_AND_ENHANCE_VCF_POLARSBIO   } from '../../../modules/local/python/filter_and_enhance_vcf_polarsbio'
include { BCFTOOLS_VIEW_AND_FILTER2          } from '../../../modules/local/bcftools/view_and_filter2'

workflow PREPARE {

    take:
    ch_input_vcf
    ch_vep_cachesubdir
    ch_all_samples

    main:
    ch_rename_chr = Channel.fromPath("${projectDir}/assets/rename_chr.txt", checkIfExists: true).first()
    python_view_and_filter2_polarsbio_script_ch = Channel.fromPath(params.view_and_filter2_polarsbio_script, checkIfExists: true)
    python_filter_and_enhance_vcf_polarsbio_script_ch = Channel.fromPath(params.filter_and_enhance_vcf_polarsbio_script, checkIfExists: true)
    ch_versions = Channel.empty()
    empty_tracking_file = file("${projectDir}/assets/empty_tracking.json", checkIfExists: true)
    trackingFirstOrEmpty = { ch -> ch.first().ifEmpty(empty_tracking_file) }

    
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

    if (params.skip_preparation && params.skip_reporting) {
        // this branch is executed by nextflow_gene_assoc_tuner - preparation is done by separate pipeline and reporting is turned off to speed up tuning

        ch_filter_polarsbio_input = ch_vcf_with_sample_names_corrected
            .join(ch_vcf_with_sample_names_corrected_tbi, by: 0)
            .join(ch_all_samples, by: 0)
        
        bcftools_filter_1_options = "--exclude 'QUAL<${params.filter_and_enhance_vcf_qual_min} || AVG(FORMAT/GQ)<${params.filter_and_enhance_vcf_avg_gq_min} || AVG(FORMAT/DP)<${params.filter_and_enhance_vcf_avg_dp_min} || AVG(FORMAT/DP)>${params.filter_and_enhance_vcf_avg_dp_max}'"
        bcftools_filter_2_options = "--exclude 'FORMAT/GQ<${params.filter_and_enhance_vcf_sample_gq_min} | FORMAT/DP < ${params.filter_and_enhance_vcf_sample_dp_min} | FORMAT/DP > ${params.filter_and_enhance_vcf_sample_dp_max}' --set-GTs '.'"
        
        BCFTOOLS_VIEW_AND_FILTER2 (
            ch_filter_polarsbio_input
                .map { meta, vcf_file, tbi_file, samples_file -> tuple(meta, vcf_file, tbi_file, [], [], samples_file, [], []) },
            Channel.value(bcftools_filter_1_options),
            Channel.value(bcftools_filter_2_options),
            Channel.value("viewfilter2")
        )
        ch_filtered_vcf  = BCFTOOLS_VIEW_AND_FILTER2.out.vcf
        ch_filtered_vcf_tbi  = BCFTOOLS_VIEW_AND_FILTER2.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_AND_FILTER2.out.versions.first())
        ch_tracking = BCFTOOLS_VIEW_AND_FILTER2.out.tracking_out.first()
    } else {
        // when reporting is not skipped (and use_dosage is false) we do the filtering after generating EDA reports

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
                    //.map { meta, vcf_file, tbi_file, fasta_path -> tuple(meta, vcf_file, tbi_file, fasta_path, []) },
                    .combine(trackingFirstOrEmpty(BCFTOOLS_VIEW_1.out.tracking_out)),
                Channel.value("norm")
            )
            ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
            ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
            ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
            ch_tracking = ch_tracking.mix(BCFTOOLS_NORM.out.tracking_out.first())
            //ch_tracking = BCFTOOLS_NORM.out.tracking_out.first()

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
            //ch_tracking = Channel.empty()
        }

        if (params.use_dosage) {
            FILTER_AND_ENHANCE_VCF (
                ch_annotated_vcf.map { meta, vcf_file -> tuple(meta, vcf_file, []) },
                Channel.value(params.filter_and_enhance_vcf_qual_min),
                Channel.value(params.filter_and_enhance_vcf_avg_gq_min),
                Channel.value(params.filter_and_enhance_vcf_avg_dp_min),
                Channel.value(params.filter_and_enhance_vcf_avg_dp_max),
                Channel.value(params.filter_and_enhance_vcf_sample_gq_min),
                Channel.value(params.filter_and_enhance_vcf_sample_dp_min),
                Channel.value(params.filter_and_enhance_vcf_sample_dp_max),
                Channel.value(params.filter_and_enhance_vcf_calc_ds_min_gq),
                Channel.value("filterhance")
            )
            ch_filtered_vcf = FILTER_AND_ENHANCE_VCF.out.vcf
            ch_versions = ch_versions.mix(FILTER_AND_ENHANCE_VCF.out.versions.first())

            BCFTOOLS_INDEX_2 (
                ch_filtered_vcf
            )
            ch_filtered_vcf_tbi = BCFTOOLS_INDEX_2.out.tbi
            ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())
        } else {
            ch_filtered_vcf = ch_annotated_vcf
            ch_filtered_vcf_tbi = ch_annotated_vcf_tbi
        }
    }


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    prepared_vcf      = ch_filtered_vcf
    prepared_vcf_tbi  = ch_filtered_vcf_tbi
    versions   = ch_versions
    tracking   = ch_tracking
}

def trackingFirstOrEmpty(ch) {
    ch.first().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}

def trackingLastOrEmpty(ch) {
    ch.last().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}
