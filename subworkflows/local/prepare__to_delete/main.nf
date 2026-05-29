
include { BCFTOOLS_REPLACE_SAMPLE_NAMES      } from '../../../modules/local/bcftools/replace_sample_names'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/local/bcftools/index'
include { BCFTOOLS_VIEW_AND_FILTER2          } from '../../../modules/local/bcftools/view_and_filter2'

workflow PREPARE {

    take:
    ch_input_vcf
    ch_all_samples

    main:
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

    bcftools_filter_1_options = "--exclude 'QUAL<${params.filter_and_enhance_vcf_qual_min} || AVG(FORMAT/GQ)<${params.filter_and_enhance_vcf_avg_gq_min} || AVG(FORMAT/DP)<${params.filter_and_enhance_vcf_avg_dp_min} || AVG(FORMAT/DP)>${params.filter_and_enhance_vcf_avg_dp_max}'"
    bcftools_filter_2_options = "--exclude 'FORMAT/GQ<${params.filter_and_enhance_vcf_sample_gq_min} | FORMAT/DP < ${params.filter_and_enhance_vcf_sample_dp_min} | FORMAT/DP > ${params.filter_and_enhance_vcf_sample_dp_max}' --set-GTs '.'"

    BCFTOOLS_VIEW_AND_FILTER2 (
        ch_vcf_with_sample_names_corrected
            .join(ch_vcf_with_sample_names_corrected_tbi, by: 0)
            .join(ch_all_samples, by: 0)
            .map { meta, vcf_file, tbi_file, samples_file -> tuple(meta, vcf_file, tbi_file, [], [], samples_file, [], []) },
        Channel.value(bcftools_filter_1_options),
        Channel.value(bcftools_filter_2_options),
        Channel.value("viewfilter2")
    )
    ch_filtered_vcf     = BCFTOOLS_VIEW_AND_FILTER2.out.vcf
    ch_filtered_vcf_tbi = BCFTOOLS_VIEW_AND_FILTER2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_AND_FILTER2.out.versions.first())
    ch_tracking = BCFTOOLS_VIEW_AND_FILTER2.out.tracking_out.first()


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    prepared_vcf     = ch_filtered_vcf
    prepared_vcf_tbi = ch_filtered_vcf_tbi
    versions         = ch_versions
    tracking         = ch_tracking
}
