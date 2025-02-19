/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RSCRIPT_BUILDREPORTS   } from '../modules/local/rscript/buildreports'
include { REGENIE_STEP2          } from '../modules/local/regenie/step2'
include { REGENIE_STEP1          } from '../modules/local/regenie/step1'
include { RSCRIPT_VCFTOAAF       } from '../modules/local/rscript/vcf2aaf'
include { RSCRIPT_ANNOTATE       } from '../modules/local/rscript/annotate'
include { BGENIX                 } from '../modules/local/bgenix'
include { QCTOOL                 } from '../modules/local/qctool'
include { PLINK2_EXPORT_BGEN     } from '../modules/local/plink2/export_bgen'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_1 } from '../modules/local/plink2/write_snplist'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_2 } from '../modules/local/plink2/write_snplist'
include { PLINK19_MAKEBED        } from '../modules/local/plink19/makebed'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rare-var-assoc-nf_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RARE_VAR_ASSOC_NF {

    take:
    ch_vcf // channel: vcf read in from --input
    ch_controls // channel: controls read in from --input
    ch_cases // channel: cases read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    PLINK19_MAKEBED (
        ch_vcf
    )
    ch_bed  = PLINK19_MAKEBED.out.bed
    ch_bim  = PLINK19_MAKEBED.out.bim
    ch_fam  = PLINK19_MAKEBED.out.fam
    ch_nosex  = PLINK19_MAKEBED.out.nosex
    ch_versions = ch_versions.mix(PLINK19_MAKEBED.out.versions.first())

    PLINK2_WRITE_SNPLIST_1 (
        ch_bed,
        ch_bim,
        ch_fam,
        Channel.of('snps_pass'),
        Channel.of(params.plink2_write_snplist_options)
    )
    ch_snplist_1  = PLINK2_WRITE_SNPLIST_1.out.snplist
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_1.out.versions.first())

    PLINK2_WRITE_SNPLIST_2 (
        ch_bed,
        ch_bim,
        ch_fam,
        Channel.of('qc_pass'),
        Channel.of(params.plink2_write_snplist_qc_options)
    )
    ch_snplist_2  = PLINK2_WRITE_SNPLIST_2.out.snplist
    ch_id_2  = PLINK2_WRITE_SNPLIST_2.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_2.out.versions.first())

    PLINK2_EXPORT_BGEN (
        ch_bed,
        ch_bim,
        ch_fam,
        Channel.of('pvcf.norm_zlib'),
        Channel.of(params.plink2_export_bgen_options)
    )
    ch_bgen  = PLINK2_EXPORT_BGEN.out.bgen
    ch_sample  = PLINK2_EXPORT_BGEN.out.sample
    ch_versions = ch_versions.mix(PLINK2_EXPORT_BGEN.out.versions.first())

    QCTOOL (
        ch_bgen,
        ch_sample,
        Channel.of('pvcf.norm'),
        Channel.of(params.qctool_options)
    )
    ch_qc_bgen  = QCTOOL.out.bgen
    ch_qc_sample  = QCTOOL.out.sample
    ch_versions = ch_versions.mix(QCTOOL.out.versions.first())

    BGENIX (
        ch_qc_bgen,
        Channel.of(params.bgenix_options)
    )
    ch_bgen_bgi  = BGENIX.out.bgen_bgi
    ch_versions = ch_versions.mix(BGENIX.out.versions.first())

    r_script_annotate_ch = Channel.fromPath(params.rscript_annotate_path, checkIfExists: true)
    // r_script_annotate_ch = Channel.fromPath("${projectDir}/modules/local/rscript/annotate/assets/test.R", checkIfExists: true)
    RSCRIPT_ANNOTATE (
        r_script_annotate_ch,
        ch_vcf,
        ch_bed,
        ch_bim,
        ch_fam,
        ch_qc_bgen,
        ch_bgen_bgi,
        ch_qc_sample,
        ch_controls,
        ch_cases,
        Channel.of(params.rscript_annotate_options)
    )
    ch_r_out_fam  = RSCRIPT_ANNOTATE.out.out_fam
    ch_r_out_sample  = RSCRIPT_ANNOTATE.out.out_sample
    ch_phenotype  = RSCRIPT_ANNOTATE.out.phenotype
    ch_annotations  = RSCRIPT_ANNOTATE.out.annotations
    ch_masks  = RSCRIPT_ANNOTATE.out.masks
    ch_setlist  = RSCRIPT_ANNOTATE.out.setlist
    ch_versions = ch_versions.mix(RSCRIPT_ANNOTATE.out.versions.first())

    r_script_vcf2aaf_ch = Channel.fromPath(params.rscript_vcf2aaf_path, checkIfExists: true)
    RSCRIPT_VCFTOAAF (
        r_script_vcf2aaf_ch,
        ch_vcf,
        Channel.of(params.rscript_vcf2aaf_options)
    )
    ch_aaf  = RSCRIPT_VCFTOAAF.out.aaf
    ch_versions = ch_versions.mix(RSCRIPT_VCFTOAAF.out.versions.first())

    REGENIE_STEP1 (
        ch_bed,
        ch_bim,
        ch_r_out_fam,
        ch_id_2,
        ch_snplist_2,
        ch_phenotype,
        Channel.of(params.regenie_step1_options)
    )
    ch_regenie_step1_loco  = REGENIE_STEP1.out.loco
    ch_regenie_step1_pred_list  = REGENIE_STEP1.out.pred_list
    ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions.first())

    REGENIE_STEP2 (
        ch_qc_bgen,
        ch_r_out_sample,
        ch_phenotype,
        ch_annotations,
        ch_setlist,
        ch_masks,
        ch_aaf,
        ch_regenie_step1_pred_list,
        Channel.of(params.regenie_step2_options)
    )
    ch_regenie_step2_masks_bed  = REGENIE_STEP2.out.masks_bed
    ch_regenie_step2_masks_bim  = REGENIE_STEP2.out.masks_bim
    ch_regenie_step2_masks_fam  = REGENIE_STEP2.out.masks_fam
    ch_regenie_step2_masks_snplist  = REGENIE_STEP2.out.masks_snplist
    ch_regenie_step2_y1_regenie  = REGENIE_STEP2.out.y1_regenie
    ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions.first())

    r_script_buildreports_ch = Channel.fromPath(params.rscript_buildreports_path, checkIfExists: true)
    RSCRIPT_BUILDREPORTS (
        r_script_buildreports_ch,
        ch_regenie_step2_masks_snplist,
        ch_regenie_step2_y1_regenie,
        ch_vcf,
        ch_phenotype,
        ch_annotations
    )
    ch_annotated_snps  = RSCRIPT_BUILDREPORTS.out.annotated_snps
    ch_res_log10p_1_annotated  = RSCRIPT_BUILDREPORTS.out.res_log10p_1_annotated
    ch_annotated_snps_with_sample_ids  = RSCRIPT_BUILDREPORTS.out.annotated_snps_with_sample_ids
    ch_versions = ch_versions.mix(RSCRIPT_BUILDREPORTS.out.versions.first())


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'rare-var-assoc-nf_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
