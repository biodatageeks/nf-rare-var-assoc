/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FILE          } from '../modules/local/cmds/download_file'
include { MERGE_RESULTS          } from '../modules/local/cmds/merge_results'
include { RENAME                 } from '../modules/local/cmds/rename'
include { RSCRIPT_BUILDREPORTS   } from '../modules/local/rscript/buildreports'
include { REGENIE_STEP2          } from '../modules/local/regenie/step2'
include { REGENIE_STEP1          } from '../modules/local/regenie/step1'
include { RSCRIPT_VCFTOAAF       } from '../modules/local/rscript/vcf2aaf'
include { RSCRIPT_BUILD_PHENOTYPES           } from '../modules/local/rscript/build_phenotypes'
include { RSCRIPT_ASSIGN_ANNOTATIONS         } from '../modules/local/rscript/assign_annotations'
include { BGENIX                 } from '../modules/local/bgenix'
include { QCTOOL                 } from '../modules/local/qctool'
include { PLINK2_EXPORT_BGEN     } from '../modules/local/plink2/export_bgen'
include { PLINK2_WRITE_SNPLIST   } from '../modules/local/plink2/write_snplist'
include { PLINK2_MAKEBED as PLINK2_MAKEBED_1 } from '../modules/local/plink2/makebed'
include { PLINK2_MAKEBED as PLINK2_MAKEBED_2 } from '../modules/local/plink2/makebed'
include { PLINK2_MAKEBED as PLINK2_MAKEBED_3 } from '../modules/local/plink2/makebed'
include { PLINK19_MAKEBED        } from '../modules/local/plink19/makebed'
include { VEP_ANNOTATE           } from '../modules/local/vep/annotate'
include { VEP_UPDATECACHE        } from '../modules/local/vep/updatecache'
include { BCFTOOLS_VCF2FRQ       } from '../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_FILTER_QUAL_DP          } from '../modules/local/bcftools/filter_qual_dp'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_1 } from '../modules/local/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_2 } from '../modules/local/bcftools/view'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../modules/local/bcftools/index'
include { BCFTOOLS_NORM          } from '../modules/local/bcftools/norm'
include { BCFTOOLS_ANNOTATE      } from '../modules/nf-core/bcftools/annotate'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { F_COEFFICIENT_FILTERING            } from '../subworkflows/local/f_coefficient_filtering'
include { PCA                    } from '../subworkflows/local/pca'
include { REPORTING              } from '../subworkflows/local/reporting'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rare-var-assoc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process JOIN_CASES_AND_CONTROLS {
    label 'process_single'

    input:
    path(cases)
    path(controls)

    output:
    path("all.samples"), emit: output_file

    script:
    """
    cut -f1 ${cases} > all.samples
    cut -f1 ${controls} >> all.samples
    """
}

process CHECK_X_CHROM_PRESENT {
    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), env(has_x), path(bed), path(bim), path(fam), emit: has_x

    script:
    """
    # Check if chromosome X (e.g., 23) exists in the .bim file
    if grep -q "^23\\s" ${bim}; then
        echo "true" > has_x.txt
    else
        echo "false" > has_x.txt
    fi
    has_x=\$(cat has_x.txt)
    """
}

workflow RARE_VAR_ASSOC {

    take:
    ch_input_vcf // channel: vcf read in from --input
    ch_controls // channel: controls read in from --input
    ch_cases // channel: cases read in from --input

    main:

    vep_cachedir = "${projectDir}/../vep_cachedir"
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_vep_cachedir = Channel.fromPath(vep_cachedir, checkIfExists: true)
    ch_masks = Channel.fromPath(params.input_masks, checkIfExists: true)
    ch_hild = Channel.fromPath(params.hild_path, checkIfExists: true)
    ch_rename_chr = Channel.fromPath("${projectDir}/assets/rename_chr.txt", checkIfExists: true)
    ch_meta = ch_input_vcf.map { t -> t[0] }


    VEP_UPDATECACHE (
        ch_meta,
        ch_vep_cachedir,
        Channel.value(params.vep_updatecache_species),
        Channel.value(params.vep_updatecache_options),
        Channel.value(params.vep_cache_url),
        Channel.value(tuple(params.ref_fasta_url, params.vep_fasta_path))
    )
    ch_vep_cachesubdir = VEP_UPDATECACHE.out.cachesubdir.first()
    ch_versions = ch_versions.mix(VEP_UPDATECACHE.out.versions.first())
    
    JOIN_CASES_AND_CONTROLS (
        ch_cases.map { t -> t[1] },
        ch_controls.map { t -> t[1] }
    )
    ch_all_samples = JOIN_CASES_AND_CONTROLS.out.output_file

    BCFTOOLS_INDEX_1 (
        ch_input_vcf
    )
    ch_input_vcf_tbi = BCFTOOLS_INDEX_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    BCFTOOLS_VIEW_1 (
        ch_input_vcf
            .join(ch_input_vcf_tbi, by: 0)                    // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) },
        Channel.of([]),                                       // No regions file
        Channel.of([]),                                       // No targets file
        ch_all_samples,                                       // Samples file
        Channel.of([]),                                       // SNPs file
        Channel.value("--output-type z --write-index=tbi")       // input args
    )
    ch_all_samples_vcf  = BCFTOOLS_VIEW_1.out.vcf
    ch_all_samples_vcf_tbi  = BCFTOOLS_VIEW_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_1.out.versions.first())


    BCFTOOLS_FILTER_QUAL_DP (
        ch_all_samples_vcf.join(ch_all_samples_vcf_tbi, by: 0),
        Channel.value(params.bcftools_qual_filter),
        Channel.value(params.bcftools_info_filter),
        Channel.value(params.bcftools_fmt_filter),
        Channel.value(params.bcftools_info_filter_ensure_field_present),
        Channel.value(params.bcftools_fmt_filter_ensure_field_present)
    )
    ch_qual_vcf = BCFTOOLS_FILTER_QUAL_DP.out.vcf
    ch_qual_vcf_tbi = BCFTOOLS_FILTER_QUAL_DP.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_QUAL_DP.out.versions.first())
    
    BCFTOOLS_NORM (
        ch_qual_vcf
            .join(ch_qual_vcf_tbi, by: 0)  // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) },
        ch_vep_cachesubdir.map { t -> "${t}/${params.vep_fasta_path}" },
        Channel.value("norm")
    )
    ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
    ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    BCFTOOLS_ANNOTATE (
        ch_normalized_vcf
            .join(ch_normalized_vcf_tbi, by: 0)  // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file, [], []) },
        Channel.of([]),
        ch_rename_chr
    )
    ch_annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf
    ch_annotated_vcf_tbi = BCFTOOLS_ANNOTATE.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    VEP_ANNOTATE (
        ch_annotated_vcf
            .join(ch_annotated_vcf_tbi, by: 0)  // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) },
        ch_vep_cachesubdir,
        Channel.value(params.vep_annotate_species),
        Channel.value(params.vep_fasta_path),
        Channel.value(params.vep_annotate_options)
    )
    ch_vep_vcf  = VEP_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions.first())

    BCFTOOLS_INDEX_2 (
        ch_vep_vcf
    )
    ch_vep_vcf_tbi = BCFTOOLS_INDEX_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    ch_vep_vcf_with_index = ch_vep_vcf
        .join(ch_vep_vcf_tbi, by: 0)
        .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) }

    BCFTOOLS_VCF2FRQ (
        ch_vep_vcf_with_index
    )
    ch_frq = BCFTOOLS_VCF2FRQ.out.frq
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2FRQ.out.versions.first())

    //PLINK2_MAKEBED_1 (
    //    ch_vep_vcf.map { t -> tuple(t[0], [], [], [], t[1], []) },
    //    Channel.of('load_vcf'),
    //    Channel.of(params.plink2_makebed_options_1)
    //)
    //ch_bed_bim_fam_1  = PLINK2_MAKEBED_1.out.out_bed_bim_fam
    //ch_versions = ch_versions.mix(PLINK2_MAKEBED_1.out.versions.first())

    PLINK19_MAKEBED (
        ch_vep_vcf
    )
    ch_bed_bim_fam_1  = PLINK19_MAKEBED.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK19_MAKEBED.out.versions.first())

    CHECK_X_CHROM_PRESENT (
        ch_bed_bim_fam_1
    )
    ch_has_x  = CHECK_X_CHROM_PRESENT.out.has_x

    // Split based on has_x
    split_data = ch_has_x.branch {
        with_x: it[1] == "true"
        without_x: it[1] == "false"
    }

    // impute sex must be a separate step as per plink2 docs
    PLINK2_MAKEBED_2 (
        split_data.with_x
            .join(ch_frq, by: 0)
            .map { meta, has_x, bed, bim, fam, frq -> tuple(meta, bed, bim, fam, [], frq, [], []) },
        Channel.value('impute_sex'),
        Channel.value(params.plink2_makebed_options_2)
    )
    ch_bed_bim_fam_2  = PLINK2_MAKEBED_2.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK2_MAKEBED_2.out.versions.first())

    ch_bed_bim_fam_before_quality_filtering = ch_bed_bim_fam_2
        .mix(split_data.without_x.map { meta, has_x, bed, bim, fam -> tuple(meta, bed, bim, fam) })

    PLINK2_MAKEBED_3 (
        ch_bed_bim_fam_before_quality_filtering
            .map { meta, bed_file, bim_file, fam_file -> tuple(meta, bed_file, bim_file, fam_file, [], [], []) }
            .merge(ch_hild),
        Channel.value('filter_pass'),
        Channel.value(params.plink2_makebed_options_3)
    )
    ch_bed_bim_fam_3  = PLINK2_MAKEBED_3.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK2_MAKEBED_3.out.versions.first())


    F_COEFFICIENT_FILTERING (
        ch_bed_bim_fam_3
    )
    ch_bed_bim_fam_4 = F_COEFFICIENT_FILTERING.out.bed_bim_fam_out
    ch_versions = ch_versions.mix(F_COEFFICIENT_FILTERING.out.versions.first())

    PCA (
        ch_bed_bim_fam_4
    )
    ch_sscore = PCA.out.sscore
    ch_versions = ch_versions.mix(PCA.out.versions.first())

    PLINK2_WRITE_SNPLIST (
        ch_bed_bim_fam_4,
        Channel.value('writesnp_pass'),
        Channel.value(params.plink2_write_snplist_qc_options)
    )
    ch_snplist  = PLINK2_WRITE_SNPLIST.out.snplist
    ch_id  = PLINK2_WRITE_SNPLIST.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST.out.versions.first())

    BCFTOOLS_VIEW_2 (
        ch_vep_vcf_with_index,
        Channel.of([]),                     // No regions file
        Channel.of([]),                     // No targets file
        ch_id.map { t -> t[1] },            // Samples file
        ch_snplist.map { t -> t[1] },       // SNPs file
        Channel.value("--output-type v")       // input args
    )
    ch_filtered_vcf  = BCFTOOLS_VIEW_2.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_2.out.versions.first())
    

    PLINK2_EXPORT_BGEN (
        ch_bed_bim_fam_before_quality_filtering,
        Channel.value('pvcf.norm_zlib'),
        Channel.value(params.plink2_export_bgen_options)
    )
    ch_bgen  = PLINK2_EXPORT_BGEN.out.bgen
    ch_sample  = PLINK2_EXPORT_BGEN.out.sample
    ch_versions = ch_versions.mix(PLINK2_EXPORT_BGEN.out.versions.first())

    QCTOOL (
        ch_bgen,
        ch_sample,
        Channel.value('pvcf.norm'),
        Channel.value(params.qctool_options)
    )
    ch_qc_bgen  = QCTOOL.out.bgen
    ch_qc_sample  = QCTOOL.out.sample
    ch_versions = ch_versions.mix(QCTOOL.out.versions.first())

    BGENIX (
        ch_qc_bgen,
        Channel.value(params.bgenix_options)
    )
    ch_bgen_bgi  = BGENIX.out.bgen_bgi
    ch_versions = ch_versions.mix(BGENIX.out.versions.first())

    r_script_build_phenotypes_ch = Channel.fromPath(params.rscript_build_phenotypes_path, checkIfExists: true)
    RSCRIPT_BUILD_PHENOTYPES (
        r_script_build_phenotypes_ch,
        ch_bed_bim_fam_4.map { meta, bed_file, bim_file, fam_file -> tuple(meta, fam_file) }
            .join(ch_controls, by: 0)
            .join(ch_cases, by: 0),
        Channel.value(params.rscript_build_phenotypes_options)
    )
    ch_r_out_fam  = RSCRIPT_BUILD_PHENOTYPES.out.out_fam
    ch_phenotype  = RSCRIPT_BUILD_PHENOTYPES.out.phenotype
    ch_versions = ch_versions.mix(RSCRIPT_BUILD_PHENOTYPES.out.versions.first())

    renamed_file_name = ch_meta.map { t -> "${t.id}_remove_inbreeding_outliers.fam" }.first()
    RENAME (
        ch_r_out_fam,
        renamed_file_name
    )
    ch_renamed_fam  = RENAME.out.output

    ch_regenie_step_1_input = ch_bed_bim_fam_4
            .join(ch_renamed_fam, by: 0)
            .map { meta, bed_file, bim_file, fam_file, new_fam_file -> tuple(meta, bed_file, bim_file, new_fam_file) }
            .join(ch_id, by: 0)
            .join(ch_snplist, by: 0)
            .join(ch_phenotype, by: 0)
            .join(ch_sscore, by: 0)

    r_script_vcf2aaf_ch = Channel.fromPath(params.rscript_vcf2aaf_path, checkIfExists: true)
    RSCRIPT_VCFTOAAF (
        r_script_vcf2aaf_ch,
        ch_vep_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) },
        Channel.value(params.rscript_vcf2aaf_options)
    )
    ch_aaf  = RSCRIPT_VCFTOAAF.out.aaf
    ch_versions = ch_versions.mix(RSCRIPT_VCFTOAAF.out.versions.first())

    REGENIE_STEP1 (
        ch_regenie_step_1_input,
        Channel.value(params.regenie_step1_options)
    )
    ch_regenie_step1_loco  = REGENIE_STEP1.out.loco
    ch_regenie_step1_pred_list  = REGENIE_STEP1.out.pred_list
    ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions.first())

    
    r_script_annotate_ch = Channel.fromPath(params.rscript_annotate_path, checkIfExists: true)
    RSCRIPT_ASSIGN_ANNOTATIONS (
        r_script_annotate_ch,
        ch_vep_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) },
        ch_qc_sample,
        ch_masks,
        Channel.value(params.rscript_annotate_options)
    )
    ch_r_out_sample  = RSCRIPT_ASSIGN_ANNOTATIONS.out.out_sample
    ch_annotations  = RSCRIPT_ASSIGN_ANNOTATIONS.out.annotations
    ch_setlist  = RSCRIPT_ASSIGN_ANNOTATIONS.out.setlist
    ch_versions = ch_versions.mix(RSCRIPT_ASSIGN_ANNOTATIONS.out.versions.first())

    REGENIE_STEP2 (
        ch_qc_bgen
            .join(ch_r_out_sample, by: 0)
            .join(ch_phenotype, by: 0)
            .join(ch_annotations, by: 0)
            .join(ch_setlist, by: 0)
            .join(ch_aaf, by: 0)
            .join(ch_regenie_step1_pred_list, by: 0)
            .join(ch_sscore, by: 0),
        ch_masks,
        Channel.value(params.regenie_step2_options)
    )
    ch_regenie_step2_masks_bed_bim_fam  = REGENIE_STEP2.out.masks_bed_bim_fam
    ch_regenie_step2_masks_snplist  = REGENIE_STEP2.out.masks_snplist
    ch_regenie_step2_regenie_out  = REGENIE_STEP2.out.regenie_out
    ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions.first())

    r_script_buildreports_ch = Channel.fromPath(params.rscript_buildreports_path, checkIfExists: true)
    RSCRIPT_BUILDREPORTS (
        r_script_buildreports_ch,
        ch_regenie_step2_masks_snplist,
        ch_regenie_step2_regenie_out,
        ch_vep_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) },
        ch_phenotype,
        ch_annotations
    )
    ch_annotated_snps  = RSCRIPT_BUILDREPORTS.out.annotated_snps
    ch_res_log10p_1_annotated  = RSCRIPT_BUILDREPORTS.out.res_log10p_1_annotated
    ch_annotated_snps_with_sample_ids  = RSCRIPT_BUILDREPORTS.out.annotated_snps_with_sample_ids
    ch_versions = ch_versions.mix(RSCRIPT_BUILDREPORTS.out.versions.first())

    ch_regenie_step2_regenie_out.map { t -> [t[0].id, t[1]] }
        .transpose()
        .map { prefix, fl -> tuple(RegenieUtil.getPhenotypeByChunk(prefix, fl), fl) }
        .set { ch_regenie_step2_by_phenotype }

    py_script_csv_concat_ch = Channel.fromPath(params.py_script_csv_concat_path, checkIfExists: true)
    MERGE_RESULTS (
        py_script_csv_concat_ch,
        ch_regenie_step2_by_phenotype.groupTuple()
    )
    ch_results_merged = MERGE_RESULTS.out.results_merged

    REPORTING (
        ch_results_merged,
        ch_masks,
        ch_phenotype
    )


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

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
