/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERGE_RESULTS          } from '../modules/local/cmds/merge_results'
include { RENAME as RENAME_1     } from '../modules/local/cmds/rename'
include { RENAME as RENAME_2     } from '../modules/local/cmds/rename'
include { CHECK_X_CHROM_PRESENT  } from '../modules/local/cmds/check_x_chrom_present'
include { GENERATE_TRACKING_REPORT           } from '../modules/local/python/generate_tracking_report'
include { EXPLORATORY_DATA_ANALYSIS          } from '../modules/local/python/eda'
include { MERGE_SEX_COVAR        } from '../modules/local/python/merge_sex_covar'
include { RSCRIPT_BUILDREPORTS   } from '../modules/local/rscript/buildreports'
include { REGENIE_STEP2          } from '../modules/local/regenie/step2'
include { REGENIE_STEP1          } from '../modules/local/regenie/step1'
include { PYTHON_VCFTOAAF        } from '../modules/local/python/vcf2aaf'
include { BCFTOOLS_ASSIGN_ANNOTATIONS        } from '../modules/local/bcftools/assign_annotations'
include { BCFTOOLS_ASSIGN_ANNOTATIONS as BCFTOOLS_ASSIGN_ANNOTATIONS_2 } from '../modules/local/bcftools/assign_annotations'
include { PLINK2_IMPORT_DOSAGE   } from '../modules/local/plink2/import_dosage'
include { PLINK2_EXPORT_OTHER    } from '../modules/local/plink2/export_other'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_1 } from '../modules/local/plink2/write_snplist'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_2 } from '../modules/local/plink2/write_snplist'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_1 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_2 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_3 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_4 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_5 } from '../modules/local/plink2/makepgen'
include { PREPARE_VCF            } from '../modules/local/nextflow_run/prepare_vcf'
include { BCFTOOLS_VCF2PSAM      } from '../modules/local/bcftools/vcf2psam'
include { BCFTOOLS_VCF2FRQ       } from '../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_INDEX         } from '../modules/local/bcftools/index'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_2   } from '../modules/local/bcftools/view'
include { BCFTOOLS_REPLACE_SAMPLE_NAMES      } from '../modules/local/bcftools/replace_sample_names'
include { BCFTOOLS_VIEW_AND_FILTER2          } from '../modules/local/bcftools/view_and_filter2'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { FILTER_MISSING_PER_PHENO           } from '../subworkflows/local/filter_missing_per_pheno'
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

workflow RARE_VAR_ASSOC {

    take:
    ch_input_vcf // channel: vcf read in from --input
    ch_phenotype
    ch_all_samples

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_masks = Channel.fromPath(params.input_masks, checkIfExists: true).first()
    ch_meta = ch_input_vcf.map { t -> t[0] }


    // call BCFTOOLS_ASSIGN_ANNOTATIONS with a dummy python script to pull the bioinf_combo image before HyperQueue provisions workers on PLGrid
    dummy_python_script_ch = Channel.fromPath("${projectDir}/modules/local/bcftools/assign_annotations/assets/dummy.py", checkIfExists: true).first()
    BCFTOOLS_ASSIGN_ANNOTATIONS_2 (
        ch_input_vcf
            .combine(ch_masks)
            .combine(dummy_python_script_ch),
        Channel.value(params.rscript_annotate_options)
    )

    if (params.skip_preparation == false) {
        ch_prep_params_file = Channel.fromPath("${projectDir}/conf/nf_prepare_params.yml", checkIfExists: true).first()

        PREPARE_VCF (
            ch_input_vcf,
            ch_prep_params_file
        )
        ch_prepared_vcf = PREPARE_VCF.out.prepared_vcf
        ch_versions = ch_versions.mix(PREPARE_VCF.out.versions)
        ch_tracking = PREPARE_VCF.out.tracking.first()
    } else {
        ch_prepared_vcf = ch_input_vcf
    }


    BCFTOOLS_REPLACE_SAMPLE_NAMES (
        ch_prepared_vcf,
        Channel.value(params.bcftools_replace_sample_names_sed_arg),
        Channel.value('replace_sample_names')
    )
    ch_vcf_with_sample_names_corrected = BCFTOOLS_REPLACE_SAMPLE_NAMES.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REPLACE_SAMPLE_NAMES.out.versions.first())


    BCFTOOLS_INDEX (
        ch_vcf_with_sample_names_corrected
    )
    ch_vcf_with_sample_names_corrected_tbi = BCFTOOLS_INDEX.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    if (!params.skip_reporting) {
        eda_script_ch = Channel.fromPath("${projectDir}/modules/local/python/eda/assets/eda_v6.py", checkIfExists: true).first()
        EXPLORATORY_DATA_ANALYSIS (
            ch_vcf_with_sample_names_corrected
                .join(ch_vcf_with_sample_names_corrected_tbi, by: 0)
                .join(ch_phenotype, by: 0)
                .combine(eda_script_ch),
            Channel.value(params.use_dosage)
        )
        ch_eda_plots = EXPLORATORY_DATA_ANALYSIS.out.plots
        ch_versions = ch_versions.mix(EXPLORATORY_DATA_ANALYSIS.out.versions.first())
    } else {
        ch_eda_plots = Channel.empty()
    }


    bcftools_filter_1_options = "--exclude 'QUAL<${params.filter_vcf_qual_min} || AVG(FORMAT/GQ)<${params.filter_vcf_avg_gq_min} || AVG(FORMAT/DP)<${params.filter_vcf_avg_dp_min} || AVG(FORMAT/DP)>${params.filter_vcf_avg_dp_max}'"
    bcftools_filter_2_options = "--exclude 'FORMAT/GQ<${params.filter_vcf_sample_gq_min} | FORMAT/DP < ${params.filter_vcf_sample_dp_min} | FORMAT/DP > ${params.filter_vcf_sample_dp_max}' --set-GTs '.'"
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
    if (params.skip_preparation) {
        ch_tracking = BCFTOOLS_VIEW_AND_FILTER2.out.tracking_out.first()
    } else {
        ch_tracking = ch_tracking.mix(BCFTOOLS_VIEW_AND_FILTER2.out.tracking_out.first())
    }
    ch_filtered_vcf_with_index = ch_filtered_vcf.join(ch_filtered_vcf_tbi, by: 0)


    BCFTOOLS_VCF2FRQ (
        ch_filtered_vcf_with_index
    )
    ch_frq = BCFTOOLS_VCF2FRQ.out.frq
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2FRQ.out.versions.first())

    BCFTOOLS_VCF2PSAM (
        ch_filtered_vcf_with_index
    )
    ch_unk_sex_psam = BCFTOOLS_VCF2PSAM.out.psam
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2PSAM.out.versions.first())

    PLINK2_MAKEPGEN_1 (
        ch_filtered_vcf_with_index
            .join(ch_unk_sex_psam, by: 0)
            .map { meta, vcf_file, tbi_file, unk_sex_psam_file -> tuple(meta, [], [], unk_sex_psam_file, vcf_file, tbi_file, [], [], [],   []) },  // no upstream tracking JSON to thread here: VIEW_AND_FILTER2 is the predecessor and its tracking is mixed into ch_tracking separately
        Channel.value(''),
        Channel.value(''),
        Channel.value(params.plink2_makepgen_1_vcf_input_options),
        Channel.value('plink2_makepgen_1'),
        Channel.value(params.plink2_makepgen_1_options)
    )
    ch_pgen_pvar_psam_1  = PLINK2_MAKEPGEN_1.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_1.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN_1.out.tracking_out.first())

    CHECK_X_CHROM_PRESENT (
        ch_pgen_pvar_psam_1
    )
    ch_has_x  = CHECK_X_CHROM_PRESENT.out.has_x

    // Split based on has_x
    split_data = ch_has_x.branch {
        with_x: it[1] == "true"
        without_x: it[1] == "false"
    }

    // impute sex must be a separate step as per plink2 docs
    PLINK2_MAKEPGEN_2 (
        split_data.with_x
            .join(ch_frq, by: 0)
            .map { meta, has_x, pgen, pvar, psam, frq -> tuple(meta, pgen, pvar, psam, [], [], frq, [], []) }
            .combine(trackingFirstOrEmpty(PLINK2_MAKEPGEN_1.out.tracking_out)),
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('plink2_makepgen_2_impute_sex'),
        Channel.value(params.plink2_makepgen_2_options)
    )
    ch_pgen_pvar_psam_2  = PLINK2_MAKEPGEN_2.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_2.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN_2.out.tracking_out.first())

    renamed_file_name = ch_pgen_pvar_psam_2.map { meta, pgen, pvar, psam -> tuple(meta, "${meta.id}_plink2_makepgen_1.psam") }
    RENAME_1 (
        ch_pgen_pvar_psam_2.map { meta, pgen, pvar, psam -> tuple(meta, psam) }
            .join(renamed_file_name, by: 0)
    )
    ch_renamed_psam  = RENAME_1.out.output

    ch_pgen_pvar_psam_before_quality_filtering = split_data.with_x
        .join(ch_renamed_psam, by: 0)
        .map { meta, has_x, pgen, pvar, psam, psam_imputesex -> tuple(meta, pgen, pvar, psam_imputesex) }
        .mix(split_data.without_x.map { meta, has_x, pgen, pvar, psam -> tuple(meta, pgen, pvar, psam) })

    FILTER_MISSING_PER_PHENO (
        ch_pgen_pvar_psam_before_quality_filtering,
        ch_phenotype,
        trackingFirstOrFallback(PLINK2_MAKEPGEN_2.out.tracking_out, PLINK2_MAKEPGEN_1.out.tracking_out)
    )
    ch_pgen_pvar_psam_filtered_per_pheno  = FILTER_MISSING_PER_PHENO.out.pgen_pvar_psam_out
    ch_versions = ch_versions.mix(FILTER_MISSING_PER_PHENO.out.versions.first())
    ch_tracking_filtered_per_pheno = ch_tracking.mix(FILTER_MISSING_PER_PHENO.out.tracking)


    PLINK2_MAKEPGEN_3 (
        ch_pgen_pvar_psam_filtered_per_pheno
            .map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], []) }
            .combine(trackingLastOrEmpty(FILTER_MISSING_PER_PHENO.out.tracking)),
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('filter_pass'),
        Channel.value(params.plink2_makepgen_3_options)
    )
    ch_pgen_pvar_psam_3  = PLINK2_MAKEPGEN_3.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_3.out.versions.first())
    ch_tracking_step1 = ch_tracking_filtered_per_pheno.mix(PLINK2_MAKEPGEN_3.out.tracking_out.first())


    F_COEFFICIENT_FILTERING (
        ch_pgen_pvar_psam_3,
        PLINK2_MAKEPGEN_3.out.tracking_out.first()
    )
    ch_pgen_pvar_psam_4 = F_COEFFICIENT_FILTERING.out.pgen_pvar_psam_out
    ch_versions = ch_versions.mix(F_COEFFICIENT_FILTERING.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(F_COEFFICIENT_FILTERING.out.tracking)

    PCA (
        ch_pgen_pvar_psam_4,
        ch_phenotype,
        ch_frq,
        F_COEFFICIENT_FILTERING.out.tracking.last()
    )
    ch_sscore = PCA.out.sscore
    ch_pca_plot_file = PCA.out.plot_file
    ch_king_cutoff_prune_in = PCA.out.king_cutoff_prune_in
    ch_versions = ch_versions.mix(PCA.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(PCA.out.tracking)

    MERGE_SEX_COVAR (
        ch_sscore.join(
            ch_pgen_pvar_psam_4.map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, psam_file) },
            by: 0
        )
    )
    ch_regenie_covar_file = MERGE_SEX_COVAR.out.covar_file
    ch_versions = ch_versions.mix(MERGE_SEX_COVAR.out.versions.first())

    if (params.regenie_step1_kinship_filtering) {
        plink2_write_snplist_1_input = ch_pgen_pvar_psam_4
            .join(ch_king_cutoff_prune_in, by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, king_cutoff_prune_in -> tuple(meta, pgen_file, pvar_file, psam_file, king_cutoff_prune_in) }
    } else {
        plink2_write_snplist_1_input = ch_pgen_pvar_psam_4
            .map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, []) }
    }

    PLINK2_WRITE_SNPLIST_1 (
        plink2_write_snplist_1_input.combine(trackingLastOrEmpty(PCA.out.tracking)),
        Channel.value('writesnp_pass'),
        Channel.value(params.plink2_write_snplist_qc_options)
    )
    ch_snplist  = PLINK2_WRITE_SNPLIST_1.out.snplist
    ch_id  = PLINK2_WRITE_SNPLIST_1.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_1.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(PLINK2_WRITE_SNPLIST_1.out.tracking_out.first())
    

    PLINK2_MAKEPGEN_4 (
        ch_pgen_pvar_psam_filtered_per_pheno
            .map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], []) }
            .combine(trackingLastOrEmpty(FILTER_MISSING_PER_PHENO.out.tracking)),
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('step2_filter'),
        Channel.value(params.plink2_makepgen_4_options)
    )
    ch_pgen_pvar_psam_5  = PLINK2_MAKEPGEN_4.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_4.out.versions.first())
    ch_tracking_step2 = PLINK2_MAKEPGEN_4.out.tracking_out.first()


    PLINK2_MAKEPGEN_5 (
        ch_pgen_pvar_psam_5
            .map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], []) }
            .combine(trackingFirstOrEmpty(PLINK2_MAKEPGEN_4.out.tracking_out)),
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('step2_input'),
        Channel.value(params.plink2_makepgen_5_options)
    )
    ch_pgen_pvar_psam_6  = PLINK2_MAKEPGEN_5.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_5.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(PLINK2_MAKEPGEN_5.out.tracking_out.first())

    
    PLINK2_WRITE_SNPLIST_2 (
        ch_pgen_pvar_psam_6.map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, []) }
            .combine(trackingFirstOrEmpty(PLINK2_MAKEPGEN_5.out.tracking_out)),
        Channel.value('writesnp_step2'),
        Channel.value(params.plink2_write_snplist_step2_options),
    )
    ch_step2_snplist = PLINK2_WRITE_SNPLIST_2.out.snplist
    ch_step2_sample_ids = PLINK2_WRITE_SNPLIST_2.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(PLINK2_WRITE_SNPLIST_2.out.tracking_out.first())

    BCFTOOLS_VIEW_2 (
        ch_filtered_vcf_with_index
            .join(ch_step2_sample_ids, by: 0)            // Samples file
            .join(ch_step2_snplist, by: 0)               // SNPs file
            .map { meta, vcf_file, tbi_file, samples_file, snplist_file -> tuple(meta, vcf_file, tbi_file, [], [], samples_file, snplist_file) }
            .combine(trackingFirstOrEmpty(PLINK2_WRITE_SNPLIST_2.out.tracking_out)),
        Channel.value("--output-type z --write-index=tbi"),         // input args
        Channel.value("view2")
    )
    ch_step2_vcf  = BCFTOOLS_VIEW_2.out.vcf
    ch_step2_vcf_tbi  = BCFTOOLS_VIEW_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(BCFTOOLS_VIEW_2.out.tracking_out.first())


    ch_regenie_step_1_input_part = ch_pgen_pvar_psam_4
            .join(ch_id, by: 0)
            .join(ch_snplist, by: 0)
            .join(ch_phenotype, by: 0)
    
    if (params.regenie_step1_options.contains("covarColList")) {
        ch_regenie_step_1_input = ch_regenie_step_1_input_part
            .join(ch_regenie_covar_file, by: 0)
    } else {
        ch_regenie_step_1_input = ch_regenie_step_1_input_part
            .map { meta, pgen_file, pvar_file, psam_file, id_file, snplist_file, phenotype_file ->
                tuple(meta, pgen_file, pvar_file, psam_file, id_file, snplist_file, phenotype_file, [])
            }
    }

    REGENIE_STEP1 (
        ch_regenie_step_1_input
            .combine(trackingFirstOrEmpty(PLINK2_WRITE_SNPLIST_1.out.tracking_out)),
        Channel.value(params.regenie_step1_options)
    )
    ch_regenie_step1_loco  = REGENIE_STEP1.out.loco
    ch_regenie_step1_pred_list  = REGENIE_STEP1.out.pred_list
    ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(REGENIE_STEP1.out.tracking_out.first())


    // Python vcf2aaf (replaces slower R version)
    python_vcf2aaf_script_ch = Channel.fromPath("${projectDir}/modules/local/python/vcf2aaf/assets/vcf2aaf.py", checkIfExists: true).first()
    // Split vcf2aaf options: if two parts use both, if one part use empty string + the value
    def vcf2aaf_opts = params.rscript_vcf2aaf_options.split(' ')
    def vcf2aaf_opt1 = vcf2aaf_opts.size() > 1 ? vcf2aaf_opts[0] : ''
    def vcf2aaf_opt2 = vcf2aaf_opts.size() > 1 ? vcf2aaf_opts[1] : vcf2aaf_opts[0]
    PYTHON_VCFTOAAF (
        ch_filtered_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) }
            .combine(python_vcf2aaf_script_ch),
        Channel.value(vcf2aaf_opt1),
        Channel.value(vcf2aaf_opt2)
    )
    ch_aaf  = PYTHON_VCFTOAAF.out.aaf
    ch_versions = ch_versions.mix(PYTHON_VCFTOAAF.out.versions.first())
    
    // bcftools/polars assign_annotations (replaces slower R version, also fixes NULL bug)
    python_assign_annotations_script_ch = Channel.fromPath("${projectDir}/modules/local/bcftools/assign_annotations/assets/assign_annotations.py", checkIfExists: true).first()
    BCFTOOLS_ASSIGN_ANNOTATIONS (
        ch_step2_vcf
            .combine(ch_masks)
            .combine(python_assign_annotations_script_ch),
        Channel.value(params.rscript_annotate_options)  // Uses same params format
    )
    ch_annotations  = BCFTOOLS_ASSIGN_ANNOTATIONS.out.annotations
    ch_setlist  = BCFTOOLS_ASSIGN_ANNOTATIONS.out.setlist
    ch_versions = ch_versions.mix(BCFTOOLS_ASSIGN_ANNOTATIONS.out.versions.first())

    if (params.use_dosage) {
        PLINK2_EXPORT_OTHER (
            ch_filtered_vcf_with_index
                .join(ch_pgen_pvar_psam_6.map { meta, pgen, pvar, psam -> tuple(meta, psam) }, by: 0),
            Channel.value('traw'),
            Channel.value(params.plink2_export_other_options)
        )
        ch_traw = PLINK2_EXPORT_OTHER.out.out_file
        ch_versions = ch_versions.mix(PLINK2_EXPORT_OTHER.out.versions.first())

        PLINK2_IMPORT_DOSAGE (
            ch_pgen_pvar_psam_6
                .join(ch_traw, by: 0)
                .map { meta, pgen, pvar, psam, traw -> tuple(meta, psam, traw) },
            Channel.value('import_dosage'),
            Channel.value(params.plink2_import_dosage_options)
        )
        ch_pgen_pvar_psam_with_dosage = PLINK2_IMPORT_DOSAGE.out.out_pgen_pvar_psam
        ch_versions = ch_versions.mix(PLINK2_IMPORT_DOSAGE.out.versions.first())
    } else {
        ch_pgen_pvar_psam_with_dosage = ch_pgen_pvar_psam_6
    }

    ch_regenie_step_2_input_part = ch_pgen_pvar_psam_with_dosage
            .join(ch_phenotype, by: 0)
            .join(ch_annotations, by: 0)
            .join(ch_setlist, by: 0)
            .join(ch_aaf, by: 0)
            .join(ch_regenie_step1_pred_list, by: 0)
    if (params.regenie_step1_options.contains("covarColList")) {
        ch_regenie_step_2_input = ch_regenie_step_2_input_part
            .join(ch_regenie_covar_file, by: 0)
    } else {
        ch_regenie_step_2_input = ch_regenie_step_2_input_part
            .map { meta, pgen, pvar, psam, pheno, anno, setlist, aaf, pred_list ->
                tuple(meta, pgen, pvar, psam, pheno, anno, setlist, aaf, pred_list, [])
            }
    }

    REGENIE_STEP2 (
        ch_regenie_step_2_input
            .combine(ch_masks)
            .combine(trackingFirstOrEmpty(BCFTOOLS_VIEW_2.out.tracking_out)),
        Channel.value(params.regenie_step2_options),
        Channel.value(!params.skip_reporting)
    )
    ch_regenie_step2_regenie_out  = REGENIE_STEP2.out.regenie_out
    ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(REGENIE_STEP2.out.tracking_out.first())

    if (!params.skip_reporting) {
        //r_script_buildreports_ch = Channel.fromPath(params.rscript_buildreports_path, checkIfExists: true).first()
        //RSCRIPT_BUILDREPORTS (
        //    REGENIE_STEP2.out.masks_snplist
        //        .join(ch_regenie_step2_regenie_out, by: 0)
        //        .join(ch_step2_vcf, by: 0)
        //        .join(ch_phenotype, by: 0)
        //        .join(ch_annotations, by: 0)
        //        .combine(r_script_buildreports_ch)
        //)
        //ch_annotated_snps  = RSCRIPT_BUILDREPORTS.out.annotated_snps
        //ch_res_log10p_1_annotated  = RSCRIPT_BUILDREPORTS.out.res_log10p_1_annotated
        //ch_annotated_snps_with_sample_ids  = RSCRIPT_BUILDREPORTS.out.annotated_snps_with_sample_ids
        //ch_versions = ch_versions.mix(RSCRIPT_BUILDREPORTS.out.versions.first())

        ch_regenie_step2_regenie_out  //.map { t -> [t[0].id, t[1]] }
            .transpose()
            .map { meta, fl -> tuple(meta, RegenieUtil.getPhenotypeByChunk(meta.id, fl), fl) }
            .groupTuple(by: [0, 1])
            .set { ch_regenie_step2_by_phenotype }

        py_script_csv_concat_ch = Channel.fromPath(params.py_script_csv_concat_path, checkIfExists: true).first()
        MERGE_RESULTS (
            ch_regenie_step2_by_phenotype
                .combine(py_script_csv_concat_ch)
        )
        ch_results_merged = MERGE_RESULTS.out.results_merged

        REPORTING (
            ch_results_merged,
            ch_masks,
            ch_phenotype,
            ch_pca_plot_file,
            ch_eda_plots,
            ch_setlist
        )
    }

    ch_tracking = ch_tracking_step1.mix(ch_tracking_step2)

    GENERATE_TRACKING_REPORT (
        ch_meta.merge(ch_tracking.collect().map { paths -> ['paths': paths] })
            .map { meta, paths -> tuple(meta, paths['paths']) }
    )
    ch_versions = ch_versions.mix(GENERATE_TRACKING_REPORT.out.versions.first())

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
    // Exposed for integration tests (IT-6/IT-7): the regenie association results and the
    // annotation artifacts the naive-LOG10P soundness check recomputes against.
    regenie_step2_out       = ch_regenie_step2_regenie_out
    regenie_step1_pred_list = ch_regenie_step1_pred_list
    regenie_step1_loco      = ch_regenie_step1_loco
    setlist                 = ch_setlist
    annotations             = ch_annotations
    vep_annotated_vcf       = ch_vcf_with_sample_names_corrected
    vep_annotated_vcf_tbi   = ch_vcf_with_sample_names_corrected_tbi
    eda_plots_out           = ch_eda_plots

}

def trackingFirstOrEmpty(ch) {
    ch.first().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}

def trackingFirstOrFallback(primaryCh, fallbackCh) {
    primaryCh
        .first()
        .concat(fallbackCh.first())
        .first()
        .ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
}

def trackingLastOrEmpty(ch) {
    ch.last().ifEmpty(file("${projectDir}/assets/empty_tracking.json", checkIfExists: true))
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
