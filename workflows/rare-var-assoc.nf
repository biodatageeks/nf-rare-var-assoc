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
include { FIX_ZERO_PL            } from '../modules/local/python/fix_zero_PL'
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
include { VEP_ANNOTATE           } from '../modules/local/vep/annotate'
include { VEP_UPDATECACHE        } from '../modules/local/vep/updatecache'
include { BCFTOOLS_VCF2PSAM      } from '../modules/local/bcftools/vcf2psam'
include { BCFTOOLS_VCF2FRQ       } from '../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_2   } from '../modules/local/bcftools/view'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_1   } from '../modules/local/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_2   } from '../modules/local/bcftools/filter'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_3 } from '../modules/local/bcftools/index'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { FILTER_MISSING_PER_PHENO           } from '../subworkflows/local/filter_missing_per_pheno'
include { F_COEFFICIENT_FILTERING            } from '../subworkflows/local/f_coefficient_filtering'
include { PCA                    } from '../subworkflows/local/pca'
include { REPORTING              } from '../subworkflows/local/reporting'
include { PREPARE                } from '../subworkflows/local/prepare'
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

    vep_cachedir = "${projectDir}/../vep_cachedir"
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_vep_cachedir = Channel.fromPath(vep_cachedir, checkIfExists: true)
    ch_masks = Channel.fromPath(params.input_masks, checkIfExists: true).first()
    ch_meta = ch_input_vcf.map { t -> t[0] }


    VEP_UPDATECACHE (
        ch_meta.first().combine(ch_vep_cachedir),
        Channel.value(params.vep_updatecache_species),
        Channel.value(params.vep_updatecache_options),
        Channel.value(params.vep_cache_url),
        Channel.value(tuple(params.ref_fasta_url, params.vep_fasta_path))
    )
    ch_vep_cachesubdir = VEP_UPDATECACHE.out.cachesubdir.first()
    ch_versions = ch_versions.mix(VEP_UPDATECACHE.out.versions.first())


    // call BCFTOOLS_ASSIGN_ANNOTATIONS with a dummy python script to pull the bioinf_combo image before HyperQueue provisions workers on PLGrid
    dummy_python_script_ch = Channel.fromPath("${projectDir}/modules/local/bcftools/assign_annotations/assets/dummy.py", checkIfExists: true).first()
    BCFTOOLS_ASSIGN_ANNOTATIONS_2 (
        ch_input_vcf
            .combine(ch_masks)
            .combine(dummy_python_script_ch),
        Channel.value(params.rscript_annotate_options)
    )


    PREPARE (
        ch_input_vcf,
        ch_vep_cachesubdir,
        ch_all_samples
    )
    ch_prepared_vcf = PREPARE.out.prepared_vcf
    ch_prepared_vcf_tbi = PREPARE.out.prepared_vcf_tbi
    ch_versions = ch_versions.mix(PREPARE.out.versions.first())
    ch_tracking = PREPARE.out.tracking

    ch_vep_vcf = Channel.empty()
    ch_eda_plots = Channel.empty()

    if (!params.skip_reporting) {
        if (!params.use_dosage) {
            // do this only to be able to produce plots
            FIX_ZERO_PL (
                ch_prepared_vcf,
                Channel.value(params.filter_and_enhance_vcf_calc_ds_min_gq),
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
            ch_vcf_with_pl_corrected = ch_prepared_vcf
            ch_vcf_with_pl_corrected_tbi = ch_prepared_vcf_tbi
        }

        EXPLORATORY_DATA_ANALYSIS (
            ch_vcf_with_pl_corrected
                .join(ch_vcf_with_pl_corrected_tbi, by: 0)
                .join(ch_phenotype, by: 0),
            Channel.value(params.use_dosage)
        )
        ch_eda_plots = EXPLORATORY_DATA_ANALYSIS.out.plots
        ch_versions = ch_versions.mix(EXPLORATORY_DATA_ANALYSIS.out.versions.first())
    }

    if ((!params.skip_preparation || !params.skip_reporting) && !params.use_dosage) {
        // turns out bcftools, even called two times, is faster than rust-htslib

        bcftools_filter_1_options = "--exclude 'QUAL<${params.filter_and_enhance_vcf_qual_min} || AVG(FORMAT/GQ)<${params.filter_and_enhance_vcf_avg_gq_min} || AVG(FORMAT/DP)<${params.filter_and_enhance_vcf_avg_dp_min} || AVG(FORMAT/DP)>${params.filter_and_enhance_vcf_avg_dp_max}' --output-type z --write-index=tbi"
        bcftools_filter_2_options = "--exclude 'FORMAT/GQ<${params.filter_and_enhance_vcf_sample_gq_min} | FORMAT/DP < ${params.filter_and_enhance_vcf_sample_dp_min} | FORMAT/DP > ${params.filter_and_enhance_vcf_sample_dp_max}' --set-GTs '.' --output-type z --write-index=tbi"
    
        BCFTOOLS_FILTER_1 (
            ch_prepared_vcf
                .join(ch_prepared_vcf_tbi, by: 0)
                .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file, [], [], [], []) }
                .combine(trackingLastOrEmpty(PREPARE.out.tracking)),
            Channel.value(bcftools_filter_1_options),
            Channel.value("filter1")
        )
        ch_filter_1_vcf  = BCFTOOLS_FILTER_1.out.vcf
        ch_filter_1_vcf_tbi  = BCFTOOLS_FILTER_1.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_1.out.versions.first())
        ch_tracking = ch_tracking.mix(BCFTOOLS_FILTER_1.out.tracking_out.first())

        BCFTOOLS_FILTER_2 (
            ch_filter_1_vcf
                .join(ch_filter_1_vcf_tbi, by: 0)
                .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file, [], [], [], []) }
                .combine(trackingFirstOrEmpty(BCFTOOLS_FILTER_1.out.tracking_out)),
            Channel.value(bcftools_filter_2_options),
            Channel.value("filter2"),
        )
        ch_filtered_vcf  = BCFTOOLS_FILTER_2.out.vcf
        ch_filtered_vcf_tbi  = BCFTOOLS_FILTER_2.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_2.out.versions.first())
        ch_tracking = ch_tracking.mix(BCFTOOLS_FILTER_2.out.tracking_out.first())
    } else {
        ch_filtered_vcf = ch_prepared_vcf
        ch_filtered_vcf_tbi = ch_prepared_vcf_tbi
    }

    if (params.skip_preparation == false) {
        VEP_ANNOTATE (
            ch_filtered_vcf
                .join(ch_filtered_vcf_tbi, by: 0)
                .combine(ch_vep_cachesubdir),
            Channel.value(params.vep_annotate_species),
            Channel.value(params.vep_fasta_path),
            Channel.value(params.vep_annotate_options)
        )
        ch_vep_vcf  = VEP_ANNOTATE.out.vcf  // also exposed as test-support emit
        ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions.first())

        BCFTOOLS_INDEX_3 (
            ch_vep_vcf
        )
        ch_vep_vcf_tbi = BCFTOOLS_INDEX_3.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_3.out.versions.first())

        ch_vep_vcf_with_index = ch_vep_vcf.join(ch_vep_vcf_tbi, by: 0)
    } else {
        ch_vep_vcf_with_index = ch_filtered_vcf.join(ch_filtered_vcf_tbi, by: 0)
    }


    BCFTOOLS_VCF2FRQ (
        ch_vep_vcf_with_index
    )
    ch_frq = BCFTOOLS_VCF2FRQ.out.frq
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2FRQ.out.versions.first())

    BCFTOOLS_VCF2PSAM (
        ch_vep_vcf_with_index
    )
    ch_unk_sex_psam = BCFTOOLS_VCF2PSAM.out.psam
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2PSAM.out.versions.first())

    PLINK2_MAKEPGEN_1 (
        ch_vep_vcf_with_index
            .join(ch_unk_sex_psam, by: 0)
            .map { meta, vcf_file, tbi_file, unk_sex_psam_file -> tuple(meta, [], [], unk_sex_psam_file, vcf_file, tbi_file, [], [], [],   []) },  // TODO add tracking generation to FILTER_AND_ENHANCE_VCF and use here the tracking out of PREPARE
            //.combine(BCFTOOLS_FILTER_2.out.tracking_out.first()),
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
        ch_vep_vcf_with_index
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
        ch_vep_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) }
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
            ch_vep_vcf_with_index
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
        r_script_buildreports_ch = Channel.fromPath(params.rscript_buildreports_path, checkIfExists: true).first()
        RSCRIPT_BUILDREPORTS (
            REGENIE_STEP2.out.masks_snplist
                .join(ch_regenie_step2_regenie_out, by: 0)
                .join(ch_step2_vcf, by: 0)
                .join(ch_phenotype, by: 0)
                .join(ch_annotations, by: 0)
                .combine(r_script_buildreports_ch)
        )
        ch_annotated_snps  = RSCRIPT_BUILDREPORTS.out.annotated_snps
        ch_res_log10p_1_annotated  = RSCRIPT_BUILDREPORTS.out.res_log10p_1_annotated
        ch_annotated_snps_with_sample_ids  = RSCRIPT_BUILDREPORTS.out.annotated_snps_with_sample_ids
        ch_versions = ch_versions.mix(RSCRIPT_BUILDREPORTS.out.versions.first())

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
    vep_annotated_vcf       = ch_vep_vcf
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
