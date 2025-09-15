/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FILE          } from '../modules/local/cmds/download_file'
include { MERGE_RESULTS          } from '../modules/local/cmds/merge_results'
include { RENAME as RENAME_1     } from '../modules/local/cmds/rename'
include { RENAME as RENAME_2     } from '../modules/local/cmds/rename'
include { GENERATE_TRACKING_REPORT           } from '../modules/local/python/generate_tracking_report'
include { EXPLORATORY_DATA_ANALYSIS          } from '../modules/local/python/eda'
include { RSCRIPT_BUILDREPORTS   } from '../modules/local/rscript/buildreports'
include { REGENIE_STEP2          } from '../modules/local/regenie/step2'
include { REGENIE_STEP1          } from '../modules/local/regenie/step1'
include { RSCRIPT_VCFTOAAF       } from '../modules/local/rscript/vcf2aaf'
include { RSCRIPT_ASSIGN_ANNOTATIONS         } from '../modules/local/rscript/assign_annotations'
include { BGENIX                 } from '../modules/local/bgenix'
include { QCTOOL                 } from '../modules/local/qctool'
include { PLINK2_IMPORT_DOSAGE   } from '../modules/local/plink2/import_dosage'
include { PLINK2_EXPORT_OTHER    } from '../modules/local/plink2/export_other'
include { PLINK2_EXPORT_BGEN     } from '../modules/local/plink2/export_bgen'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_1 } from '../modules/local/plink2/write_snplist'
include { PLINK2_WRITE_SNPLIST as PLINK2_WRITE_SNPLIST_2 } from '../modules/local/plink2/write_snplist'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_1 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_2 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_3 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_4 } from '../modules/local/plink2/makepgen'
include { PLINK2_MAKEPGEN as PLINK2_MAKEPGEN_5 } from '../modules/local/plink2/makepgen'
include { VEP_ANNOTATE           } from '../modules/local/vep/annotate'
include { VEP_UPDATECACHE        } from '../modules/local/vep/updatecache'
include { BCFTOOLS_REPLACE_SAMPLE_NAMES      } from '../modules/local/bcftools/replace_sample_names'
include { BCFTOOLS_VCF2PSAM      } from '../modules/local/bcftools/vcf2psam'
include { BCFTOOLS_VCF2FRQ       } from '../modules/local/bcftools/vcf2frq'
include { BCFTOOLS_TAG2TAG       } from '../modules/local/bcftools/tag2tag'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_1   } from '../modules/local/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_2   } from '../modules/local/bcftools/view'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_1   } from '../modules/local/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_2   } from '../modules/local/bcftools/filter'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../modules/local/bcftools/index'
include { BCFTOOLS_NORM          } from '../modules/local/bcftools/norm'
include { BCFTOOLS_ANNOTATE      } from '../modules/nf-core/bcftools/annotate'
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

process CHECK_X_CHROM_PRESENT {
    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), env(has_x), path(pgen), path(pvar), path(psam), emit: has_x

    script:
    """
    # Check if chromosome X or chrX exists in the .pvar file
    if grep -E "^(X|chrX)\\s" ${pvar}; then
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
    ch_phenotype
    ch_all_samples

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
    

    BCFTOOLS_INDEX_1 (
        ch_input_vcf
    )
    ch_input_vcf_tbi = BCFTOOLS_INDEX_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    BCFTOOLS_VIEW_1 (
        ch_input_vcf.join(ch_input_vcf_tbi, by: 0),
        Channel.of([]),                                       // No regions file
        Channel.of([]),                                       // No targets file
        ch_all_samples,                                       // Samples file
        Channel.of([]),                                       // SNPs file
        Channel.value("--output-type z --write-index=tbi"),       // input args
        Channel.value("view1"),
        Channel.of([])
    )
    ch_all_samples_vcf  = BCFTOOLS_VIEW_1.out.vcf
    ch_all_samples_vcf_tbi  = BCFTOOLS_VIEW_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_1.out.versions.first())
    ch_tracking = BCFTOOLS_VIEW_1.out.tracking_out.first()

    EXPLORATORY_DATA_ANALYSIS (
        ch_all_samples_vcf
            .join(ch_all_samples_vcf_tbi, by: 0)
            .join(ch_phenotype, by: 0)
    )
    ch_eda_plots = EXPLORATORY_DATA_ANALYSIS.out.plots
    ch_versions = ch_versions.mix(EXPLORATORY_DATA_ANALYSIS.out.versions.first())

    BCFTOOLS_FILTER_1 (
        ch_all_samples_vcf.join(ch_all_samples_vcf_tbi, by: 0),
        Channel.of([]),                                       // No regions file
        Channel.of([]),                                       // No targets file
        Channel.of([]),                                       // No samples file
        Channel.of([]),                                       // SNPs file
        Channel.value(params.bcftools_filter_1_options),        // input args
        Channel.value("view2"),
        BCFTOOLS_VIEW_1.out.tracking_out.first()
    )
    ch_filter_1_vcf  = BCFTOOLS_FILTER_1.out.vcf
    ch_filter_1_vcf_tbi  = BCFTOOLS_FILTER_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_1.out.versions.first())
    ch_tracking = BCFTOOLS_FILTER_1.out.tracking_out.first()
    // shouldn't this be: ch_tracking = ch_tracking.mix(BCFTOOLS_FILTER_1.out.tracking_out.first())

    BCFTOOLS_FILTER_2 (
        ch_filter_1_vcf.join(ch_filter_1_vcf_tbi, by: 0),
        Channel.of([]),                                       // No regions file
        Channel.of([]),                                       // No targets file
        Channel.of([]),                                       // No samples file
        Channel.of([]),                                       // SNPs file
        Channel.value(params.bcftools_filter_2_options),        // input args
        Channel.value("view3"),
        BCFTOOLS_FILTER_1.out.tracking_out.first()
    )
    ch_filter_2_vcf  = BCFTOOLS_FILTER_2.out.vcf
    ch_filter_2_vcf_tbi  = BCFTOOLS_FILTER_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_2.out.versions.first())
    ch_tracking = BCFTOOLS_FILTER_2.out.tracking_out.first()
    // shouldn't this be: ch_tracking = ch_tracking.mix(BCFTOOLS_FILTER_2.out.tracking_out.first())
    
    BCFTOOLS_NORM (
        ch_filter_2_vcf.join(ch_filter_2_vcf_tbi, by: 0),
        ch_vep_cachesubdir.map { t -> "${t}/${params.vep_fasta_path}" },
        Channel.value("norm"),
        BCFTOOLS_FILTER_2.out.tracking_out.first()
    )
    ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
    ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
    ch_tracking = ch_tracking.mix(BCFTOOLS_NORM.out.tracking_out.first())

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
        ch_annotated_vcf.join(ch_annotated_vcf_tbi, by: 0),
        ch_vep_cachesubdir,
        Channel.value(params.vep_annotate_species),
        Channel.value(params.vep_fasta_path),
        Channel.value(params.vep_annotate_options)
    )
    ch_vep_vcf  = VEP_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions.first())

    BCFTOOLS_REPLACE_SAMPLE_NAMES (
        ch_vep_vcf,
        Channel.value(params.bcftools_replace_sample_names_sed_arg),
        Channel.value('replace_sample_names')
    )
    ch_vcf_with_sample_names_corrected = BCFTOOLS_REPLACE_SAMPLE_NAMES.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REPLACE_SAMPLE_NAMES.out.versions.first())

    BCFTOOLS_INDEX_2 (
        ch_vcf_with_sample_names_corrected
    )
    ch_vep_vcf_tbi = BCFTOOLS_INDEX_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    ch_vep_vcf_corr_sample_names_with_index = ch_vcf_with_sample_names_corrected
        .join(ch_vep_vcf_tbi, by: 0)
        .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) }


    BCFTOOLS_VCF2FRQ (
        ch_vep_vcf_corr_sample_names_with_index
    )
    ch_frq = BCFTOOLS_VCF2FRQ.out.frq
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2FRQ.out.versions.first())

    BCFTOOLS_TAG2TAG (
        ch_vep_vcf_corr_sample_names_with_index,
        Channel.value(params.bcftools_tag2tag_tag_from),
        Channel.value(params.bcftools_tag2tag_tag_to)
    )
    ch_vcf_with_dosage_tag = BCFTOOLS_TAG2TAG.out.vcf
    ch_vcf_with_dosage_tag_tbi = BCFTOOLS_TAG2TAG.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_TAG2TAG.out.versions.first())

    ch_vcf_with_dosage_tag_with_index = ch_vcf_with_dosage_tag
        .join(ch_vcf_with_dosage_tag_tbi, by: 0)
        .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file) }

    BCFTOOLS_VCF2PSAM (
        ch_vep_vcf_corr_sample_names_with_index
    )
    ch_unk_sex_psam = BCFTOOLS_VCF2PSAM.out.psam
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2PSAM.out.versions.first())

    PLINK2_MAKEPGEN_1 (
        ch_vcf_with_dosage_tag_with_index
            .join(ch_unk_sex_psam, by: 0)
            .map { meta, vcf_file, tbi_file, unk_sex_psam_file -> tuple(meta, [], [], unk_sex_psam_file, vcf_file, tbi_file, [], [], []) },
        Channel.value(''),
        Channel.value(''),
        Channel.value(params.plink2_makepgen_1_vcf_input_options),
        Channel.value('plink2_makepgen_1'),
        Channel.value(params.plink2_makepgen_1_options),
        BCFTOOLS_NORM.out.tracking_out.first()
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
            .map { meta, has_x, pgen, pvar, psam, frq -> tuple(meta, pgen, pvar, psam, [], [], frq, [], []) },
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('plink2_makepgen_2_impute_sex'),
        Channel.value(params.plink2_makepgen_2_options),
        PLINK2_MAKEPGEN_1.out.tracking_out.first()
    )
    ch_pgen_pvar_psam_2  = PLINK2_MAKEPGEN_2.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_2.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN_2.out.tracking_out.first())

    renamed_file_name = ch_meta.map { t -> "${t.id}_plink2_makepgen_1.psam" }.first()
    RENAME_1 (
        ch_pgen_pvar_psam_2.map { meta, pgen, pvar, psam -> tuple(meta, psam) },
        renamed_file_name
    )
    ch_renamed_psam  = RENAME_1.out.output

    ch_pgen_pvar_psam_before_quality_filtering = split_data.with_x
        .join(ch_renamed_psam, by: 0)
        .map { meta, has_x, pgen, pvar, psam, psam_imputesex -> tuple(meta, pgen, pvar, psam_imputesex) }
        .mix(split_data.without_x.map { meta, has_x, pgen, pvar, psam -> tuple(meta, pgen, pvar, psam) })

    FILTER_MISSING_PER_PHENO (
        ch_pgen_pvar_psam_before_quality_filtering,
        ch_phenotype,
        PLINK2_MAKEPGEN_2.out.tracking_out.first().mix(PLINK2_MAKEPGEN_1.out.tracking_out.first()).collect()
    )
    ch_pgen_pvar_psam_filtered_per_pheno  = FILTER_MISSING_PER_PHENO.out.pgen_pvar_psam_out
    ch_versions = ch_versions.mix(FILTER_MISSING_PER_PHENO.out.versions.first())
    ch_tracking_filtered_per_pheno = ch_tracking.mix(FILTER_MISSING_PER_PHENO.out.tracking)


    PLINK2_MAKEPGEN_3 (
        ch_pgen_pvar_psam_filtered_per_pheno
            .join(ch_meta.merge(ch_hild), by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, hild_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], hild_file) },
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('filter_pass'),
        Channel.value(params.plink2_makepgen_3_options),
        FILTER_MISSING_PER_PHENO.out.tracking.last()
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
    ch_versions = ch_versions.mix(PCA.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(PCA.out.tracking)

    PLINK2_WRITE_SNPLIST_1 (
        ch_pgen_pvar_psam_4.map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, []) },
        Channel.value('writesnp_pass'),
        Channel.value(params.plink2_write_snplist_qc_options),
        PCA.out.tracking.last()
    )
    ch_snplist  = PLINK2_WRITE_SNPLIST_1.out.snplist
    ch_id  = PLINK2_WRITE_SNPLIST_1.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_1.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(PLINK2_WRITE_SNPLIST_1.out.tracking_out.first())
    

    PLINK2_MAKEPGEN_4 (
        ch_pgen_pvar_psam_filtered_per_pheno
            .join(ch_meta.merge(ch_hild), by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, hild_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], hild_file) },
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('step2_filter'),
        Channel.value(params.plink2_makepgen_4_options),
        FILTER_MISSING_PER_PHENO.out.tracking.last()
    )
    ch_pgen_pvar_psam_5  = PLINK2_MAKEPGEN_4.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_4.out.versions.first())
    ch_tracking_step2 = PLINK2_MAKEPGEN_4.out.tracking_out.first()


    PLINK2_MAKEPGEN_5 (
        ch_pgen_pvar_psam_5
            .map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], []) },
        Channel.value('--remove'),
        Channel.value('--exclude'),
        Channel.value(''),
        Channel.value('step2_input'),
        Channel.value(params.plink2_makepgen_5_options),
        PLINK2_MAKEPGEN_4.out.tracking_out.first()
    )
    ch_pgen_pvar_psam_6  = PLINK2_MAKEPGEN_5.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN_5.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(PLINK2_MAKEPGEN_5.out.tracking_out.first())

    
    PLINK2_WRITE_SNPLIST_2 (
        ch_pgen_pvar_psam_6.map { meta, pgen_file, pvar_file, psam_file -> tuple(meta, pgen_file, pvar_file, psam_file, []) },
        Channel.value('writesnp_step2'),
        Channel.value(params.plink2_write_snplist_step2_options),
        PLINK2_MAKEPGEN_5.out.tracking_out.first()
    )
    ch_step2_snplist = PLINK2_WRITE_SNPLIST_2.out.snplist
    ch_step2_sample_ids = PLINK2_WRITE_SNPLIST_2.out.id
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST_2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(PLINK2_WRITE_SNPLIST_2.out.tracking_out.first())

    BCFTOOLS_VIEW_2 (
        ch_vcf_with_dosage_tag_with_index,
        Channel.of([]),                     // No regions file
        Channel.of([]),                     // No targets file
        ch_step2_sample_ids.map { t -> t[1] },    // Samples file
        ch_step2_snplist.map { t -> t[1] },       // SNPs file
        Channel.value("--output-type z --write-index=tbi"),         // input args
        Channel.value("view2"),
        PLINK2_WRITE_SNPLIST_2.out.tracking_out.first()
    )
    ch_step2_vcf  = BCFTOOLS_VIEW_2.out.vcf
    ch_step2_vcf_tbi  = BCFTOOLS_VIEW_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(BCFTOOLS_VIEW_2.out.tracking_out.first())


    ch_regenie_step_1_input_part = ch_pgen_pvar_psam_4
            .join(ch_id, by: 0)
            .join(ch_snplist, by: 0)
            .join(ch_phenotype, by: 0)
    
    if (params.regenie_step1_options.contains("covarColList")) {
        ch_regenie_step_1_input = ch_regenie_step_1_input_part
            .join(ch_sscore, by: 0)
    } else {
        ch_regenie_step_1_input = ch_regenie_step_1_input_part
            .map { meta, pgen_file, pvar_file, psam_file, id_file, snplist_file, phenotype_file ->
                tuple(meta, pgen_file, pvar_file, psam_file, id_file, snplist_file, phenotype_file, [])
            }
    }

    REGENIE_STEP1 (
        ch_regenie_step_1_input,
        Channel.value(params.regenie_step1_options),
        PLINK2_WRITE_SNPLIST_1.out.tracking_out.first()
    )
    ch_regenie_step1_loco  = REGENIE_STEP1.out.loco
    ch_regenie_step1_pred_list  = REGENIE_STEP1.out.pred_list
    ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions.first())
    ch_tracking_step1 = ch_tracking_step1.mix(REGENIE_STEP1.out.tracking_out.first())


    r_script_vcf2aaf_ch = Channel.fromPath(params.rscript_vcf2aaf_path, checkIfExists: true)
    RSCRIPT_VCFTOAAF (
        r_script_vcf2aaf_ch,
        ch_vcf_with_dosage_tag_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) },
        Channel.value(params.rscript_vcf2aaf_options)
    )
    ch_aaf  = RSCRIPT_VCFTOAAF.out.aaf
    ch_versions = ch_versions.mix(RSCRIPT_VCFTOAAF.out.versions.first())
    
    r_script_annotate_ch = Channel.fromPath(params.rscript_annotate_path, checkIfExists: true)
    RSCRIPT_ASSIGN_ANNOTATIONS (
        r_script_annotate_ch,
        ch_step2_vcf,
        ch_masks,
        Channel.value(params.rscript_annotate_options)
    )
    ch_annotations  = RSCRIPT_ASSIGN_ANNOTATIONS.out.annotations
    ch_setlist  = RSCRIPT_ASSIGN_ANNOTATIONS.out.setlist
    ch_versions = ch_versions.mix(RSCRIPT_ASSIGN_ANNOTATIONS.out.versions.first())


    PLINK2_EXPORT_OTHER (
        ch_vcf_with_dosage_tag_with_index
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


    ch_regenie_step_2_input_part = ch_pgen_pvar_psam_with_dosage
            .join(ch_phenotype, by: 0)
            .join(ch_annotations, by: 0)
            .join(ch_setlist, by: 0)
            .join(ch_aaf, by: 0)
            .join(ch_regenie_step1_pred_list, by: 0)
    if (params.regenie_step1_options.contains("covarColList")) {
        ch_regenie_step_2_input = ch_regenie_step_2_input_part
            .join(ch_sscore, by: 0)
    } else {
        ch_regenie_step_2_input = ch_regenie_step_2_input_part
            .map { meta, pgen, pvar, psam, pheno, anno, setlist, aaf, pred_list ->
                tuple(meta, pgen, pvar, psam, pheno, anno, setlist, aaf, pred_list, [])
            }
    }

    REGENIE_STEP2 (
        ch_regenie_step_2_input,
        ch_masks,
        Channel.value(params.regenie_step2_options),
        BCFTOOLS_VIEW_2.out.tracking_out.first()
    )
    ch_regenie_step2_masks_bed_bim_fam  = REGENIE_STEP2.out.masks_bed_bim_fam
    ch_regenie_step2_masks_snplist  = REGENIE_STEP2.out.masks_snplist
    ch_regenie_step2_regenie_out  = REGENIE_STEP2.out.regenie_out
    ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions.first())
    ch_tracking_step2 = ch_tracking_step2.mix(REGENIE_STEP2.out.tracking_out.first())

    r_script_buildreports_ch = Channel.fromPath(params.rscript_buildreports_path, checkIfExists: true)
    RSCRIPT_BUILDREPORTS (
        r_script_buildreports_ch,
        ch_regenie_step2_masks_snplist,
        ch_regenie_step2_regenie_out,
        ch_step2_vcf,
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
        ch_phenotype,
        ch_pca_plot_file,
        ch_eda_plots.collect(),
        ch_setlist
    )

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
