include { PLINK2_WRITE_SNPLIST            } from '../../../modules/local/plink2/write_snplist'
include { PLINK2_MAKEPGEN                 } from '../../../modules/local/plink2/makepgen'
include { EXTRACT_PHENOTYPES_AND_SAMPLES  } from '../../../modules/local/cmds/extract_phenotypes_and_samples'


workflow FILTER_MISSING_PER_PHENO {

    take:
    ch_pgen_pvar_psam
    ch_phenotype
    ch_tracking_in

    main:

    ch_versions = Channel.empty()


    EXTRACT_PHENOTYPES_AND_SAMPLES(
        ch_phenotype
    )
    ch_sample_files = EXTRACT_PHENOTYPES_AND_SAMPLES.out.samples
    //    .map { meta, files -> files }.flatten().view()


    PLINK2_WRITE_SNPLIST (
        ch_pgen_pvar_psam
            .join(ch_sample_files, by: 0)
            .flatMap { meta, pgen_file, pvar_file, psam_file, sample_files ->
                sample_files.collect { v -> tuple(meta, pgen_file, pvar_file, psam_file, v) } 
            }
            .map { meta, pgen_file, pvar_file, psam_file, sample_file ->
                tuple(meta + ['orig_id': meta.id, 'id': sample_file.getBaseName()], pgen_file, pvar_file, psam_file, sample_file)
            }
            .combine(ch_tracking_in.collect()),
        Channel.value('identify_acceptable_variants'),
        Channel.value(params.plink2_missing_per_pheno_options)
    )
    ch_snplist = PLINK2_WRITE_SNPLIST.out.snplist
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST.out.versions.first())
    ch_tracking = PLINK2_WRITE_SNPLIST.out.tracking_out


    PLINK2_MAKEPGEN (
        ch_pgen_pvar_psam
            .join(ch_snplist.map { meta, snplist_files ->
                tuple(meta + ['id': meta.orig_id] - ['orig_id': meta.orig_id], snplist_files)
            }.groupTuple(), by: 0)
            .map { meta, pgen_file, pvar_file, psam_file, snplist_files ->
                tuple(meta, pgen_file, pvar_file, psam_file, [], [], [], [], snplist_files)
            }
            .combine(PLINK2_WRITE_SNPLIST.out.tracking_out.first()),
        Channel.value(''),
        Channel.value('--extract-intersect'),
        Channel.value(''),
        Channel.value('intersect_variants_to_keep'),
        Channel.value('')
    )
    ch_pgen_pvar_psam_out = PLINK2_MAKEPGEN.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN.out.tracking_out.first())


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    pgen_pvar_psam_out = ch_pgen_pvar_psam_out
    versions        = ch_versions
    tracking        = ch_tracking
}