include { PLINK2_WRITE_SNPLIST   } from '../../../modules/local/plink2/write_snplist'
include { PLINK2_MAKEBED         } from '../../../modules/local/plink2/makebed'

// Process to extract unique phenotypes and corresponding sample lists
process EXTRACT_PHENOTYPES_AND_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(pheno_file)

    output:
    tuple val(meta), path("*_samples_*.txt"), emit: samples

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract unique phenotype values (excluding header)
    cut -f 3 ${pheno_file} | sort | uniq | grep -v "Y1" > pheno_values.txt

    # Loop through phenotypes to create sample lists
    while read pheno; do
        awk -F'\\t' -v p="\$pheno" '\$3 == p {print \$1 "\\t" \$2}' ${pheno_file} > ${prefix}_samples_\${pheno}.txt
    done < pheno_values.txt
    """
}


workflow FILTER_MISSING_PER_PHENO {

    take:
    ch_bed_bim_fam
    ch_phenotype

    main:

    ch_versions = Channel.empty()


    EXTRACT_PHENOTYPES_AND_SAMPLES(
        ch_phenotype
    )
    ch_sample_files = EXTRACT_PHENOTYPES_AND_SAMPLES.out.samples
        .map { meta, files -> files }.flatten()
    
    //ch_sample_files = ch_sample_files_together.map { meta, files -> meta }
    //    .combine(ch_sample_files_together.map { meta, files -> files }.flatten())


    PLINK2_WRITE_SNPLIST (
        ch_bed_bim_fam
            .combine(ch_sample_files)
            .map { meta, bed_file, bim_file, fam_file, sample_file ->
                tuple(meta + ['orig_id': meta.id, 'id': sample_file.getBaseName()], bed_file, bim_file, fam_file, sample_file)
            },
        Channel.value('identify_acceptable_variants'),
        Channel.value(params.plink2_missing_per_pheno_options)
    )
    ch_snplist = PLINK2_WRITE_SNPLIST.out.snplist
    ch_versions = ch_versions.mix(PLINK2_WRITE_SNPLIST.out.versions.first())


    PLINK2_MAKEBED (
        ch_bed_bim_fam
            .join(ch_snplist.map { meta, snplist_files ->
                tuple(['id': meta.orig_id], snplist_files)
            }.groupTuple(), by: 0)
            .map { meta, bed_file, bim_file, fam_file, snplist_files ->
                tuple(meta, bed_file, bim_file, fam_file, [], [], [], snplist_files)
            },
        Channel.value(''),
        Channel.value('--extract-intersect'),
        Channel.value('intersect_variants_to_keep'),
        Channel.value('')
    )
    ch_bed_bim_fam_out = PLINK2_MAKEBED.out.out_bed_bim_fam
    ch_versions = ch_versions.mix(PLINK2_MAKEBED.out.versions.first())


    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }

    emit:
    bed_bim_fam_out = ch_bed_bim_fam_out
    versions        = ch_versions
}