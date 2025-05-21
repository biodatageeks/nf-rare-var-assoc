process REGENIE_STEP2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/rgcgithub/regenie/regenie:v4.1.gz'

    input:
    tuple val(meta), path(bgen), path(sample), path(phenotype), path(annotations), path(setlist), path(aaf), path(step1_pred_list), path(covar_file)
    path(masks)
    val(input_args)

    output:
    tuple val(meta), path("*_masks.bed"), path("*_masks.bim"), path("*_masks.fam"), emit: masks_bed_bim_fam
    tuple val(meta), path("*_masks.snplist"), emit: masks_snplist
    tuple val(meta), path("*.regenie"), emit: regenie_out
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def covar_file_corr_format = covar_file + "_correct_format.txt"
    def covar_file_input = covar_file ? "--covarFile " + covar_file_corr_format : ""
    """
    if [ "${covar_file}" != ""]; then
        sed '1s/^#//' ${covar_file} > ${covar_file_corr_format}
    fi
    regenie \\
       --step 2 \\
       --threads ${task.cpus} \\
       --bgen ${bgen} \\
       --sample ${sample} \\
       --phenoFile ${phenotype} \\
       ${covar_file_input} \\
       --anno-file ${annotations} \\
       --set-list ${setlist} \\
       --mask-def ${masks} \\
       --pred ${step1_pred_list} \\
       --aaf-file ${aaf} \\
       --out ${prefix}_step2 \\
       $args $input_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_masks.bed
    touch ${prefix}_masks.bim
    touch ${prefix}_masks.fam
    touch ${prefix}_masks.snplist
    touch ${prefix}_Y1.regenie

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """
}
