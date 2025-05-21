process REGENIE_STEP1 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/rgcgithub/regenie/regenie:v4.1.gz'

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(qc_pass_id), path(qc_pass_snplist), path(phenotype), path(covar_file)
    val(input_args)

    output:
    tuple val(meta), path("*.loco"), emit: loco
    tuple val(meta), path("*_pred.list"), emit: pred_list
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
       --step 1 \\
       --threads ${task.cpus} \\
       --bed ${bed.baseName} \\
       --keep ${qc_pass_id} \\
       --extract ${qc_pass_snplist} \\
       --phenoFile ${phenotype} \\
       ${covar_file_input} \\
       --out ${prefix}_step1 \\
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
    touch ${prefix}.loco
    touch ${prefix}_pred.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """
}
