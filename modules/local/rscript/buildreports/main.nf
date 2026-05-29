process RSCRIPT_BUILDREPORTS {
    tag "$meta.id"
    label 'process_5'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    tuple val(meta), path(regenie_step2_masks_snplist), path(regenie_step2_Y1_regenie), path(vcf), path(phenotype), path(annotations), path(r_script_ch)

    output:
    tuple val(meta), path("*_annotated_snps.csv"), emit: annotated_snps
    tuple val(meta), path("*_res_log10p_1_annotated.csv"), emit: res_log10p_1_annotated
    tuple val(meta), path("*_annotated_snps_with_sample_ids.csv"), emit: annotated_snps_with_sample_ids
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.4.2-0.1'
    """
    Rscript \\
        ${r_script_ch} \\
        ${regenie_step2_masks_snplist} \\
        ${regenie_step2_Y1_regenie} \\
        ${vcf} \\
        ${phenotype} \\
        ${annotations} \\
        ${prefix}_annotated_snps.csv \\
        ${prefix}_res_log10p_1_annotated.csv \\
        ${prefix}_annotated_snps_with_sample_ids.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptvcftoaaf: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated_snps.csv
    touch ${prefix}_res_log10p_1_annotated.csv
    touch ${prefix}_annotated_snps_with_sample_ids.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptvcftoaaf: ${VERSION}
    END_VERSIONS
    """
}
