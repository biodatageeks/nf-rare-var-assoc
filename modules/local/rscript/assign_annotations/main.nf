process RSCRIPT_ASSIGN_ANNOTATIONS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    path(r_script_ch)
    tuple val(meta), path(vcf)
    tuple val(meta), path(sample)
    path(masks)
    val(input_args)

    output:
    tuple val(meta), path("*.annotations"), emit: annotations
    tuple val(meta), path("*.setlist"), emit: setlist
    tuple val(meta), path("*_r_out.sample"), emit: out_sample
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.4.2-0.2'
    """
    Rscript \\
        ${r_script_ch} \\
        --vcf-path ${vcf} \\
        --sample-path ${sample} \\
        --masks-path ${masks} \\
        --out-sample-path ${prefix}_r_out.sample \\
        --out-anno-path ${prefix}.annotations \\
        --out-setlist-path ${prefix}.setlist \\
        $input_args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotations
    touch ${prefix}.masks
    touch ${prefix}.setlist
    touch ${prefix}_r_out.sample

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """
}
