process BCFTOOLS_ASSIGN_ANNOTATIONS_SHELL {
    tag "$meta.id"
    label 'process_1'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(masks), path(shell_script)
    val(input_args)

    output:
    tuple val(meta), path("*.annotations"), emit: annotations
    tuple val(meta), path("*.setlist")    , emit: setlist
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = input_args ?: ''
    """
    echo 'task.cpus: ${task.cpus}'

    bash ${shell_script} \\
        --vcf-path ${vcf} \\
        --masks-path ${masks} \\
        --out-anno-path ${prefix}.annotations \\
        --out-setlist-path ${prefix}.setlist \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
        bash: \$(bash --version | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotations
    touch ${prefix}.setlist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
        bash: \$(bash --version | head -n1)
    END_VERSIONS
    """
}
