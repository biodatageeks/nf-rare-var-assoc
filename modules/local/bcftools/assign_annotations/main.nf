process BCFTOOLS_ASSIGN_ANNOTATIONS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21 conda-forge::polars=1.0.0"
    container 'docker.io/psuszynski/bioinf_combo:1.3.0'

    input:
    tuple val(meta), path(vcf), path(masks), path(python_script)
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

    python3 ${python_script} \\
        --vcf-path ${vcf} \\
        --masks-path ${masks} \\
        --out-anno-path ${prefix}.annotations \\
        --out-setlist-path ${prefix}.setlist \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        polars: \$(python3 -c "import polars; print(polars.__version__)")
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
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        polars: \$(python3 -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
