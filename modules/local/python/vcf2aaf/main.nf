process PYTHON_VCFTOAAF {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::polars=1.0.0 bioconda::polars-bio=0.5.0"
    container 'docker.io/psuszynski/bioinf_combo:1.3.0'

    input:
    tuple val(meta), path(vcf), path(python_script)
    val(tag_name)           // e.g., "AF_nfe" - primary tag to extract
    val(default_tag_name)   // e.g., "AF" - fallback tag when primary is missing

    output:
    tuple val(meta), path("*_aaf.tsv"), emit: aaf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${python_script} \\
        ${vcf} \\
        ${prefix}_aaf.tsv \\
        ${tag_name} \\
        ${default_tag_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        polars-bio: \$(python3 -c "import polars_bio; print(polars_bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_aaf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        polars-bio: \$(python3 -c "import polars_bio; print(polars_bio.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
