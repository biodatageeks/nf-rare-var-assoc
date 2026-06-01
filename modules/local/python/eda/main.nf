process EXPLORATORY_DATA_ANALYSIS {

    tag "$meta.id"
    label 'process_2'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/python_tools:1.0.11'

    input:
    tuple val(meta), path(vcf), path(tbi), path(phenotype_file), path(python_script)
    val(use_dosage)

    output:
    tuple val(meta), path("plots/*.png"), emit: plots
    tuple val(meta), path("plots/*.svg"), emit: plots_svg
    path "versions.yml", emit: versions
    path "eda_stats", emit: stats, optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${python_script} \\
        --vcf ${vcf} \\
        --phenotype ${phenotype_file} \\
        --use-dosage ${use_dosage} \\
        --process-name ${task.process}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir plots
    touch plots/heterozygosity_by_phenotype_samples.png
    touch plots/heterozygosity_by_phenotype_samples.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}
