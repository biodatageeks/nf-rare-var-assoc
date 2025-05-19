process PLINK19_MAKESET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h7b50bb2_6':
        'biocontainers/plink:1.90b6.21--h7b50bb2_6' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(makeset_file)

    output:
    tuple val(meta), path("*.set"), emit: out_set
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()

    def bed_input = bed ? "--bed ${bed}" : ""
    def bim_input = bim ? "--bim ${bim}" : ""
    def fam_input = fam ? "--fam ${fam}" : ""
    """
    plink \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        $args \\
        ${bed_input} \\
        ${bim_input} \\
        ${fam_input} \\
        --make-set ${makeset_file} \\
        --write-set \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink19: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.set
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink19: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
