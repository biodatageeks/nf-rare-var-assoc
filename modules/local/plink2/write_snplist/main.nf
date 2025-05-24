
process PLINK2_WRITE_SNPLIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0':
        params.cpu_support_avx2 ? 'docker.io/psuszynski/plink:2.0-alpha.6.9': 'docker.io/psuszynski/plink:2.0-alpha.6.9.noavx2' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(keep)
    val(out_name_part)
    val(input_args)

    output:
    tuple val(meta), path("*.snplist"), emit: snplist
    tuple val(meta), path("*${out_name_part}.id"), emit: id, optional: true
    tuple val(meta), path("*${out_name_part}.mindrem.id"), emit: mindremid, optional: true
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    def keep_samples = keep ? "--keep ${keep}" : ""
    """
    plink2 \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        $args \\
        $input_args \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        ${keep_samples} \\
        --write-samples --write-snplist \\
        --out ${prefix}_${out_name_part}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.snplist
    touch ${prefix}_${out_name_part}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
