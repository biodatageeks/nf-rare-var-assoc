process PLINK2_IMPORT_DOSAGE {
    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.0.0a.6.9--h9948957_0':
        params.cpu_support_avx2 ? 'docker.io/psuszynski/plink:2.0-alpha.6.9': 'docker.io/psuszynski/plink:2.0-alpha.6.9.noavx2' }"

    input:
    tuple val(meta), path(psam), path(traw)
    val(out_name_part)
    val(input_args)

    output:
    tuple val(meta), path("*_${out_name_part}.pgen"), path("*_${out_name_part}.pvar"), path("*_${out_name_part}.psam"), emit: out_pgen_pvar_psam
    tuple val(meta), path("*_${out_name_part}.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()

    def psam_input = psam ? "--psam ${psam}" : ""

    """
    plink2 \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        ${psam_input} \\
        --import-dosage ${traw} \\
        $args $input_args \\
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
    touch ${prefix}_${out_name_part}.pgen
    touch ${prefix}_${out_name_part}.pvar
    touch ${prefix}_${out_name_part}.psam
    touch ${prefix}_${out_name_part}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
