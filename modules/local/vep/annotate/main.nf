
process VEP_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:release_113.4':
        'docker.io/psuszynski/ensembl-vep:113.4.3' }"

    input:
    tuple val(meta), path(input_vcf), path(input_vcf_tbi), path(vep_cache)
    val(species)
    val(fasta_path)
    val(input_args)

    output:
    tuple val(meta), path("*_vep.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vep \\
        --dir ${vep_cache}/.. \\
        --fasta ${vep_cache}/${fasta_path} \\
        -i ${input_vcf} \\
        -o ${prefix}_vep.vcf.gz \\
        --fork ${task.cpus} \\
        --species $species \\
        $args $input_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
