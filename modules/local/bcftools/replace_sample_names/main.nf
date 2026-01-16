process BCFTOOLS_REPLACE_SAMPLE_NAMES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf_in)
    val(sed_arg)
    val(out_name_part)

    output:
    tuple val(meta), path("*_${out_name_part}.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_gz = vcf_in.extension.endsWith("gz") ? vcf_in : "${vcf_in.name}.gz"

    """
    # Check if the input is uncompressed VCF
    if [[ "${vcf_in}" =~ \\.vcf\$ ]]; then
        bgzip -c ${vcf_in} > ${vcf_gz}
    fi

    # Extract sample names
    bcftools query -l ${vcf_in} > samples.txt
    
    # Replace underscores with specified character
    sed '${sed_arg}' samples.txt > new_samples.txt
    
    # Update VCF with new sample names
    bcftools reheader -s new_samples.txt ${vcf_gz} \\
       $args \\
       --threads $task.cpus \\
       -o ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
