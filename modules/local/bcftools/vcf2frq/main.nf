process BCFTOOLS_VCF2FRQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.frq"), emit: frq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools query -f '%CHROM\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/AC\\n' $vcf | \\
    awk 'BEGIN {print "CHR\\tSNP\\tA1\\tA2\\tMAF\\tNCHROBS"} \\
        {maf=\$5; a1=\$3; a2=\$4; if (\$5 > 0.5) {maf=1-\$5; a1=\$4; a2=\$3} \\
        nchrom=(maf > 0 ? int(2*\$6/maf + 0.5) : 0); if (\$2==".") \$2=\$1"_"\$3"_"\$4; \\
        print \$1"\\t"\$2"\\t"a1"\\t"a2"\\t"maf"\\t"nchrom}' > ${prefix}.frq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.frq \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
