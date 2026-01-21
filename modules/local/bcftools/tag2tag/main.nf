process BCFTOOLS_TAG2TAG {
    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(index)
    val(tag_from)
    val(tag_to)

    output:
    tuple val(meta), path("*_tag2tag.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*_tag2tag*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*_tag2tag*.csi")                    , emit: csi, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    # Check if FMT/${tag_from} and FMT/${tag_to} are present in the VCF header

    if [ -n "${tag_from}" ]; then
        FMT_TAG_FROM_FIELD=\$(bcftools view -h "${vcf}" | grep '^##FORMAT=<ID=' | grep -c -e '${tag_from}' || true)
    else
        FMT_TAG_FROM_FIELD=0
    fi

    if [ -n "${tag_to}" ]; then
        FMT_TAG_TO_FIELD=\$(bcftools view -h "${vcf}" | grep '^##FORMAT=<ID=' | grep -c -e '${tag_to}' || true)
    else
        FMT_TAG_TO_FIELD=0
    fi

    # Ensure single integer output (handle empty or invalid output)
    FMT_TAG_FROM_FIELD=\${FMT_TAG_FROM_FIELD:-0}
    FMT_TAG_TO_FIELD=\${FMT_TAG_TO_FIELD:-0}

    echo "FMT/${tag_from}, present: \$FMT_TAG_FROM_FIELD"
    echo "FMT/${tag_to}, present: \$FMT_TAG_TO_FIELD"
    
    # Apply filters with bcftools if tag_from is present and tag_to is NOT present
    if [ "\$FMT_TAG_FROM_FIELD" -gt 0 ] && [ "\$FMT_TAG_TO_FIELD" -eq 0 ]; then
        bcftools +tag2tag \\
            ${vcf} \\
            --output ${prefix}_tag2tag.${extension} \\
            --threads $task.cpus \\
            $args \\
            -- --${tag_from}-to-${tag_to}
    else
        echo "Error: FMT/${tag_from} is not present in the VCF header or FMT/${tag_to} is present in the VCF header"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    def stub_index = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
                     args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
                     args.contains("--write-index")     || args.contains("-W") ? "csi" :
                     ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && stub_index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${stub_index}" : ""
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
