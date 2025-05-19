process BCFTOOLS_FILTER_QUAL_DP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(index)
    val(qual_filter)
    val(info_filter)
    val(fmt_filter)
    val(info_filter_ensure_field_present)
    val(fmt_filter_ensure_field_present)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
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
    # Check if INFO/${info_filter_ensure_field_present} and FMT/${fmt_filter_ensure_field_present} are present in the VCF header
    INFO_FIELD=\$(bcftools view -h "${vcf}" | grep -c '^##INFO=<ID=${info_filter_ensure_field_present},' || true)
    FMT_FIELD=\$(bcftools view -h "${vcf}" | grep -c '^##FORMAT=<ID=${fmt_filter_ensure_field_present},' || true)

    # Ensure single integer output (handle empty or invalid output)
    INFO_FIELD=\${INFO_FIELD:-0}
    FMT_FIELD=\${FMT_FIELD:-0}

    echo "INFO/${info_filter_ensure_field_present}, present: \$INFO_FIELD"
    echo "FMT/${fmt_filter_ensure_field_present} present: \$FMT_FIELD"

    # Initialize filter expression with QUAL
    FILTER="${qual_filter}"

    # Add INFO_FIELD filters if present
    if [ "\$INFO_FIELD" -gt 0 ]; then
        if [ "\$FILTER" = "" ];then
            FILTER="${info_filter}"
        else
            FILTER="\$FILTER && ${info_filter}"
        fi
    fi

    # Add FMT_FIELD filters if present
    if [ "\$FMT_FIELD" -gt 0 ]; then
        if [ "\$FILTER" = "" ];then
            FILTER="${fmt_filter}"
        else
            FILTER="\$FILTER && ${fmt_filter}"
        fi
    fi

    # Debug: Print final filter expression
    echo "Applying filter: \$FILTER"
    
    # Apply filters with bcftools
    bcftools view \\
        --output ${prefix}_filter_qual_dp.${extension} \\
        $args \\
        --include "\$FILTER" \\
        --threads $task.cpus \\
        ${vcf}

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
