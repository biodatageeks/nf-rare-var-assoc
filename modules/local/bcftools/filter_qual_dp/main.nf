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
    path(tracking_in)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
    path "versions.yml"                               , emit: versions
    path "*_filter_qual_dp_tracking.json"             , emit: tracking_out

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
    if [ -n "${info_filter_ensure_field_present}" ]; then
        INFO_FIELD=\$(bcftools view -h "${vcf}" | grep '^##INFO=<ID=' | grep -c -e '${info_filter_ensure_field_present}' || true)
    else
        INFO_FIELD=0
    fi
    if [ -n "${fmt_filter_ensure_field_present}" ]; then
        FMT_FIELD=\$(bcftools view -h "${vcf}" | grep '^##FORMAT=<ID=' | grep -c -e '${fmt_filter_ensure_field_present}' || true)
    else
        FMT_FIELD=0
    fi

    # Ensure single integer output (handle empty or invalid output)
    INFO_FIELD=\${INFO_FIELD:-0}
    FMT_FIELD=\${FMT_FIELD:-0}

    echo "INFO/${info_filter_ensure_field_present}, present: \$INFO_FIELD"
    echo "FMT/${fmt_filter_ensure_field_present} present: \$FMT_FIELD"

    # Initialize filter expression with QUAL
    FILTER="${qual_filter}"

    # Add INFO_FIELD filters if present
    if [ "\$INFO_FIELD" -gt 0 ]; then
        if [ -n "\$FILTER" ]; then
            FILTER="\$FILTER && ${info_filter}"
        else
            FILTER="${info_filter}"
        fi
    fi

    # Add FMT_FIELD filters if present
    if [ "\$FMT_FIELD" -gt 0 ]; then
        if [ -n "\$FILTER" ]; then
            FILTER="\$FILTER && ${fmt_filter}"
        else
            FILTER="${fmt_filter}"
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

    bcftools stats ${vcf} > ${prefix}_filter_qual_dp_${extension}_in_stats.txt
    bcftools stats ${prefix}_filter_qual_dp.${extension} > ${prefix}_filter_qual_dp_${extension}_out_stats.txt


    # Extract counts from stats files
    samples_in=\$(grep "number of samples:" ${prefix}_filter_qual_dp_${extension}_in_stats.txt | cut -f4)
    variants_in=\$(grep "number of records:" ${prefix}_filter_qual_dp_${extension}_in_stats.txt | cut -f4)
    samples_out=\$(grep "number of samples:" ${prefix}_filter_qual_dp_${extension}_out_stats.txt | cut -f4)
    variants_out=\$(grep "number of records:" ${prefix}_filter_qual_dp_${extension}_out_stats.txt | cut -f4)
    echo "samples_in: \$samples_in"
    echo "variants_in: \$variants_in"
    echo "samples_out: \$samples_out"
    echo "variants_out: \$variants_out"

    echo "tracking_in: ${tracking_in}"
    predecessor="none"
    if [ -s "${tracking_in}" ]; then
        predecessor=\$(grep '"process_name"' ${tracking_in} | sed 's/.*"process_name": "\\([^"]*\\)".*/\\1/')
    fi
    echo "predecessor: \$predecessor"

    workflow_name=\$(echo "${task.process}" | awk -F: '{print \$(NF-1)}')
    echo "workflow_name: \$workflow_name"

    out_tracking_file_name=\$(echo "${task.process}_${prefix}_filter_qual_dp_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
    echo "out_tracking_file_name: \$out_tracking_file_name"

    # Create tracking JSON
    cat <<-END_TRACKING_JSON > \$out_tracking_file_name
    {
        "process_name": "${task.process}_filter_qual_dp_${prefix}",
        "workflow_name": "\$workflow_name",
        "inputs": {
            "variants": \$variants_in,
            "samples": \$samples_in
        },
        "outputs": {
            "variants": \$variants_out,
            "samples": \$samples_out
        },
        "parameters": "$args --include '\$FILTER'",
        "predecessor": "\$predecessor"
    }
    END_TRACKING_JSON

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
