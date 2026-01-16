process PLINK19_MAKESET {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h7b50bb2_6':
        'biocontainers/plink:1.90b6.21--h7b50bb2_6' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(makeset_file), path(tracking_in)

    output:
    tuple val(meta), path("*.set"), emit: out_set
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions
    path "*_plink19_makeset_tracking.json"        , emit: tracking_out

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


    # Extract counts from PLINK log
    samples_in=\$(grep "people.*loaded from" ${prefix}.log | head -1 | awk '{print \$1}' || echo "-1")
    variants_in=\$(grep "variants loaded from" ${prefix}.log | head -1 | awk '{print \$1}' || echo "-1")
    samples_out=\$(grep "variants and .* people pass filters and QC" ${prefix}.log | head -1 | awk '{print \$4}' || echo "-1")
    variants_out=\$(grep "variants and .* people pass filters and QC" ${prefix}.log | head -1 | awk '{print \$1}' || echo "-1")

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

    out_tracking_file_name=\$(echo "${task.process}_plink19_makeset_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
    echo "out_tracking_file_name: \$out_tracking_file_name"

    # Create tracking JSON
    cat <<-END_TRACKING_JSON > \$out_tracking_file_name
    {
        "process_name": "${task.process}_${prefix}",
        "workflow_name": "\$workflow_name",
        "inputs": {
            "variants": \$variants_in,
            "samples": \$samples_in
        },
        "outputs": {
            "variants": \$variants_out,
            "samples": \$samples_out
        },
        "parameters": "$args",
        "predecessor": "\$predecessor"
    }
    END_TRACKING_JSON


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
