process PLINK2_PROJECTION_SCORE {
    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.0.0a.6.9--h9948957_0':
        params.cpu_support_avx2 ? 'docker.io/psuszynski/plink:2.0-alpha.6.9': 'docker.io/psuszynski/plink:2.0-alpha.6.9.noavx2' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam), path(frq), path(eigenvec_allele), path(tracking_in)
    val(score_options)
    val(out_name_part)
    val(input_args)

    output:
    tuple val(meta), path("*.sscore"), emit: sscore
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions
    path "*_${out_name_part}_tracking.json"        , emit: tracking_out

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()

    def pgen_input = pgen ? "--pgen ${pgen}" : ""
    def pvar_input = pvar ? "--pvar ${pvar}" : ""
    def psam_input = psam ? "--psam ${psam}" : ""

    """
    plink2 \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        $args $input_args \\
        ${pgen_input} \\
        ${pvar_input} \\
        ${psam_input} \\
        --read-freq ${frq} \\
        --score ${eigenvec_allele} ${score_options} \\
        --out ${prefix}_${out_name_part}


    # Extract counts from PLINK log
    samples_in=\$(grep "samples.*loaded from" ${prefix}_${out_name_part}.log | head -1 | awk '{print \$1}' || echo "-1")
    variants_in=\$(grep "variants loaded from" ${prefix}_${out_name_part}.log | head -1 | awk '{print \$1}' || echo "-1")
    samples_out=\$(grep "samples (.*) remaining" ${prefix}_${out_name_part}.log | tail -1 | awk '{print \$1}' || echo "-1")
    variants_out=\$(grep "variants remaining after" ${prefix}_${out_name_part}.log | tail -1 | awk '{print \$1}' || echo "-1")
   
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

    out_tracking_file_name=\$(echo "${task.process}_${prefix}_${out_name_part}_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
    echo "out_tracking_file_name: \$out_tracking_file_name"

    # Create tracking JSON
    cat <<-END_TRACKING_JSON > \$out_tracking_file_name
    {
        "process_name": "${task.process}_${prefix}_${out_name_part}",
        "workflow_name": "\$workflow_name",
        "inputs": {
            "variants": \$variants_in,
            "samples": \$samples_in
        },
        "outputs": {
            "variants": \$variants_out,
            "samples": \$samples_out
        },
        "parameters": "$args $input_args --score ${eigenvec_allele} ${score_options}",
        "predecessor": "\$predecessor"
    }
    END_TRACKING_JSON


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.sscore
    touch ${prefix}_${out_name_part}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
