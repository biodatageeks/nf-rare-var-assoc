process REGENIE_STEP1 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/rgcgithub/regenie/regenie:v4.1.gz'

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(qc_pass_id), path(qc_pass_snplist), path(phenotype), path(covar_file)
    val(input_args)
    path(tracking_in)

    output:
    tuple val(meta), path("*.loco"), emit: loco
    tuple val(meta), path("*_pred.list"), emit: pred_list
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions
    path "*_regenie_step1_tracking.json"        , emit: tracking_out

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def covar_file_corr_format = covar_file + "_correct_format.txt"
    def covar_file_input = covar_file ? "--covarFile " + covar_file_corr_format : ""
    """
    if [ -n "${covar_file}" ]; then
        sed '1s/^#//' ${covar_file} > ${covar_file_corr_format}
    fi
    regenie \\
       --step 1 \\
       --threads ${task.cpus} \\
       --bed ${bed.baseName} \\
       --keep ${qc_pass_id} \\
       --extract ${qc_pass_snplist} \\
       --phenoFile ${phenotype} \\
       ${covar_file_input} \\
       --out ${prefix}_step1 \\
       $args $input_args


    # Extract counts from PLINK log
    samples_in=\$(grep "n_samples =" ${prefix}_step1.log | head -1 | awk '{print \$NF}' || echo "-1")
    variants_in=\$(grep "n_snps =" ${prefix}_step1.log | head -1 | awk '{print \$NF}' || echo "-1")
   
    echo "samples_in: \$samples_in"
    echo "variants_in: \$variants_in"

    echo "tracking_in: ${tracking_in}"
    predecessor="none"
    if [ -s "${tracking_in}" ]; then
        predecessor=\$(grep '"process_name"' ${tracking_in} | sed 's/.*"process_name": "\\([^"]*\\)".*/\\1/')
    fi
    echo "predecessor: \$predecessor"

    workflow_name=\$(echo "${task.process}" | awk -F: '{print \$(NF-1)}')
    echo "workflow_name: \$workflow_name"

    out_tracking_file_name=\$(echo "${task.process}_regenie_step1_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
    echo "out_tracking_file_name: \$out_tracking_file_name"

    # Create tracking JSON
    cat <<-END_TRACKING_JSON > \$out_tracking_file_name
    {
        "process_name": "${task.process}",
        "workflow_name": "\$workflow_name",
        "inputs": {
            "variants": \$variants_in,
            "samples": \$samples_in
        },
        "parameters": "$args $input_args --step 1 --keep ${qc_pass_id} --extract ${qc_pass_snplist} ${covar_file_input}",
        "predecessor": "\$predecessor"
    }
    END_TRACKING_JSON


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.loco
    touch ${prefix}_pred.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """
}
