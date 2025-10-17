process REGENIE_STEP2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/rgcgithub/regenie/regenie:v4.1.gz'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam), path(phenotype), path(annotations), path(setlist), path(aaf), path(step1_pred_list), path(covar_file), path(masks), path(tracking_in)
    val(input_args)

    output:
    tuple val(meta), path("*_masks.bed"), path("*_masks.bim"), path("*_masks.fam"), emit: masks_bed_bim_fam
    tuple val(meta), path("*_masks.snplist"), emit: masks_snplist
    tuple val(meta), path("*.regenie"), emit: regenie_out
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions
    path "*_regenie_step2_tracking.json"        , emit: tracking_out

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
       --step 2 \\
       --threads ${task.cpus} \\
       --pgen ${pgen.baseName} \\
       --phenoFile ${phenotype} \\
       ${covar_file_input} \\
       --anno-file ${annotations} \\
       --set-list ${setlist} \\
       --mask-def ${masks} \\
       --pred ${step1_pred_list} \\
       --aaf-file ${aaf} \\
       --out ${prefix}_step2 \\
       $args $input_args


    # Extract counts from PLINK log
    samples_in=\$(grep "number of individuals used in analysis =" ${prefix}_step2.log | head -1 | awk '{print \$NF}' || echo "-1")
    variants_in=\$(grep "bgen file .* with .* samples and .* variants" ${prefix}_step2.log | head -1 | grep -oP '\\d+(?=\\s*variants)' || echo "-1")
   
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

    out_tracking_file_name=\$(echo "${task.process}_regenie_step2_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
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
        "parameters": "$args $input_args --step 2 ${covar_file_input}",
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
    touch ${prefix}_masks.bed
    touch ${prefix}_masks.bim
    touch ${prefix}_masks.fam
    touch ${prefix}_masks.snplist
    touch ${prefix}_Y1.regenie

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regenie: \$(regenie --version)
    END_VERSIONS
    """
}
