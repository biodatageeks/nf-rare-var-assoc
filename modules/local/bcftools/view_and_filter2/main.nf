process BCFTOOLS_VIEW_AND_FILTER2 {
    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(index), path(regions), path(targets), path(samples), path(snplist), path(tracking_in)
    val(input_args_filter_1)
    val(input_args_filter_2)
    val(out_name_part)

    output:
    tuple val(meta), path("*_${out_name_part}.{vcf,vcf.gz,bcf,bcf.gz}"),     emit: vcf
    tuple val(meta), path("*_${out_name_part}.{vcf,vcf.gz,bcf,bcf.gz}.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*_${out_name_part}.{vcf,vcf.gz,bcf,bcf.gz}.csi"), emit: csi, optional: true
    path "versions.yml"                               , emit: versions
    path "*_${out_name_part}_tracking.json"           , emit: tracking_out

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file  = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file samples_transformed.txt" : ""
    def snplist_file =  snplist ? "-i 'ID=@snplist_transformed.txt'" : ""

    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    echo "args:"
    echo ${args}

    if [ -s "${samples}" ]; then
        if [ \$(head -1 ${samples} | cut -c 1) == '#' ]; then
            cut -f2 ${samples} | tail -n +2 > samples_transformed.txt
        else
            ln ${samples} samples_transformed.txt
        fi
    fi
    if [ -s "${snplist}" ]; then
        sed 's/^chr//' ${snplist} > snplist_transformed.txt
    fi

    bcftools view \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        ${snplist_file} \\
        --threads $task.cpus \\
        ${vcf} | \\
    bcftools filter \\
        $input_args_filter_1 \\
        --threads $task.cpus | \\
    bcftools filter \\
        --output ${prefix}_${out_name_part}.${extension} \\
        $args $input_args_filter_2 \\
        --threads $task.cpus
    

    bcftools stats ${vcf} > ${prefix}_${out_name_part}_${extension}_in_stats.txt
    bcftools stats ${prefix}_${out_name_part}.${extension} > ${prefix}_${out_name_part}_${extension}_out_stats.txt


    # Extract counts from stats files
    samples_in=\$(grep "number of samples:" ${prefix}_${out_name_part}_${extension}_in_stats.txt | cut -f4)
    variants_in=\$(grep "number of records:" ${prefix}_${out_name_part}_${extension}_in_stats.txt | cut -f4)
    samples_out=\$(grep "number of samples:" ${prefix}_${out_name_part}_${extension}_out_stats.txt | cut -f4)
    variants_out=\$(grep "number of records:" ${prefix}_${out_name_part}_${extension}_out_stats.txt | cut -f4)
    echo "samples_in: \$samples_in"
    echo "variants_in: \$variants_in"
    echo "samples_out: \$samples_out"
    echo "variants_out: \$variants_out"

    echo "tracking_in: ${tracking_in}"
    predecessor="none"
    if [ -s "${tracking_in}" ]; then
        predecessor=\$(grep '"process_name"' ${tracking_in} | sed 's/.*"process_name": "\\([^"]*\\)".*/\\1/' || true)
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
        "parameters": "${regions_file} ${targets_file} ${samples_file} ${snplist_file} $args",
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
