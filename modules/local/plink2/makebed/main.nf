process PLINK2_MAKEBED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0':
        params.cpu_support_avx2 ? 'docker.io/psuszynski/plink:2.0-alpha.6.9': 'docker.io/psuszynski/plink:2.0-alpha.6.9.noavx2' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(vcf), path(frq), path(samples_filtering_file), path(variants_filtering_file)
    val(samples_filtering_type)  // for example: '--remove', '--keep'
    val(variants_filtering_type) // for example: '--extract', '--extract-intersect', '--exclude'
    val(out_name_part)
    val(input_args)

    output:
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam"), emit: out_bed_bim_fam
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()

    def bed_input = bed ? "--bed ${bed}" : ""
    def bim_input = bim ? "--bim ${bim}" : ""
    def fam_input = fam ? "--fam ${fam}" : ""
    def vcf_input = vcf ? "--vcf ${vcf}" : ""
    def frq_input = frq ? "--read-freq ${frq}" : ""
    def samples_filtering_input = samples_filtering_file ? "${samples_filtering_type} ${samples_filtering_file}" : ""
    def variants_filtering_input = variants_filtering_file ? "${variants_filtering_type} ${variants_filtering_file}" : ""

    """
    plink2 \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        $args \\
        $input_args \\
        ${bed_input} \\
        ${bim_input} \\
        ${fam_input} \\
        ${vcf_input} \\
        ${frq_input} \\
        ${samples_filtering_input} \\
        ${variants_filtering_input} \\
        --make-bed \\
        --out ${prefix}_${out_name_part}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.bed
    touch ${prefix}_${out_name_part}.bim
    touch ${prefix}_${out_name_part}.fam
    touch ${prefix}_${out_name_part}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
