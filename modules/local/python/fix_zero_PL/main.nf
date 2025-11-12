process FIX_ZERO_PL {

    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/bioinf_combo:1.1.1'

    input:
    tuple val(meta), path(vcf_file)
    val(min_gq)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # comp_ds_htslib is a Rust project that has been compiled and built into the bioinf_combo container
    
    comp_ds_htslib \\
        --input-vcf-path ${vcf_file} \\
        --output-vcf-path ${prefix}_${out_name_part}.vcf.gz \\
        --min-gq ${min_gq} \\
        --threads-num $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        python: \$(python3 --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}