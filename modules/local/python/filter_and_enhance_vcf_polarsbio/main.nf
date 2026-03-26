process FILTER_AND_ENHANCE_VCF_POLARSBIO {

    tag "$meta.id"
    label 'process_1'

    container 'docker.io/psuszynski/python_tools:1.0.10'

    input:
    tuple val(meta), path(vcf_file), path(tbi_file), path(samples_file), path(python_script)
    val(qual_min)
    val(avg_gq_min)
    val(avg_dp_min)
    val(avg_dp_max)
    val(sample_gq_min)
    val(sample_dp_min)
    val(sample_dp_max)
    val(calc_ds_min_gq)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samples_arg =  samples_file ? "--samples-file samples_transformed.txt" : ""
    """
    if [ -s "${samples_file}" ]; then
        if [ \$(head -1 ${samples_file} | cut -c 1) == '#' ]; then
            cut -f2 ${samples_file} | tail -n +2 > samples_transformed.txt
        else
            ln ${samples_file} samples_transformed.txt
        fi
    fi

    python3 ${python_script} \\
        --input-vcf-path ${vcf_file} \\
        --output-vcf-path ${prefix}_${out_name_part}.vcf.gz \\
        ${samples_arg} \\
        --qual-min ${qual_min} \\
        --avg-gq-min ${avg_gq_min} \\
        --avg-dp-min ${avg_dp_min} \\
        --avg-dp-max ${avg_dp_max} \\
        --sample-gq-min ${sample_gq_min} \\
        --sample-dp-min ${sample_dp_min} \\
        --sample-dp-max ${sample_dp_max} \\
        --calc-ds-min-gq ${calc_ds_min_gq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars-bio: 0.26.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars-bio: 0.26.1
    END_VERSIONS
    """
}