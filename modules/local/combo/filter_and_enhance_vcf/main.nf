process FILTER_AND_ENHANCE_VCF {

    tag "$meta.id"
    label 'process_low'

    container 'docker.io/psuszynski/bioinf_combo:1.3.0'

    input:
    tuple val(meta), path(vcf_file), path(samples_file)
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

    # comp_ds_htslib is a Rust project that has been compiled and built into the bioinf_combo container
    #
    # Usage: comp_ds_htslib [OPTIONS] --input-vcf-path <INPUT_VCF_PATH> --output-vcf-path <OUTPUT_VCF_PATH>
    #
    # Options:
    #     --input-vcf-path <INPUT_VCF_PATH>    
    #     --output-vcf-path <OUTPUT_VCF_PATH>  
    #     --samples-file <SAMPLES_FILE>        
    #     --qual-min <QUAL_MIN>                [default: 25]
    #     --avg-gq-min <AVG_GQ_MIN>            [default: 25]
    #     --avg-dp-min <AVG_DP_MIN>            [default: 25]
    #     --avg-dp-max <AVG_DP_MAX>            [default: 200]
    #     --sample-gq-min <SAMPLE_GQ_MIN>      [default: 20]
    #     --sample-dp-min <SAMPLE_DP_MIN>      [default: 20]
    #     --sample-dp-max <SAMPLE_DP_MAX>      [default: 250]
    #     --calc-ds-min-gq <CALC_DS_MIN_GQ>    [default: 1]
    #     --threads-num <THREADS_NUM>          [default: 1]
    # -h, --help                               Print help
    # -V, --version                            Print version
    
    comp_ds_htslib \\
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
        --calc-ds-min-gq ${calc_ds_min_gq} \\
        --threads-num $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comp_ds_htslib: 1.2.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comp_ds_htslib: 1.2.0
    END_VERSIONS
    """
}