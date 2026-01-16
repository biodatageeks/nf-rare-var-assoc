process PLINK2_EXPORT_OTHER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.0.0a.6.9--h9948957_0':
        params.cpu_support_avx2 ? 'docker.io/psuszynski/plink:2.0-alpha.6.9': 'docker.io/psuszynski/plink:2.0-alpha.6.9.noavx2' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(psam)
    val(extension)
    val(input_args)

    output:
    tuple val(meta), path("*.${extension}"), emit: out_file
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()

    //def extract_ids_from_psam = psam ? "tail -n +2 ${psam} | cut -f 2 | tr '\\n' ' ' | sed 's/ \$//' > sample_ids.txt": ""
    def extract_ids_from_psam = psam ? "(echo '#IID'; tail -n +2 ${psam} | cut -f 2) > sample_ids.txt": ""
    def vcf_input = vcf ? "--vcf ${vcf}" : ""
    def psam_input = psam ? "--update-sex ${psam} --keep ${psam}" : ""  // sample_ids.txt" : ""

    
    """
    ${extract_ids_from_psam}

    plink2 \\
        --threads ${task.cpus} \\
        --memory $mem_mb \\
        ${psam_input} \\
        ${vcf_input} \\
        $args $input_args \\
        --out ${prefix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(echo \$(plink2 --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """
}
