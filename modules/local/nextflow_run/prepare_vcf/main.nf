process PREPARE_VCF {
    tag "$meta.id"
    label 'process_1'

    // Invokes nf-prepare-vcf as a child nextflow run on the same node (no container).
    // nextflow must be on PATH in the host environment.
    // Profile forwarding: workflow.profile (e.g. "low_resources,podman") is forwarded
    // verbatim to the child so container settings and resource limits are consistent.
    // Resume: no -resume flag; child always runs fresh when parent cache misses.
    // publish_intermediate=true in prep.yml is required so tracking JSONs are published
    // from BCFTOOLS_NORM and PLINK2_MAKEPGEN to results/bcftools_norm/ and
    // results/plink2_makepgen/ respectively.

    input:
    tuple val(meta), path(vcf)
    path(params_file)
    path(add_config)

    output:
    tuple val(meta), path("results/bcftools_reheader/*_reheader.vcf.gz"),    emit: prepared_vcf
    tuple val(meta), path("results/bcftools_index/*_reheader.vcf.gz.tbi"),   emit: prepared_vcf_tbi
    path("results/**/*tracking*.json"),                                        emit: tracking
    path("results/pipeline_info/nf-prepare-vcf_software_versions.yml"),       emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def child_pipeline = "${projectDir}/../nf-prepare-vcf/main.nf"
    def profile_arg    = workflow.profile ? "-profile ${workflow.profile}" : ''
    def add_config_arg = add_config ? "-c ${add_config}" : ''
    """
    nextflow run ${child_pipeline} \\
        ${profile_arg} \\
        -params-file ${params_file} \\
        ${add_config_arg} \\
        --input_vcf ${vcf} \\
        --outdir results \\
        --cpu_support_avx2 ${params.cpu_support_avx2} \\
        -work-dir \${PWD}/child_work \\
        -ansi-log false
    """
}
