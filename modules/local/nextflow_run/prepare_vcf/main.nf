process PREPARE_VCF {
    tag "$meta.id"
    label 'process_9'

    // Invokes nf-prepare-vcf as a child nextflow run on the same node (no container).
    // nextflow must be on PATH in the host environment.
    // Profile forwarding: workflow.profile (e.g. "low_resources,podman") is forwarded
    // verbatim to the child so container settings and resource limits are consistent.
    // Resume: no -resume flag; child always runs fresh when parent cache misses.
    // publish_intermediate=true in prep.yml is required so tracking JSONs are published
    // from BCFTOOLS_NORM and PLINK2_MAKEPGEN to child_results/bcftools_norm/ and
    // child_results/plink2_makepgen/ respectively.

    input:
    tuple val(meta), path(vcf), path(ref_fasta), path(ref_fasta_fai)
    path(params_file)

    output:
    tuple val(meta), path("child_results/bcftools_reheader/*_reheader.vcf.gz"),    emit: prepared_vcf
    tuple val(meta), path("child_results/bcftools_index/*_reheader.vcf.gz.tbi"),   emit: prepared_vcf_tbi
    path("child_results/**/*tracking*.json"),                                      emit: tracking, optional: true
    path("child_results/pipeline_info/nf-prepare-vcf_software_versions.yml"),      emit: versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def child_pipeline = "${projectDir}/../nf-prepare-vcf/main.nf"
    def ref_fasta_arg  = ref_fasta ? "--input_ref_fasta ${ref_fasta}" : ''
    def profile_arg    = workflow.profile ? "-profile ${workflow.profile}" : ''
    // Values given on a nextflow command line are always strings, and the string "false"
    // is truthy in Groovy -- so --cpu_support_avx2 false would leave the child on the AVX2
    // images. Hand it over in a config file instead, where it stays a real boolean.
    """
    echo 'params.cpu_support_avx2 = ${params.cpu_support_avx2}' > child_avx2.config

    nextflow run ${child_pipeline} \\
        ${profile_arg} \\
        -c child_avx2.config \\
        -params-file ${params_file} \\
        --input_vcf ${vcf} \\
        ${ref_fasta_arg} \\
        --outdir child_results \\
        -work-dir \${PWD}/child_work \\
        -ansi-log false
    """
}
