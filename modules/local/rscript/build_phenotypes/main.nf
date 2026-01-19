process RSCRIPT_BUILD_PHENOTYPES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    tuple val(meta), path(controls), path(cases), path(r_script_ch)
    val(input_args)

    output:
    tuple val(meta), path("*_phenotype.txt"), emit: phenotype
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.4.2-0.2'
    """
    Rscript \\
        ${r_script_ch} \\
        --controls-path ${controls} \\
        --cases-paths ${cases} \\
        --out-pheno-path ${prefix}_phenotype.txt \\
        $input_args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_phenotype.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """
}
