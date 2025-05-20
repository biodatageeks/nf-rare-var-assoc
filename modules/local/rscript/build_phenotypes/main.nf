process RSCRIPT_BUILD_PHENOTYPES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    path(r_script_ch)
    tuple val(meta), path(fam), path(controls), path(cases)
    val(input_args)

    output:
    tuple val(meta), path("*_phenotype.txt"), emit: phenotype
    tuple val(meta), path("*_r_out.fam"), emit: out_fam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.4.2-0.2'
    """
    Rscript \\
        ${r_script_ch} \\
        --fam-path ${fam} \\
        --controls-path ${controls} \\
        --cases-paths ${cases} \\
        --out-fam-path ${prefix}_r_out.fam \\
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
    touch ${prefix}_r_out.fam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """
}
