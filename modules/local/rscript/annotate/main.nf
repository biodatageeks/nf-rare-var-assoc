process RSCRIPT_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    path(r_script_ch)
    tuple val(meta), path(vcf)
    tuple val(meta), path(bed)
    tuple val(meta), path(bim)
    tuple val(meta), path(fam)
    tuple val(meta), path(bgen)
    tuple val(meta), path(bgen_bgi)
    tuple val(meta), path(sample)
    tuple val(meta), path(controls)
    tuple val(meta), path(cases)
    val(input_args)

    output:
    tuple val(meta), path("*_phenotype.txt"), emit: phenotype
    tuple val(meta), path("*.annotations"), emit: annotations
    tuple val(meta), path("*.setlist"), emit: setlist
    tuple val(meta), path("*_r_out.fam"), emit: out_fam
    tuple val(meta), path("*_r_out.sample"), emit: out_sample
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
        --vcf-path ${vcf} \\
        --sample-path ${sample} \\
        --out-fam-path ${prefix}_r_out.fam \\
        --out-sample-path ${prefix}_r_out.sample \\
        --out-pheno-path ${prefix}_phenotype.txt \\
        --out-anno-path ${prefix}.annotations \\
        --out-setlist-path ${prefix}.setlist \\
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
    touch ${prefix}.annotations
    touch ${prefix}.masks
    touch ${prefix}.setlist
    touch ${prefix}_r_out.fam
    touch ${prefix}_r_out.sample

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptannotate: ${VERSION}
    END_VERSIONS
    """
}
