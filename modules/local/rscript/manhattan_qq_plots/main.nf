
process RSCRIPT_MANHATTAN_QQ_PLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.8'

    input:
    tuple val(phenotype), path(regenie_merged)
    path(gwas_report_template)
    path(r_functions_file)
    tuple val(meta), path(mask_file)
    tuple val(meta), path(phenotype_file)
    path(rmd_pheno_stats_file)

    output:
    path("*.html"), emit: plots_report
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.4.2.5-0.1'
    def annotation_as_string = params.manhattan_annotation_enabled.toString().toUpperCase()

    """
    Rscript -e "require( 'rmarkdown' ); render('${gwas_report_template}', \\
        params = list( \\
          project = '${params.project_name}', \\
          date = '${params.project_date}', \\
          version = '${workflow.manifest.version}', \\
          regenie_out='${regenie_merged}', \\
          regenie_filename='${regenie_merged.baseName}', \\
          phenotype='${phenotype}', \\
          phenotype_file='${phenotype_file}', \\
          plot_ylimit=${params.plot_ylimit}, \\
          manhattan_annotation_enabled = $annotation_as_string, \\
          annotation_min_log10p = ${params.annotation_min_log10p}, \\
          mask_file='${mask_file}', \\
          r_functions='${r_functions_file}', \\
          rmd_pheno_stats='${rmd_pheno_stats_file}' \\
        ), \\
        intermediates_dir='\$PWD', \\
        knit_root_dir='\$PWD', \\
        output_file='${prefix}_manhattan_qq_plots.html' \\
    )"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptmanhattanqqplots: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rscriptmanhattanqqplots: ${VERSION}
    END_VERSIONS
    """
}
