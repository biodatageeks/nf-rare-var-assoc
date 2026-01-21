
process RSCRIPT_MANHATTAN_QQ_PLOTS {
    tag "$meta.id"
    label 'process_4'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    tuple val(meta), path(phenotype_file), path(pc_plot_file), path(eda_plots), path(setlist_file), val(phenotype), path(regenie_merged), path(mask_file), path(gwas_report_template), path(r_functions_file), path(rmd_pheno_stats_file)

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
    def all_params = params.toString().replaceAll("'", "\\\\'")

    """
    gene_count=\$(cat ${setlist_file} | wc -l)
    echo "gene count: \$gene_count"
    echo "PWD: \$PWD"

    Rscript -e "require( 'rmarkdown' ); render('${gwas_report_template}', \\
        params = list( \\
          project = '${params.project_name}', \\
          date = '${params.project_date}', \\
          version = '${workflow.manifest.version}', \\
          regenie_out='${regenie_merged}', \\
          regenie_filename='${regenie_merged.baseName}', \\
          phenotype='${phenotype}', \\
          phenotype_file='${phenotype_file}', \\
          pc_plot_file='${pc_plot_file}', \\
          eda_plots='${eda_plots}', \\
          plot_ylimit=${params.plot_ylimit}, \\
          genes_num=\$gene_count, \\
          manhattan_annotation_enabled = ${annotation_as_string}, \\
          annotation_min_log10p = ${params.annotation_min_log10p}, \\
          mask_file='${mask_file}', \\
          r_functions='${r_functions_file}', \\
          rmd_pheno_stats='${rmd_pheno_stats_file}', \\
          all_params='${all_params}' \\
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
