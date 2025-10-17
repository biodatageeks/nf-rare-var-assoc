include { RSCRIPT_MANHATTAN_QQ_PLOTS } from '../../../modules/local/rscript/manhattan_qq_plots'

workflow REPORTING {

    take:
    results_merged
    masks_file
    phenotype_file
    pca_plot_file
    eda_plots
    setlist_file

    main:
    gwas_report_template = Channel.fromPath("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/gene_level_report_template.Rmd", checkIfExists: true)
    r_functions_file = Channel.fromPath("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/functions.R", checkIfExists: true)

    if(!params.phenotypes_apply_rint) {
        rmd_pheno_stats_file = Channel.fromPath("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/child_phenostatistics.Rmd", checkIfExists: true)
    } else {
        rmd_pheno_stats_file = Channel.fromPath("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/child_phenostatistics_rint.Rmd", checkIfExists: true)
    }

    RSCRIPT_MANHATTAN_QQ_PLOTS (
        phenotype_file
            .join(pca_plot_file, by: 0)
            .join(eda_plots, by: 0)
            .join(setlist_file, by: 0)
            .join(results_merged, by: 0)
            .combine(masks_file)
            .combine(gwas_report_template)
            .combine(r_functions_file)
            .combine(rmd_pheno_stats_file)
    )
}


