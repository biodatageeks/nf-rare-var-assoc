include { RSCRIPT_MANHATTAN_QQ_PLOTS } from '../../../modules/local/rscript/manhattan_qq_plots'

workflow REPORTING {

    take:
    results_merged
    masks_file
    phenotype_file
    ch_plot_file

    main:
    gwas_report_template = file("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/gene_level_report_template.Rmd", checkIfExists: true)
    r_functions_file = file("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/functions.R", checkIfExists: true)

    if(!params.phenotypes_apply_rint) {
        rmd_pheno_stats_file = file("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/child_phenostatistics.Rmd", checkIfExists: true)
    } else {
        rmd_pheno_stats_file = file("$baseDir/modules/local/rscript/manhattan_qq_plots/assets/child_phenostatistics_rint.Rmd", checkIfExists: true)
    }

    RSCRIPT_MANHATTAN_QQ_PLOTS (
        results_merged,
        gwas_report_template,
        r_functions_file,
        masks_file,
        phenotype_file,
        ch_plot_file,
        rmd_pheno_stats_file
    )
}


