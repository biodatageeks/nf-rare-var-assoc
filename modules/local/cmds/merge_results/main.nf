process MERGE_RESULTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.9'

    input:
    path(csv_concat_py_script)
    tuple val(phenotype), path(regenie_chromosomes)

    output:
    tuple val(phenotype), path ("${phenotype}.regenie.gz"), emit: results_merged

    script:
    """
    python3 ${csv_concat_py_script} --input_sep ' ' --output_sep '\t' --output ${phenotype}.regenie.tmp.gz --inputs ${regenie_chromosomes}
    zcat ${phenotype}.regenie.tmp.gz | awk 'NR<=1{print \$0;next}{print \$0| "sort -n -k1 -k2 -T \$PWD"}' | bgzip -c > ${phenotype}.regenie.gz
    """

}
