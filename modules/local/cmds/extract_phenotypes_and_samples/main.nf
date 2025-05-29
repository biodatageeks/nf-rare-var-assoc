// Process to extract unique phenotypes and corresponding sample lists
process EXTRACT_PHENOTYPES_AND_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(pheno_file)

    output:
    tuple val(meta), path("*_samples_*.txt"), emit: samples

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract unique phenotype values (excluding header)
    cut -f 3 ${pheno_file} | sort | uniq | grep -v "Y1" > pheno_values.txt

    # Loop through phenotypes to create sample lists
    while read pheno; do
        awk -F'\\t' -v p="\$pheno" '\$3 == p {print \$1 "\\t" \$2}' ${pheno_file} > ${prefix}_samples_\${pheno}.txt
    done < pheno_values.txt
    """
}