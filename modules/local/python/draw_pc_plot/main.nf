process DRAW_PC_PLOT {

    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seaborn:0.13.2':
        'biocontainers/seaborn:0.13.2' }"

    input:
    tuple val(meta), path(sscore_file), path(pheno_file)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.png"), emit: plot_file
    tuple val(meta), path("*_${out_name_part}.svg"), emit: plot_file_svg
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python
    import sys
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.rcParams.update({'font.size': 16})

    pca = pd.read_table("${sscore_file}", sep="\t")
    pheno = pd.read_table("${pheno_file}", sep="\t")
    pheno_cols = [c for c in pheno.columns if c not in ['FID', 'IID']]
    first_pheno_col = pheno_cols[0]

    merged_data = pd.merge(pca, pheno, on='IID', how='left')

    plt.figure(figsize=(9, 6))
    plt.title(f'Plot of first two principal components')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    sns.scatterplot(x='PC1_AVG', y='PC2_AVG', hue=first_pheno_col, data=merged_data, s=65)
    plt.savefig('${prefix}_${out_name_part}.png')
    plt.savefig('${prefix}_${out_name_part}.svg', format="svg")
    plt.close()

    # Write versions.yml
    with open('versions.yml', 'w') as f:
        f.write('${task.process}:\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    seaborn: {sns.__version__}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}