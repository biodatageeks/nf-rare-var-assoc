process CALCULATE_F_OUTLIERS {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(het_file)
    val(range_stds)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.txt"), emit: outliers
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python
    import sys
    import numpy as np
    import pandas as pd

    # Read the first line to extract column names
    with open('${het_file}', 'r') as file:
        header_line = file.readline().strip()  # Read the first line and remove trailing newline
        # Remove the '#' and split into column names based on tabs
        columns = header_line.lstrip('#').strip().split('\t')

    data = pd.read_csv('${het_file}', sep='\t', skiprows=1, names=columns)

    f_values = data['F']
    mean_f = np.mean(f_values)
    std_f = np.std(f_values)
    lower_bound = mean_f - ${range_stds} * std_f
    upper_bound = mean_f + ${range_stds} * std_f
    print(f"mean F: {mean_f}  std F: {std_f}  lower bound: {lower_bound}  upper bound: {upper_bound}")

    if 'FID' in data.columns:
        ids_cols = ['FID', 'IID']
    else:
        ids_cols = ['IID', 'IID']

    outliers = data[(data['F'] < lower_bound) | (data['F'] > upper_bound)][ids_cols]
    outliers.to_csv('${prefix}_${out_name_part}.txt', sep='\\t', index=False, header=False)

    # Write versions.yml
    with open('versions.yml', 'w') as f:
        f.write('${task.process}:\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    numpy: {np.__version__}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}