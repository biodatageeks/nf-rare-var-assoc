process MERGE_SEX_COVAR {
    tag "$meta.id"
    label 'process_1'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(covar_file), path(psam_file)

    output:
    tuple val(meta), path("*_sex_covar.txt"), emit: covar_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python
    from pathlib import Path
    import sys
    import csv

    covar_path = Path('${covar_file}')
    psam_path = Path('${psam_file}')
    output_path = Path('${prefix}_sex_covar.txt')

    sex_by_sample = {}
    with psam_path.open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\\t')
        header = next(reader, None)
        if header is None:
            raise SystemExit('Empty PSAM file')
        try:
            fid_index = header.index('#FID')
        except ValueError:
            fid_index = header.index('FID')
        iid_index = header.index('IID')
        sex_index = header.index('SEX')

        for row in reader:
            if not row:
                continue
            sex_by_sample[(row[fid_index], row[iid_index])] = row[sex_index]

    with covar_path.open(newline='') as handle, output_path.open('w', newline='') as output_handle:
        reader = csv.reader(handle, delimiter='\\t')
        writer = csv.writer(output_handle, delimiter='\\t', lineterminator='\\n')
        header = next(reader, None)
        if header is None:
            raise SystemExit('Empty covariate file')
        writer.writerow(header + ['SEX'])

        for row in reader:
            if not row:
                continue
            key = (row[0], row[1])
            writer.writerow(row + [sex_by_sample.get(key, '')])


    # Write versions.yml
    with open('versions.yml', 'w') as f:
        f.write('${task.process}:\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sex_covar.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS    
    """
}
