process RENAME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/r-ver:4.4.2.1'

    input:
    tuple val(meta), path(input_path)
    val(output_path)

    output:
    tuple val(meta), path(output_path), emit: output

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cp ${input_path} ${output_path}
    """
}
