process DOWNLOAD_FILE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(url)
    val(extension)

    output:
    tuple val(meta), path("*.${extension}"), emit: output_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wget -O ${prefix}.${extension} ${url}
    """
}
