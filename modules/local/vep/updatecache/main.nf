process VEP_UPDATECACHE {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:release_113.4':
        'docker.io/psuszynski/ensembl-vep:113.4.3' }"

    input:
    tuple val(meta), path(vep_cache)
    val(species)
    val(input_args)
    val(vep_cache_url)
    tuple val(ref_fasta_url), val(vep_fasta_path)

    output:
    path("${vep_cache}/${species}"), emit: cachesubdir
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -d "${vep_cache}/$species" ]; then
        if [ -z "${vep_cache_url}" ]; then
            perl /opt/vep/src/ensembl-vep/INSTALL.pl --CACHEDIR ${vep_cache} --SPECIES $species $args $input_args
            perl /opt/vep/src/ensembl-vep/convert_cache.pl --dir ${vep_cache} --species $species --version all
        else
            wget -O cache_file.tar.gz ${vep_cache_url}
            tar -xzf cache_file.tar.gz -C ${vep_cache}/
        fi
    fi
    if [ ! -z "${ref_fasta_url}" ]; then
        if [ ! -f "${vep_cache}/${species}/${vep_fasta_path}" ]; then
            wget -O ${vep_cache}/${species}/${vep_fasta_path} ${ref_fasta_url}
            gunzip -c ${vep_cache}/${species}/${vep_fasta_path} | bgzip --index --index-name ${vep_cache}/${species}/${vep_fasta_path}.gzi > ${vep_cache}/${species}/${vep_fasta_path}_2
            mv ${vep_cache}/${species}/${vep_fasta_path}_2 ${vep_cache}/${species}/${vep_fasta_path}
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
