process CHECK_X_CHROM_PRESENT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), env(has_x), path(pgen), path(pvar), path(psam), emit: has_x

    script:
    """
    # Check if chromosome X or chrX exists in the .pvar file
    if grep -E "^(X|chrX)\\s" ${pvar}; then
        echo "true" > has_x.txt
    else
        echo "false" > has_x.txt
    fi
    has_x=\$(cat has_x.txt)
    """
}