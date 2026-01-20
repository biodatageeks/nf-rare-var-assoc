process BCFTOOLS_VCFTOAAF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    val(tag_name)           // e.g., "AF_nfe" - primary tag to extract
    val(default_tag_name)   // e.g., "AF" - fallback tag when primary is missing

    output:
    tuple val(meta), path("*_aaf.tsv"), emit: aaf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract allele frequencies from VCF using bcftools query and awk
    # This replaces the slower R-based vcf2aaf.R script
    #
    # Logic:
    # 1. Extract CHROM, POS, REF, ALT, and the two AF tags from INFO field
    # 2. Create position identifier: chr_pos_ref_alt (stripping 'chr' prefix)
    # 3. Use coalesce logic: prefer tag_name, fallback to default_tag_name, then "0"

    # Create awk script that extracts AF from INFO field
    # This handles the case where tag_name may not exist in the VCF
    # Compatible with BusyBox awk (no capturing groups in match)
    cat > process_aaf.awk << 'AWKEOF'
BEGIN { 
    OFS="\\t"
}
function extract_tag(info, tag,    i, n, parts, kv) {
    # Split INFO by semicolons and find exact tag match
    n = split(info, parts, ";")
    for (i = 1; i <= n; i++) {
        # Check if this part starts with "TAG="
        if (index(parts[i], tag "=") == 1) {
            # Extract value after the "="
            kv = parts[i]
            sub("^" tag "=", "", kv)
            # Only return if it looks like a number
            if (kv ~ /^[0-9.]+\$/) {
                return kv
            }
        }
    }
    return ""
}
{
    chrom = \$1
    gsub(/^[Cc][Hh][Rr]/, "", chrom)
    pos = chrom "_" \$2 "_" \$3 "_" \$4
    info = \$5
    
    af = ""
    
    # Try primary tag first, then fallback
    af = extract_tag(info, tag1)
    if (af == "") {
        af = extract_tag(info, tag2)
    }
    
    # Default to 0 if no AF found
    if (af == "" || af == ".") {
        af = "0"
    }
    
    print pos, af
}
AWKEOF
    
    # Extract CHROM, POS, REF, ALT, INFO and process with awk
    # Pass tag names as awk variables using -v
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO\\n' ${vcf} | \\
        awk -v tag1="${tag_name}" -v tag2="${default_tag_name}" -f process_aaf.awk > ${prefix}_aaf.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_aaf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
