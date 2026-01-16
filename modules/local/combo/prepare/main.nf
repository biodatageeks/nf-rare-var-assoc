process PREPARE {

    tag "$meta.id"
    label 'process_low'

    container 'docker.io/psuszynski/bioinf_combo:1.1.1'

    input:
    tuple val(meta), path(vcf_in), path(samples), path(cachesubdir), path(rename_chrs), path(fix_zero_pl_script_path), path(tracking_in)
    val(sed_arg)
    val(min_gq)
    val(vep_fasta_subpath)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*_${out_name_part}*.tbi"), emit: tbi
    path "versions.yml", emit: versions
    path "*_${out_name_part}_tracking.json", emit: tracking_out

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_gz = vcf_in.extension.endsWith("gz") ? vcf_in : "${vcf_in.name}.gz"
    def rename_chrs_file = rename_chrs ? "--rename-chrs ${rename_chrs}" : ''
    def samples_file =  samples ? "--samples-file samples_transformed.txt" : ""

    """
    echo "args: ${args}"
    echo "meta.id: ${meta.id}"
    date
    echo ""


    echo "BCFTOOLS_REPLACE_SAMPLE_NAMES"
    date
    echo ""

    # Check if the input is uncompressed VCF
    if [[ "${vcf_in}" =~ \\.vcf\$ ]]; then
        bgzip -c ${vcf_in} > ${vcf_gz}
    fi

    # Extract sample names
    bcftools query -l ${vcf_in} > samples.txt
    
    # For sed_arg="s/_/-/g"  replace underscores with specified character
    sed '${sed_arg}' samples.txt > new_samples.txt
    
    # Update VCF with new sample names
    bcftools reheader -s new_samples.txt ${vcf_gz} \\
       --threads $task.cpus \\
       -o vcf_with_sample_names_corrected.vcf.gz


    echo "BCFTOOLS_INDEX"
    date
    echo ""

    bcftools index -t \\
        --threads $task.cpus \\
        vcf_with_sample_names_corrected.vcf.gz
    

    echo "BCFTOOLS_NORM"
    date
    echo ""

    bcftools norm \\
        --fasta-ref ${cachesubdir}/${vep_fasta_subpath} \\
        --output vcf_norm.vcf.gz \\
        --rm-dup exact -m -any --do-not-normalize --output-type z --write-index=tbi \\
        --threads $task.cpus \\
        vcf_with_sample_names_corrected.vcf.gz


    echo "BCFTOOLS_ANNOTATE"
    date
    echo ""

    bcftools annotate \\
        --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' --output-type z --write-index=tbi \\
        $rename_chrs_file \\
        --output vcf_annotate.vcf.gz \\
        --threads $task.cpus \\
        vcf_norm.vcf.gz
    

    if [[ "${params.use_dosage}" = true\$ ]]; then
        echo "FIX_ZERO_PL"
        date
        echo ""

        python3 ${fix_zero_pl_script_path} \\
            --input-vcf-path vcf_annotate.vcf.gz \\
            --output-vcf-path vcf_fixed_zero_pl.vcf.gz \\
            --min-gq ${min_gq} \\
            --threads-num $task.cpus
    else
        cp vcf_annotate.vcf.gz vcf_fixed_zero_pl.vcf.gz
    fi
    

    echo "BCFTOOLS_INDEX"
    date
    echo ""

    bcftools index -t \\
        --threads $task.cpus \\
        vcf_fixed_zero_pl.vcf.gz
    

    echo "BCFTOOLS_VIEW"
    date
    echo ""

    if [ -s "${samples}" ]; then
        if [ \$(head -1 ${samples} | cut -c 1) == '#' ]; then
            cut -f2 ${samples} | tail -n +2 > samples_transformed.txt
        else
            ln ${samples} samples_transformed.txt
        fi
    fi

    bcftools view \\
        --output ${prefix}_${out_name_part}.vcf.gz \\
        ${samples_file} \\
        --output-type z --write-index=tbi \\
        --threads $task.cpus \\
        vcf_fixed_zero_pl.vcf.gz
    
    if [[ "${vcf_in}" =~ \\.vcf\$ ]]; then
        rm ${vcf_gz}
    fi
    rm vcf_with_sample_names_corrected.vcf.gz
    rm vcf_with_sample_names_corrected.vcf.gz.tbi
    rm vcf_annotate.vcf.gz
    rm vcf_annotate.vcf.gz.tbi
    rm vcf_norm.vcf.gz
    rm vcf_norm.vcf.gz.tbi
    rm vcf_fixed_zero_pl.vcf.gz
    rm vcf_fixed_zero_pl.vcf.gz.tbi

    echo "computing tracking json"
    date
    echo ""

    bcftools stats ${vcf_in} > ${prefix}_${out_name_part}_vcf_gz_in_stats.txt
    bcftools stats ${prefix}_${out_name_part}.vcf.gz > ${prefix}_${out_name_part}_vcf_gz_out_stats.txt


    # Extract counts from stats files
    samples_in=\$(grep "number of samples:" ${prefix}_${out_name_part}_vcf_gz_in_stats.txt | cut -f4)
    variants_in=\$(grep "number of records:" ${prefix}_${out_name_part}_vcf_gz_in_stats.txt | cut -f4)
    samples_out=\$(grep "number of samples:" ${prefix}_${out_name_part}_vcf_gz_out_stats.txt | cut -f4)
    variants_out=\$(grep "number of records:" ${prefix}_${out_name_part}_vcf_gz_out_stats.txt | cut -f4)
    echo "samples_in: \$samples_in"
    echo "variants_in: \$variants_in"
    echo "samples_out: \$samples_out"
    echo "variants_out: \$variants_out"

    echo "tracking_in: ${tracking_in}"
    predecessor="none"
    if [ -s "${tracking_in}" ]; then
        predecessor=\$(grep '"process_name"' ${tracking_in} | sed 's/.*"process_name": "\\([^"]*\\)".*/\\1/')
    fi
    echo "predecessor: \$predecessor"

    workflow_name=\$(echo "${task.process}" | awk -F: '{print \$(NF-1)}')
    echo "workflow_name: \$workflow_name"

    out_tracking_file_name=\$(echo "${task.process}_${prefix}_${out_name_part}_tracking.json" | sed 's/[^:]*://' | sed 's/:/_/g')
    echo "out_tracking_file_name: \$out_tracking_file_name"

    # Create tracking JSON
    cat <<-END_TRACKING_JSON > \$out_tracking_file_name
    {
        "process_name": "${task.process}_${prefix}_${out_name_part}",
        "workflow_name": "\$workflow_name",
        "inputs": {
            "variants": \$variants_in,
            "samples": \$samples_in
        },
        "outputs": {
            "variants": \$variants_out,
            "samples": \$samples_out
        },
        "parameters": "reheader -s new_samples.txt  |  --rm-dup exact -m -any --do-not-normalize  |  --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${rename_chrs_file}  |  ${fix_zero_pl_script_path} --min-gq ${min_gq}  |  ${samples_file}",
        "predecessor": "\$predecessor"
    }
    END_TRACKING_JSON


    date
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        python: \$(python3 --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}
