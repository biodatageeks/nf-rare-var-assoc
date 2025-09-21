process FIX_ZERO_PL {

    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cyvcf2:0.31.1-3':
        'biocontainers/cyvcf2:0.31.1--py39h13a86c0_3' }"

    input:
    tuple val(meta), path(vcf_file)
    val(min_gq)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python
    from cyvcf2 import VCF, Writer
    import sys
    import math
    import numpy as np

    min_gq = float(${min_gq})  # Minimum GQ threshold
    input_vcf = '${vcf_file}'
    output_vcf = '${prefix}_${out_name_part}.vcf.gz'

    vcf = VCF(input_vcf)

    # Add DS field to header
    vcf.add_format_to_header({
        'ID': 'DS',
        'Number': '1',
        'Type': 'Float',
        'Description': 'Genotype dosage (expected number of non-reference alleles)'
    })

    w = Writer(output_vcf, vcf)

    for variant in vcf:
        pl = variant.format('PL', int)  # 2D array: samples x genotypes
        gt = variant.gt_types  # 0 = 0/0, 1 = 0/1, 2 = 1/1, 3 = ./. (missing)
        gq = variant.format('GQ', int)
        
        # Initialize new PL array with existing values
        new_pl = pl.copy() if pl is not None else np.zeros((len(variant.genotypes), 3), dtype=int)
        new_gt = variant.genotypes.copy()  # Copy genotypes to modify if needed
        new_ds = np.array([np.nan for _ in range(len(variant.genotypes))], dtype=float)  # Default to missing

        for i in range(len(variant.genotypes)):
            # Step 1: Fix PL for GT:0/0 and PL:0,0,0
            if pl[i].tolist() == [0, 0, 0] and gt[i] == 0:  # GT:0/0 and PL:0,0,0
                if gq[i][0] >= min_gq:
                    # Calculate likelihoods with L_1/1 = (L_0/1)^2
                    p_wrong = 10 ** (-gq[i][0] / 10.0)
                    
                    # Solve x^2 + x - p_wrong = 0
                    discriminant = 1 + 4 * p_wrong
                    x = (-1 + math.sqrt(discriminant)) / 2  # Positive root
                    l_het = x
                    l_alt = x * x
                    l_ref = 1 - (l_het + l_alt)
                    
                    # Convert to Phred-scaled PL
                    pl_ref = -10 * math.log10(l_ref) if l_ref > 0 else 0
                    pl_het = -10 * math.log10(l_het) if l_het > 0 else 255
                    pl_alt = -10 * math.log10(l_alt) if l_alt > 0 else 255
                    
                    # Update PL for this sample
                    new_pl[i] = [int(round(pl_ref)), int(round(pl_het)), int(round(pl_alt))]
                else:
                    # Set GT to missing (./.)
                    new_gt[i] = [0, 0, False]  # [allele1, allele2, is_phased]; 0,0,False = ./. (unphased missing)
            
            # Step 2: Calculate DS for all samples
            if gq[i][0] < min_gq:
                new_ds[i] = np.nan  #'.'  # Missing genotype
                # print(f"i={i}  gt[i] = {gt[i]}  new_pl[i] = {new_pl[i]}  new_ds[i] = {new_ds[i]}")
            else:
                # Use current PL (corrected or original)
                pl_values = new_pl[i]
                l_ref = 10 ** (-pl_values[0] / 10.0) if pl_values[0] < 255 else 0.0
                l_het = 10 ** (-pl_values[1] / 10.0) if pl_values[1] < 255 else 0.0
                l_alt = 10 ** (-pl_values[2] / 10.0) if pl_values[2] < 255 else 0.0
                l_sum = l_ref + l_het + l_alt
                if l_sum > 0:
                    p_ref = l_ref / l_sum
                    p_het = l_het / l_sum
                    p_alt = l_alt / l_sum
                dosage = p_het + 2 * p_alt
                new_ds[i] = dosage  # np.str_("{:.5f}".format(dosage))
                if new_ds[i] < 0.0001 and new_ds[i] > 0:
                    new_ds[i] = 0.0
                # if p_alt > 0.95:
                #     print(f"i={i}  gt[i] = {gt[i]}  pl_values = {pl_values}  l_ref = {l_ref}  l_het = {l_het}  l_alt = {l_alt}  p_ref = {p_ref}  p_het = {p_het}  p_alt = {p_alt}  dosage = {dosage}  new_ds[i] = {new_ds[i]}")

        # Update PL and GT fields
        if np.any(new_pl != pl):  # Only update if PL changed
            variant.set_format('PL', new_pl)
        if new_gt != variant.genotypes:  # Only update if GT changed
            variant.genotypes = new_gt
        variant.set_format('DS', new_ds)  # Always set DS

        w.write_record(variant)

    w.close()
    vcf.close()

    # Write versions.yml
    with open('versions.yml', 'w') as f:
        f.write('${task.process}:\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    numpy: {np.__version__}\\n')
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}