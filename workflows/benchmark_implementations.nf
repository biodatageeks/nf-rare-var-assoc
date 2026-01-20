/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BENCHMARK WORKFLOW - Compare R vs bcftools vs Python implementations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow benchmarks the performance of:
    - vcf2aaf: R (original) vs bcftools vs polars-bio
    - assign_annotations: R (original) vs bcftools+polars vs pure shell
    
    Performance metrics (time, memory, CPU) are captured by Nextflow's trace feature.
    See benchmark_results/pipeline_info/benchmark_trace.txt for detailed metrics.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// Import all implementations
include { RSCRIPT_VCFTOAAF                    } from '../modules/local/rscript/vcf2aaf'
include { BCFTOOLS_VCFTOAAF                   } from '../modules/local/bcftools/vcf2aaf'
include { PYTHON_VCFTOAAF                     } from '../modules/local/python/vcf2aaf'
include { RSCRIPT_ASSIGN_ANNOTATIONS          } from '../modules/local/rscript/assign_annotations'
include { BCFTOOLS_ASSIGN_ANNOTATIONS         } from '../modules/local/bcftools/assign_annotations'
include { BCFTOOLS_ASSIGN_ANNOTATIONS_SHELL   } from '../modules/local/bcftools/assign_annotations/main_shell'
include { BCFTOOLS_INDEX                      } from '../modules/local/bcftools/index'


process COMPARE_OUTPUTS {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/benchmark", mode: 'copy'
    
    conda "conda-forge::python=3.11"
    container 'docker.io/psuszynski/bioinf_combo:1.3.0'

    input:
    tuple val(meta), 
          path(r_aaf, stageAs: 'r_aaf.tsv'), 
          path(bcftools_aaf, stageAs: 'bcftools_aaf.tsv'), 
          path(python_aaf, stageAs: 'python_aaf.tsv'),
          path(r_anno, stageAs: 'r_anno.txt'), 
          path(bcftools_anno, stageAs: 'bcftools_anno.txt'), 
          path(shell_anno, stageAs: 'shell_anno.txt'),
          path(r_setlist, stageAs: 'r_setlist.txt'), 
          path(bcftools_setlist, stageAs: 'bcftools_setlist.txt'), 
          path(shell_setlist, stageAs: 'shell_setlist.txt')

    output:
    tuple val(meta), path("*.comparison_report.txt"), emit: report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 << 'PYEOF'
import os
import re
from pathlib import Path
from collections import Counter

def load_tsv(filepath):
    # Load TSV file and return list of rows (each row is a list of columns).
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\\n')
            if line:
                rows.append(line.split('\\t'))
    return rows

def compare_files_detailed(file1, file2, name1="File1", name2="File2", tolerance=0.0001):
    # Comprehensive comparison of two TSV files.
    # Returns a detailed report as a list of strings.
    report = []
    
    try:
        rows1 = load_tsv(file1)
        rows2 = load_tsv(file2)
    except Exception as e:
        return [f"  ERROR loading files: {e}"], False
    
    # 1. Compare number of lines
    report.append(f"  Line count: {name1}={len(rows1)}, {name2}={len(rows2)}")
    if len(rows1) != len(rows2):
        report.append(f"    ⚠ Line counts differ by {abs(len(rows1) - len(rows2))}")
    
    if not rows1 or not rows2:
        if not rows1 and not rows2:
            return report + ["  ✓ Both files are empty"], True
        return report + ["  ✗ One file is empty"], False
    
    # 2. Compare number of columns (check first few rows)
    cols1 = [len(r) for r in rows1[:min(10, len(rows1))]]
    cols2 = [len(r) for r in rows2[:min(10, len(rows2))]]
    col_count1 = Counter(cols1).most_common(1)[0][0]
    col_count2 = Counter(cols2).most_common(1)[0][0]
    
    report.append(f"  Column count (typical): {name1}={col_count1}, {name2}={col_count2}")
    if col_count1 != col_count2:
        report.append(f"    ⚠ Column counts differ")
    
    # 3. Get first column (key column) from both files
    keys1 = set(r[0] for r in rows1 if r)
    keys2 = set(r[0] for r in rows2 if r)
    
    common_keys = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    report.append(f"  Keys (first column): {name1}={len(keys1)}, {name2}={len(keys2)}, common={len(common_keys)}")
    
    if only_in_1:
        examples = list(only_in_1)[:5]
        report.append(f"    ⚠ {len(only_in_1)} keys only in {name1}: {examples}{'...' if len(only_in_1) > 5 else ''}")
    if only_in_2:
        examples = list(only_in_2)[:5]
        report.append(f"    ⚠ {len(only_in_2)} keys only in {name2}: {examples}{'...' if len(only_in_2) > 5 else ''}")
    
    # 4. Check if rows match when sorted by first column
    dict1 = {r[0]: r for r in rows1 if r}
    dict2 = {r[0]: r for r in rows2 if r}
    
    # 5. Compare values for common keys
    value_diffs = []
    exact_matches = 0
    numeric_matches = 0
    
    for key in sorted(common_keys):
        row1 = dict1[key]
        row2 = dict2[key]
        
        if row1 == row2:
            exact_matches += 1
            continue
        
        # Check column by column
        max_cols = max(len(row1), len(row2))
        row_diffs = []
        all_numeric_match = True
        
        for col_idx in range(max_cols):
            val1 = row1[col_idx] if col_idx < len(row1) else "<missing>"
            val2 = row2[col_idx] if col_idx < len(row2) else "<missing>"
            
            if val1 == val2:
                continue
            
            # Try numeric comparison
            try:
                num1, num2 = float(val1), float(val2)
                if abs(num1 - num2) <= tolerance:
                    continue  # Close enough
                row_diffs.append((col_idx, val1, val2, "numeric"))
                all_numeric_match = False
            except ValueError:
                row_diffs.append((col_idx, val1, val2, "string"))
                all_numeric_match = False
        
        if not row_diffs:
            numeric_matches += 1
        else:
            value_diffs.append((key, row_diffs))
    
    report.append(f"  Row comparison (common keys):")
    report.append(f"    Exact matches: {exact_matches}")
    report.append(f"    Numeric matches (within tolerance): {numeric_matches}")
    report.append(f"    Differences: {len(value_diffs)}")
    
    # 6. Show sample differences
    if value_diffs:
        report.append(f"  Sample differences (first 10):")
        for key, diffs in value_diffs[:10]:
            report.append(f"    Key '{key}':")
            for col_idx, val1, val2, diff_type in diffs[:3]:
                report.append(f"      Col {col_idx}: '{val1}' vs '{val2}' ({diff_type})")
            if len(diffs) > 3:
                report.append(f"      ... and {len(diffs) - 3} more column differences")
    
    # 7. Overall verdict
    is_match = (len(only_in_1) == 0 and len(only_in_2) == 0 and len(value_diffs) == 0)
    
    if is_match:
        report.append(f"  ✓ FILES MATCH")
    else:
        report.append(f"  ✗ FILES DIFFER")
        if only_in_1 or only_in_2:
            report.append(f"    - Key differences: {len(only_in_1)} only in {name1}, {len(only_in_2)} only in {name2}")
        if value_diffs:
            report.append(f"    - Value differences in {len(value_diffs)} rows")
    
    return report, is_match

# File mappings - using staged names
aaf_files = {
    'r_script': ('r_aaf.tsv', 'R'),
    'bcftools': ('bcftools_aaf.tsv', 'bcftools'),
    'python_polars': ('python_aaf.tsv', 'python_polars')
}

anno_files = {
    'r_script': ('r_anno.txt', 'R'),
    'bcftools_polars': ('bcftools_anno.txt', 'bcftools_polars'),
    'shell': ('shell_anno.txt', 'shell')
}

setlist_files = {
    'r_script': ('r_setlist.txt', 'R'),
    'bcftools_polars': ('bcftools_setlist.txt', 'bcftools_polars'),
    'shell': ('shell_setlist.txt', 'shell')
}

# Generate report
report_lines = []
report_lines.append("=" * 80)
report_lines.append(f"BENCHMARK COMPARISON REPORT: ${meta.id}")
report_lines.append("=" * 80)
report_lines.append("")
report_lines.append("Performance metrics are available in:")
report_lines.append("  - benchmark_results/pipeline_info/benchmark_trace.txt")
report_lines.append("  - benchmark_results/pipeline_info/benchmark_report.html")
report_lines.append("  - benchmark_results/pipeline_info/benchmark_timeline.html")
report_lines.append("")

# Compare AAF outputs
report_lines.append("=" * 80)
report_lines.append("VCF2AAF OUTPUT COMPARISON")
report_lines.append("=" * 80)

for impl in ['bcftools', 'python_polars']:
    report_lines.append("")
    report_lines.append(f"--- R vs {impl} ---")
    comparison, is_match = compare_files_detailed(
        aaf_files['r_script'][0], 
        aaf_files[impl][0],
        "R", impl
    )
    report_lines.extend(comparison)

# Compare Annotations outputs
report_lines.append("")
report_lines.append("=" * 80)
report_lines.append("ASSIGN ANNOTATIONS OUTPUT COMPARISON")
report_lines.append("=" * 80)

for impl in ['bcftools_polars', 'shell']:
    report_lines.append("")
    report_lines.append(f"--- R vs {impl} (annotations) ---")
    comparison, is_match = compare_files_detailed(
        anno_files['r_script'][0], 
        anno_files[impl][0],
        "R", impl
    )
    report_lines.extend(comparison)
    
    report_lines.append("")
    report_lines.append(f"--- R vs {impl} (setlist) ---")
    comparison, is_match = compare_files_detailed(
        setlist_files['r_script'][0], 
        setlist_files[impl][0],
        "R", impl
    )
    report_lines.extend(comparison)

report_lines.append("")
report_lines.append("=" * 80)

# Write report
with open('${prefix}.comparison_report.txt', 'w') as f:
    f.write('\\n'.join(report_lines))

print('\\n'.join(report_lines))
PYEOF
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BENCHMARK {
    take:
    ch_vcf_files  // channel: [ meta, vcf ]

    main:
    ch_versions = Channel.empty()
    
    // Script paths
    ch_r_vcf2aaf_script = Channel.fromPath("${projectDir}/../modules/local/rscript/vcf2aaf/assets/vcf2aaf.R", checkIfExists: true)
    ch_python_vcf2aaf_script = Channel.fromPath("${projectDir}/../modules/local/python/vcf2aaf/assets/vcf2aaf.py", checkIfExists: true)
    ch_r_annotate_script = Channel.fromPath("${projectDir}/../modules/local/rscript/assign_annotations/assets/annotate.R", checkIfExists: true)
    ch_python_annotate_script = Channel.fromPath("${projectDir}/../modules/local/bcftools/assign_annotations/assets/assign_annotations.py", checkIfExists: true)
    ch_shell_annotate_script = Channel.fromPath("${projectDir}/../modules/local/bcftools/assign_annotations/assets/assign_annotations.sh", checkIfExists: true)
    ch_masks = Channel.fromPath("${projectDir}/../assets/default.masks", checkIfExists: true)

    // Index the VCF files first (not counted in benchmark time)
    BCFTOOLS_INDEX(ch_vcf_files)
    ch_vcf_with_index = ch_vcf_files
        .join(BCFTOOLS_INDEX.out.tbi, by: 0)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    
    // For processes that don't need the index, just use vcf
    ch_vcf_only = ch_vcf_with_index.map { meta, vcf, tbi -> tuple(meta, vcf) }

    // VCF2AAF benchmarks - run all 3 implementations
    // R version
    RSCRIPT_VCFTOAAF(
        ch_vcf_only.combine(ch_r_vcf2aaf_script),
        params.vcf2aaf_tag_name + " " + params.vcf2aaf_default_tag
    )
    
    // bcftools version
    BCFTOOLS_VCFTOAAF(
        ch_vcf_only,
        params.vcf2aaf_tag_name,
        params.vcf2aaf_default_tag
    )
    
    // Python/polars-bio version
    PYTHON_VCFTOAAF(
        ch_vcf_only.combine(ch_python_vcf2aaf_script),
        params.vcf2aaf_tag_name,
        params.vcf2aaf_default_tag
    )

    // Assign annotations benchmarks - run all 3 implementations
    // R version
    RSCRIPT_ASSIGN_ANNOTATIONS(
        ch_vcf_only.combine(ch_masks).combine(ch_r_annotate_script),
        params.annotate_options
    )
    
    // bcftools+polars version
    BCFTOOLS_ASSIGN_ANNOTATIONS(
        ch_vcf_only.combine(ch_masks).combine(ch_python_annotate_script),
        params.annotate_options
    )
    
    // Pure shell version
    BCFTOOLS_ASSIGN_ANNOTATIONS_SHELL(
        ch_vcf_only.combine(ch_masks).combine(ch_shell_annotate_script),
        params.annotate_options
    )

    // Collect all outputs for comparison
    ch_comparison_input = RSCRIPT_VCFTOAAF.out.aaf
        .join(BCFTOOLS_VCFTOAAF.out.aaf, by: 0)
        .join(PYTHON_VCFTOAAF.out.aaf, by: 0)
        .join(RSCRIPT_ASSIGN_ANNOTATIONS.out.annotations, by: 0)
        .join(BCFTOOLS_ASSIGN_ANNOTATIONS.out.annotations, by: 0)
        .join(BCFTOOLS_ASSIGN_ANNOTATIONS_SHELL.out.annotations, by: 0)
        .join(RSCRIPT_ASSIGN_ANNOTATIONS.out.setlist, by: 0)
        .join(BCFTOOLS_ASSIGN_ANNOTATIONS.out.setlist, by: 0)
        .join(BCFTOOLS_ASSIGN_ANNOTATIONS_SHELL.out.setlist, by: 0)
        .map { meta, r_aaf, bc_aaf, py_aaf, r_anno, bc_anno, sh_anno, r_set, bc_set, sh_set ->
            tuple(meta, r_aaf, bc_aaf, py_aaf, r_anno, bc_anno, sh_anno, r_set, bc_set, sh_set)
        }

    COMPARE_OUTPUTS(ch_comparison_input)

    emit:
    reports = COMPARE_OUTPUTS.out.report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY POINT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create channel from input files
    ch_input_vcf = Channel.fromPath(params.input_vcf, checkIfExists: true)
        .map { vcf -> 
            def basename = vcf.baseName.replaceAll(/\.vcf$/, '').replaceAll(/\.gz$/, '')
            tuple([id: basename], vcf)
        }

    BENCHMARK(ch_input_vcf)
}
