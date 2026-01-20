#!/usr/bin/env bash
#
# Benchmark Runner Script
# =======================
# Runs benchmark comparisons between R, bcftools, and Python implementations
# of vcf2aaf and assign_annotations processes.
#
# Usage:
#   ./run_benchmark.sh small      # Test with 17 MB VEP-annotated file
#   ./run_benchmark.sh large      # Test with 1.8 GB VEP-annotated file  
#   ./run_benchmark.sh both       # Test with both files
#   ./run_benchmark.sh small podman   # Use podman containers
#   ./run_benchmark.sh small conda    # Use conda environments
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Parse arguments
TEST_SIZE="${1:-small}"
RUNTIME="${2:-docker}"

# Validate arguments
if [[ ! "$TEST_SIZE" =~ ^(small|large|both)$ ]]; then
    echo "Error: Invalid test size '$TEST_SIZE'. Use: small, large, or both"
    exit 1
fi

if [[ ! "$RUNTIME" =~ ^(docker|podman|singularity|conda|local)$ ]]; then
    echo "Error: Invalid runtime '$RUNTIME'. Use: docker, podman, singularity, conda, or local"
    exit 1
fi

echo "========================================"
echo "  Benchmark Implementation Comparison"
echo "========================================"
echo ""
echo "Test size: $TEST_SIZE"
echo "Runtime:   $RUNTIME"
echo ""

# Check if required files exist
echo "Checking input files..."
if [[ "$TEST_SIZE" == "small" ]] || [[ "$TEST_SIZE" == "both" ]]; then
    SMALL_FILE="assets/test_sim_chr22_vep.vcf.gz"
    if [[ ! -f "$SMALL_FILE" ]]; then
        echo "Error: Small test file not found: $SMALL_FILE"
        exit 1
    fi
    echo "  ✓ Small file found ($(du -h "$SMALL_FILE" | cut -f1))"
fi

if [[ "$TEST_SIZE" == "large" ]] || [[ "$TEST_SIZE" == "both" ]]; then
    LARGE_FILE="assets/prepare_vep_medium.vcf.gz"
    if [[ ! -f "$LARGE_FILE" ]]; then
        echo "Error: Large test file not found: $LARGE_FILE"
        exit 1
    fi
    echo "  ✓ Large file found ($(du -h "$LARGE_FILE" | cut -f1))"
fi

# Check required scripts exist
echo ""
echo "Checking implementation scripts..."
SCRIPTS=(
    "modules/local/rscript/vcf2aaf/assets/vcf2aaf.R"
    "modules/local/python/vcf2aaf/assets/vcf2aaf.py"
    "modules/local/rscript/assign_annotations/assets/annotate.R"
    "modules/local/bcftools/assign_annotations/assets/assign_annotations.py"
    "modules/local/bcftools/assign_annotations/assets/assign_annotations.sh"
)

for script in "${SCRIPTS[@]}"; do
    if [[ -f "$script" ]]; then
        echo "  ✓ $(basename "$script")"
    else
        echo "  ✗ Missing: $script"
        exit 1
    fi
done

echo ""
echo "Starting benchmark..."
echo "========================================"
echo ""

# Build the nextflow command
NF_CMD="nextflow run workflows/benchmark_implementations.nf"
NF_CMD="$NF_CMD -c conf/benchmark.config"
NF_CMD="$NF_CMD -profile test_${TEST_SIZE},${RUNTIME}"
NF_CMD="$NF_CMD -resume"

# Add additional options for local execution
if [[ "$RUNTIME" == "local" ]]; then
    NF_CMD="$NF_CMD -profile local"
fi

echo "Running: $NF_CMD"
echo ""

# Execute
eval "$NF_CMD"

# Show results summary
OUTDIR="benchmark_results_${TEST_SIZE}"
if [[ -f "${OUTDIR}/final_benchmark_report.txt" ]]; then
    echo ""
    echo "========================================"
    echo "  Benchmark Complete!"
    echo "========================================"
    echo ""
    echo "Results saved to: ${OUTDIR}/"
    echo ""
    echo "--- Summary ---"
    cat "${OUTDIR}/final_benchmark_report.txt"
fi
