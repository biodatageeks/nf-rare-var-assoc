#!/usr/bin/env bash
set -euo pipefail

export DATA=/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison
export RVA_REPO=/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc
export REFDIR=/data/doktorat/biodatageeks/genome_in_a_bottle/reference
export CHRS="12 22"
export VERSIONS="full"
export CLEANUP=false
export SCORE=false
export THREADS=1

SCRIPT="$RVA_REPO/docs/tool-comparison/chained-benchmark/run_chain_ricopili_staar.sh"

# dataset pairs, in the same order as your original list
PAIRS=(
  "0 2 30 20 17"
  "3 4 31 21 25"
  "5 7 28 23 14"
  "8 9 29 18 15"
  "10 11 26 19 24"
  "12 13 27 22 16"
)

pids=()

for pair in "${PAIRS[@]}"; do
  tag=$(echo "$pair" | tr ' ' '_')
  LOG=~/ricopili_staar_chr12_22__ds_${tag}.log
  echo "Launching DATASET_IDXS=\"$pair\" -> $LOG"
  DATASET_IDXS="$pair" bash "$SCRIPT" > "$LOG" 2>&1 &
  pids+=($!)
done

echo "Launched ${#pids[@]} jobs: ${pids[*]}"
echo "Waiting for all to finish..."
wait "${pids[@]}"
echo "All jobs completed."

