#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_BASE="${ROOT_DIR}/tests/results/initial"
TRUE_OUTDIR="${OUT_BASE}/use_dosage_true"
FALSE_OUTDIR="${OUT_BASE}/use_dosage_false"

mkdir -p "${OUT_BASE}"
rm -rf "${TRUE_OUTDIR}" "${FALSE_OUTDIR}" \
       "${OUT_BASE}/work_true" "${OUT_BASE}/work_false" \
       "${OUT_BASE}/regenie_hashes_true.txt" "${OUT_BASE}/regenie_hashes_false.txt"

run_pipeline() {
    local use_dosage="$1"
    local outdir="$2"
    local workdir="$3"

    echo "Running pipeline with --use_dosage=${use_dosage}"
    nextflow run "${ROOT_DIR}/main.nf" \
        -profile podman,test_skip_preparation_and_reporting \
        -c "${ROOT_DIR}/tests/nextflow.config" \
        -w "${workdir}" \
        --use_dosage "${use_dosage}" \
        --outdir "${outdir}" \
        --project_name "initial_${use_dosage}"
}

assert_pipeline_info_exists() {
    local outdir="$1"

    compgen -G "${outdir}/pipeline_info/execution_report_*" > /dev/null
    compgen -G "${outdir}/pipeline_info/execution_trace_*" > /dev/null
    compgen -G "${outdir}/pipeline_info/execution_timeline_*" > /dev/null
}

latest_trace_file() {
    local outdir="$1"
    ls -1t "${outdir}"/pipeline_info/execution_trace_* 2>/dev/null | head -n 1
}

assert_required_processes_completed() {
    local outdir="$1"
    local trace_file
    trace_file="$(latest_trace_file "${outdir}")"

    if [[ -z "${trace_file}" || ! -f "${trace_file}" ]]; then
        echo "Missing execution trace in ${outdir}/pipeline_info" >&2
        exit 1
    fi

    # Require a minimal set of core processes to complete for the skip-prep/skip-report branch.
    local required=(
        "REGENIE_STEP1"
        "REGENIE_STEP2"
        "GENERATE_TRACKING_REPORT"
    )

    local p
    for p in "${required[@]}"; do
        if ! awk -F '\t' -v proc="${p}" '
            NR==1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "name") {
                        name_col = i
                    } else if ($i == "status") {
                        status_col = i
                    }
                }
                next
            }
            !name_col || !status_col { next }
            index($name_col, proc) && $status_col == "COMPLETED" { found = 1 }
            END {
                if (!name_col || !status_col) {
                    exit 2
                }
                exit(found ? 0 : 1)
            }
        ' "${trace_file}"; then
            echo "Required process ${p} was not COMPLETED according to ${trace_file}" >&2
            exit 1
        fi
    done

    # MERGE_RESULTS should not run when skip_reporting=true.
    if awk -F '\t' '
        NR==1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "name") {
                    name_col = i
                }
            }
            next
        }
        !name_col { next }
        index($name_col, "MERGE_RESULTS") { found = 1 }
        END {
            if (!name_col) {
                exit 2
            }
            exit(found ? 0 : 1)
        }
    ' "${trace_file}"; then
        echo "MERGE_RESULTS appears in trace despite skip_reporting=true (${trace_file})" >&2
        exit 1
    fi
}

semantic_check_regenie_outputs() {
    local outdir="$1"

    local regenie_dir="${outdir}/regenie_step2"
    if [[ ! -d "${regenie_dir}" ]]; then
        echo "Missing directory: ${regenie_dir}" >&2
        exit 1
    fi

    find "${regenie_dir}" -type f -name "*.regenie" | sort > "${regenie_dir}/.paths"
    if [[ ! -s "${regenie_dir}/.paths" ]]; then
        echo "No REGENIE output files found in ${regenie_dir}" >&2
        exit 1
    fi

    # Validate table-level semantics and p-value bounds for every REGENIE file.
    python3 - "${regenie_dir}/.paths" << 'PY'
import pathlib
import sys

paths_file = pathlib.Path(sys.argv[1])
paths = [pathlib.Path(p.strip()) for p in paths_file.read_text().splitlines() if p.strip()]

def parse_table(path):
    header = None
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if header is None:
                header = cols
                continue
            rows.append(cols)
    if header is None:
        raise RuntimeError(f"{path}: missing header")
    if len(rows) == 0:
        raise RuntimeError(f"{path}: no data rows")
    return header, rows

for p in paths:
    header, rows = parse_table(p)
    if "LOG10P" not in header:
        raise RuntimeError(f"{p}: column LOG10P not found in header: {header}")
    p_idx = header.index("LOG10P")
    for i, row in enumerate(rows, start=2):
        if len(row) <= p_idx:
            raise RuntimeError(f"{p}: row {i} shorter than header")
        try:
            pval = float(row[p_idx])
        except ValueError as exc:
            raise RuntimeError(f"{p}: row {i} has non-numeric LOG10P value: {row[p_idx]}") from exc
        if not (0.0 <= pval):
            raise RuntimeError(f"{p}: row {i} has LOG10P below 0: {pval}")

print(f"Validated {len(paths)} REGENIE files with numeric non-negative LOG10P values.")
PY
}

compare_dosage_modes_semantically() {
    local true_dir="$1"
    local false_dir="$2"

    python3 - "${true_dir}" "${false_dir}" << 'PY'
import math
import pathlib
import re
import sys

true_dir = pathlib.Path(sys.argv[1])
false_dir = pathlib.Path(sys.argv[2])

def parse_regenie(path):
    def parse_optional_float(value):
        if value is None:
            return None
        value = value.strip()
        if value == "" or value.upper() in {"NA", "NAN", ".", "NULL"}:
            return None
        try:
            return float(value)
        except ValueError:
            return None

    header = None
    data = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if header is None:
                header = cols
                continue
            if len(cols) != len(header):
                continue
            rec = dict(zip(header, cols))
            key = rec.get("ID") or rec.get("MASK") or rec.get("GENE") or rec.get("SET") or cols[0]
            p_raw = rec["LOG10P"]
            p = parse_optional_float(p_raw)
            if p is None:
                raise RuntimeError(f"{path}: LOG10P is missing/non-numeric for key '{key}': {p_raw}")
            beta = None
            for b in ("BETA", "Beta"):
                if b in rec:
                    beta = parse_optional_float(rec[b])
                    break
            data[key] = (p, beta)
    if header is None:
        raise RuntimeError(f"Missing header in {path}")
    if "LOG10P" not in header:
        raise RuntimeError(f"Column LOG10P not found in {path}")
    return data

true_files = {p.name: p for p in true_dir.glob("*.regenie")}
false_files = {p.name: p for p in false_dir.glob("*.regenie")}
print(f"true_files: {true_files}")
print(f"false_files: {false_files}")

def canonical_name(name):
    # Strip run-mode prefixes, e.g. initial_true_ / initial_false_
    # so files like initial_true_step2_Y1.regenie and
    # initial_false_step2_Y1.regenie can be paired.
    return re.sub(r"^initial_(true|false)_", "", name)

true_by_canon = {canonical_name(name): path for name, path in true_files.items()}
false_by_canon = {canonical_name(name): path for name, path in false_files.items()}

common = sorted(set(true_by_canon).intersection(false_by_canon))
if not common:
    # Fallback for single-file runs where names may differ for other reasons.
    if len(true_files) == 1 and len(false_files) == 1:
        common = ["__single_file_pair__"]
        t_single = next(iter(true_files.values()))
        f_single = next(iter(false_files.values()))
    else:
        raise RuntimeError(
            "No matching .regenie files between dosage modes after canonicalization"
        )

changes = []
for name in common:
    if name == "__single_file_pair__":
        t_path = t_single
        f_path = f_single
        label = f"{t_path.name} <-> {f_path.name}"
    else:
        t_path = true_by_canon[name]
        f_path = false_by_canon[name]
        label = name

    t_data = parse_regenie(t_path)
    f_data = parse_regenie(f_path)
    keys = sorted(set(t_data).intersection(f_data))
    for key in keys:
        tp, tb = t_data[key]
        fp, fb = f_data[key]
        dp = abs(tp - fp)
        db = abs(tb - fb) if (tb is not None and fb is not None) else 0.0
        if dp > 1e-12 or db > 1e-12:
            changes.append((dp, db, label, key, tp, fp, tb, fb))

if not changes:
    raise RuntimeError(
        "No semantic differences detected between use_dosage=true and use_dosage=false. "
        "Expected at least one change in REGENIE statistics."
    )

changes.sort(reverse=True)
print("Top dosage-mode differences (up to 10):")
for dp, db, name, key, tp, fp, tb, fb in changes[:10]:
    print(
        f"{name}\t{key}\tdelta_p={dp:.3e}\tdelta_beta={db:.3e}\t"
        f"p_true={tp:.6g}\tp_false={fp:.6g}"
    )

print(f"Detected {len(changes)} changed rows across {len(common)} REGENIE files.")
PY

}

run_pipeline "true" "${TRUE_OUTDIR}" "${OUT_BASE}/work_true"
run_pipeline "false" "${FALSE_OUTDIR}" "${OUT_BASE}/work_false"

assert_pipeline_info_exists "${TRUE_OUTDIR}"
assert_pipeline_info_exists "${FALSE_OUTDIR}"
assert_required_processes_completed "${TRUE_OUTDIR}"
assert_required_processes_completed "${FALSE_OUTDIR}"

semantic_check_regenie_outputs "${TRUE_OUTDIR}"
semantic_check_regenie_outputs "${FALSE_OUTDIR}"
compare_dosage_modes_semantically "${TRUE_OUTDIR}/regenie_step2" "${FALSE_OUTDIR}/regenie_step2"

echo "PASS: Initial runnable tests completed."
echo "- Smoke checks passed for both runs."
echo "- REGENIE files passed semantic validation (shape + numeric p-values in range)."
echo "- At least one REGENIE result row changed between dosage modes."
