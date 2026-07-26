#!/usr/bin/env bash
# =============================================================================
# run_benchmark_production.sh — FlexAIDdS Production Benchmark Runner
#
# MacBook Pro M3 Pro 18 GB · macOS 14+ · --workers 6 --omp-threads 2
#
# Execution order:
#   0. Pre-flight: build freshness, constants, ulimits, disk space
#   1. Thermodynamic validation gate (C-1 through C-5 + GPF)
#   2. Pilot calibration: 5 complexes → empirical wall-clock estimate
#   3. Full Astex Diverse (85 complexes) production run
#   4. Summary CSV aggregation
#
# Usage:
#   ./scripts/run_benchmark_production.sh [--dry-run] [--pilot-only]
#        [--skip-pilot] [--phase2] [--seed N] [--threads N]
#
#   caffeinate -i ./scripts/run_benchmark_production.sh
#
# Apache-2.0 · Le Bonhomme Pharma / NRGlab, Université de Montréal
# =============================================================================

set -euo pipefail

# ─── Cleanup trap (define before use) ─────────────────────────────────────────

ACTIVE_PID=""
on_exit() {
    local exit_code="$1"
    if [[ -n "${ACTIVE_PID}" ]] && kill -0 "${ACTIVE_PID}" 2>/dev/null; then
        warn "Killing active docking process PID=${ACTIVE_PID}"
        kill -TERM "${ACTIVE_PID}" 2>/dev/null || true
        sleep 2
        kill -KILL "${ACTIVE_PID}" 2>/dev/null || true
    fi
    if [[ "${exit_code}" -ne 0 ]]; then
        fail "Script exited with code ${exit_code}"
        fail "Check logs in: ${LOG_DIR}"
    fi
}

trap 'on_exit $?' EXIT

# ─── Constants ────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${FLEXAIDDS_BUILD_DIR:-${REPO_ROOT}/build}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
# Prefer iCloud via FLEXAIDDS_ICLOUD (or FLEXAIDDS_RESULTS*) for benchmark outputs/logs.
# Active under working/ subdirs when iCloud used.
if [[ -z "${FLEXAIDDS_RESULTS_ROOT:-}" ]]; then
    if [[ -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
        FLEXAIDDS_RESULTS_ROOT="${FLEXAIDDS_ICLOUD}/results"
    else
        FLEXAIDDS_RESULTS_ROOT="${HOME}/flexaidds_benchmark_results"
    fi
fi
RESULTS_ROOT="${FLEXAIDDS_RESULTS_ROOT}"
if [[ -z "${FLEXAIDDS_LOG_ROOT:-}" ]]; then
    if [[ -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
        FLEXAIDDS_LOG_ROOT="${FLEXAIDDS_ICLOUD}/logs"
    else
        FLEXAIDDS_LOG_ROOT="${HOME}/flexaidds_logs"
    fi
fi
LOG_ROOT="${FLEXAIDDS_LOG_ROOT}"
# For iCloud-based results, route active benchmark outputs through working/ (timestamp keeps unique)
if [[ "${RESULTS_ROOT}" == *"/results" && -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
    RESULTS_ROOT="${RESULTS_ROOT}/working"
fi
LOG_DIR="${LOG_ROOT}/benchmark/${TIMESTAMP}"
RESULTS_DIR="${RESULTS_ROOT}/${TIMESTAMP}"
FIGURES_DIR="${RESULTS_DIR}/figures"
ASTEX_CSV="${REPO_ROOT}/benchmarks/astex_diverse/astex_diverse_set.csv"
ASTEX_STRUCT="${REPO_ROOT}/benchmarks/astex_diverse/structures"
ASTEX_DATA="${REPO_ROOT}/benchmarks/astex_diverse/data"
SUMMARY_CSV="${RESULTS_DIR}/summary.csv"
PILOT_WALL_MIN=300   # 5 min — lower bound acceptable pilot wall clock
PILOT_WALL_MAX=600   # 10 min — upper bound acceptable pilot wall clock
# Sol #9 / 18 GB box: 20 GiB free floor (fail closed). Override only with LP approval.
MIN_FREE_GB="${FLEXAIDDS_MIN_FREE_GB:-20}"

# ─── Colours ─────────────────────────────────────────────────────────────────

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
info()  { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { printf "${GREEN}[ OK ]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
fail()  { printf "${RED}[FAIL]${NC}  %s\n" "$*" >&2; }
phase() { printf "\n${BOLD}══════════════════════════════════════════════════${NC}\n"
          printf "${BOLD}  %s${NC}\n" "$*"
          printf "${BOLD}══════════════════════════════════════════════════${NC}\n\n"; }

# ─── Parse arguments ─────────────────────────────────────────────────────────

DRY_RUN=false
PILOT_ONLY=false
SKIP_PILOT=false
RUN_PHASE2=false
RESUME=false
TWO_PASS=false
BENCHMARK=""
OUTPUT_OVERRIDE=""
SEED=42
N_THREADS=4          # concurrent FlexAIDdS workers (18 GB cap: WORKERS≤4)
OMP_PER_WORKER=2     # OMP threads per worker (4 workers × 2 = 8 threads)

# Two-pass parameters
PASS1_NCHROM=250
PASS1_NGEN=200
PASS1_GRID=0.5
PASS2_NCHROM=1000
PASS2_NGEN=500
PASS2_GRID=0.375
# H_final threshold for flagging (= ln2 = kHSC_hard_nats)
H_FINAL_FLAG_THRESHOLD=0.693

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)      DRY_RUN=true  ;;
        --pilot-only)   PILOT_ONLY=true ;;
        --skip-pilot)   SKIP_PILOT=true ;;
        --phase2)       RUN_PHASE2=true ;;
        --resume)       RESUME=true ;;
        --two-pass)     TWO_PASS=true ;;
        --benchmark)    BENCHMARK="$2"; shift ;;
        --out)          OUTPUT_OVERRIDE="$2"; shift ;;
        --seed)         SEED="$2";     shift ;;
        --threads)      N_THREADS="$2"; shift ;;
        --workers)      N_THREADS="$2"; shift ;;       # alias for --threads
        --omp-threads)  OMP_PER_WORKER="$2"; shift ;;  # override OMP threads/worker
        -h|--help)
            grep '^#' "$0" | head -30 | sed 's/^# \?//'
            exit 0 ;;
        *) fail "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

# Apply output override if provided
if [[ -n "${OUTPUT_OVERRIDE}" ]]; then
    RESULTS_DIR="${OUTPUT_OVERRIDE}"
    LOG_DIR="${OUTPUT_OVERRIDE}/logs"
    FIGURES_DIR="${RESULTS_DIR}/figures"
    SUMMARY_CSV="${RESULTS_DIR}/summary.csv"
fi

# ─── Dry-run wrapper ──────────────────────────────────────────────────────────

run() {
    if [[ "${DRY_RUN}" == true ]]; then
        printf "${YELLOW}[DRY]${NC}   %s\n" "$*"
    else
        "$@"
    fi
}

# ─── OMP / environment ────────────────────────────────────────────────────────

setup_env() {
    export OMP_NUM_THREADS="${OMP_PER_WORKER}"   # 2 per worker; DatasetRunner also sets per-subprocess
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    export OMP_WAIT_POLICY=passive
    export FLEXAID_MAX_MEM_MB=16384
    export FLEXAID_SIMD=NEON          # M3 Pro: Neon replaces AVX2
    export FLEXAID_SEED="${SEED}"
    export SHANNON_TRACE_LEVEL=2      # per-step CSV for Astex runs
    export FLEXAIDS_METAL_DEVICE=default  # Enable Metal GPU acceleration
    export FLEXAIDDS_REPO="${REPO_ROOT}"
    ulimit -n 65536 2>/dev/null || warn "Could not raise file descriptor limit"
    ulimit -s unlimited 2>/dev/null || warn "Could not raise stack limit"

    info "WORKERS=${N_THREADS}  OMP/worker=${OMP_NUM_THREADS}  SEED=${SEED}  DRY_RUN=${DRY_RUN}  METAL=ON"
}

ensure_runtime_data_files() {
    local dir="$1"
    local sources=(
        "${REPO_ROOT}/build"
        "${REPO_ROOT}"
        "/private/tmp/flexaidds-build-tests"
        "/private/tmp/flexaidds-build-prod"
    )

    mkdir -p "${dir}"
    for data_file in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def; do
        if [[ -s "${dir}/${data_file}" ]]; then
            continue
        fi
        for src in "${sources[@]}"; do
            [[ "${src}" == "${dir}" ]] && continue
            if [[ -s "${src}/${data_file}" ]]; then
                cp -p "${src}/${data_file}" "${dir}/${data_file}"
                break
            fi
        done
    done

    for data_file in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def; do
        [[ -s "${dir}/${data_file}" ]] || return 1
    done
}

benchmark_manifest_path() {
    case "${BENCHMARK:-astex}" in
        casf|casf2016)
            printf '%s\n' "${REPO_ROOT}/benchmarks/datasets/casf2016.yaml"
            ;;
        astex_non_native|astex_nonnative)
            printf '%s\n' "${REPO_ROOT}/benchmarks/datasets/astex_nonnative.yaml"
            ;;
        hap2)
            printf '%s\n' "${REPO_ROOT}/benchmarks/datasets/hap2.yaml"
            ;;
        *)
            printf '%s\n' "${REPO_ROOT}/benchmarks/datasets/astex_diverse.yaml"
            ;;
    esac
}

# ─── Locate binary ────────────────────────────────────────────────────────────

locate_binary() {
    if [[ -n "${FLEXAIDDS_BINARY:-}" && -x "${FLEXAIDDS_BINARY}" ]]; then
        FLEXAIDDS_BIN="${FLEXAIDDS_BINARY}"
        BUILD_DIR="$(cd "$(dirname "${FLEXAIDDS_BIN}")" && pwd)"
        ensure_runtime_data_files "${BUILD_DIR}" || return 1
        ok "Binary: ${FLEXAIDDS_BIN}"
        export FLEXAIDDS_BUILD_DIR="${BUILD_DIR}"
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        export FLEXAIDDS_REPO="${REPO_ROOT}"
        return 0
    fi
    if [[ -n "${FLEXAID_BINARY:-}" && -x "${FLEXAID_BINARY}" ]]; then
        FLEXAIDDS_BIN="${FLEXAID_BINARY}"
        BUILD_DIR="$(cd "$(dirname "${FLEXAIDDS_BIN}")" && pwd)"
        ensure_runtime_data_files "${BUILD_DIR}" || return 1
        ok "Binary: ${FLEXAIDDS_BIN}"
        export FLEXAIDDS_BUILD_DIR="${BUILD_DIR}"
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        export FLEXAIDDS_REPO="${REPO_ROOT}"
        return 0
    fi
    if [[ -n "${FLEXAIDDS_BUILD_DIR:-}" ]]; then
        local build_candidates=(
            "${BUILD_DIR}/FlexAIDdS"
            "${BUILD_DIR}/FlexAID"
        )
        for b in "${build_candidates[@]}"; do
            if [[ -x "${b}" ]]; then
                FLEXAIDDS_BIN="${b}"
                BUILD_DIR="$(cd "$(dirname "${FLEXAIDDS_BIN}")" && pwd)"
                ensure_runtime_data_files "${BUILD_DIR}" || return 1
                ok "Binary: ${FLEXAIDDS_BIN}"
                export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
                export FLEXAIDDS_BUILD_DIR="${BUILD_DIR}"
                export FLEXAIDDS_REPO="${REPO_ROOT}"
                return 0
            fi
        done
    fi

    local candidates=(
        "/private/tmp/flexaidds-build-prod/FlexAIDdS"
        "/private/tmp/flexaidds-build-tests/FlexAIDdS"
        "${BUILD_DIR}/FlexAIDdS"
        "${REPO_ROOT}/WRK/FlexAID"
        "${BUILD_DIR}/FlexAID"
        "${REPO_ROOT}/BIN/FlexAIDdS"
        "${REPO_ROOT}/BIN/FlexAID"
    )
    for b in "${candidates[@]}"; do
        if [[ -x "${b}" ]]; then
            FLEXAIDDS_BIN="${b}"
            BUILD_DIR="$(cd "$(dirname "${FLEXAIDDS_BIN}")" && pwd)"
            ensure_runtime_data_files "${BUILD_DIR}" || return 1
            ok "Binary: ${FLEXAIDDS_BIN}"
            export FLEXAIDDS_BUILD_DIR="${BUILD_DIR}"
            export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
            export FLEXAIDDS_REPO="${REPO_ROOT}"
            return 0
        fi
    done
    fail "FlexAIDdS binary not found in: ${candidates[*]}"
    fail "Run: cmake --build ${BUILD_DIR} --target FlexAIDdS -j ${N_THREADS}"
    return 1
}

locate_dataset_runner() {
    DATASET_BIN="${BUILD_DIR}/benchmark_datasets"
    if [[ -x "${DATASET_BIN}" ]]; then
        ok "benchmark_datasets: ${DATASET_BIN}"
        return 0
    fi
    warn "benchmark_datasets not found — will use individual FlexAIDdS calls"
    DATASET_BIN=""
    return 0
}

# ─── Locate structure files ───────────────────────────────────────────────────

# Returns receptor and ligand paths for a given PDB code.
# Searches multiple candidate layouts produced by download.sh / prepare-only.
find_structure_files() {
    local pdb="$1"
    local pdb_lower="$(printf '%s\n' "$pdb" | tr '[:upper:]' '[:lower:]')"
    local pdb_upper="$(printf '%s\n' "$pdb" | tr '[:lower:]' '[:upper:]')"

    local receptor="" ligand=""

    # Layout 1: structures/{PDB}/{pdb}_protein.pdb  (plan Section 9 canonical)
    if [[ -f "${ASTEX_STRUCT}/${pdb_upper}/${pdb_lower}_protein.pdb" ]]; then
        receptor="${ASTEX_STRUCT}/${pdb_upper}/${pdb_lower}_protein.pdb"
        ligand="${ASTEX_STRUCT}/${pdb_upper}/${pdb_lower}_ligand.mol2"
    elif [[ -f "${ASTEX_STRUCT}/${pdb_lower}/${pdb_lower}_protein.pdb" ]]; then
        receptor="${ASTEX_STRUCT}/${pdb_lower}/${pdb_lower}_protein.pdb"
        ligand="${ASTEX_STRUCT}/${pdb_lower}/${pdb_lower}_ligand.mol2"
    # Layout 2: data/{pdb}.pdb  (download.sh flat layout)
    elif [[ -f "${ASTEX_DATA}/${pdb_lower}.pdb" ]]; then
        receptor="${ASTEX_DATA}/${pdb_lower}.pdb"
        ligand="${ASTEX_DATA}/${pdb_lower}_ligand.mol2"
    fi

    if [[ -z "${receptor}" ]] || [[ ! -f "${receptor}" ]]; then
        echo "MISSING_RECEPTOR"
        return 1
    fi
    if [[ -z "${ligand}" ]] || [[ ! -f "${ligand}" ]]; then
        echo "MISSING_LIGAND"
        return 1
    fi

    echo "${receptor}:${ligand}"
}

# ─── Pre-flight checks ────────────────────────────────────────────────────────

preflight() {
    phase "PRE-FLIGHT CHECKLIST"
    local n_fail=0

    # 0. Sol #9 multi-session guard (hold + mkdir lock + disk + workers + binary pin)
    info "0. Dock session guard (Sol #9)"
    local guard="${REPO_ROOT}/scripts/dock_session_guard.py"
    if [[ ! -f "${guard}" ]]; then
        fail "   dock_session_guard.py missing at ${guard} (Sol #9 fail-closed)"
        exit 77
    fi
    local binary_candidate="${FLEXAIDDS_BINARY:-${BUILD_DIR}/FlexAIDdS}"
    if [[ ! -x "${binary_candidate}" ]]; then
        binary_candidate="${BUILD_DIR}/FlexAIDdS"
    fi
    local guard_args=(
        preflight
        --out-dir "${RESULTS_DIR}"
        --workers "${N_THREADS}"
        --repo-root "${REPO_ROOT}"
        --max-workers 4
        --min-free-gb "${MIN_FREE_GB}"
        --owner "run_benchmark_production"
        --note "production benchmark preflight"
    )
    if [[ -x "${binary_candidate}" ]]; then
        guard_args+=(--binary "${binary_candidate}")
    else
        # Still require guard; fail if engine cannot be pinned for a real run.
        if [[ "${DRY_RUN}" != true ]]; then
            fail "   no executable binary to pin at ${binary_candidate}"
            exit 77
        fi
        guard_args+=(--no-copy-binary)
    fi
    if [[ "${DRY_RUN}" == true ]]; then
        guard_args+=(--no-lock)
    fi
    if ! python3 "${guard}" "${guard_args[@]}"; then
        fail "   dock_session_guard refused launch (hold/lock/disk/workers/binary)"
        # Hold/lock are hard stops even in dry-run.
        exit 78
    fi
    ok "   hold/lock/disk/workers preflight accepted"
    # Rebind to run-namespace pin so rebuild of shared build/ cannot hit live exec.
    local pinned_engine="${RESULTS_DIR}/bin/$(basename "${binary_candidate}")"
    if [[ -x "${pinned_engine}" ]]; then
        FLEXAIDDS_BIN="${pinned_engine}"
        export FLEXAIDDS_BINARY="${pinned_engine}"
        export FLEXAIDDS_BIN="${pinned_engine}"
        ok "   rebound FLEXAIDDS_BIN -> ${pinned_engine}"
    elif [[ "${DRY_RUN}" != true ]]; then
        fail "   pinned engine missing after preflight: ${pinned_engine}"
        exit 77
    fi
    if [[ "${DRY_RUN}" != true ]]; then
        # Foreground production waits for docks; release lock only after they finish.
        # Do NOT force-release while a recorded dock_pid is live.
        trap 'python3 "'"${guard}"'" release-lock 2>/dev/null || true; on_exit $?' EXIT
    fi

    # 1. Git state
    info "1. Git state"
    local git_commit git_dirty
    git_commit="$(cd "${REPO_ROOT}" && git rev-parse --short HEAD 2>/dev/null || echo 'UNKNOWN')"
    git_dirty="$(cd "${REPO_ROOT}" && git status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
    if [[ "${git_dirty}" -eq 0 ]]; then
        ok "   HEAD=${git_commit} (clean)"
    else
        warn "   HEAD=${git_commit} — ${git_dirty} uncommitted file(s)"
    fi

    # 2. Binary exists
    info "2. Binary"
    if locate_binary; then
        # Check freshness: find any source newer than binary using POSIX -newer
        # Guard: only pass directories that exist (set -o pipefail + nonexistent
        # dir causes find to exit 1, poisoning the pipeline even with 2>/dev/null)
        local stale_count=0
        local _search_dirs=()
        for _d in "${REPO_ROOT}/LIB" "${REPO_ROOT}/src"; do
            [[ -d "${_d}" ]] && _search_dirs+=("${_d}")
        done
        if [[ ${#_search_dirs[@]} -gt 0 ]]; then
            stale_count="$(find "${_search_dirs[@]}" \
                            \( -name '*.cpp' -o -name '*.h' \) \
                            -newer "${FLEXAIDDS_BIN}" 2>/dev/null \
                          | wc -l | tr -d ' ')" || stale_count=0
        fi
        if [[ "${stale_count}" -gt 0 ]]; then
            warn "   Binary is STALE — ${stale_count} source file(s) newer than binary"
            warn "   Run: cmake --build ${BUILD_DIR} --target FlexAIDdS -j ${N_THREADS}"
            n_fail=$((n_fail + 1))
        else
            ok "   Binary is fresh"
        fi
    else
        n_fail=$((n_fail + 1))
    fi
    locate_dataset_runner

    # 3. Unit tests (last ctest run)
    info "3. Unit test gate"
    local failed_log="${BUILD_DIR}/Testing/Temporary/LastTestsFailed.log"
    if [[ -f "${failed_log}" ]] && [[ -s "${failed_log}" ]]; then
        local n_failed
        n_failed="$(wc -l < "${failed_log}" | tr -d ' ')"
        fail "   ${n_failed} test(s) failed in last ctest run:"
        cat "${failed_log}" | sed 's/^/       /'
        fail "   Run: cd ${BUILD_DIR} && ctest --output-on-failure -j${N_THREADS}"
        n_fail=$((n_fail + 1))
    else
        ok "   No failed tests in last ctest run"
        warn "   (Verify with: cd ${BUILD_DIR} && ctest -j${N_THREADS})"
    fi

    # 4. Disk space
    info "4. Disk space"
    local free_gb
    free_gb="$(df -BG "${REPO_ROOT}" 2>/dev/null \
               | awk 'NR==2{gsub(/G/,"",$4); print $4}' \
               || df -g "${REPO_ROOT}" 2>/dev/null \
               | awk 'NR==2{print $4}' \
               || echo 0)"
    if [[ "${free_gb}" -ge "${MIN_FREE_GB}" ]] 2>/dev/null; then
        ok "   ${free_gb} GB free (≥ ${MIN_FREE_GB} GB minimum)"
    else
        warn "   Only ${free_gb} GB free — output is configured outside the repo, but space is still tight"
        [[ "${free_gb}" -lt "${MIN_FREE_GB}" ]] && n_fail=$((n_fail + 1))
    fi

    # 5. Astex CSV
    info "5. Astex Diverse CSV"
    if [[ -f "${ASTEX_CSV}" ]]; then
        local n_complexes
        n_complexes="$(tail -n +2 "${ASTEX_CSV}" | wc -l | tr -d ' ')"
        ok "   ${ASTEX_CSV} (${n_complexes} complexes)"
        if [[ "${n_complexes}" -lt 85 ]]; then
            warn "   Expected 85 complexes, found ${n_complexes}"
        fi
    else
        fail "   Astex CSV not found: ${ASTEX_CSV}"
        n_fail=$((n_fail + 1))
    fi

    # 6. Structure files (sample check)
    info "6. Astex structure files"
    local pilot_pdb
    pilot_pdb="$(tail -n +2 "${ASTEX_CSV}" | head -1 | cut -d',' -f1)"
    local struct_info
    struct_info="$(find_structure_files "${pilot_pdb}" 2>/dev/null || echo 'MISSING')"
    if [[ "${struct_info}" == MISSING* ]]; then
        warn "   Structures NOT downloaded yet — run before pilot:"
        warn "   bash benchmarks/astex_diverse/download.sh"
        warn "   OR: ${BUILD_DIR}/benchmark_datasets --benchmark astex --prepare-only"
        n_fail=$((n_fail + 1))
    else
        ok "   Sample structure found: ${pilot_pdb}"
    fi

    # 7. kB_kcal constant sanity (grep source)
    info "7. kB_kcal constant"
    local kb_count
    kb_count="$(grep -r 'kB_kcal.*0\.001987206\|0\.001987206.*kB_kcal' \
                "${REPO_ROOT}/LIB" --include='*.h' --include='*.cpp' 2>/dev/null \
                | wc -l | tr -d ' ')"
    if [[ "${kb_count}" -ge 2 ]]; then
        ok "   kB_kcal = 0.001987206 confirmed in ${kb_count} source location(s)"
    else
        warn "   kB_kcal constant not verified (grep matched ${kb_count} lines)"
    fi

    # 8. Shannon thresholds
    info "8. Shannon HSC thresholds"
    local hsc_soft hsc_hard
    hsc_soft="$(grep -r 'kHSC_soft_nats\|1\.3863\|kHSC_soft_bits.*2\.0' \
                "${REPO_ROOT}/LIB" --include='*.h' 2>/dev/null | head -1 | xargs echo)"
    hsc_hard="$(grep -r 'kHSC_hard_nats\|0\.6931\|kHSC_hard_bits.*1\.0' \
                "${REPO_ROOT}/LIB" --include='*.h' 2>/dev/null | head -1 | xargs echo)"
    if [[ -n "${hsc_soft}" ]] && [[ -n "${hsc_hard}" ]]; then
        ok "   HSC thresholds: soft=2·ln2 nats (1.3863), hard=1·ln2 nats (0.6931)"
    else
        warn "   HSC thresholds not verified in source"
    fi

    # 9. Git tag
    info "9. Git tag (pre-benchmark snapshot)"
    local tag_name="benchmark-$(date +%Y%m%d)"
    local existing_tag
    existing_tag="$(cd "${REPO_ROOT}" && git tag -l "${tag_name}" 2>/dev/null)"
    if [[ -z "${existing_tag}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            info "   [DRY] git tag ${tag_name} -m 'Pre-benchmark snapshot'"
        else
            cd "${REPO_ROOT}" && git tag "${tag_name}" \
                -m "Pre-benchmark snapshot (${git_commit})" 2>/dev/null \
                && ok "   Tagged: ${tag_name}" \
                || warn "   Could not create tag ${tag_name} (already exists?)"
        fi
    else
        ok "   Tag already exists: ${tag_name}"
    fi

    # Summary
    echo ""
    if [[ "${n_fail}" -gt 0 ]]; then
        fail "Pre-flight: ${n_fail} BLOCKING issue(s) — address before running full benchmark"
        if [[ "${DRY_RUN}" == false ]] && [[ "${SKIP_PILOT}" == false ]]; then
            exit 1
        else
            warn "Continuing despite ${n_fail} failure(s) (dry-run or --skip-pilot mode)"
        fi
    else
        ok "Pre-flight: ALL CLEAR"
    fi
}

# ─── Thermodynamic validation gate ────────────────────────────────────────────

thermo_gate() {
    phase "THERMODYNAMIC VALIDATION GATE (C-1 … C-5 + GPF)"
    local thermo_script="${REPO_ROOT}/tests/thermodynamic_validation.py"
    local thermo_out="${RESULTS_DIR}/thermo_validation.json"

    if [[ ! -f "${thermo_script}" ]]; then
        warn "Thermodynamic validation script not found: ${thermo_script}"
        warn "Skipping C-1…C-5 gate — build with Python bindings to enable"
        return 0
    fi

    run python3 "${thermo_script}" \
        --build-dir "${BUILD_DIR}" \
        --output "${thermo_out}" \
        || { fail "THERMODYNAMIC VALIDATION FAILED — abort"; exit 1; }
    ok "Thermodynamic gate passed"
}

# ─── Dock one complex ─────────────────────────────────────────────────────────

# dock_complex <pdb_id> <receptor.pdb> <ligand.mol2> <out_dir> <log_file>
# Writes wall_time_s to <out_dir>/wall_time_s
dock_complex() {
    local pdb="$1" receptor="$2" ligand="$3" out_dir="$4" log_file="$5"

    run mkdir -p "${out_dir}"

    if [[ "${DRY_RUN}" == true ]]; then
        printf "[DRY]   FlexAIDdS --input %s --ligand %s --output %s --seed %d --nthreads %d\n" \
               "${receptor}" "${ligand}" "${out_dir}" "${SEED}" "${N_THREADS}"
        echo "0" > "${out_dir}/wall_time_s" 2>/dev/null || true
        return 0
    fi

    local t_start t_end wall_s
    t_start="$(date +%s%N)"

    ACTIVE_PID=""
    # Note: seed and threads are typically set via environment or config
    # Binary uses positional args: receptor ligand [options]
    "${FLEXAIDDS_BIN}" \
        "${receptor}" \
        "${ligand}" \
        -o "${out_dir}" \
        >> "${log_file}" 2>&1 &
    ACTIVE_PID=$!
    wait "${ACTIVE_PID}" || {
        fail "Docking failed for ${pdb} (exit code $?)"
        ACTIVE_PID=""
        return 1
    }
    ACTIVE_PID=""

    t_end="$(date +%s%N)"
    wall_s=$(( (t_end - t_start) / 1000000000 ))
    echo "${wall_s}" > "${out_dir}/wall_time_s"
    return 0
}

# ─── Extract top-1 score and RMSD from output ─────────────────────────────────

extract_result() {
    local out_dir="$1" pdb="$2" rmsd_threshold="$3"
    local wall_s top1_score rmsd success

    wall_s="$(cat "${out_dir}/wall_time_s" 2>/dev/null || echo 0)"

    # Try JSON binding modes first
    if [[ -f "${out_dir}/binding_modes.json" ]]; then
        top1_score="$(python3 -c "
import json, sys
try:
    d = json.load(open('${out_dir}/binding_modes.json'))
    modes = d.get('binding_modes', [])
    if modes:
        print(modes[0].get('best_score', modes[0].get('score', 'N/A')))
    else:
        print('N/A')
except Exception as e:
    print('N/A')
" 2>/dev/null || echo 'N/A')"

        rmsd="$(python3 -c "
import json, sys
try:
    d = json.load(open('${out_dir}/binding_modes.json'))
    modes = d.get('binding_modes', [])
    if modes:
        val = modes[0].get('best_pose_rmsd', modes[0].get('rmsd_to_crystal', 'N/A'))
        print(val)
    else:
        print('N/A')
except Exception as e:
    print('N/A')
" 2>/dev/null || echo 'N/A')"
    else
        top1_score="N/A"
        rmsd="N/A"
    fi

    # Determine success
    if [[ "${rmsd}" != "N/A" ]] && python3 -c \
        "import sys; sys.exit(0 if float('${rmsd}') < ${rmsd_threshold} else 1)" \
        2>/dev/null; then
        success=1
    else
        success=0
    fi

    echo "${pdb},${wall_s},${top1_score},${rmsd},${success}"
}

# ─── Two-pass helpers ─────────────────────────────────────────────────────────

# Write a minimal JSON config override for FlexAIDdS
# Usage: write_flexaid_config <file> <nchrom> <ngen> <grid>
write_flexaid_config() {
    local cfg_file="$1" nchrom="$2" ngen="$3" grid="$4"
    cat > "${cfg_file}" <<EOF
{
  "optimization": { "grid_spacing": ${grid} },
  "ga": {
    "num_chromosomes": ${nchrom},
    "num_generations": ${ngen}
  }
}
EOF
}

# Two-pass docking for a single complex.
# Usage: dock_complex_twopass <pdb> <receptor> <ligand> <out_dir> <log_file>
# Runs Pass 1 (coarse); flags for Pass 2 if H_final >= threshold or RMSD > 3Å.
# Copies the final-pass binding_modes.json to out_dir for extract_result.
dock_complex_twopass() {
    local pdb="$1" receptor="$2" ligand="$3" out_dir="$4" log_file="$5"
    local pass1_dir="${out_dir}/pass1"
    local pass2_dir="${out_dir}/pass2"

    run mkdir -p "${pass1_dir}"

    # ── Pass 1 ──
    local pass1_cfg pass1_exit=0
    pass1_cfg="$(mktemp /tmp/flexaid_p1_XXXXXX.json)"
    write_flexaid_config "${pass1_cfg}" ${PASS1_NCHROM} ${PASS1_NGEN} ${PASS1_GRID}

    local t0 t1
    t0="$(date +%s%N)"

    if [[ "${DRY_RUN}" == false ]]; then
        ACTIVE_PID=""
        "${FLEXAIDDS_BIN}" "${receptor}" "${ligand}" \
            -o "${pass1_dir}" -c "${pass1_cfg}" \
            >> "${log_file}" 2>&1 &
        ACTIVE_PID=$!
        wait "${ACTIVE_PID}" || pass1_exit=$?
        ACTIVE_PID=""
    else
        printf "${YELLOW}[DRY]${NC}   Pass 1: %s × %d chrom × %d gen, grid=%.3f\n" \
               "${pdb}" ${PASS1_NCHROM} ${PASS1_NGEN} ${PASS1_GRID}
    fi
    rm -f "${pass1_cfg}"

    t1="$(date +%s%N)"
    local pass1_wall=$(( (t1 - t0) / 1000000000 ))
    echo "${pass1_wall}" > "${pass1_dir}/wall_time_s"

    # Parse H_final from the log (matches "H_final = <float>" from gaboom)
    local H_FINAL
    H_FINAL="$(grep -oE 'H_final[[:space:]]*=[[:space:]]*[0-9]+\.[0-9]+' \
                   "${log_file}" 2>/dev/null | tail -1 | grep -oE '[0-9]+\.[0-9]+' || echo 999)"

    # Parse RMSD from Pass 1 binding_modes.json
    local pass1_rmsd="N/A"
    if [[ -f "${pass1_dir}/binding_modes.json" ]]; then
        pass1_rmsd="$(python3 -c "
import json
try:
    d=json.load(open('${pass1_dir}/binding_modes.json'))
    m=d.get('binding_modes',[])
    print(m[0].get('best_pose_rmsd', m[0].get('rmsd_to_crystal','N/A')) if m else 'N/A')
except: print('N/A')
" 2>/dev/null || echo 'N/A')"
    fi

    # ── Flagging ──
    local needs_pass2=false flag_reason=""
    if [[ "${pass1_exit}" -ne 0 ]]; then
        needs_pass2=true; flag_reason="pass1_exit=${pass1_exit}"
    elif python3 -c "import sys; sys.exit(0 if float('${H_FINAL}') >= ${H_FINAL_FLAG_THRESHOLD} else 1)" 2>/dev/null; then
        needs_pass2=true; flag_reason="H_final=${H_FINAL}"
    elif [[ "${pass1_rmsd}" != "N/A" ]] && \
         python3 -c "import sys; sys.exit(0 if float('${pass1_rmsd}') > 3.0 else 1)" 2>/dev/null; then
        needs_pass2=true; flag_reason="rmsd=${pass1_rmsd}"
    fi

    local final_pass="pass1"
    local total_wall="${pass1_wall}"

    # ── Pass 2 (if flagged) ──
    if [[ "${needs_pass2}" == true ]]; then
        info "  ${pdb}: → Pass 2 (${flag_reason})"
        run mkdir -p "${pass2_dir}"
        t0="$(date +%s%N)"
        if [[ "${DRY_RUN}" == false ]]; then
            ACTIVE_PID=""
            "${FLEXAIDDS_BIN}" "${receptor}" "${ligand}" \
                -o "${pass2_dir}" \
                >> "${log_file}" 2>&1 &
            ACTIVE_PID=$!
            wait "${ACTIVE_PID}" || true
            ACTIVE_PID=""
        else
            printf "${YELLOW}[DRY]${NC}   Pass 2: %s × %d chrom × %d gen, grid=%.3f\n" \
                   "${pdb}" ${PASS2_NCHROM} ${PASS2_NGEN} ${PASS2_GRID}
        fi
        t1="$(date +%s%N)"
        local pass2_wall=$(( (t1 - t0) / 1000000000 ))
        echo "${pass2_wall}" > "${pass2_dir}/wall_time_s"
        total_wall=$(( pass1_wall + pass2_wall ))
        final_pass="pass2"
    fi

    # Write totals and promote final binding_modes.json to out_dir
    echo "${total_wall}" > "${out_dir}/wall_time_s"
    echo "${final_pass}"  > "${out_dir}/final_pass"
    local final_dir="${out_dir}/${final_pass}"
    if [[ -f "${final_dir}/binding_modes.json" ]] && [[ "${DRY_RUN}" == false ]]; then
        cp "${final_dir}/binding_modes.json" "${out_dir}/binding_modes.json" 2>/dev/null || true
    fi
}

# Two-pass run using benchmark_datasets binary (CASF / Astex Non-Native / etc.)
# Usage: run_twopass_dataset <dataset_name> <output_base_dir>
run_twopass_dataset() {
    local dataset="$1" out_base="$2"
    local pass1_out="${out_base}/pass1"
    local pass2_out="${out_base}/pass2"
    local flagged_file="${out_base}/flagged_for_pass2.txt"

    if [[ -z "${DATASET_BIN:-}" ]]; then
        fail "benchmark_datasets binary not found — required for --benchmark ${dataset} --two-pass"
        exit 1
    fi

    # ── Pass 1: coarse run ──
    phase "TWO-PASS PASS 1: ${dataset} (${PASS1_NCHROM}×${PASS1_NGEN}, grid=${PASS1_GRID}Å)"
    run mkdir -p "${pass1_out}" "${LOG_DIR}"
    export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"

    run "${DATASET_BIN}" \
        --benchmark "${dataset}" \
        --ga-population  ${PASS1_NCHROM} \
        --ga-generations ${PASS1_NGEN} \
        --grid-spacing   ${PASS1_GRID} \
        --output     "${pass1_out}" \
        --threads    "${N_THREADS}" \
        --omp-threads "${OMP_PER_WORKER}" \
        --job-timeout-seconds 1800 \
        2>&1 | tee "${LOG_DIR}/${dataset}_pass1.log"

    # ── Flag complexes for Pass 2 ──
    info "Analysing Pass 1 results for Pass 2 flagging..."
    python3 - "${pass1_out}" "${flagged_file}" "${H_FINAL_FLAG_THRESHOLD}" <<'PYEOF'
import csv, os, sys, json, re, glob

out_dir      = sys.argv[1]
flagged_file = sys.argv[2]
threshold    = float(sys.argv[3])
flagged      = []

for cdir in sorted(glob.glob(os.path.join(out_dir, '*/'))):
    pdb = os.path.basename(cdir.rstrip('/'))
    # Parse H_final from any .log in the complex dir
    h_final = 999.0
    for lf in glob.glob(os.path.join(cdir, '*.log')) + [os.path.join(cdir, 'stdout.log')]:
        try:
            for line in open(lf):
                m = re.search(r'H_final\s*=\s*([\d.]+)', line)
                if m:
                    h_final = float(m.group(1))
        except Exception:
            pass
    # Parse RMSD from the benchmark_datasets result.csv for this complex
    rmsd = None
    rcsv = os.path.join(cdir, 'result.csv')
    if os.path.exists(rcsv):
        try:
            with open(rcsv) as fh:
                rows = list(csv.DictReader(fh))
            if rows:
                v = rows[0].get('rmsd_to_crystal')
                if v not in (None, '', 'N/A'):
                    rmsd = float(v)
        except Exception:
            pass
    flag = (h_final >= threshold) or (rmsd is not None and rmsd > 3.0)
    if flag:
        flagged.append(pdb)

with open(flagged_file, 'w') as f:
    for code in flagged:
        f.write(code + '\n')
print(f"Flagged {len(flagged)} complexes for Pass 2 (threshold H>={threshold})")
PYEOF

    local n_flagged=0
    [[ -f "${flagged_file}" ]] && n_flagged="$(wc -l < "${flagged_file}" | tr -d ' ')"
    ok "Pass 1 done. ${n_flagged} complexes flagged for Pass 2."

    if [[ "${n_flagged}" -gt 0 ]]; then
        # ── Pass 2: full resolution on flagged only ──
        phase "TWO-PASS PASS 2: ${n_flagged} flagged complexes (${PASS2_NCHROM}×${PASS2_NGEN})"
        run mkdir -p "${pass2_out}"

        run "${DATASET_BIN}" \
            --benchmark "pdb_list:${flagged_file}" \
            --ga-population  ${PASS2_NCHROM} \
            --ga-generations ${PASS2_NGEN} \
            --grid-spacing   ${PASS2_GRID} \
            --output     "${pass2_out}" \
            --threads    "${N_THREADS}" \
            --omp-threads "${OMP_PER_WORKER}" \
            --job-timeout-seconds 7200 \
            2>&1 | tee "${LOG_DIR}/${dataset}_pass2.log"
    fi

    # ── Merge: Pass 1 results + Pass 2 overrides → summary CSV ──
    info "Merging Pass 1 + Pass 2 into summary CSV..."
    python3 - "${pass1_out}" "${pass2_out}" "${out_base}/summary.csv" <<'PYEOF'
import csv, os, sys, json, glob

pass1_dir  = sys.argv[1]
pass2_dir  = sys.argv[2]
out_csv    = sys.argv[3]

def best_result(cdir):
    wall = 0
    try: wall = int(open(os.path.join(cdir, 'wall_time_s')).read().strip())
    except: pass
    score, rmsd, success = 'N/A', 'N/A', 0
    rcsv = os.path.join(cdir, 'result.csv')
    if os.path.exists(rcsv):
        try:
            with open(rcsv) as fh:
                rows = list(csv.DictReader(fh))
            if rows:
                row = rows[0]
                score = row.get('best_score', 'N/A')
                r = row.get('rmsd_to_crystal', 'N/A')
                if r not in (None, '', 'N/A'):
                    rmsd = str(r)
                    success = 1 if float(r) < 2.0 else 0
        except: pass
    return wall, score, rmsd, success

rows = []
for cdir in sorted(glob.glob(os.path.join(pass1_dir, '*/'))):
    pdb = os.path.basename(cdir.rstrip('/'))
    p2_cdir = os.path.join(pass2_dir, pdb)
    if os.path.isdir(p2_cdir):
        w, s, r, ok = best_result(p2_cdir)
        rows.append((pdb, w, s, r, ok, 'pass2'))
    else:
        w, s, r, ok = best_result(cdir)
        rows.append((pdb, w, s, r, ok, 'pass1'))

with open(out_csv, 'w') as f:
    f.write('complex_id,wall_time_s,top1_score,rmsd_to_crystal,success,pass\n')
    for pdb, w, s, r, ok, p in rows:
        f.write(f'{pdb},{w},{s},{r},{ok},{p}\n')

n_p1 = sum(1 for *_, p in rows if p == 'pass1')
n_p2 = sum(1 for *_, p in rows if p == 'pass2')
print(f'Summary: {len(rows)} complexes → Pass1={n_p1} ({100*n_p1//max(len(rows),1)}%), Pass2={n_p2} ({100*n_p2//max(len(rows),1)}%)')
PYEOF

    ok "Two-pass run complete for ${dataset}. Summary: ${out_base}/summary.csv"
}

# ─── Pilot calibration (5 complexes) ─────────────────────────────────────────

run_pilot() {
    phase "PILOT CALIBRATION (5 complexes)"

    # Take lines 2-6 from CSV (skip header), use first column (pdb_id)
    local pilot_ids=()
    while IFS=',' read -r pdb_id lig_id res rmsd_thr ref; do
        pilot_ids+=("${pdb_id}")
    done < <(tail -n +2 "${ASTEX_CSV}" | head -5)

    info "Pilot complexes: ${pilot_ids[*]}"
    info "Expected wall-clock: 5–8 min/complex (${PILOT_WALL_MIN}–${PILOT_WALL_MAX} s)"
    echo ""

    local pilot_dir="${RESULTS_DIR}/pilot"
    run mkdir -p "${pilot_dir}"

    local total_wall=0 n_docked=0

    for pdb in "${pilot_ids[@]}"; do
        local struct_info
        struct_info="$(find_structure_files "${pdb}" 2>/dev/null | head -1 || echo 'MISSING')"

        if [[ "${struct_info}" == MISSING* ]] || [[ -z "${struct_info}" ]]; then
            warn "  ${pdb}: structures not found — skipping (run download.sh first)"
            continue
        fi

        local receptor="${struct_info%%:*}"
        local ligand="${struct_info##*:}"
        local out_dir="${pilot_dir}/${pdb}"
        local log_file="${LOG_DIR}/pilot_${pdb}.log"

        run mkdir -p "${LOG_DIR}"
        info "  Docking ${pdb} ..."

        if dock_complex "${pdb}" "${receptor}" "${ligand}" "${out_dir}" "${log_file}"; then
            local wall_s
            wall_s="$(cat "${out_dir}/wall_time_s" 2>/dev/null || echo 0)"
            ok "  ${pdb}: ${wall_s}s"
            total_wall=$((total_wall + wall_s))
            n_docked=$((n_docked + 1))
        else
            warn "  ${pdb}: docking failed"
        fi
    done

    if [[ "${n_docked}" -eq 0 ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            info "Dry-run: skipping pilot timing check"
            return 0
        fi
        fail "No complexes docked in pilot — check binary and structure files"
        exit 1
    fi

    local mean_wall=$(( total_wall / n_docked ))
    local est_full_h=$(( mean_wall * 85 / 3600 ))

    echo ""
    info "Pilot results:"
    info "  Complexes docked: ${n_docked}/5"
    info "  Mean wall-clock:  ${mean_wall}s/complex"
    info "  Estimated full Astex (85): ~${est_full_h}h"
    echo ""

    if [[ "${mean_wall}" -lt "${PILOT_WALL_MIN}" ]] && [[ "${DRY_RUN}" == false ]]; then
        warn "Mean wall-clock ${mean_wall}s < ${PILOT_WALL_MIN}s — unusually fast"
        warn "Possible GA misconfiguration (check num_chromosomes=1000, num_generations=500)"
    elif [[ "${mean_wall}" -gt "${PILOT_WALL_MAX}" ]] && [[ "${DRY_RUN}" == false ]]; then
        warn "Mean wall-clock ${mean_wall}s > ${PILOT_WALL_MAX}s — slower than expected"
        warn "Verify OMP_NUM_THREADS=${N_THREADS} and FLEXAID_SIMD=NEON are active"
    else
        ok "Wall-clock within expected range (${PILOT_WALL_MIN}–${PILOT_WALL_MAX}s/complex)"
    fi

    # Write pilot wall-clock to results dir for downstream use
    echo "${mean_wall}" > "${RESULTS_DIR}/pilot_mean_wall_s"
}

# ─── Full Astex Diverse run (85 complexes) ────────────────────────────────────

run_full_astex() {
    phase "FULL ASTEX DIVERSE RUN (85 complexes)"

    local ast_out="${RESULTS_DIR}/astex_diverse"
    # Always create log dir (infrastructure, not docking — no dry-run guard here)
    mkdir -p "${LOG_DIR}" "${ast_out}" 2>/dev/null || true

    # Write summary CSV header
    if [[ "${DRY_RUN}" == false ]]; then
        echo "complex_id,wall_time_s,top1_score,rmsd_to_crystal,success" > "${SUMMARY_CSV}"
    else
        printf "[DRY]   echo 'complex_id,wall_time_s,top1_score,rmsd_to_crystal,success' > %s\n" \
               "${SUMMARY_CSV}"
    fi

    # In dry-run mode, tee to /dev/null so the pipeline doesn't fail on missing log dirs
    local _astex_log="${LOG_DIR}/astex_full.log"
    [[ "${DRY_RUN}" == true ]] && _astex_log="/dev/null"

    # Prefer benchmark_datasets binary if available (handles download + dock atomically)
    if [[ -n "${DATASET_BIN:-}" ]] && [[ "${TWO_PASS}" == false ]]; then
        info "Using benchmark_datasets binary for full Astex run"
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        run "${DATASET_BIN}" \
            --benchmark   astex \
            --output      "${ast_out}" \
            --threads     "${N_THREADS}" \
            --omp-threads "${OMP_PER_WORKER}" \
            --cache       "${REPO_ROOT}/benchmarks/astex_diverse" \
            --job-timeout-seconds 7200 \
            2>&1 | tee "${_astex_log}"

        # benchmark_datasets already emits the authoritative Astex CSV.
        # Use that file directly so validation sees the actual benchmark rows
        # instead of a reconstructed summary with missing binding_modes.json.
        if [[ "${DRY_RUN}" == false ]]; then
            local dataset_results_csv="${ast_out}/astex_diverse_results.csv"
            if [[ -f "${dataset_results_csv}" ]]; then
                cp "${dataset_results_csv}" "${SUMMARY_CSV}"
            else
                warn "Expected benchmark results CSV not found: ${dataset_results_csv}"
            fi
        fi
        return 0
    fi

    # Fallback: individual FlexAIDdS calls per complex
    info "Falling back to per-complex FlexAIDdS invocations"
    local n_total n_success=0 n_fail=0 n_skip=0

    # Read all complexes from CSV
    local all_pdbs=() all_rmsd_thrs=()
    while IFS=',' read -r pdb_id lig_id res rmsd_thr ref; do
        all_pdbs+=("${pdb_id}")
        all_rmsd_thrs+=("${rmsd_thr}")
    done < <(tail -n +2 "${ASTEX_CSV}")
    n_total="${#all_pdbs[@]}"

    info "Total complexes: ${n_total}"

    local idx=0
    for pdb in "${all_pdbs[@]}"; do
        local rmsd_thr="${all_rmsd_thrs[$idx]}"
        idx=$((idx + 1))

        local struct_info
        struct_info="$(find_structure_files "${pdb}" 2>/dev/null | head -1 || echo 'MISSING')"

        if [[ "${struct_info}" == MISSING* ]] || [[ -z "${struct_info}" ]]; then
            warn "[${idx}/${n_total}] ${pdb}: structures not found (run download.sh)"
            n_skip=$((n_skip + 1))
            echo "${pdb},0,N/A,N/A,0" >> "${SUMMARY_CSV}"
            continue
        fi

        local receptor="${struct_info%%:*}"
        local ligand="${struct_info##*:}"
        local out_dir="${ast_out}/${pdb}"
        local log_file="${LOG_DIR}/astex_${pdb}.log"

        info "[${idx}/${n_total}] Docking ${pdb} ..."

        local dock_ok=true
        if [[ "${TWO_PASS}" == true ]]; then
            dock_complex_twopass "${pdb}" "${receptor}" "${ligand}" "${out_dir}" "${log_file}" || dock_ok=false
        else
            dock_complex "${pdb}" "${receptor}" "${ligand}" "${out_dir}" "${log_file}" || dock_ok=false
        fi

        if [[ "${dock_ok}" == true ]]; then
            local result_line
            result_line="$(extract_result "${out_dir}" "${pdb}" "${rmsd_thr}")"
            echo "${result_line}" >> "${SUMMARY_CSV}"
            local success="${result_line##*,}"
            if [[ "${success}" -eq 1 ]]; then
                n_success=$((n_success + 1))
            fi
            local wall_s
            wall_s="$(cat "${out_dir}/wall_time_s" 2>/dev/null || echo 0)"
            local pass_used="single"
            [[ -f "${out_dir}/final_pass" ]] && pass_used="$(cat "${out_dir}/final_pass")"
            ok "[${idx}/${n_total}] ${pdb}: ${wall_s}s (${pass_used})"
        else
            fail "[${idx}/${n_total}] ${pdb}: docking failed"
            echo "${pdb},0,N/A,N/A,0" >> "${SUMMARY_CSV}"
            n_fail=$((n_fail + 1))
        fi
    done

    local sr=0
    [[ $((n_success + n_fail)) -gt 0 ]] && \
        sr=$(( n_success * 100 / (n_success + n_fail) ))

    echo ""
    ok "Astex Diverse complete"
    info "  Success: ${n_success}/${n_total} (${sr}%)"
    info "  Failed:  ${n_fail}"
    info "  Skipped: ${n_skip} (structures missing)"
    info "  Summary: ${SUMMARY_CSV}"
}

# ─── Phase 2: Astex Non-Native + CASF-2016 ────────────────────────────────────

run_phase2() {
    phase "PHASE 2: ASTEX NON-NATIVE + CASF-2016"
    export SHANNON_TRACE_LEVEL=1   # final H only for non-native cross-docking

    if [[ -n "${DATASET_BIN:-}" ]]; then
        info "Running Astex Non-Native..."
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        run "${DATASET_BIN}" \
            --benchmark   astex_nonnative \
            --output      "${RESULTS_DIR}/astex_nonnative" \
            --threads     "${N_THREADS}" \
            --omp-threads "${OMP_PER_WORKER}" \
            2>&1 | tee "${LOG_DIR}/astex_nonnative.log"

        if [[ -d "${REPO_ROOT}/benchmarks/casf2016" ]]; then
            export SHANNON_TRACE_LEVEL=2
            info "Running CASF-2016..."
            export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
            run "${DATASET_BIN}" \
                --benchmark   casf2016 \
                --output      "${RESULTS_DIR}/casf2016" \
                --threads     "${N_THREADS}" \
                --omp-threads "${OMP_PER_WORKER}" \
                2>&1 | tee "${LOG_DIR}/casf2016.log"
        else
            warn "CASF-2016 data not found — skipping (download to benchmarks/casf2016/)"
        fi
    else
        warn "benchmark_datasets binary not available — Phase 2 requires it"
        warn "Build with: cmake -DENABLE_BENCHMARK_DATASETS=ON ..."
    fi
}

# ─── Final report ─────────────────────────────────────────────────────────────

print_final_report() {
    phase "BENCHMARK COMPLETE"
    info "Timestamp:    ${TIMESTAMP}"
    info "Results dir:  ${RESULTS_DIR}"
    info "Logs dir:     ${LOG_DIR}"
    info "Summary CSV:  ${SUMMARY_CSV}"
    echo ""

    if [[ -f "${SUMMARY_CSV}" ]] && [[ "${DRY_RUN}" == false ]]; then
        python3 - "${SUMMARY_CSV}" <<'PYEOF'
import sys, csv, math

csv_file = sys.argv[1]
rows = []
with open(csv_file) as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

n = len(rows)
if n == 0:
    print("No results in summary CSV")
    sys.exit(0)

wall_times = [float(r['wall_time_s']) for r in rows if r['wall_time_s'] not in ('0','N/A','')]
successes  = [int(r['success']) for r in rows]
n_success  = sum(successes)
sr         = n_success / n * 100 if n > 0 else 0.0

mean_wall = sum(wall_times) / len(wall_times) if wall_times else 0.0
var_wall  = sum((x - mean_wall)**2 for x in wall_times) / max(len(wall_times)-1, 1)
std_wall  = math.sqrt(var_wall)

print(f"  Complexes:     {n}")
print(f"  Success rate:  {n_success}/{n}  ({sr:.1f}%)")
print("  Target:        see manifest validation log")
print("  Status:        see manifest validation log")
print(f"  Mean wall:     {mean_wall:.0f}s ± {std_wall:.0f}s/complex")
PYEOF

        manifest_path="$(benchmark_manifest_path)"
        if [[ -f "${manifest_path}" ]]; then
            info "Validating against manifest: ${manifest_path}"
            python3 "${REPO_ROOT}/scripts/validate_benchmark_results.py" \
                "${SUMMARY_CSV}" \
                --manifest "${manifest_path}" \
                --out-dir "${FIGURES_DIR}" \
                2>&1 | tee "${RESULTS_DIR}/validation.log"
        else
            warn "Manifest not found: ${manifest_path}"
        fi
    fi

    echo ""
    info "Next steps:"
    info "  1. Inspect validation log: ${RESULTS_DIR}/validation.log"
    info "  2. Deposit raw results to Zenodo (DOI required for thesis)"
    echo ""
}

# ─── Entry point ─────────────────────────────────────────────────────────────

main() {
    echo ""
    printf "${BOLD}FlexAIDdS Production Benchmark Runner${NC}\n"
    printf "  SEED=${SEED}  THREADS=${N_THREADS}  DRY_RUN=${DRY_RUN}\n"
    printf "  TIMESTAMP=${TIMESTAMP}\n\n"

    setup_env
    run mkdir -p "${LOG_DIR}" "${RESULTS_DIR}" "${FIGURES_DIR}"
    run mkdir -p "${RESULTS_DIR}"
    run xattr -w com.apple.fileprovider.ignore#P 1 "${RESULTS_DIR}" 2>/dev/null || true
    run touch "${RESULTS_DIR}/.metadata_never_index" 2>/dev/null || true

    # Record run metadata
    if [[ "${DRY_RUN}" == false ]]; then
        {
            echo "timestamp=${TIMESTAMP}"
            echo "seed=${SEED}"
            echo "threads=${N_THREADS}"
            echo "two_pass=${TWO_PASS}"
            echo "benchmark=${BENCHMARK:-astex}"
            echo "git_commit=$(cd "${REPO_ROOT}" && git rev-parse HEAD 2>/dev/null || echo UNKNOWN)"
            echo "hostname=$(hostname)"
            echo "os=$(uname -srm)"
        } > "${RESULTS_DIR}/run_metadata.txt"
    fi

    # ── Route benchmark-specific runs (casf, astex_non_native) ──
    if [[ "${BENCHMARK}" == "casf" || "${BENCHMARK}" == "casf2016" ]]; then
        locate_binary
        locate_dataset_runner
        info "Running CASF-2016 benchmark (two_pass=${TWO_PASS})"
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        if [[ "${TWO_PASS}" == true ]]; then
            run_twopass_dataset "casf2016" "${RESULTS_DIR}"
        else
            run "${DATASET_BIN}" \
                --benchmark   casf2016 \
                --output      "${RESULTS_DIR}" \
                --threads     "${N_THREADS}" \
                --omp-threads "${OMP_PER_WORKER}" \
                --job-timeout-seconds 7200 \
                2>&1 | tee "${LOG_DIR}/casf2016.log" || true
        fi
        print_final_report
        return
    fi

    if [[ "${BENCHMARK}" == "astex_non_native" || "${BENCHMARK}" == "astex_nonnative" ]]; then
        locate_binary
        locate_dataset_runner
        info "Running Astex Non-Native benchmark (two_pass=${TWO_PASS})"
        export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
        if [[ "${TWO_PASS}" == true ]]; then
            run_twopass_dataset "astex_nonnative" "${RESULTS_DIR}"
        else
            run "${DATASET_BIN}" \
                --benchmark   astex_nonnative \
                --output      "${RESULTS_DIR}" \
                --threads     "${N_THREADS}" \
                --omp-threads "${OMP_PER_WORKER}" \
                --job-timeout-seconds 7200 \
                2>&1 | tee "${LOG_DIR}/astex_nonnative.log" || true
        fi
        print_final_report
        return
    fi

    # ── Default: Astex Diverse ──
    preflight

    if [[ "${SKIP_PILOT}" == false ]]; then
        run_pilot
    else
        warn "Pilot calibration skipped (--skip-pilot)"
    fi

    if [[ "${PILOT_ONLY}" == false ]]; then
        thermo_gate
        run_full_astex
        if [[ "${RUN_PHASE2}" == true ]]; then
            run_phase2
        fi
    else
        info "Pilot-only mode — stopping after calibration"
    fi

    print_final_report
}

main "$@"
