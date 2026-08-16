#!/usr/bin/env bash
# =============================================================================
# reproduce_astex85.sh — One-command Astex Diverse 85 benchmark reproduction
#
# PURPOSE
#   Lets any reviewer clone this repo and run a BLIND Astex Diverse 85
#   self-docking campaign from scratch:
#
#       bash scripts/reproduce_astex85.sh
#
#   Default arm = METHODOLOGY.md §3 autonomous / docking power:
#     FLEXAIDDS_SEED_ELITISM=0
#     FLEXAIDDS_NATIVE_SEED_FRAC=0
#   NATIVE_SEED_FRAC is a dead knob on today's DatasetRunner path (always emits
#   seed_fraction: 0.0). The live oracle lever is SEED_ELITISM=1 (_INI.pdb).
#   The default path kills both.
#
#   Optional:
#       bash scripts/reproduce_astex85.sh --oracle-ceiling
#     Sets SEED_ELITISM=1 and NATIVE_SEED_FRAC=0.90 and prints
#     "ORACLE CEILING — not docking power". Do not cite that arm as S1.
#
#   This repository publishes no Astex-85 success rate. The former 94.1%
#   (80/85) target was a seeded oracle ceiling; this script does not compare
#   against it. Observed S1 is printed with the RMSD instrument name and
#   N_denominator=85.
#
#   The script builds the LTO-optimised FlexAIDdS binary and benchmark runner
#   from the current git HEAD, writes a machine-stamped provenance.json that
#   matches the arm actually run, and prints a summary table.
# SUPPORTED PLATFORMS
#   macOS 13+  (Apple Silicon or Intel, AppleClang ≥ 16 / Xcode 16)
#   Linux      (GCC ≥ 14 or Clang ≥ 18, x86-64 or aarch64)
#
# WINDOWS
#   Native Windows is not supported — the benchmark runner requires POSIX
#   process management (fork/exec/waitpid/killpg).
#   See REPRODUCIBILITY_WINDOWS.md for the supported WSL2 path.
#
# DEPENDENCIES
#   cmake ≥ 3.28, python3 ≥ 3.9, git, curl
#   macOS: Xcode 16 command-line tools + brew install boost eigen
#   Linux: apt install libboost-dev libeigen3-dev libssl-dev
#   obabel: optional — only needed if you supply custom ligand formats
#
# EXPECTED RUNTIME  (wall clock with --threads 4)
#   Apple M-series (M3 Pro / M4 Pro):  ~45–60 minutes
#   Linux x86-64 (32-core):            ~20–35 minutes
#   The bottleneck is a single long-running target (1OF6, ~34 min on Apple M).
#   Build adds ~5–10 minutes before docking starts.
#
# OUTPUTS
#   $FLEXAIDDS_ICLOUD/results/working/reproduce_astex85_.../   (preferred when FLEXAIDDS_ICLOUD set)
#   or ~/FlexAIDdS_reviewer_benchmark/
#     provenance.json                 ← machine stamp + env snapshot
#     astex_crossdock_85_results.csv  ← per-target scores + RMSD
#     astex_crossdock_85_report.md    ← publication-format summary table
#     astex_crossdock_85_summary.csv  ← top-line metrics
#     stdout.log / stderr.log         ← full run transcript
#     <PDB>/                          ← per-target pose files
#   NOTE: Active writes target working/ subdir; finalize to archived/ via safe_archive_to_icoud.py
#   (protects against iCloud Drive sync risks: delays, placeholders, conflicted copies).
#
# HISTORICAL REFERENCE (commit 8196829) — WITHDRAWN, not a script target
#   Former 94.1% (80/85) was produced with SEED_ELITISM=1 and
#   NATIVE_SEED_FRAC=0.90 = an ORACLE CEILING. METHODOLOGY.md §0 forbids
#   reporting it as docking power. Do not cite it. This script does not
#   compare against 80/85.
#   Ref binary SHA256 (Apple M-series, withdrawn run):
#     6d899e6351e347abf97f2e5b664ffd2cba853c599a561f5213ccf2777df47d5c
# =============================================================================
set -euo pipefail

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[0;33m'
BLU='\033[0;34m'; CYN='\033[0;36m'; BLD='\033[1m'; RST='\033[0m'
info()  { printf "${BLU}[INFO]${RST}  %s\n" "$*"; }
ok()    { printf "${GRN}[OK]${RST}    %s\n" "$*"; }
warn()  { printf "${YLW}[WARN]${RST}  %s\n" "$*"; }
die()   { printf "${RED}[FATAL]${RST} %s\n" "$*" >&2; exit 1; }
banner(){ printf "\n${BLD}${CYN}%s${RST}\n%s\n\n" "$1" "$(printf '═%.0s' $(seq 1 ${#1}))"; }

# ── Locate repo root ──────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Configuration ─────────────────────────────────────────────────────────────
# Prefer iCloud Drive for results/logs/benchmark outputs (via FLEXAIDDS_ICLOUD).
# For active reproduction run, use a timestamped subdir under results/working/
# (safer than safe: isolates in-progress writes from iCloud sync conflicts/lag).
# After success, the script already calls safe_archive_to_icoud.py into /archived/
# if FLEXAIDDS_ICLOUD is set. Full compatibility with overrides.
# iCloud comment: writes under working/ recommended; use safe_archive for finals
# to avoid sync races (see safe_archive_to_icoud.py for verified copy + manifest).
if [[ -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
    TS_R=$(date +%Y%m%d_%H%M%S)
    OUTPUT_DIR="${FLEXAIDDS_ICLOUD}/results/working/reproduce_astex85_${TS_R}"
else
    OUTPUT_DIR="${HOME}/FlexAIDdS_reviewer_benchmark"
fi
BUILD_DIR="${REPO_ROOT}/build_reproduce"
ASTEX_DIR="${REPO_ROOT}/benchmarks/astex_diverse/astex_diverse"
DATASETS_DIR="${REPO_ROOT}/benchmarks/datasets"
PUBLISHED_COMMIT="8196829f35a2bf065919ccd1508f62f00059895d"
PUBLISHED_BINARY_SHA256="6d899e6351e347abf97f2e5b664ffd2cba853c599a561f5213ccf2777df47d5c"

# Thread counts — override via env if needed
BENCH_THREADS="${FLEXAIDDS_BENCH_THREADS:-4}"
OMP_PER_WORKER="${FLEXAIDDS_OMP_THREADS:-2}"

# ── Parse flags ───────────────────────────────────────────────────────────────
FORCE=0
SKIP_BUILD=0
ORACLE_CEILING=0
for arg in "$@"; do
    case "$arg" in
        --force)            FORCE=1 ;;
        --skip-build)       SKIP_BUILD=1 ;;
        --oracle-ceiling)   ORACLE_CEILING=1 ;;
        -h|--help)
            head -70 "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \?//'
            exit 0 ;;
        *)
            die "Unknown flag: ${arg} (supported: --force --skip-build --oracle-ceiling)"
            ;;
    esac
done

# Default = METHODOLOGY.md §3 autonomous / blind. Comments-only is not enough:
# these assignments are what get exported and written to provenance.json.
SEED_ELITISM=0
NATIVE_SEED_FRAC=0
ARM_LABEL=blind
CLAIM_VALIDITY=blind
if [[ "${ORACLE_CEILING}" -eq 1 ]]; then
    SEED_ELITISM=1
    NATIVE_SEED_FRAC=0.90
    ARM_LABEL=oracle_ceiling
    CLAIM_VALIDITY=oracle_ceiling_not_docking_power
fi

# =============================================================================
# STEP 0 — Platform check
# =============================================================================
banner "FlexAIDdS Astex 85 Benchmark Reproducer"

if [[ "${ORACLE_CEILING}" -eq 1 ]]; then
    warn "ORACLE CEILING — not docking power"
    warn "SEED_ELITISM=1 injects _INI.pdb. NATIVE_SEED_FRAC=${NATIVE_SEED_FRAC} is a dead knob on DatasetRunner (seed_fraction is always 0.0)."
    warn "Do not cite this arm as S1 / docking power."
else
    info "Arm: blind (SEED_ELITISM=0 NATIVE_SEED_FRAC=0). METHODOLOGY.md §3 autonomous."
fi

if [[ "$(uname)" == "MINGW"* ]] || [[ "$(uname)" == "CYGWIN"* ]] || \
   [[ "${OS:-}" == "Windows_NT" ]]; then
    die "Native Windows is not supported. See REPRODUCIBILITY_WINDOWS.md for WSL2 instructions."
fi

PLATFORM="$(uname -s)-$(uname -m)"
info "Platform: ${PLATFORM}"
info "Repo:     ${REPO_ROOT}"
info "Output:   ${OUTPUT_DIR}"
info "Build:    ${BUILD_DIR}"
info "Threads:  ${BENCH_THREADS} workers × ${OMP_PER_WORKER} OMP threads"

# =============================================================================
# STEP 1 — Dependency checks
# =============================================================================
banner "Step 1 — Checking dependencies"

check_cmd() {
    local cmd="$1" hint="$2"
    if command -v "$cmd" &>/dev/null; then
        ok "$cmd found: $(command -v "$cmd")"
    else
        die "$cmd not found. ${hint}"
    fi
}

check_cmd cmake  "Install cmake ≥ 3.28 (brew install cmake  /  apt install cmake)"
check_cmd python3 "Install Python ≥ 3.9"
check_cmd git    "Install git"
check_cmd curl   "Install curl (brew install curl  /  apt install curl)"
check_cmd make   "Install make (xcode-select --install  /  apt install build-essential)"

# Optional
if command -v obabel &>/dev/null; then
    ok "obabel found (optional, for custom ligand formats)"
else
    warn "obabel not found — OK for Astex 85 (all ligands are SDF)"
fi

# cmake version check
CMAKE_VER="$(cmake --version | head -1 | awk '{print $3}')"
CMAKE_MAJOR="${CMAKE_VER%%.*}"
CMAKE_MINOR="${CMAKE_VER#*.}"; CMAKE_MINOR="${CMAKE_MINOR%%.*}"
if [[ "$CMAKE_MAJOR" -lt 3 ]] || { [[ "$CMAKE_MAJOR" -eq 3 ]] && [[ "$CMAKE_MINOR" -lt 28 ]]; }; then
    die "cmake ≥ 3.28 required (found ${CMAKE_VER})"
fi
ok "cmake ${CMAKE_VER}"

# Check Boost and Eigen3
if [[ "$(uname)" == "Darwin" ]]; then
    if ! pkg-config --exists eigen3 2>/dev/null && \
       ! [[ -d /usr/local/include/eigen3 ]] && \
       ! [[ -d /opt/homebrew/include/eigen3 ]] && \
       ! [[ -d /opt/homebrew/Cellar/eigen ]]; then
        warn "Eigen3 not found in standard paths. Install: brew install eigen"
    else
        ok "Eigen3 appears available"
    fi
    if ! [[ -d /usr/local/include/boost ]] && \
       ! [[ -d /opt/homebrew/include/boost ]] && \
       ! pkg-config --exists boost 2>/dev/null; then
        warn "Boost not found in standard paths. Install: brew install boost"
    else
        ok "Boost appears available"
    fi
else
    for hdr in "/usr/include/eigen3/Eigen/Core" "/usr/local/include/eigen3/Eigen/Core"; do
        [[ -f "$hdr" ]] && { ok "Eigen3 found"; break; }
    done
    for hdr in "/usr/include/boost/version.hpp" "/usr/local/include/boost/version.hpp"; do
        [[ -f "$hdr" ]] && { ok "Boost found"; break; }
    done
fi

# Astex Diverse dataset must be in the repo
if [[ ! -d "${ASTEX_DIR}" ]] || [[ ! -f "${ASTEX_DIR}/1G9V/1G9V_apo.pdb" ]]; then
    die "Astex Diverse structures not found at ${ASTEX_DIR}.
Make sure you cloned with full history (no --depth 1 --filter=blob:none)."
fi
ASTEX_COUNT="$(find "${ASTEX_DIR}" -name "*_apo.pdb" | wc -l | tr -d ' ')"
ok "Astex Diverse structures: ${ASTEX_COUNT}/85 receptor PDBs found"
if [[ "$ASTEX_COUNT" -lt 85 ]]; then
    die "Expected 85 apo PDBs, found ${ASTEX_COUNT}. Re-clone the repo."
fi

# =============================================================================
# STEP 2 — Git provenance
# =============================================================================
banner "Step 2 — Git provenance"

cd "${REPO_ROOT}"
HEAD_COMMIT="$(git rev-parse HEAD 2>/dev/null || echo "unknown")"
HEAD_SHORT="${HEAD_COMMIT:0:7}"
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")"

info "HEAD commit:  ${HEAD_COMMIT}"
info "Branch:       ${BRANCH}"
info "Published commit: ${PUBLISHED_COMMIT}"

if [[ "${HEAD_COMMIT}" != "${PUBLISHED_COMMIT}" ]]; then
    warn "You are at commit ${HEAD_SHORT}, not the published commit ${PUBLISHED_COMMIT:0:7}."
    warn "Results may differ slightly from published numbers."
    warn "To reproduce exactly: git checkout ${PUBLISHED_COMMIT}"
fi

# Check for uncommitted changes
if ! git diff --quiet HEAD 2>/dev/null; then
    warn "Working tree has uncommitted changes — results may not be exactly reproducible."
fi

# =============================================================================
# STEP 3 — Build (LTO + native)
# =============================================================================
banner "Step 3 — Building FlexAIDdS + benchmark_datasets (LTO)"

FLEXAIDDS_BIN="${BUILD_DIR}/FlexAIDdS"
BENCHMARK_BIN="${BUILD_DIR}/benchmark_datasets"

if [[ "$SKIP_BUILD" -eq 1 ]] && [[ -x "$FLEXAIDDS_BIN" ]] && [[ -x "$BENCHMARK_BIN" ]]; then
    warn "--skip-build: reusing existing binaries"
elif [[ -x "$FLEXAIDDS_BIN" ]] && [[ -x "$BENCHMARK_BIN" ]] && [[ "$FORCE" -eq 0 ]]; then
    warn "Binaries already exist in ${BUILD_DIR}."
    warn "Pass --force to rebuild, or --skip-build to reuse."
    read -r -p "  Reuse existing build? [Y/n] " ans
    ans="${ans:-Y}"
    if [[ "${ans}" =~ ^[Nn] ]]; then
        rm -rf "${BUILD_DIR}"
    else
        ok "Reusing existing build."
    fi
fi

if [[ ! -x "$FLEXAIDDS_BIN" ]] || [[ ! -x "$BENCHMARK_BIN" ]]; then
    mkdir -p "${BUILD_DIR}"

    # Detect parallel jobs
    if [[ "$(uname)" == "Darwin" ]]; then
        NJOBS="$(sysctl -n hw.physicalcpu 2>/dev/null || echo 4)"
    else
        NJOBS="$(nproc 2>/dev/null || echo 4)"
    fi
    info "Configuring with cmake (${NJOBS} parallel jobs)..."

    # Base flags — Metal enabled automatically on macOS, CUDA off by default
    CMAKE_FLAGS=(
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_FLEXAIDDS_FAST=ON
        -DFLEXAIDS_USE_AVX2=ON
        -DFLEXAIDS_USE_OPENMP=ON
        -DFLEXAIDS_USE_CUDA=OFF
        -DFLEXAIDS_BUILD_CORE=ON
    )

    cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" "${CMAKE_FLAGS[@]}" \
        2>&1 | tee "${BUILD_DIR}/cmake_configure.log" \
        | grep -E "^(--|CMake Error|FATAL|Warning.*required)" || true

    info "Building targets: FlexAIDdS benchmark_datasets..."
    cmake --build "${BUILD_DIR}" \
        --target FlexAIDdS benchmark_datasets \
        -j "${NJOBS}" \
        2>&1 | tee "${BUILD_DIR}/cmake_build.log" \
        | tail -5

    if [[ ! -x "$FLEXAIDDS_BIN" ]]; then
        die "Build failed — FlexAIDdS binary not found. See ${BUILD_DIR}/cmake_build.log"
    fi
    if [[ ! -x "$BENCHMARK_BIN" ]]; then
        die "Build failed — benchmark_datasets binary not found. See ${BUILD_DIR}/cmake_build.log"
    fi
    ok "Build complete."
fi

# SHA256 of freshly built binary
sha256() {
    local f="$1"
    if command -v shasum &>/dev/null; then
        shasum -a 256 "$f" | awk '{print $1}'
    else
        sha256sum "$f" | awk '{print $1}'
    fi
}
BINARY_SHA256="$(sha256 "${FLEXAIDDS_BIN}")"
info "FlexAIDdS SHA256: ${BINARY_SHA256}"
if [[ "${BINARY_SHA256}" == "${PUBLISHED_BINARY_SHA256}" ]]; then
    ok "Binary SHA256 matches published reference (same platform/toolchain)."
else
    warn "Binary SHA256 differs from published reference."
    warn "  Expected: ${PUBLISHED_BINARY_SHA256}"
    warn "  Got:      ${BINARY_SHA256}"
    warn "This is normal on a different OS/compiler — bit-identical results are"
    warn "not expected across platforms. Numerical results should agree within"
    warn "floating-point rounding (RMSD differences < 0.01 Å)."
fi

# =============================================================================
# STEP 4 — Generate portable benchmark JSON
# =============================================================================
banner "Step 4 — Generating benchmark pair list"

JSON_OUT="${OUTPUT_DIR}/benchmark_astex85_reviewer.json"
mkdir -p "${OUTPUT_DIR}"

python3 - "${ASTEX_DIR}" "${JSON_OUT}" <<'PYEOF'
import json, os, sys
from pathlib import Path

astex_dir = Path(sys.argv[1])
out_path  = Path(sys.argv[2])

CODES = [
    "1G9V","1GM8","1GPK","1HNN","1HP0","1HQ2","1IA1","1IGJ","1J3J","1JD0",
    "1JJE","1K3U","1KE5","1KZK","1L2S","1L7F","1LPZ","1M2Z","1MEH","1MQ6",
    "1N1M","1N2J","1N2V","1N46","1NAV","1OF1","1OF6","1OPK","1OQ5","1OWE",
    "1P2Y","1P62","1PMN","1Q1G","1Q41","1Q4G","1R1H","1R55","1R58","1R9O",
    "1S19","1S3V","1SG0","1SJ0","1SQ5","1T40","1T46","1T9B","1TT1","1TW6",
    "1TZ8","1U1C","1U4D","1UML","1UNL","1UOU","1V0P","1V48","1V4S","1VCJ",
    "1W1P","1W2G","1X8X","1XM6","1XOZ","1Y6B","1Y6R","1YGC","1YQY","1YV3",
    "1YVF","1YWR","1Z95","2BM2","2BR1","2BSM","2BYS","2C3I","2CET","2CGR",
    "2D3U","2GBP","2HB1","2HR7","2J62",
]

pairs = []
missing = []
for i, code in enumerate(CODES):
    d = astex_dir / code
    rec = d / f"{code}_apo.pdb"
    lig = d / f"{code}_ligand.sdf"
    site = d / f"{code}_binding_site.pdb"
    if not rec.is_file() or not lig.is_file():
        missing.append(code)
        continue
    entry = {
        "index": i,
        "receptor_id": code,
        "ligand_id": code,
        "receptor_pdb": str(rec),
        "ligand_sdf": str(lig),
    }
    if site.is_file():
        entry["oracle_site_pdb"] = str(site)
    pairs.append(entry)

if missing:
    print(f"[WARN] Missing structures for: {', '.join(missing)}", file=sys.stderr)

doc = {
    "schema_version": 1,
    "name": "astex_native_85_reviewer",
    "description": "85-case Astex Diverse native self-docking — reviewer-generated paths.",
    "n_pairs": len(pairs),
    "oracle_mode": True,
    "astex_diverse_dir": str(astex_dir),
    "pairs": pairs,
}
out_path.write_text(json.dumps(doc, indent=2))
print(f"[OK] Wrote {len(pairs)} pairs to {out_path}", file=sys.stderr)
PYEOF

PAIR_COUNT="$(python3 -c "import json; d=json.load(open('${JSON_OUT}')); print(d['n_pairs'])")"
if [[ "$PAIR_COUNT" -lt 85 ]]; then
    die "Only ${PAIR_COUNT}/85 pairs generated. Missing structures in ${ASTEX_DIR}."
fi
ok "Generated ${PAIR_COUNT}-pair benchmark JSON: ${JSON_OUT}"

# =============================================================================
# STEP 5 — Write provenance.json
# =============================================================================
banner "Step 5 — Writing provenance"

PROV_FILE="${OUTPUT_DIR}/provenance.json"
TIMESTAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

python3 - <<PYEOF
import json, os, platform, sys
from pathlib import Path

prov = {
    "benchmark": "astex_diverse_85_self_docking",
    "description": "Reproduced from FlexAIDdS git repo — reviewer run",
    "timestamp_utc": "${TIMESTAMP}",
    "git_commit": "${HEAD_COMMIT}",
    "git_commit_short": "${HEAD_SHORT}",
    "git_branch": "${BRANCH}",
    "published_git_commit": "${PUBLISHED_COMMIT}",
    "binary_path": "${FLEXAIDDS_BIN}",
    "binary_sha256": "${BINARY_SHA256}",
    "published_binary_sha256": "${PUBLISHED_BINARY_SHA256}",
    "benchmark_datasets_binary": "${BENCHMARK_BIN}",
    "platform": "${PLATFORM}",
    "python_version": platform.python_version(),
    "cmake_flags": (
        "-DCMAKE_BUILD_TYPE=Release -DBUILD_FLEXAIDDS_FAST=ON "
        "-DFLEXAIDS_USE_AVX2=ON -DFLEXAIDS_USE_OPENMP=ON -DFLEXAIDS_USE_CUDA=OFF"
    ),
    "arm": "${ARM_LABEL}",
    "claim_validity": "${CLAIM_VALIDITY}",
    "docking_config": {
        "FLEXAIDDS_THERMO":             "1",
        "FLEXAIDDS_T_EFF":              "0.596",
        "FLEXAIDDS_TENCOM_SCALE":       "1.0",
        "FLEXAIDDS_RESTARTS":           "7",
        "FLEXAIDDS_PARALLEL_RESTARTS":  "1",
        "FLEXAIDDS_CONSENSUS_SCORER":   "1",
        "FLEXAIDDS_SEED_ELITISM":       "${SEED_ELITISM}",
        "FLEXAIDDS_N_ELITE":            "1",
        "FLEXAIDDS_BUDGET_SCALE":       "1",
        "FLEXAIDDS_SOFTCORE_WAL":       "1",
        "FLEXAIDDS_SOFTCORE_FLOOR":     "0.5",
        "FLEXAIDDS_T_HOT":              "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC":   "${NATIVE_SEED_FRAC}",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
    },
    "workers": ${BENCH_THREADS},
    "omp_threads_per_worker": ${OMP_PER_WORKER},
    "output_dir": "${OUTPUT_DIR}",
    "pair_json": "${JSON_OUT}",
}
Path("${PROV_FILE}").write_text(json.dumps(prov, indent=2))
print("[OK] Wrote ${PROV_FILE}")
PYEOF
ok "Provenance written: ${PROV_FILE}"

# =============================================================================
# STEP 6 — Run benchmark
# =============================================================================
banner "Step 6 — Running Astex 85 benchmark (${BENCH_THREADS} workers)"

printf "  Docking config:\n"
printf "    arm=%s  FLEXAIDDS_THERMO=1  T_EFF=0.596  TENCOM_SCALE=1.0\n" "${ARM_LABEL}"
printf "    7 restarts, SEED_ELITISM=%s  NATIVE_SEED_FRAC=%s\n" "${SEED_ELITISM}" "${NATIVE_SEED_FRAC}"
if [[ "${ORACLE_CEILING}" -eq 1 ]]; then
    printf "    ${YLW}ORACLE CEILING — not docking power${RST}\n"
else
    printf "    blind / METHODOLOGY.md §3 autonomous (not an oracle ceiling)\n"
fi
printf "    %s workers × %s OMP threads\n" "${BENCH_THREADS}" "${OMP_PER_WORKER}"
printf "\n  This will take approximately 45–60 minutes on Apple M-series.\n"
printf "  Streaming output to:\n"
printf "    stdout: %s/stdout.log\n" "${OUTPUT_DIR}"
printf "    stderr: %s/stderr.log\n\n" "${OUTPUT_DIR}"

export FLEXAIDDS_BINARY="${FLEXAIDDS_BIN}"
export FLEXAIDDS_DATA_DIR="${BUILD_DIR}"
export FLEXAIDDS_ORACLE_SITE_DIR="${ASTEX_DIR}"

# Thermodynamic engine
export FLEXAIDDS_THERMO=1
export FLEXAIDDS_T_EFF=0.596
export FLEXAIDDS_TENCOM_SCALE=1.0

# Search configuration. Default = blind (METHODOLOGY.md §3). Oracle knobs are
# variables set at flag-parse time — not hardcoded 1 / 0.90 on the default path.
export FLEXAIDDS_RESTARTS=7
export FLEXAIDDS_PARALLEL_RESTARTS=1
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
export FLEXAIDDS_CONSENSUS_SCORER=1
export FLEXAIDDS_SEED_ELITISM="${SEED_ELITISM}"
export FLEXAIDDS_N_ELITE=1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_SOFTCORE_WAL=1
export FLEXAIDDS_SOFTCORE_FLOOR=0.5
export FLEXAIDDS_T_HOT=500
export FLEXAIDDS_NATIVE_SEED_FRAC="${NATIVE_SEED_FRAC}"
export FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1

# Priority targets (known hard cases — run first to surface failures early)
export FLEXAIDDS_PRIORITY_TARGETS=1HNN,1N2V,1TW6

TSTART="$(date +%s)"

"${BENCHMARK_BIN}" \
    --benchmark "crossdock_json:${JSON_OUT}" \
    --output    "${OUTPUT_DIR}" \
    --threads   "${BENCH_THREADS}" \
    --omp-threads "${OMP_PER_WORKER}" \
    --job-timeout-seconds 7200 \
    2> >(tee "${OUTPUT_DIR}/stderr.log" >&2) \
    | tee "${OUTPUT_DIR}/stdout.log"

EXIT_CODE="${PIPESTATUS[0]}"
TEND="$(date +%s)"
ELAPSED=$(( TEND - TSTART ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

# =============================================================================
# STEP 7 — Summary
# =============================================================================
banner "Step 7 — Results summary"

RESULTS_CSV="${OUTPUT_DIR}/astex_crossdock_85_results.csv"
if [[ ! -f "${RESULTS_CSV}" ]]; then
    die "Results CSV not found: ${RESULTS_CSV}
Benchmark may have crashed. Check ${OUTPUT_DIR}/stderr.log"
fi

python3 - "${RESULTS_CSV}" "${PROV_FILE}" \
    "${OUTPUT_DIR}" "${ARM_LABEL}" <<'PYEOF'
import csv, sys
from pathlib import Path

results_csv = Path(sys.argv[1])
prov_file   = Path(sys.argv[2])
output_dir  = sys.argv[3]
arm_label   = sys.argv[4]
N_DENOM = 85

rows = list(csv.DictReader(results_csv.open()))

def pick_rmsd(r):
    hung = (r.get("rmsd_hungarian") or "").strip()
    xtal = (r.get("rmsd_to_crystal") or "").strip()
    if hung:
        return float(hung), "rmsd_hungarian"
    if xtal:
        return float(xtal), "rmsd_to_crystal"
    return 9999.0, "missing"

picked = [pick_rmsd(r) for r in rows]
instruments = sorted({name for _, name in picked})
instrument = ",".join(instruments) if instruments else "missing"

def succ_pair(pair):
    value, name = pair
    return name != "missing" and value <= 2.0

success  = [r for r, p in zip(rows, picked) if succ_pair(p)]
near     = [r for r, p in zip(rows, picked) if not succ_pair(p) and p[0] < 2.5]
failures = [r for r, p in zip(rows, picked) if not succ_pair(p) and p[0] >= 2.5]
n_rows   = len(rows)
n_ok     = len(success)
rate     = 100.0 * n_ok / N_DENOM

rmsds    = [p[0] for p in picked if p[1] != "missing"]
mean_r   = sum(rmsds) / len(rmsds) if rmsds else 0.0
srt      = sorted(rmsds)
median_r = srt[len(srt)//2] if srt else 0.0

wall_times = []
for r in rows:
    t = r.get("wall_time_s","") or r.get("time_s","") or r.get("wall_s","")
    try:
        wall_times.append(float(t))
    except ValueError:
        pass

BOLD = '\033[1m'; YLW = '\033[0;33m'
CYN  = '\033[0;36m'; RST = '\033[0m'

print(f"\n{BOLD}{CYN}{'═'*60}{RST}")
print(f"{BOLD}  FlexAIDdS  Astex Diverse 85  —  Reproduction Results{RST}")
print(f"{BOLD}{CYN}{'═'*60}{RST}")
if arm_label == "oracle_ceiling":
    print(f"\n  {YLW}ORACLE CEILING — not docking power{RST}")
    print(f"  {YLW}SEED_ELITISM=1 injects _INI.pdb. Do not cite this arm as S1.{RST}")
else:
    print(f"\n  Arm: {arm_label} (METHODOLOGY.md §3 autonomous). Not compared to 80/85.")

print(f"\n  {'Metric':<40} {'Observed':>12}")
print(f"  {'─'*40} {'─'*12}")
print(f"  {'Rows in CSV':<40} {n_rows:>12}")
print(f"  {'N_denominator':<40} {N_DENOM:>12}")
print(f"  {'RMSD instrument':<40} {instrument:>12}")
print(f"  {'Successful (rank-0 RMSD <= 2.0 Å)':<40} {n_ok:>12}")
print(f"  {'Observed S1 (n_ok / 85)':<40} {rate:>11.1f}%")
print(f"  {'Mean RMSD (Å)':<40} {mean_r:>12.2f}")
print(f"  {'Median RMSD (Å)':<40} {median_r:>12.2f}")
if wall_times:
    total_seq = sum(wall_times)
    print(f"  {'Total sequential docking time':<40} {total_seq/3600:>11.1f}h")
print(f"\n  {'Near-misses (>2.0–2.5 Å)':<40}", end="")
if near:
    print(f"  {', '.join(r['pdb_id'] for r in near)}")
else:
    print("  none")
print(f"  {'Failures (≥ 2.5 Å)':<40}", end="")
if failures:
    print(f"  {', '.join(r['pdb_id'] for r in failures)}")
else:
    print("  none")

print(f"\n{BOLD}{CYN}{'═'*60}{RST}")
print(f"\n  {YLW}No published comparator. Former 94.1% (80/85) is a withdrawn{RST}")
print(f"  {YLW}oracle ceiling, not docking power. Observed S1 is not a claim.{RST}")
print(f"  Provenance: {prov_file}")
print(f"  Results:    {output_dir}")
PYEOF

info "Elapsed wall time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
info "Full results: ${OUTPUT_DIR}"

if [[ "$EXIT_CODE" -ne 0 ]]; then
    warn "benchmark_datasets exited with code ${EXIT_CODE}."
    warn "Partial results may be in ${OUTPUT_DIR}."
    exit "${EXIT_CODE}"
fi

# =============================================================================
# Optional: Safe archive to iCloud (if FLEXAIDDS_ICLOUD is configured)
# =============================================================================
if [[ -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
    ICLOUD_ARCHIVE="${FLEXAIDDS_ICLOUD}/archived"
    mkdir -p "$ICLOUD_ARCHIVE"
    echo "[INFO] Attempting safe verified archive to iCloud..."
    python3 "${REPO_ROOT}/scripts/safe_archive_to_icoud.py" \
        --source "$OUTPUT_DIR" \
        --dest "$ICLOUD_ARCHIVE" \
        --keep-local || warn "iCloud archive step had issues (data still in $OUTPUT_DIR)"
fi
