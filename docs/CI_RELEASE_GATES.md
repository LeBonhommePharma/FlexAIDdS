# CI Release Gates

This document defines **what must be green** before claiming a platform or
backend is release-validated. It addresses the audit finding that macOS was
allowed to fail, C++ tests were disabled on macOS, and Metal had no blocking
gate.

## Blocking gates (required on every PR / push)

| Gate | Workflow / job | Runner | Fail behavior |
|:--|:--|:--|:--|
| Pure Python results | `ci.yml` → `pure_python_results` | `ubuntu-latest` | Blocking |
| C++ core Linux GCC | `ci.yml` → `cxx_core_build` (`linux-gcc`) | `ubuntu-latest` | Blocking + **ctest** |
| C++ core Linux Clang | `ci.yml` → `cxx_core_build` (`linux-clang`) | `ubuntu-latest` | Blocking + **ctest** |
| C++ MPI / ASan Linux | `ci.yml` → `cxx_core_build` | `ubuntu-latest` | Blocking + **ctest** |
| **macOS CPU C++ tests** | `ci.yml` → `macos_cpu_tests` | `macos-15` | **Blocking + ctest** (no `continue-on-error`) |
| **Metal shader compile smoke** | `ci.yml` → `macos_metal_compile_smoke` | `macos-15` | Blocking if Metal toolchain present; soft-skip (job success with notice) only if `xcrun metal` is absent |
| Python bindings smoke (Linux) | `ci.yml` → `python_bindings_smoke` | `ubuntu-latest` | Blocking |
| Repo hygiene + skill tests | `ci.yml` → `repo_hygiene` | `ubuntu-latest` | Blocking |
| DatasetRunner dry-run | `ci.yml` → `flexaid_docking_datasetrunner` | `ubuntu-latest` | Blocking |

### Before / after (macOS audit)

| Behavior | Before | After |
|:--|:--|:--|
| macOS in `cxx_core_build` | Matrix entry with `allow_failure: true` | Removed from soft-fail matrix |
| macOS `BUILD_TESTING` | Forced `OFF` in configure | Dedicated job forces **`ON`** |
| macOS `ctest` | Not run | **`macos_cpu_tests` runs ctest** |
| Metal | No CI gate | Hosted **shader compile smoke** + self-hosted full gate |

## Soft / experimental (non-blocking)

| Gate | Why soft |
|:--|:--|
| `linux-gcc-avx512` | Architecture may be unavailable on GHA hardware; `allow_failure: true` |
| Windows Python bindings | Known flaky / deep MSVC issues; `allow_failure: true` |
| Coverage, sanitizers, TSAN, perf workflows | Extra signal, not merge gates unless required by a release checklist |

## Metal: hosted smoke vs self-hosted full gate

### Hosted (every PR): shader compile only

- Job: `macos_metal_compile_smoke` in `.github/workflows/ci.yml`
- Script: `scripts/ci/metal_compile_smoke.sh`
- Compiles tracked `*.metal` sources to `.air` (and `.metallib` when `metallib` exists)
- Does **not** full-link Metal OBJCXX bridges into `FlexAID` / `FlexAIDdS`
- Runtime: seconds to a few minutes — not a docking campaign

If GitHub-hosted `macos-15` lacks the Metal compiler component, the job prints a
clear skip notice and exits 0 so Linux/macOS CPU gates stay green.

**Soft-skip is not a Metal release claim.** A green or soft-skipped
`macos_metal_compile_smoke` only means “hosted shader compile was attempted;
when metalc was missing the job skipped cleanly.” It does **not** validate
Metal OBJCXX link, `metal_eval_*` / `metal_rmsd` symbols in production
binaries, GPU runtime, or docking performance.

### Self-hosted full gate (Metal release claims)

- Workflow: `.github/workflows/metal-self-hosted.yml`
- Trigger: **`workflow_dispatch` only** (manual or release checklist)
  - GitHub UI: **Actions → “Metal full gate (self-hosted M3)” → Run workflow**
  - Inputs: `run_ctest` (default true), `nm_symbol_smoke` (default true)
- Runner labels (both required): **`self-hosted`**, **`self-hosted-m3`**
- Builds: **`FlexAID` and `FlexAIDdS`** with `FLEXAIDS_USE_METAL=ON`
- Optional **`nm` symbol smoke** for `metal_eval_get_capabilities`, other
  `metal_eval_*` entry points, and `metal_rmsd` (C++ namespace) on **both**
  binaries
- Behavior: **fail-closed**
  - Missing Metal compiler → job fails
  - CMake configure with `FLEXAIDS_USE_METAL=ON` fails or resolves Metal OFF → job fails
  - Link of either binary fails → job fails
  - Missing required Metal symbols → job fails (when `nm_symbol_smoke=true`)
  - Optional `ctest` (default on) failure → job fails
  - No `continue-on-error`

#### Runner registration

```bash
# On the Apple Silicon host (after installing the GitHub Actions runner):
# Labels must include BOTH of:
#   self-hosted
#   self-hosted-m3
#
# Verify for this repo (empty list ⇒ no M3 runner is registered; do not claim one exists):
gh api repos/LeBonhommePharma/FlexAIDdS/actions/runners \
  --jq '.runners[] | {name,status,labels:[.labels[].name]}'
```

Without a matching runner, a dispatched workflow **stays queued** until a runner
appears or the run is cancelled — it never soft-passes. **Do not claim an M3
self-hosted runner is online unless the API above lists one with label
`self-hosted-m3`.**

## Allowed release language vs forbidden claims

Use this table for GitHub releases, Homebrew caveats, README badges, papers,
and verbal product claims. “Green” means the named gate passed for the
**same commit/tag** being shipped.

### Allowed without self-hosted Metal green

| You may say | Only if |
|:--|:--|
| “macOS CPU build and unit tests are CI-gated.” | `macos_cpu_tests` green on the release commit |
| “Hosted CI compiles Metal shaders when the Metal toolchain is present.” | `macos_metal_compile_smoke` status is `ok` (not soft-skip) on that commit |
| “Metal acceleration is **experimental / best-effort** on Apple Silicon; not a production guarantee.” | Always allowed; preferred default when self-hosted gate is not green |
| “Metal full link/runtime validation requires a self-hosted M3 runner (`metal-self-hosted.yml`).” | Always allowed (factual process note) |
| “Linux GCC/Clang `ctest` is green.” | Corresponding `cxx_core_build` matrix cells green |

### Forbidden without green `Metal full gate (self-hosted M3)`

Do **not** use any of the following (or equivalent) unless
`metal-self-hosted.yml` is green for that commit/tag:

| Forbidden claim | Why |
|:--|:--|
| “Metal is release-validated / production-ready.” | Needs full self-hosted link (+ optional ctest) |
| “Metal-accelerated docking is CI-proven on Apple Silicon.” | Hosted jobs do not link/run Metal docking |
| “`FlexAID` / `FlexAIDdS` ship with verified Metal GPU kernels.” | Needs OBJCXX link + symbol smoke (and preferably runtime) |
| “Homebrew `--with-metal` is fully validated by GitHub Actions.” | Hosted GHA is shader-smoke only; full gate is self-hosted |
| “Soft-skipped Metal smoke means Metal works.” | Soft-skip = metalc missing; **not** a Metal pass |
| “We have a self-hosted M3 runner online.” | Only if `gh api …/actions/runners` lists `self-hosted-m3` |

### Allowed **with** green self-hosted Metal full gate

| You may say | Scope of the green run |
|:--|:--|
| “Metal OBJCXX bridges **link** into `FlexAID` and `FlexAIDdS` on Apple Silicon (self-hosted M3 gate).” | Configure + build both targets with `FLEXAIDS_USE_METAL=ON` |
| “Linked binaries expose `metal_eval_get_capabilities` / `metal_rmsd` symbols (`nm` smoke).” | `nm_symbol_smoke=true` (default) |
| “Metal-enabled unit tests passed on the self-hosted M3 runner.” | Only if `run_ctest=true` and ctest green |
| “Metal full gate green for tag/commit X.” | Cite the Actions run URL for that ref |

Still **not** automatic from the full gate alone (unless you add separate
evidence): full Astex / campaign docking throughput claims, numerical parity
vs CPU for every kernel, or third-party hardware beyond the registered runner.

## What still needs a self-hosted M3

| Claim | Hosted GHA sufficient? | Required gate |
|:--|:--|:--|
| macOS CPU build + unit tests | Yes (`macos_cpu_tests`) | PR CI |
| Metal `.metal` → `.air` compile | Usually yes (when Xcode Metal toolchain present) | PR CI smoke |
| Metal OBJCXX bridge **link** into `FlexAID` **and** `FlexAIDdS` | **No** (treat as self-hosted) | `metal-self-hosted.yml` |
| `metal_eval_*` / `metal_rmsd` symbols in both binaries | **No** | self-hosted + `nm` smoke |
| Metal **runtime** correctness / GPU kernels | **No** | self-hosted + local validation |
| Full docking campaign on Apple Silicon | **No** | Local / fleet benchmarks (not this gate) |

## Release checklist (LP release managers)

Copy this into the release issue or PR checklist.

### Always (every release)

- [ ] Green blocking jobs on `ci.yml` for the commit/tag, including **`macos_cpu_tests`** (blocking; no `continue-on-error`).
- [ ] Linux `ctest` clean (`linux-gcc` and `linux-clang` at minimum).
- [ ] `python3 scripts/check_repo_hygiene.py` clean.
- [ ] Release notes use only **allowed** language from the table above for platforms not fully gated.

### macOS CPU claims

- [ ] Cite green `macos_cpu_tests` for the release SHA.
- [ ] Do **not** imply Metal from macOS CPU green alone.

### Metal claims

- [ ] Decide claim level:
  - **Experimental only** → no self-hosted run required; use experimental wording.
  - **Link / production Metal** → dispatch `Metal full gate (self-hosted M3)` on the release SHA.
- [ ] Confirm a runner exists: `gh api repos/LeBonhommePharma/FlexAIDdS/actions/runners --jq '.runners[] | {name,status,labels:[.labels[].name]}'` shows `self-hosted-m3`.
- [ ] Green self-hosted run for that SHA builds **both** `FlexAID` and `FlexAIDdS`.
- [ ] Prefer leave `nm_symbol_smoke=true` and `run_ctest=true` (defaults).
- [ ] Paste Actions run URL into the release notes / checklist.
- [ ] Never treat hosted `macos_metal_compile_smoke` soft-skip as Metal validation.

### Explicit non-claims (unless separate evidence)

- [ ] No RMSD-only “success” without PoseBusters when citing docking benchmarks.
- [ ] No “true ΔG” language for CF-only scoring (see `AGENTS.md` scientific guardrails).

## Local reproduction

```bash
# macOS CPU (matches macos_cpu_tests)
cmake -B build-macos-cpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_CUDA=OFF \
  -DFLEXAIDS_USE_METAL=OFF \
  -DBUILD_TESTING=ON \
  -DCMAKE_PREFIX_PATH="$(brew --prefix)" \
  -DOpenMP_ROOT="$(brew --prefix libomp)"
cmake --build build-macos-cpu --parallel
ctest --test-dir build-macos-cpu --output-on-failure

# Metal shader smoke (hosted-equivalent)
bash scripts/ci/metal_compile_smoke.sh /tmp/metal-smoke

# Metal full link (self-hosted-equivalent)
cmake -B build-metal -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_METAL=ON \
  -DBUILD_TESTING=ON \
  -DBUILD_FLEXAIDDS_FAST=ON
cmake --build build-metal --parallel --target FlexAID FlexAIDdS
# Optional symbol smoke (matches metal-self-hosted.yml nm step)
nm -gU build-metal/FlexAID build-metal/FlexAIDdS | grep -E 'metal_eval_get_capabilities|metal_rmsd'
```
