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
- Does **not** full-link Metal OBJCXX bridges into `FlexAID` (that path needs a stable Apple Silicon + Xcode Metal toolchain environment)
- Runtime: seconds to a few minutes — not a docking campaign

If GitHub-hosted `macos-15` lacks the Metal compiler component, the job prints a clear skip notice and exits 0 so Linux/macOS CPU gates stay green. That **does not** authorize Metal release claims.

### Self-hosted full gate (Metal release claims)

- Workflow: `.github/workflows/metal-self-hosted.yml`
- Trigger: **`workflow_dispatch` only** (manual or release checklist)
- Runner labels: **`self-hosted`**, **`self-hosted-m3`**
- Behavior: **fail-closed**
  - Missing Metal compiler → job fails
  - CMake `FLEXAIDS_USE_METAL=ON` configure/build/link failure → job fails
  - Optional `ctest` (default on) failure → job fails
  - No `continue-on-error`

Register an Apple Silicon (M-series) machine as a GitHub Actions self-hosted runner with both labels. Without a matching runner, the workflow stays queued until a runner appears or the run is cancelled — it never soft-passes.

## What still needs a self-hosted M3

| Claim | Hosted GHA sufficient? | Required gate |
|:--|:--|:--|
| macOS CPU build + unit tests | Yes (`macos_cpu_tests`) | PR CI |
| Metal `.metal` → `.air` compile | Usually yes (when Xcode Metal toolchain present) | PR CI smoke |
| Metal OBJCXX bridge **link** into `FlexAID` | **No** (treat as self-hosted) | `metal-self-hosted.yml` |
| Metal **runtime** correctness / GPU kernels | **No** | self-hosted + local validation |
| Full docking campaign on Apple Silicon | **No** | Local / fleet benchmarks (not this gate) |

## Release checklist (minimum)

1. Green blocking jobs on `ci.yml` for the commit/tag (including **`macos_cpu_tests`**).
2. Linux `ctest` clean (`linux-gcc` and `linux-clang` at minimum).
3. If advertising Metal acceleration in release notes:
   - Green `Metal full gate (self-hosted M3)` for that commit/tag, **or**
   - Explicitly mark Metal as experimental with no production guarantee (see `docs/SUPPORT_MATRIX.md`).
4. No secrets / machine-absolute paths: `python3 scripts/check_repo_hygiene.py`.

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
  -DBUILD_TESTING=ON
cmake --build build-metal --parallel --target FlexAID
```
