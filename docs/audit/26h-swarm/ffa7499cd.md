# Audit: `ffa7499cd` — macOS CPU ctest blocking + Metal CI release gates

| Field | Value |
|-------|--------|
| **Commit** | `ffa7499cdf1464709b7d34760fc1c3c3be588564` |
| **Short** | `ffa7499cd` |
| **Subject** | Fix: Make macOS CPU ctest blocking and add Metal CI release gates |
| **Parent** | `964bec0a2bb5ea6518868a9487288d0ddb81ac5c` (Merge PR #258 Bonhomme Fleet dataset-runner) |
| **AuthorDate** | 2026-07-14 23:29:24 −0400 |
| **Merged via** | PR **#261** → `e2b799495` (“Blocking macOS CPU ctest + Metal CI release gates”) |
| **Files** | 6 changed, **+419 / −19** |
| **Scope** | CI workflows + release-gate docs only (**no `LIB/` / scoring / ranking / engine source**) |
| **Audit focus** | **CI gate strength for science releases** — macOS soft-fail removal, blocking ctest, Metal hosted smoke vs self-hosted full gate, claim hygiene |
| **Verdict** | **CONDITIONAL PASS** for **macOS CPU** science-merge gates; **FAIL-CLOSED intent** for Metal **link/runtime** is documented but **not enforceable on PR merge** at this commit (dispatch-only + incomplete FlexAID-only full gate). Follow-ups (`8c42517bd` / PR #271) harden Metal; infrastructure (`self-hosted-m3` runner) remains the real Metal release blocker. |

---

## 1. Executive summary

This commit is the **first concrete repair** of a real audit hole: on parent, macOS sat in `cxx_core_build` with `allow_failure: true`, **forced `BUILD_TESTING=OFF`**, and was **excluded from the ctest step** (Linux-only `if:`). That meant Apple Silicon / AppleClang regressions could never block a PR, while Metal had **zero** CI surface.

**What this commit correctly does**

1. **Removes** soft-fail `macos-clang` from the Linux-oriented matrix.
2. Adds dedicated job **`macos_cpu_tests`**: `BUILD_TESTING=ON`, Metal **OFF**, **`ctest --output-on-failure --timeout 120`**, and **no** `continue-on-error` / `allow_failure` on the job (verified by job-block parse).
3. Adds hosted **`macos_metal_compile_smoke`**: compiles the five tracked `*.metal` sources to `.air` (optional `.metallib`); soft-skips only when `xcrun metal` is absent (`exit 2` → job success + notice).
4. Adds **`metal-self-hosted.yml`**: fail-closed Metal **link** path on labels `self-hosted` + `self-hosted-m3`, **`workflow_dispatch` only**.
5. Documents the before/after contract in **`docs/CI_RELEASE_GATES.md`**, with pointers from `CONTRIBUTING.md` and `docs/SUPPORT_MATRIX.md`.

**What remains weak for science-release claims at this exact SHA**

| Gap | Severity for science release |
|-----|------------------------------|
| Metal full gate is **manual dispatch**, not a PR required check | **High** if advertising Metal; **OK** if Metal stays “experimental” |
| Hosted Metal job is “blocking” only when metalc exists; otherwise **soft-success** | **Medium** (docs mostly honest; table language can oversell) |
| Self-hosted builds **`FlexAID` only**, not **`FlexAIDdS`** (primary LTO binary) | **High** for production Metal binary claims |
| No `CMakeCache` assert that `FLEXAIDS_USE_METAL` stayed **ON** | **Medium** (CMake can FORCE Metal OFF without OBJCXX) |
| No nm / runtime / GPU kernel proof | **High** for any “Metal-accelerated docking validated” claim |
| Tag `release.yml` still builds macOS with **Metal OFF + `BUILD_TESTING=OFF`** | **Medium** (release artifacts ≠ gated test path) |
| Branch-protection “required checks” are **outside** the repo | **Medium** (YAML blocking ≠ merge-blocking without admin config) |
| Smoke script **skips missing shaders with WARN** if any remain | **Low–Medium** (deletion of one shader can go green) |

**Downstream:** PR #271 (`8c42517bd` / merge `125f76cc8`) hardens the self-hosted gate (dual binary + nm smoke + claim language). That is **out of scope as code of this commit** but is evidence the residuals below were real.

---

## 2. Change inventory (exact commit)

| Path | Status | Role |
|------|--------|------|
| `.github/workflows/ci.yml` | M (+93/−~19 net in hunks) | Drop soft macOS matrix; add `macos_cpu_tests`; add `macos_metal_compile_smoke` |
| `.github/workflows/metal-self-hosted.yml` | **A** (+117) | Dispatch-only full Metal configure/build/link (+ optional ctest) |
| `scripts/ci/metal_compile_smoke.sh` | **A** (+93, mode 100755) | Hosted/self-hosted shared `.metal` → `.air`/`.metallib` smoke |
| `docs/CI_RELEASE_GATES.md` | **A** (+106) | Canonical gate table + Metal hosted vs self-hosted contract |
| `CONTRIBUTING.md` | M (+17) | §2b CI release gates; forbid re-softening `macos_cpu_tests` |
| `docs/SUPPORT_MATRIX.md` | M (+12/−1) | Metal status row + platform CI summary table |

**Not touched:** `LIB/**`, CMake Metal wiring, ranking/CF, GA, StatMech, datasets, claim launchers, `release.yml`.

---

## 3. Parent defect (what was broken)

Parent `ci.yml` `cxx_core_build`:

```yaml
continue-on-error: ${{ matrix.allow_failure || false }}
# matrix include:
- name: macos-clang
  os: macos-15
  allow_failure: true
# configure:
if [ "$RUNNER_OS" = "macOS" ]; then
  CMAKE_ARGS+=(-DBUILD_TESTING=OFF)   # overrides global BUILD_TESTING=ON
fi
# test:
if: matrix.name == 'linux-gcc' || matrix.name == 'linux-clang' || ...
```

**Triple soft-fail on macOS:**

1. Job could red without failing the workflow (`continue-on-error`).
2. Even on green configure/build, **tests were not built**.
3. Even if someone flipped testing on, the ctest step **never selected** `macos-clang`.

Metal: no workflow, no script, no release checklist language.

For a science codebase that ships Homebrew/macOS paths and optional Metal acceleration, parent state was **inadequate for “platform validated” language**.

---

## 4. Gate-by-gate strength analysis

### 4.1 `macos_cpu_tests` — **STRONG (PR CI)**

| Property | At `ffa7499cd` |
|----------|----------------|
| Runner | `macos-15` (GitHub-hosted) |
| `continue-on-error` | **Absent** on job |
| Configure | `BUILD_TESTING=ON`, `FLEXAIDS_USE_METAL=OFF`, CUDA OFF, strict source validation |
| Build | full `cmake --build build --parallel` (all default targets, not FlexAID-only) |
| Test | `ctest --test-dir build --output-on-failure --timeout 120` |
| Intentional exclusion | Metal OFF (correct separation) |

**Strength for science:** Catches AppleClang / libc++ / OpenMP-prefix / macOS-specific compile and unit-test regressions that Linux never sees. Aligns with AGENTS.md “verify with actual execution” for the **macOS CPU** product surface.

**Residuals (CPU gate):**

| ID | Finding | Severity |
|----|---------|----------|
| C1 | `sudo xcode-select -s /Applications/Xcode_16.app \|\| true` — missing Xcode 16 soft-continues | **Low–Medium** (may use older default; C++26 gate may then fail closed at configure — acceptable if fatal) |
| C2 | No job-level `timeout-minutes` (unlike self-hosted 90) | **Low** |
| C3 | Does not build/test with `FLEXAIDS_USE_METAL=ON` (by design) | **Info** |
| C4 | `release.yml` tag builds still use `BUILD_TESTING=OFF` on macOS | **Medium** for “release artifact was tested” claims (orthogonal workflow) |
| C5 | Required-status-check enrollment is not in-repo | **Medium** process |

**Verdict C:** **PASS** as a real blocking workflow job for macOS CPU unit coverage, subject to C5.

---

### 4.2 `macos_metal_compile_smoke` — **WEAK–MODERATE (shader syntax only)**

**Script contract (`scripts/ci/metal_compile_smoke.sh`):**

| Exit | Meaning | Hosted job treatment |
|------|---------|----------------------|
| 0 | ≥1 shader compiled; optional metallib | success, `status=ok` |
| 2 | no `xcrun -f metal` | **success**, `status=skipped` + notice |
| 1 | compile failure | job fails |

**Shaders listed (complete set of repo `*.metal` at this commit):**

1. `LIB/CavityDetect/CavityDetect.metal`
2. `LIB/ShannonThermoStack/shannon_metal.metal`
3. `LIB/TurboQuant.metal`
4. `LIB/MetalRMSD.metal`
5. `LIB/gpu_fast_optics_metal.metal`

**What it proves:** Metal **shader frontend** accepts those sources (and optionally packs `.metallib`).

**What it does not prove:**

- OBJCXX bridge compile (`metal_eval.mm`, `*MetalBridge.mm`, `gpu_fast_optics_metal.mm`, `tencm_metal.mm`)
- Framework link (`Metal`, `MetalKit`, `Foundation`)
- CMake custom commands / `FLEXAIDS_USE_METAL=ON` path
- Runtime GPU dispatch, numerical parity vs CPU, docking CF path
- That **all** listed shaders still exist (missing file → WARN + skip; success if any remain)

| ID | Finding | Severity |
|----|---------|----------|
| M1 | Soft-skip exit 2 makes the job **non-informative** when hosted metalc is absent — still “green” | **Medium** for “Metal is CI-gated on every PR” marketing |
| M2 | Docs table labels smoke as a **blocking** PR gate “if toolchain present” — correct in prose, easy to misread in summary tables | **Low–Medium** (SUPPORT_MATRIX row is clearer: “shaders” only) |
| M3 | Missing-shader WARN continues; deleting Shannon metal alone can leave 4/5 green | **Medium** for shader inventory integrity |
| M4 | Smoke compiles `gpu_fast_optics_metal.metal`, but CMake at this commit **does not** register a metallib custom target for it (only `.mm` bridge) — smoke is **stricter than CMake** for that file (good), yet production may not ship that metallib via the same path | **Low–Medium** (CMake asymmetry, pre-existing) |
| M5 | Shannon `.metallib` CMake path is under **`FlexAIDdS`**, while FlexAID Metal block links Cavity/Turbo/RMSD metallibs — dual product surface | **Medium** for which binary is “Metal complete” |

**Verdict M (hosted):** **PASS as smoke**, **FAIL as Metal release evidence**. Docs mostly state this; residual overclaim risk if operators equate green PR CI with Metal validation.

---

### 4.3 `metal-self-hosted.yml` — **MODERATE intent, INCOMPLETE product coverage at this SHA**

| Property | At `ffa7499cd` |
|----------|----------------|
| Trigger | **`workflow_dispatch` only** (not push/PR/tag) |
| Runner | `[self-hosted, self-hosted-m3]` |
| Timeout | 90 minutes |
| Soft-fail | None (`continue-on-error` absent) |
| Preflight | `xcrun` + metalc executable + `ctypes` load of `Metal.framework` |
| Configure | `FLEXAIDS_USE_METAL=ON`, `BUILD_TESTING=ON` |
| Build | **`--target FlexAID` only** |
| Extra | Re-runs `metal_compile_smoke.sh`; optional ctest (default `true`) |
| Deps | `brew install … \|\| true` then require `cmake`/`ninja` |

**Fail-closed properties that work:**

- Missing metalc → job fails (unlike hosted).
- Configure/build/link failure → job fails.
- Missing runner → job **queues**, never green-soft-passes (documented).

**Science / product holes at this commit:**

| ID | Finding | Severity |
|----|---------|----------|
| S1 | **Not a PR gate** — Metal regressions can merge forever without this job ever running | **High** for Metal release claims; **acceptable** for “Metal experimental” |
| S2 | Builds **FlexAID only**; **`FlexAIDdS`** (LTO production binary with its **own** Metal source block + Shannon metallib dependency) is **unproven** | **High** for Homebrew/campaign “Metal ON” binaries that are FlexAIDdS |
| S3 | No post-configure check that `FLEXAIDS_USE_METAL:BOOL=ON` in `CMakeCache.txt` (OBJCXX absence **FORCE**s Metal OFF at top of `CMakeLists.txt`) — a CPU-only build could theoretically be mislabeled if CMake did not fatal on metalc (here metalc preflight reduces risk) | **Medium** |
| S4 | Optional ctest is CPU unit tests under Metal-enabled **link**, not Metal kernel correctness tests | **Medium** (honest as “build+link”, not “GPU validated”) |
| S5 | No `nm` / symbol smoke for `metal_eval_*` (added later in #271) | **Medium** |
| S6 | `brew install … \|\| true` can hide package failure; later steps may fail closed or use stale tools | **Low** |
| S7 | Without a registered runner, checklist item “run Metal full gate” is **impossible** — gate is paper until ops registers `self-hosted-m3` | **Critical infrastructure** (not a YAML logic bug) |

**Verdict S:** **Architecture correct (dispatch + fail-closed + labeled runner)**; **coverage incomplete** for the science/production binary set at this SHA. Follow-up hardening is justified.

---

### 4.4 Documentation contract — **STRONG intent, process-dependent enforcement**

`docs/CI_RELEASE_GATES.md` (this commit) correctly:

- Tables blocking vs soft jobs.
- Before/after macOS audit matrix.
- Separates hosted shader smoke vs self-hosted full gate.
- States soft-skip **does not authorize Metal release claims**.
- Gives local reproduction commands.
- Release checklist: green `macos_cpu_tests` + self-hosted Metal **or** mark Metal experimental.

`CONTRIBUTING.md` forbids re-adding soft-fail to `macos_cpu_tests` without issue + waiver.

`SUPPORT_MATRIX.md` keeps Metal **Experimental** and points at gates doc.

| ID | Finding | Severity |
|----|---------|----------|
| D1 | Checklist is human process, not an automated release job | **Medium** |
| D2 | No cross-link from `release.yml` / tag publish path to this checklist | **Medium** |
| D3 | Language is careful enough that residual risk is **operator non-compliance**, not silent overclaim in the gate doc itself | **Info** |

---

## 5. Science-release claim matrix (what may be said after this commit)

Assume PR CI is green including `macos_cpu_tests`, and Metal full gate was **not** run (default).

| Claim | Allowed at this commit? | Why |
|-------|-------------------------|-----|
| “macOS CPU build + unit tests are CI-gated (blocking).” | **Yes**, if `macos_cpu_tests` is required in branch protection **or** observed green on the SHA | Job design is sound |
| “Linux and macOS CPU ctest both block merges.” | **Yes** (Linux matrix + new macOS job), subject to required-check config | Real improvement vs parent |
| “Metal shaders are compile-checked on every PR.” | **Only if** hosted metalc present; else soft-skip | Do not claim without smoke `status=ok` |
| “Metal acceleration is CI release-validated.” | **No** | No PR link gate; no runtime |
| “FlexAIDdS Metal link validated.” | **No** at this SHA (FlexAID-only self-hosted target) | Fixed later in #271 for workflow; still needs runner |
| “Tag `v*` macOS artifact is unit-tested.” | **No** via `release.yml` (still `BUILD_TESTING=OFF`) | Separate gap |
| CF/contact-function scoring or ΔG / ranking changed | **N/A — no** | CI/docs only |

This commit **does not** alter ranking, clustering, CF scoring proxy, StatMech, or ensemble thermodynamics. Science impact is **platform trust**, not docking scores.

---

## 6. Interaction with later history (context only)

| SHA | Relation |
|-----|----------|
| `e2b799495` | Merge of this commit via PR #261 |
| `e752df7aa` / `8c42517bd` | Harden self-hosted gate (FlexAID+FlexAIDdS, nm smoke, claim language); see audits `e752df7aa.md`, `125f76cc8.md` |
| `4d51e413e` | Revert of premature land on a **docs** branch — not undoing #261 macOS CPU gate |
| `125f76cc8` | PR #271 merge of hardened Metal gate onto main |

**Do not score this commit as if it already contains #271.** Residuals S2/S5 are **valid findings of `ffa7499cd`**, later partially addressed.

---

## 7. Security / hygiene

- No secrets, no machine-absolute paths in new scripts (ROOT derived from `BASH_SOURCE`).
- Pin `actions/checkout@9c091bb…` (v6) consistent with existing CI.
- Self-hosted workflow is **dispatch-only** — reduces unsolicited runner load (good); also means no ambient Metal security surface from PR forks (good for self-hosted safety).
- `metal_compile_smoke.sh` uses fixed relative paths under repo ROOT — no user-controlled path injection beyond `OUT_DIR` arg (CI-controlled).

---

## 8. Ranking / reproducibility / tests

| Axis | Impact |
|------|--------|
| **Ranking/scoring** | **None** (no engine changes) |
| **Reproducibility of science numbers** | **Indirect positive**: macOS unit-test regressions less likely to ship silently |
| **Reproducibility of Metal performance claims** | **Not yet** — full gate unenforced without runner + dispatch |
| **Tests adequate for stated goal** | **Yes for macOS CPU blocking**; **No for Metal production validation** |

---

## 9. Findings summary (prioritized)

### P0 / Critical (for Metal production claims only)

1. **S7 + S1:** No automated Metal link gate on PR/tag; full gate needs live `self-hosted-m3`. Without runner + green dispatch URL, **forbid** “Metal release-validated” language.
2. **S2:** Self-hosted builds FlexAID only — **FlexAIDdS Metal link untested** at this SHA.

### P1 / High–Medium

3. **M1:** Hosted smoke soft-success when metalc missing — PR green ≠ shader compile proven.
4. **S3/S5:** No CMakeCache Metal-ON assert / nm symbol smoke (later partially fixed).
5. **C4:** `release.yml` still ships untested macOS artifacts (`BUILD_TESTING=OFF`, Metal OFF).
6. **C5:** Ensure GitHub required checks include exact job name `macOS CPU tests (blocking)`.

### P2 / Medium–Low

7. **M3:** Fail closed if any listed shader path is missing (not only if all missing).
8. **M4/M5:** Align CMake metallib targets with smoke list and FlexAID vs FlexAIDdS Metal source parity (engine/build follow-up, not this commit’s job).
9. **C1/C2:** Harden Xcode select fail-closed; optional job timeout.

### Intentional / acceptable

10. Metal OFF on CPU job; soft AVX-512 / Windows bindings remain experimental.
11. Metal remains **Experimental** in SUPPORT_MATRIX for 1.0 CPU guarantees.

---

## 10. Recommendations (no source edits in this audit)

1. **Keep** `macos_cpu_tests` blocking forever; treat any reintroduction of soft-fail as a process incident (already written in CONTRIBUTING).
2. **Register** durable `self-hosted` + `self-hosted-m3` and store green run URLs next to any Metal claim in release notes.
3. **Adopt #271-class** dual-binary + nm smoke (if not already on main when reading this).
4. **Tighten smoke:** exit 1 on any missing path in `SHADERS` array; optionally upload `status` as a check annotation consumers can query.
5. **Wire release checklist:** either document that tags must cite green `macos_cpu_tests` SHA, or add a release job that fails if that check is red.
6. **Never** equate hosted Metal smoke green/skip with thermodynamic or docking performance claims (AGENTS.md CF vs thermo discipline is orthogonal but same claim hygiene).

---

## 11. Verdict

| Dimension | Result |
|-----------|--------|
| Addresses parent macOS soft-fail + no-ctest audit? | **YES — fixed in workflow** |
| Metal “blocking gate” for science Metal releases? | **PARTIAL — process + infra, not PR-enforced** |
| Safe for docking science kernel? | **YES** (no ranking/CF/ΔG code touch) |
| Overall | **CONDITIONAL PASS** |

**Ship/merge assessment:** Correct and necessary CI/docs commit for platform honesty. Treat **macOS CPU** as materially strengthened. Treat **Metal** as **documented experimental with a fail-closed *optional* full gate**, not as production-validated, until a labeled runner produces a green full-gate run for the release SHA (and preferably until FlexAIDdS + symbol smoke land — later commits).

---

## 12. Evidence anchors (repo paths at commit)

| Artifact | Path |
|----------|------|
| Blocking CPU job | `.github/workflows/ci.yml` → job `macos_cpu_tests` |
| Hosted Metal smoke | `.github/workflows/ci.yml` → job `macos_metal_compile_smoke` |
| Smoke script | `scripts/ci/metal_compile_smoke.sh` |
| Full Metal gate | `.github/workflows/metal-self-hosted.yml` |
| Gate contract | `docs/CI_RELEASE_GATES.md` |
| Contributing waiver rule | `CONTRIBUTING.md` §2b |
| Support matrix | `docs/SUPPORT_MATRIX.md` |
| Parent soft-fail pattern | `ffa7499cd^:.github/workflows/ci.yml` matrix `macos-clang` + `BUILD_TESTING=OFF` |

**Audit method:** `git show` full patch + parent workflow + CMake Metal blocks at `ffa7499cd`; job-block parse for absence of `continue-on-error`; cross-check later PR #271 history for residual confirmation. **No source edits. No build/ctest run required for this CI-only audit scope.**
