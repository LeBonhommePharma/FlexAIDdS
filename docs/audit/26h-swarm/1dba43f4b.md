# Audit: `1dba43f4b` — Merge PR #260 Homebrew Metal link (Metal bridges on flexaid_core)

| Field | Value |
|-------|-------|
| **Short** | `1dba43f4b` |
| **Full** | `1dba43f4b3191ffbd079a13611564a26f192a390` |
| **Subject** | Merge pull request #260 from LeBonhommePharma/fix/homebrew-metal-link |
| **PR title** | Fix: Homebrew `--with-metal` link (Metal bridges on `flexaid_core`) |
| **Author / committer** | Louis-Philippe Morency (`LeBonhommePharma`) / GitHub merge bot |
| **AuthorDate** | 2026-07-15 00:31:56 -0400 |
| **Parents** | `d9bee5a38` (main: CI default branch master→main) · `99e5b35eb` (PR tip) |
| **PR commits (first-parent..second-parent)** | `5eb1be79e` → `283d9ee19` → `3e059594b` → `5a3b95430` → `12f3a6a08` → `99e5b35eb` |
| **PR URL** | https://github.com/LeBonhommePharma/FlexAIDdS/pull/260 |
| **Delta vs first parent** | **+246 / −294** across **7 files** (net −48 lines; CMake de-duplication dominates) |
| **Scope of this audit** | Merge quality + CMake/Homebrew Metal linkage correctness. **No source edits** in this session. |
| **Verdict (merge hygiene)** | **CLEAN merge** — no conflict markers; first-parent tree delta equals PR tip for all 7 paths. |
| **Verdict (fix quality @ merge)** | **ACCEPT with follow-ups** — link-time root cause is correctly fixed; residual packaging/docs gaps remain (metallib install, doc drift, weak attach-helper guard). |

---

## Summary (2–4 sentences)

PR #260 fixes a real Homebrew macOS link failure: under `FLEXAIDS_USE_METAL=ON`, `flexaid_core` OBJECT TUs (`gaboom.cpp`, `FOPTICS.cpp`, `hardware_detect.cpp`, `UnifiedHardwareDispatch.cpp`) call `metal_eval_*` / `metal_rmsd::*` while the OBJCXX implementations lived only on `FlexAID` / `FlexAIDdS` (and partially on `benchmark_datasets`). Moving all Metal `.mm` bridges + frameworks onto `flexaid_core` **PUBLIC**, compiling metallibs once, and stripping per-executable duplicate `target_sources` is the correct membership model for multi-target + LTO links. Formula/docs changes correctly refuse stable `v2.0.2 --with-metal`, force `head` onto `main`, and document Homebrew 6 tap + `install -s --HEAD` semantics. Residual risks are **runtime metallib packaging** (absolute build-dir paths, formula does not install `.metallib`), a fragile `flexaids_attach_metal_dispatch_runtime` duplicate-guard, and a leftover “`head` tracks `master`” line in `docs/INSTALLATION.md`.

## Severity: MEDIUM

---

## 1. Executive summary

| Concern | Assessment |
|---------|------------|
| Link-time undefined Metal symbols on brew | **Fixed** by attaching bridges to `flexaid_core` |
| Architecture (single membership + shared metallib build) | **Sound** for OBJECT-library consumers |
| Homebrew 6 install UX (tap, HEAD, reinstall) | **Improved** |
| Stable bottle path at this merge tip | **Intentionally blocked** for Metal (`odie` unless `--HEAD`) until a later tag |
| Ranking / CF / thermo / GA budget | **Unchanged** (build + packaging only) |
| Runtime GPU shader availability after `brew install` | **Still incomplete** (pre-existing metallib path model; not solved here) |

Immediate follow-up in history (out of this merge’s tree but important for operators): `600e95ba9` / v2.0.3 removes the temporary refuse-stable Metal gate and bumps formula/Python versions to ship this fix as stable.

---

## 2. Merge topology & hygiene

### 2.1 Graph

```
d9bee5a38  Update: CI default branch master → main          (first parent)
    \
     \     5eb1be79e  Fix: Homebrew install docs (Homebrew 6+)
      \       |
       \     283d9ee19  Fix: Align Homebrew/PyPI install docs + pypi-release
        \       |
         \     3e059594b  Fix: Link Metal OBJCXX bridges via flexaid_core   ← core fix
          \       |
           \     5a3b95430  Fix: Route --with-metal through HEAD until fix ships
            \       |
             \     12f3a6a08  Fix: Document install -s --HEAD (not reinstall)
              \       |
               \     99e5b35eb  Fix: Point Formula head at main
                \       /
                 1dba43f4b  Merge PR #260  ← audited
```

| Role | SHA | Subject |
|------|-----|---------|
| Merge | `1dba43f4b` | Merge PR #260 homebrew-metal-link |
| Parent1 (ours / main) | `d9bee5a38` | CI default branch master → main |
| Parent2 (theirs / PR tip) | `99e5b35eb` | Formula head → main after Metal fix |
| Core engineering commit | `3e059594b` | Metal bridges on `flexaid_core` |

### 2.2 Files changed (first-parent `d9bee5a38...1dba43f4b`)

| Path | Δ | Role |
|------|---|------|
| `CMakeLists.txt` | +94 / −255 | Centralize Metal bridges + metallibs on `flexaid_core`; drop FlexAID/FlexAIDdS/cavity/benchmark duplicate lists |
| `LIB/CMakeLists.txt` | +7 / −4 | Comment + keep PUBLIC `FLEXAIDS_USE_METAL` on core; remove wrong “.mm only on FlexAIDdS” note |
| `Formula/flexaidds.rb` | +59 / −13 | `head`→`main`; refuse stable `--with-metal`; shrink brew CMake surface; caveats |
| `docs/INSTALLATION.md` | +49 / −14 | Tap-based install; Homebrew vs pip matrix; HEAD semantics |
| `README.md` | +8 / −5 | Tap install; native vs Python split |
| `python/README.md` | +8 / −1 | Analysis package vs native tools |
| `.github/workflows/pypi-release.yml` | +21 / −2 | sdist version must match `pyproject.toml` (no hardcoded `2.0.0`) |

**No C++ docking/scoring/thermo source TUs changed.** No ranking, GA, PoseBusters, DatasetRunner, or claim-path edits.

### 2.3 Merge quality scorecard

| Check | Result |
|-------|--------|
| Conflict markers in merge tree | **None** (`git grep '<<<<<<'` clean on merge SHA for touched paths) |
| First-parent path set equals PR tip | **Yes** (7 files only) |
| Secrets / `.env` / absolute user paths | **None** |
| GPL / forbidden deps | **None** |
| License line on formula | **Apache-2.0** retained |
| Science / ranking surface | **Untouched** |

---

## 3. Root cause (pre-merge) — verified against code

`flexaid_core` is an **OBJECT** library (`LIB/CMakeLists.txt` `add_library(flexaid_core OBJECT …)`). Core TUs under `FLEXAIDS_USE_METAL` reference Metal entry points:

| TU (in `flexaid_core`) | Symbols / headers |
|------------------------|-------------------|
| `LIB/gaboom.cpp` | `metal_eval.h` → `metal_eval_init` / `metal_eval_batch` |
| `LIB/hardware_detect.cpp` | `metal_eval_get_capabilities` |
| `LIB/UnifiedHardwareDispatch.cpp` | `metal_eval_get_capabilities` |
| `LIB/FOPTICS.cpp` | `MetalRMSDBridge.h` → `metal_rmsd::*` |

Pre-parent (`d9bee5a38` / first parent tree) attached implementations only via:

- `target_sources(FlexAID PRIVATE LIB/metal_eval.mm … MetalRMSDBridge.mm …)`
- Duplicate full list on `FlexAIDdS`
- Partial list on `benchmark_datasets`
- `cavity_detect_cli` only got CavityDetect bridge pieces

Any target that **links `flexaid_core`** and pulls those Metal call sites (notably `cavity_detect_cli`, and multi-target LTO links where the linker does not pull another executable’s private `.mm` objects) failed with:

```text
ld: symbol(s) not found for architecture arm64
# metal_eval_* / metal_rmsd::*
```

This matches the PR description and formula `odie` text. **Diagnosis is correct.**

---

## 4. Fix design (post-merge) — deep read

### 4.1 Single attachment point on `flexaid_core`

Root `CMakeLists.txt` (merge tree ~L271–411):

```cmake
set(_FLEXAIDS_METAL_BRIDGES
    LIB/metal_eval.mm
    LIB/MetalRMSDBridge.mm
    LIB/ShannonThermoStack/ShannonMetalBridge.mm
    LIB/tENCoM/tencm_metal.mm
    LIB/gpu_fast_optics_metal.mm
    LIB/CavityDetect/CavityDetectMetalBridge.mm
    LIB/TurboQuantMetalBridge.mm
)
target_sources(flexaid_core PRIVATE ${_FLEXAIDS_METAL_BRIDGES})
target_link_libraries(flexaid_core PUBLIC
    ${METAL_LIBRARY} ${FOUNDATION_LIBRARY} ${METALKIT_LIBRARY})
target_compile_definitions(flexaid_core PUBLIC
    FLEXAIDS_USE_METAL FLEXAIDS_HAS_METAL_SHANNON FLEXAIDS_HAS_METAL_TENCM
    CAVITY_METALLIB_PATH=… TURBOQUANT_METALLIB_PATH=…
    METALRMSD_METALLIB_PATH=… SHANNON_METALLIB_PATH=…)
add_dependencies(flexaid_core CavityDetectMetal TurboQuantMetal MetalRMSDMetal ShannonEntropyMetal)
```

**Why this works for OBJECT libs (CMake ≥ 3.12 / project requires 3.28):**

1. `target_sources(flexaid_core PRIVATE …)` adds `.mm` objects into the same OBJECT set as `gaboom.o` / `FOPTICS.o`.
2. Consumers `target_link_libraries(… PRIVATE flexaid_core)` pull those objects into the final link.
3. `PUBLIC` framework linkage propagates Metal/Foundation/MetalKit to dependents.
4. `PUBLIC` `FLEXAIDS_USE_METAL` keeps call-site `#ifdef`s consistent across TUs that inherit core usage requirements.
5. Metallib custom commands run **once**; FlexAIDdS no longer builds parallel `*_ds.air` / `*_ds.metallib` trees (large duplication removed — the −255 lines in root CMake).

OBJCXX is enabled at `project()` on Apple when `xcrun -f clang++` succeeds; `set_source_files_properties(… LANGUAGE OBJCXX … -fno-objc-arc)` is set on the bridge list. This matches prior per-exe settings.

### 4.2 Consumer cleanups

| Target | Change |
|--------|--------|
| `FlexAID` | No longer privately lists Metal `.mm`; inherits via `flexaid_core` |
| `FlexAIDdS` | Same; status message only |
| `cavity_detect_cli` | Drops private Cavity bridge + frameworks; optional `add_dependencies(… CavityDetectMetal)` only |
| `benchmark_datasets` | Drops private bridge list; depends on shared metallib custom targets |
| `flexaids_attach_metal_dispatch_runtime` | Still used for **non-core** targets (`benchmark_dispatch`, `test_hardware_detect_dispatch`, `test_unified_dispatch`); early-out if `LINK_LIBRARIES` matches `flexaid_core` |

### 4.3 Formula surface (merge tip)

- `head "…", branch: "main"` (was `master`) — aligned with parent1’s default-branch rename.
- **If** `build.with?("metal") && !build.head?` → `odie` with install `-s --HEAD --with-metal` instructions. Correct for **v2.0.2** tarball still lacking the membership fix.
- Extra CMake OFF flags for brew: `ENABLE_CAVITY_DETECT_CLI`, `ENABLE_BENCHMARK_DATASETS`, `ENABLE_DUAL_ASSEMBLY_TOOL`, `ENABLE_DIFT_TOOL` — shrinks build surface / disk; also **avoids** building the auxiliary targets that historically exposed the bug under brew (see F4).
- Caveats: native vs Python, tap path, Homebrew 6 `--HEAD` is install-only.

### 4.4 Docs / CI side changes

- README / INSTALLATION / python README: replace raw formula URL installs with `brew tap lebonhommepharma/flexaidds …` (Homebrew 6+).
- `pypi-release.yml`: parse `[project] version` from `python/pyproject.toml` and assert equality with `flexaidds.__version__` instead of hardcoding `"2.0.0"`.

---

## 5. Findings

### F1. Metal bridge membership fix is correct and necessary (INFO — good)

- **Evidence:** Pre-merge `target_sources` only on FlexAID/FlexAIDdS; core TUs reference `metal_eval` / `metal_rmsd` under `FLEXAIDS_USE_METAL`; post-merge `_FLEXAIDS_METAL_BRIDGES` on `flexaid_core` PRIVATE + frameworks PUBLIC; `LIB/CMakeLists.txt` comment updated to forbid “implementations only on executables.”
- **Why it matters:** Without this, `brew … --with-metal` (and any multi-target Metal build that enables cavity/benchmark) is link-broken. LTO on `FlexAIDdS` (`-flto`) made undefined-ref failures more visible.
- **Fix recommendation:** None for membership model. Keep future Metal bridges on `flexaid_core` (or a dedicated `flexaid_metal` OBJECT/STATIC that core PUBLIC-links), never reintroduce per-exe-only lists.

### F2. Metallib runtime paths remain build-dir absolute; Homebrew does not install them (MEDIUM)

- **Evidence:** Bridges load compile-time macros, e.g.:
  - `MetalRMSDBridge.mm`: `@METALRMSD_METALLIB_PATH` → `newLibraryWithURL`
  - `ShannonMetalBridge.mm`: `@SHANNON_METALLIB_PATH`
  - `CavityDetectMetalBridge.mm`: `@CAVITY_METALLIB_PATH` (fallback: fail → CPU)
  - `TurboQuantMetalBridge.mm`: `@TURBOQUANT_METALLIB_PATH`
  Paths are `"${CMAKE_CURRENT_BINARY_DIR}/….metallib"`.
- Formula install copies only named **executables** into `libexec/bin` (`FlexAIDdS`, `tENCoM`, `FlexAID`, …) and MC matrices / `*.def` — **no** `*.metallib`.
- `metal_eval.mm` is an exception: shaders are compiled from embedded source at runtime (no metallib file).
- `gpu_fast_optics_metal.mm` / `tencm_metal.mm` look next to the executable / main bundle for named metallibs that CMake in this PR **does not** even define as installable custom products under those names.
- **Why it matters for science/repro:** After a successful `brew install -s --HEAD --with-metal`, the binary may **link** and advertise Metal, while GPU paths for RMSD clustering, Shannon histograms, CavityDetect, TurboQuant silently fall back or fail when the Cellar build directory is gone. Operators can over-claim “Metal-accelerated Homebrew install” for cluster/entropy GPU paths. This packaging gap is **pre-existing** relative to the membership bug, but PR #260 markets Metal as the fixed brew option without installing shaders.
- **Fix recommendation:**
  1. `install(FILES ${…_METALLIB} DESTINATION …)` next to binaries **or** under `libexec/share/metal/`.
  2. Prefer runtime discovery: env `FLEXAIDDS_METALLIB_DIR`, then `dirname(argv0)`, then compile-time path.
  3. Formula: install metallibs beside `libexec/bin` (or set env in wrappers).
  4. Optional brew `test do` that touches Metal capability + one metallib load when `--with-metal`.

### F3. `flexaids_attach_metal_dispatch_runtime` duplicate guard is fragile (LOW–MEDIUM)

- **Evidence:**
  ```cmake
  get_target_property(_flexaids_metal_link_libs ${target_name} LINK_LIBRARIES)
  if(_flexaids_metal_link_libs MATCHES "flexaid_core")
      return()
  endif()
  ```
  Issues: (a) only inspects `LINK_LIBRARIES`, not `INTERFACE_LINK_LIBRARIES` / full link interface; (b) property may be `…-NOTFOUND`; (c) targets linking `FlexAID::Core` alias would **not** match `flexaid_core`; (d) does not detect targets that already have `metal_eval.mm` via other means.
- **Current callers** (`benchmark_dispatch`, `test_hardware_detect_dispatch`, `test_unified_dispatch`) do **not** link `flexaid_core`, so behavior is OK today.
- **Why it matters:** Future target that links `FlexAID::Core` and also calls the helper could **double-compile** `metal_eval.mm` → duplicate-symbol link errors (exactly the class of bug this PR is trying to eliminate).
- **Fix recommendation:** Prefer `if(TARGET …)` + `$<TARGET_PROPERTY:…,LINK_LIBRARIES>` generator-aware check, or document “never call this on core consumers” and assert with a CMake check for both `flexaid_core` and `FlexAID::Core`. Stronger: make helper a no-op if any metal bridge source is already on the target.

### F4. Brew formula disables the multi-consumer targets that motivated the fix (LOW)

- **Evidence:** Formula always passes `-DENABLE_CAVITY_DETECT_CLI=OFF -DENABLE_BENCHMARK_DATASETS=OFF …` even for `--HEAD --with-metal`.
- **Why it matters:** User-facing brew Metal validation only builds `FlexAID` / `FlexAIDdS` / tENCoM tools — the historical failure mode on `cavity_detect_cli` is **not** exercised in the formula path. In-repo full builds with those options ON are still the real regression test.
- **Fix recommendation:** Keep OFF for bottle size/disk by default, but add a maintainer CI job (or optional formula option) that builds `cavity_detect_cli` + `benchmark_datasets` with `FLEXAIDS_USE_METAL=ON` to lock the membership fix. Document that brew OFF is intentional, not a remaining link failure.

### F5. Temporary stable Metal refuse is correct at this tip; supersession is external (INFO)

- **Evidence:** Merge-tip formula `odie` when `--with-metal` without `--HEAD` because stable URL is still **v2.0.2**. Commit message / PR body acknowledge next release must carry the CMake fix.
- **Follow-up (not in this merge):** `600e95ba9` Release v2.0.3 removes the `odie` and points stable at `v2.0.3` tarball (sha later fixed in `00a7b6eae`).
- **Fix recommendation:** None for merge tip. Auditors reading only `1dba43f4b` must not claim “stable brew --with-metal works” until a post-260 tag is installed.

### F6. Doc drift: INSTALLATION still says head tracks `master` (LOW)

- **Evidence:** `docs/INSTALLATION.md` release bullet (merge tree ~L290): *“`head` tracks `master`”* while `Formula/flexaidds.rb` and the same file’s install commands use / imply **`main`**.
- **Why it matters:** Maintainer checklist conflicts with formula; reintroduces the exact branch foot-gun fixed in `99e5b35eb` / parent CI rename.
- **Fix recommendation:** One-line edit: `head tracks main`. Grep remaining `master` install URLs in docs.

### F7. PyPI sdist version check is better but parser is minimal (LOW)

- **Evidence:** Workflow Python snippet walks `pyproject.toml` lines for first `version =` under `[project]`. No TOML library; fails closed if missing/mismatch.
- **Why it matters:** Stops the historical `assert == "2.0.0"` drift. Can false-fail if version is multiline / annotated differently (unlikely with current PEP 621 layout).
- **Fix recommendation:** Optional: `tomllib` (Py≥3.11) for robust parse. Acceptable as-is for current packaging.

### F8. Verification claimed for this PR is incomplete vs AGENTS.md ideal (MEDIUM process)

- **Evidence (PR body):** Configure + build `FlexAIDdS`/`FlexAID` with `-DFLEXAIDS_USE_METAL=ON`; `nm` for Metal symbols. Explicitly **did not** re-run full `brew reinstall --with-metal` (disk ~14 GiB free). No `ctest` Metal matrix in this PR. Formula `test do` only checks paths + `--help`.
- **Why it matters:** AGENTS.md requires verify-with-execution before “done.” Link fix is highly plausible and architecture-correct; residual risk is (F2) runtime metallib + (F4) unbuilt brew targets.
- **Fix recommendation:** Add a self-hosted or macOS CI job: `FLEXAIDS_USE_METAL=ON` configure, build FlexAID + FlexAIDdS + cavity_detect_cli, link check, optional smoke that metallibs exist under build dir. Separate job for formula packaging of metallibs (F2).

### F9. Duplicate `FLEXAIDS_USE_METAL` PUBLIC defs (INFO)

- Set in both `LIB/CMakeLists.txt` and root Metal block. Harmless redundancy; slightly noisy compile lines.
- **Fix recommendation:** Optional single source of truth in root after bridges attach, or only in LIB with root adding metallib path macros.

### F10. Security / hygiene (INFO)

- No secrets, no machine-local absolute paths in formula/scripts, Apache-2.0 retained, tap remote is the public GitHub monorepo.
- Formula wrappers set `FLEXAIDDS_DATA_DIR` — pre-existing, appropriate for matrix lookup.

### F11. Science / ranking language (INFO)

- Docs carefully split **native docking tools** vs **Python analysis package**. No CF/ΔG overclaim introduced. GA budget / clustering / scoring untouched.

---

## 6. Ranking/scoring impact: **NO**

No changes to CF contact scoring, GA operators, BindingMode clustering, StatMech partition functions, PoseBusters gates, or claim CSV columns. Any performance difference appears only when Metal code paths successfully initialize at runtime (and only if metallibs/capabilities allow — see F2).

## 7. Reproducibility impact: **YES (mostly positive)**

| Axis | Impact |
|------|--------|
| Install commands (Homebrew 6 tap) | **Positive** — raw URL installs were already broken |
| Metal link of multi-target trees | **Positive** — source of truth for membership |
| Stable Metal brew at this exact tag | **Blocked by design** until next release (F5) |
| Runtime GPU shader repro after brew install | **Still weak** (F2) |
| HEAD branch pin | **Positive** (`main` not feature branch / not stale `master`) |

## 8. Tests adequate: **PARTIAL**

| Layer | At merge |
|-------|----------|
| In-repo CMake Metal link of FlexAID/FlexAIDdS | Claimed in PR; not re-executed in this audit session |
| `ctest` / Metal-specific tests | Not part of PR verification |
| Homebrew formula `test do` | Path + `--help` only |
| Full brew `--with-metal` | Explicitly skipped (disk) |
| metallib install / load after install | **Missing** |
| Regression: cavity + benchmark with Metal | Formula disables them; needs CI |

---

## 9. Follow-ups already in history (context for auditors)

| Later commit | Relation to this merge |
|--------------|------------------------|
| `600e95ba9` Release v2.0.3 | Ships membership fix as stable; removes temporary `odie` |
| `00a7b6eae` Formula sha256 for v2.0.3 | Completes bottle/tarball pin |
| `d280a8fc5` Pin Homebrew HEAD to main only | Reinforces `99e5b35eb` head policy |
| Metal CI gate PRs (#261 / later) | Separate hardening of macOS ctest / release language |

Do **not** credit those to `1dba43f4b` itself; they close F5 packaging versioning, not F2 metallib install.

---

## 10. Verdict: **ACCEPT_WITH_FOLLOWUPS** (merge hygiene: **CLEAN**)

**Keep this merge.** The CMake membership fix is the right repair for Homebrew `--with-metal` undefined Metal symbols and removes a large class of per-executable duplication. Formula/docs correctly bridge the gap until a post-v2.0.2 tag.

**Before treating “Homebrew Metal” as production-complete**, address at least:

1. **F2** — install + discover metallibs for installed binaries (or document CPU fallback as the supported brew Metal subset: e.g. `metal_eval` only).
2. **F6** — fix residual `head tracks master` in INSTALLATION.
3. **F3/F8** — harden attach helper; add CI that builds Metal + at least one non-FlexAID core consumer (`cavity_detect_cli`) with `FLEXAIDS_USE_METAL=ON`.

No ranking/scoring/thermo regressions introduced. No source edits performed for this audit.
)
