# PoseBust — science and code audit (2026-08-18)

**Scope:** `LIB/PoseBust/` as shipped on this branch, plus DatasetRunner claim wiring and the scoring-time `pb_clash` term (not in `LIB/PoseBust/`, but named so chemists will confuse it).
**Mode:** Diagnostic. No engine, ranking, or thermodynamic behaviour was changed.
**Chemist-facing introduction:** [`LIB/PoseBust/README.md`](../../LIB/PoseBust/README.md)
**Auditor:** Cursor Grok 4.6 cloud agent; every factual claim below was read from source in this session.

Literature reference for the *official* gate (not incorporated as source):
Buttenschoen, Morris, Deane, *Chem. Sci.* **15**, 3130–3139 (2024), doi:10.1039/D3SC04185A.

---

## 1. Executive verdict

PoseBust is a **two-backend post-election physical-validity stack**. The claim-grade backend is an argv bridge to PoseBusters 0.6.5 `bust`. The in-tree NativePoseQC is a serious, fail-closed C++26 diagnostic that **must not be sold as PoseBusters**.

| Question | Answer on this tip |
|----------|--------------------|
| Can you report “PoseBusters pass” from NativePoseQC? | **No.** Shared *key names*, different physics. |
| Is official `pb_pass` the 27-boolean 0.6.5 redock set, water included? | **Yes.** Schema pin in `BustCli.cpp` is set-equality, not substring. |
| Does `success_pb` require RMSD ≤ 2 Å on the same pose? | **Yes**, via `ElectedPoseBustOutcome::finalize_success_pb`. |
| Is `claim_ready` reachable without `bust_cli`? | **No.** DatasetRunner requires `pb_backend == "bust_cli"` plus tENCoM/Eigen on the pose SHA-256. |
| Does missing `bust` silently become a native pass that still looks like `pb_pass`? | **Yes, unless** `FLEXAIDDS_POSEBUSTERS_REQUIRE_CLI` is set. Fallback backend string is `native_pose_qc_fallback`. A banner is printed. `claim_ready` still cannot fire. Intermediate `success_pb` **can**. |
| Does NativePoseQC test waters / HEM / Zn? | **No.** Protein loader keeps standard-AA ATOM heavy atoms only. Cofactor/water keys are **vacuous True**. |
| Does PoseBust rank poses or compute ΔS? | **No.** |
| Is `CF.pb_clash` the same object? | **No.** Optional GA penalty, default weight 0, different vdW table. |

**One-line verdict:** BustCli is a careful, pin-locked referee. NativePoseQC is a useful intramolecular/protein-clash smoke test that will green-light poses official PoseBusters would reject (waters, cofactors, RDKit sanitization, UFF energy, CIP stereo) and that can fail-open on InChI.

---

## 2. Architecture (as implemented)

```
elected_pose.pdb ──► load_pdb_flexaid_ligand (CONECT → REMARK → HETATM fallback)
                 ──► assign_topology_from_reference (crystal SDF, fail-closed)
                 ──► write_sdf (V2000)                ──► run_upstream_bust  ──► pb_pass (27 booleans)
                 ──► evaluate_paths NativePoseQC      ──► native_qc_* (diagnostic)
                 ──► SHA-256 elected bytes == validator input
                 ──► finalize_success_pb(success_rmsd)
```

Locked API: `Engine.h` (`evaluate`, `evaluate_paths`, `validate_elected_pose`, `resolve_backend_from_env`).
Default backend: `Backend::BustCli` (`Engine.cpp` `resolve_backend_from_env`).

DatasetRunner (`LIB/DatasetRunner.cpp` ~7624–7770) always calls `validate_elected_pose` after a completed dock, with `force_native_when_off=true` and `native_fallback_if_bust_missing=true`.

---

## 3. What NativePoseQC gets right

These are load-bearing and tested (`tests/test_posebust.cpp`).

1. **Ligand identity on FlexAID complexes.** `load_pdb_flexaid_ligand` prefers CONECT serials, skips standard AA and the cofactor blacklist (HEM, NAD, …). `Extract1G9VNotHEM` asserts 25 heavies, not heme-scale.
2. **Topology is not guessed for the claim path.** `evaluate_paths` refuses missing crystal SDF. `assign_topology_from_reference` requires identical atom counts (including H/Du) and element sequence; permuted order fails closed.
3. **Shared PDB XYZ decoder.** `PdbCoords.h` handles FlexAID compact negatives (`-80.275-146.614`) and rejects NaN/junk. Same parser as DatasetRunner RMSD.
4. **Fail-closed election.** Empty/missing elected path → `pb_ran=false`, `pb_pass=false`. `finalize_success_pb` cannot invent a pass.
5. **Provenance.** Elected file is hashed before and after; mismatch appends `validator_input_provenance` and clears `pb_pass`. Bust receipts record binary path, SHA-256, argv, raw CSV.
6. **Schema pin (BustCli).** Canonical 27 columns in emission order. Duplicate headers, column-count mismatch, missing columns, extra scored columns (including extra `rmsd_*` names) fail closed. Blank/NaN/non-boolean values count as fails. Metadata exemption is an exact-name list of four columns, not `find("rmsd")`.
7. **Water stays in the official gate.** Comments in `BustCli.cpp` document the 2026-07-31 pin bump that removed the phantom `no_protein_clashes` column (upstream 0.6.5 never emitted it, so every real run schema-failed and fell back to native).
8. **JSON report labels native success as diagnostic.** `native_qc_diagnostic_pass` / `success_pb_campaign` are explicitly not DatasetRunner `success_pb`.

---

## 4. Science gaps (NativePoseQC vs published PoseBusters)

Anchored to Buttenschoen et al. 2024 Table 4 / §2.2 and to `Checks*.cpp`.

### 4.1 Intermolecular physics

| Quantity | Paper / `bust` 0.6.5 | NativePoseQC |
|----------|----------------------|--------------|
| Min lig–protein distance | \(d/(r_i+r_j)\) vs 0.75; RDKit vdW; covalent radii for inorganic cofactors | Absolute 1.5 Å **and** 0.75 × **Bondi** (`ChecksProtein.cpp`) |
| Volume overlap | RDKit ShapeTverskyIndex; vdW **×0.8**; threshold 7.5% of ligand volume | 0.5 Å voxels of **unscaled** Bondi spheres; ≤7.5% (`kMaxVolumeOverlap=0.075`) |
| Waters / organic / inorganic cofactors | Classified from the condition molecule; both distance and volume tests selected in `redock.yml` | Six keys emitted as vacuous pass (`ChecksProtein.cpp` ~429–452) |
| Protein atoms seen | Whatever PDB/SDF was passed to `-p` | `load_pdb_protein_heavy`: ATOM + standard AA only; **no HETATM metals, no HOH, no HEM** |
| Spatial support | Full condition molecule | Crop to 10 Å of ligand heavy COM (`Engine.cpp`) |

Consequence: NativePoseQC cannot fail a water clash or a Zn coordination-distance check. On holo receptors that still contain crystal waters, official `bust` will fail many poses NativePoseQC reports green. That matches the campaign note in `workorders/PHASE4_NEAR_MISS_NULL_STACK.md` (water dominates observed official failures).

### 4.2 Intramolecular physics

| Quantity | Paper / `bust` | NativePoseQC |
|----------|----------------|--------------|
| Bond lengths / angles | RDKit DistanceGeometry bounds, 25% | Cordero radii ±25%; angles from a hybridization heuristic + 25° absolute floor |
| Internal clash | DG 1–5+ , 30% | Graph distance ≥4, 0.70 × Bondi |
| Internal energy | UFF(pose) / mean UFF(50 ETKDGv3 conformers) ≤ 100 | Mean squared relative bond strain vs covalent radii |
| Aromatic flatness | RDKit aromatic rings | 5–6 cycles, majority degree-3 vote. **Dead zone:** MDL aromatic rings that fail the degree-3 majority (typical heavy-only phenyl) are skipped here **and** skipped by the aliphatic test (order-4 bond). 1G9V crystal self-dock this session: `aromatic_ring_flatness n_checked=0` yet PASS. |
| Sanitization | RDKit SanitizeMol | Finite XYZ + known Z + bond index/order ∈ {1,2,3,4} |
| Radicals | RDKit radical electrons | Over-valence only |
| Stereo | RDKit CIP / E-Z | Geometric signs; **vacuous True without crystal** |
| InChI | RDKit | `inchi-1` or **soft pass if missing** (`ChecksChemistry.cpp` ~473–479) |

`Suite` in `EvaluateOptions` is dead: `evaluate()` never reads `opt.suite`. Dock vs redock vs mol is approximated by “is the crystal pointer null?”.

### 4.3 Split-brain van der Waals tables

Three tables exist in one repository:

| Location | N | O | Cl | P | I | Claimed role |
|----------|---|---|----|---|---|--------------|
| `LIB/soft_wall.h` `posebusters_vdw_radius` | 1.60 | 1.55 | 1.80 | 1.95 | 2.10 | “RDKit periodic-table vdW used by PoseBusters” for **GA `pb_clash`** |
| `LIB/PoseBust/ChecksGeometry.cpp` `vdw_radius` | 1.55 | 1.52 | 1.75 | 1.80 | 1.98 | NativePoseQC clash + volume |
| `ChecksProtein.h` comment | default 1.70 if unknown | | | | | **Wrong vs implementation default 2.00** in `ChecksGeometry.cpp` |

A chemist comparing Native `min_lig_prot_dist` to a PoseBusters waterfall, or to a `CF.pb_clash` debug dump, is not comparing one physical model.

`CF.pb_clash` is also **off by default** (`FLEXAIDDS_PB_CLASH_WEIGHT=0`). Root README architecture art still draws it as if it were in the live scoring loop.

---

## 5. Code-quality findings (non-ranking)

Severity is “correctness of the *validator*,” not “moves docked coordinates.”

| ID | Finding | Evidence | Severity |
|----|---------|----------|----------|
| C1 | `inchi_convertible` fail-open when `inchi-1` absent | `ChecksChemistry.cpp` 473–479; crystal self-dock test expects pass even without InChI | High for native-as-chemistry |
| C2 | Vacuous cofactor/water keys included in `success_pb_full()` | `Types.h` `all_passed()` ANDs every check; `validate_elected_pose` maps native pass from `success_pb_full()` | High if `pb_backend!=bust_cli` is reported as PB |
| C3 | `popen("command -v inchi-1")` | `ChecksChemistry.cpp` 450–461. Constant string, but a shell. Rest of BustCli uses argv `shell_exec`. | Medium hygiene |
| C4 | Homebrew absolute paths as InChI candidates | `/opt/homebrew/bin/inchi-1`, `/usr/local/bin/inchi-1` | Low (not user-home; still host-specific) |
| C5 | `Suite` unused | `Engine.cpp` `evaluate()` ignores `opt.suite` | Low / API lie |
| C6 | Cell-list path reports `min_dist = 5 Å` when no neighbour pair is found | `ChecksProtein.cpp` 224–228 | Low (no clash, metric is a sentinel) |
| C7 | `write_sdf` V2000, 3-digit counts | `Loaders.cpp` 509–512 | Low for drug-like ligands |
| C8 | CONECT bonds always order 1 until crystal topology is applied | `Loaders.cpp` 691 | Mitigated: `evaluate_paths` requires topology assign |
| C9 | `evaluate()` intermolecular keys omitted when protein is empty (not cropped-empty) | `Engine.cpp` 268–288 vs empty protein | Dock-without-protein looks like ligand-only; `mol_cond_loaded` fails |
| C10 | Sibling GitHub README is a 353-byte placeholder claiming NEON/Metal/CUDA and “ΔS modeling” | `gh api` README decode this session | Docs hazard on the standalone repo; **this tree’s code has none of that** |
| C11 | `test_posebust_upstream_parity.py` compiles a throwaway `clang++ -std=c++26` tool | Skips/fails on GCC-13 clouds without clang++ 18 | Test portability |
| C12 | `catch (...) {}` around bust receipt write | `Engine.cpp` 692–693 | Receipt can vanish without failing `pb_pass` |
| C13 | Aromatic-ring dead zone on heavy-only graphs | `ChecksGeometry.cpp` `majority_aromatic_or_sp2` (degree==3) vs aliphatic skip on order==4. 1G9V crystal smoke this session: both flatness checks `n_checked=0`, still PASS | High for native-as-chemistry |

None of C1–C13 change GA ranking. C1, C2, and C13 change what a native `pb_pass` *means*.

---

## 5b. Runtime smoke (this session, GCC 13, C++23)

Compiled `LIB/PoseBust/*.cpp` and ran `evaluate(crystal_sdf, apo_pdb, &crystal)` on Astex 1G9V.

| Observation | Value |
|-------------|--------|
| Ligand atoms / bonds | 25 / 26 (all heavy; formula C20NO4) |
| Protein ATOM heavies loaded | 4384 |
| After 10 Å COM crop | 135 |
| Native checks | 27 pass / 0 fail |
| `inchi_convertible` | `soft=true inchi-1_missing` |
| Aromatic / aliphatic ring checks | **0 rings scored** |
| Min lig–protein distance | 3.22 Å, 3375 pairs, 0 relative vdW clashes |
| Volume overlap fraction | 0.0 (0.5 Å voxels, unscaled Bondi) |
| Six cofactor/water keys | vacuous pass |
| `resolve_backend_from_env()` | `BustCli` (0) |

Official `bust` was **not** installed in this environment, so native-vs-upstream boolean parity was not re-run here. The GoogleTest `PoseBustParity.CrystalSelfDockAgreesWithUpstreamBust` is the pin for that comparison when `FLEXAIDDS_POSEBUSTERS_BIN` is set.

---

## 6. Claim-contract wiring (DatasetRunner)

From `LIB/DatasetRunner.cpp`:

- `success_pb := success_rmsd && pb_pass` after `finalize_success_pb`.
- `claim_ready` additionally requires `pb_backend == "bust_cli"`, matching RMSD/PB/tENCoM pose SHA-256, non-empty `posebusters_input_sha256`, `tencom_status==ok`, `eigen_status==ok`, protocol eligibility, and `|score_pose_delta| ≤ 1e-4`.
- Native suite is always logged as `[NATIVE-POSE-QC]` “parity diagnostic.”
- Official result is logged as `[POSEBUSTERS] backend=...`.

So: **STRICT headline cannot be faked by NativePoseQC.** Intermediate column `success_pb` **can** be native-backed. Aggregators that headline `success_pb` without filtering `pb_backend` will over-call “PoseBusters.” `scripts/aggregate_claim_metrics.py` is specified to use `claim_ready` as headline (`admission_metrics_contract.md` §1).

---

## 7. Tests that pin the contract

`tests/test_posebust.cpp` (target `PoseBustTests`):

- PDB compact-negative XYZ
- BustCli schema: 27 True → pass; extra scored column → pin bump; extra `rmsd_*` → pin bump; chemistry False ≠ schema error
- Loaders: Cl recovery from misaligned PDB, Du hydrogen, topology permute fail, 1G9V 25 heavy, not HEM
- Engine: required upstream key names emitted; crystal self-dock core checks
- Optional parity: native vs `bust` on 1G9V crystal self-dock, all 27 keys
- DatasetRunner mapping: native `pb_pass` := `success_pb_full()`
- Default backend is `BustCli`
- Elected-path fail-closed; Off → Native floor; `success_pb` algebra

`tests/p0_claim_contract/test_posebusters_fixtures.py`: real PoseBusters Python API discriminates a clean MMFF pose from a coincident-atom clash (skips if rdkit/posebusters absent).

`tests/test_soft_wall.cpp`: search-time `posebusters_vdw_radius` table (not NativePoseQC).

---

## 8. Licensing

- NativePoseQC + BustCli: Apache-2.0, first-party, clean-room (`docs/licensing/clean-room-policy.md` §3.1).
- `bust` / RDKit: BSD, subprocess only, not vendored.
- Check *names* match published PoseBusters columns for report parity; algorithms are independent.
- No GPL.

---

## 9. What this audit did not do

- Did not run a live Astex-85 campaign or quote a success rate.
- Did not byte-compare NativePoseQC vs every PoseBusters 0.6.5 CSV in an OUT tree (that evidence lives in prior workorders; this session re-read the *code* that would produce those CSVs).
- Did run a native `evaluate()` smoke on Astex 1G9V crystal vs apo (27/27 native pass, including vacuous water/cofactor and zero ring-flatness checks).
- Did not modify `LIB/` sources. Documentation only.

---

## 10. Recommended follow-ups (not done here)

1. Treat native `inchi_convertible` missing-binary as **fail** (or a distinct `soft` flag that cannot enter `all_passed()`).
2. Stop counting vacuous cofactor/water keys inside native `success_pb_full()`; keep emitting them as `n_checked=0` diagnostics.
3. Score aromatic rings from MDL order-4 bonds (or RDKit-like SSSR aromaticity), not a heavy-only degree-3 majority — otherwise phenyls silently skip both flatness tests.
4. Align NativePoseQC vdW with `soft_wall.h` *or* document three tables in one chemist-visible place (now done in `LIB/PoseBust/README.md`).
5. Either implement `opt.suite` or delete the enum.
6. Default claim campaigns: `FLEXAIDDS_POSEBUSTERS_REQUIRE_CLI=1` in DatasetRunner claim protocols so `success_pb` cannot be native-backed.
7. Replace the standalone PoseBust GitHub README placeholder (hardware/ΔS fiction) with this chemist text.
