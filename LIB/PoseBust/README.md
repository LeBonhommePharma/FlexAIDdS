# PoseBust — pose physical-validity for computational chemists

**In-tree module:** `LIB/PoseBust/` (namespace `flexaids::posebust`)
**License:** Apache-2.0 (clean-room). Official PoseBusters `bust` is an optional BSD subprocess, not vendored.
**Sibling package:** [LeBonhommePharma/PoseBust](https://github.com/LeBonhommePharma/PoseBust) (same lineage, no docking-engine dependency)

This note is written for people who dock ligands, not for people who want another C++ API dump. It says what PoseBust *is*, what a “pass” *means*, and which sentences you must not write in a paper, SI, or internal claim table.

---

## Why this exists

A pose with crystallographic RMSD ≤ 2 Å can still be chemically garbage: stretched bonds, inverted stereo, a phenyl ring folded like a boat, or a ligand occupying the same van der Waals volume as a protein side chain. Buttenschoen, Morris and Deane showed that this is not a theoretical quibble. On their PoseBusters Benchmark, several deep-learning dockers that looked competitive on RMSD produced structures that fail ordinary physical-plausibility tests. Their prescription is now the field’s minimum bar:

> Report a pose as successful only if it is **near-native (RMSD ≤ 2 Å) and physically valid (PoseBusters)**.

That is the S2 / `success_pb` contract in FlexAIDdS. RMSD-only is S1, a diagnostic. It is not docking success.

**Citation (the method we gate against, not a dependency we copy):**
M. Buttenschoen, G. M. Morris, C. M. Deane, *Chem. Sci.* **15**, 3130–3139 (2024). [doi:10.1039/D3SC04185A](https://doi.org/10.1039/D3SC04185A)

---

## Two layers. Do not conflate them.

PoseBust is a **post-election validator**. It does not search, does not rank BindingModes, and does not compute a free energy. The genetic algorithm still optimizes the Voronoi **contact-function (CF) scoring proxy**. After a pose is elected, PoseBust asks whether that pose is a chemically and sterically admissible molecule in the pocket.

| Layer | What it actually is | What you may call a pass |
|-------|---------------------|---------------------------|
| **BustCli** (`BustCli.cpp`) | Argv bridge to the installed PoseBusters 0.6.5 CLI (`bust`). RDKit chemistry, UFF internal energy, RDKit distance-geometry bounds, ShapeTversky volume overlap. BSD tool, user-installed. | **`pb_pass` for claims.** `pb_backend` must be `bust_cli`. |
| **NativePoseQC** (`Checks*.cpp`, `Engine.cpp`) | Apache-2.0 clean-room C++26 suite that *reuses PoseBusters column names* so reports can be compared. Algorithms are original. No RDKit, no posebusters source. | **Diagnostic / fallback only.** Never “PoseBusters passed.” |

Default backend is official `bust`. NativePoseQC always still runs as a parity diagnostic. If `bust` is missing, NativePoseQC still runs (`native_qc_*`) and `pb_backend` becomes `native_pose_qc_fallback`, but **`pb_ran` and `pb_pass` stay false**. `success_pb` cannot be minted from the in-house suite. `claim_ready` stays unreachable because DatasetRunner requires `pb_backend == "bust_cli"`.

```
success_rmsd  =  elected_pose RMSD ≤ 2.0 Å          (S1)
pb_pass       =  all 27 PoseBusters 0.6.5 redock booleans True
success_pb    =  success_rmsd  ∧  pb_pass            (S2)
claim_ready   =  success_pb
              ∧  pb_backend == bust_cli
              ∧  tENCoM/Eigen on the same pose SHA-256
              ∧  protocol eligibility + score–pose consistency
```

RMSD is **not** one of the 27 booleans. The CSV column `rmsd_≤_2å` is recorded and then ignored for `pb_pass`. Mixing the two would double-count geometry and let a 3 Å pose “pass PoseBusters” while failing the docking criterion, or vice versa.

Normative aggregation: [`benchmarks/protocols/admission_metrics_contract.md`](../../benchmarks/protocols/admission_metrics_contract.md).

---

## What “physically valid” means (the 27 redock checks)

PoseBusters 0.6.5 `redock.yml` emits 31 CSV columns: 4 metadata (`file`, `molecule`, `position`, `rmsd_≤_2å`) and **exactly 27 scored booleans**. BustCli pins that set by **set equality**. A future PoseBusters that adds or drops a scored column fails closed and demands a deliberate pin bump. Gate membership is never decided by substring matching.

Water is **in**. Upstream redock.yml selects `minimum_distance_to_waters` and `volume_overlap_with_waters`. Dropping them would raise the apparent pass rate (water often dominates failures) and would no longer be “PoseBusters pass.” Report a non-water diagnostic beside `pb_pass` if you need it; do not mutate `pb_pass`.

| # | Check key | Official PoseBusters 0.6.5 (RDKit `bust`) | NativePoseQC (this tree) |
|---|-----------|-------------------------------------------|---------------------------|
| 1–3 | `mol_pred_loaded`, `mol_true_loaded`, `mol_cond_loaded` | RDKit load of predicted ligand, crystal ligand, protein | Non-empty atom arrays |
| 4 | `sanitization` | RDKit `SanitizeMol` (valence, aromaticity, kekulize) | Finite coordinates, known Z, valid bond indices |
| 5 | `inchi_convertible` | RDKit → InChI | Real `inchi-1` (`FLEXAIDDS_INCHI_BIN` / PATH). **Missing binary → skipped** (`passed=false`, ignored by `all_passed()`) |
| 6 | `all_atoms_connected` | Single RDKit fragment | Single connected component on the **heavy-atom** bond graph |
| 7 | `no_radicals` | RDKit radical-electron count | Over-valence only (`bos > max valence + 1.05`). Under-valence allowed (heavy-only FlexAID poses omit H) |
| 8 | `molecular_formula` | Identity vs crystal | Heavy-atom element **multiset** equality |
| 9 | `molecular_bonds` | Bond-table identity vs crystal | Bond **count** within 20% of crystal — not graph isomorphism |
| 10 | `double_bond_stereochemistry` | RDKit E/Z vs crystal | Geometric torsion sign vs crystal; **vacuous True if no crystal or atom-count mismatch** |
| 11 | `tetrahedral_chirality` | RDKit CIP vs crystal | Signed tetrahedral volume vs crystal, neighbor-index order, **not CIP**; vacuous True without crystal |
| 12 | `bond_lengths` | Within 25% of RDKit distance-geometry bounds | Within 25% of Cordero covalent-radius sum |
| 13 | `bond_angles` | Same DG 25% window | Hybridization heuristic (109.5° / 120° / 180°) with 25% relative **or** 25° absolute floor |
| 14 | `internal_steric_clash` | DG 1–5+ distances, 30% tolerance | Heavy pairs with graph distance ≥ 4: \(d \ge 0.70\,(r_i+r_j)\) using the `soft_wall.h` RDKit-like vdW table |
| 15 | `aromatic_ring_flatness` | RDKit aromatic rings, out-of-plane | Size 5–6 cycles; aromatic if **MDL order-4 bond present or** majority degree-3; max OOP ≤ 0.25 Å |
| 16 | `non-aromatic_ring_non-flatness` | Aliphatic rings must pucker | Skips cycles that have an MDL order-4 bond (those go to row 15) |
| 17 | `double_bond_flatness` | DG / RDKit | \(\lvert\sin\phi\rvert \le 0.25\) on first substituents |
| 18 | `internal_energy` | UFF energy ≤ 100 × mean of 50 ETKDGv3+UFF conformers | Mean **squared relative bond strain** vs covalent radii; pass if ≤ 0.0625, or ≤ 2× crystal |
| 19 | `protein-ligand_maximum_distance` | Ligand is in the pocket | Some protein heavy atom within 5 Å of some ligand heavy atom |
| 20 | `minimum_distance_to_protein` | \(d/(r_i+r_j) \ge 0.75\) (RDKit vdW; covalent radii for inorganic cofactors) | \(d \ge 1.5\) Å **and** \(d \ge 0.75\,(r_i+r_j)\) using the same RDKit-like table as `CF.pb_clash` |
| 21–23 | min distance to organic / inorganic cofactors / waters | Real classification from the condition molecule | **Skipped** (`n_checked=0`) on apo protein crop — keys are emitted but do not inflate `all_passed()` |
| 24 | `volume_overlap_with_protein` | RDKit `ShapeTverskyIndex`; vdW scaled **0.8**; overlap **< 7.5%** of ligand volume | 0.5 Å voxel occupancy of **unscaled** spheres; overlap **≤ 7.5%** |
| 25–27 | volume overlap with cofactors / waters | Same Tversky, scale 0.8 organic / 0.5 inorganic | **Skipped** (same as 21–23) |

NativePoseQC also crops the protein to heavy atoms within **10 Å of the ligand heavy-atom centre of mass** before intermolecular checks. Official `bust` sees the condition molecule the CLI was given (typically the full receptor PDB).

---

## How FlexAIDdS uses this (and what it does *not* do)

1. DatasetRunner elects a BindingMode pose (`elected_pose.pdb`).
2. Ligand atoms are taken from FlexAID **CONECT** (serials typically 90001+), not “all HETATM”. That is the 1G9V HEM-contamination fix: 25 ligand heavies, not 128 heme atoms.
3. Bond orders come from the **crystal SDF** (`assign_topology_from_reference`). Element sequence and full atom counts including explicit H/Du must match. Positional remapping of repeated elements is refused.
4. NativePoseQC runs. Official `bust` runs on the extracted ligand SDF vs the receptor PDB vs the crystal SDF.
5. `success_pb` is AND-ed with RMSD. `claim_ready` further requires `bust_cli` + tENCoM/Eigen on the **same** pose SHA-256.

PoseBust never feeds back into CF ranking. A physically invalid pose can still win the GA.

### Do not confuse this with `CF.pb_clash`

`FLEXAIDDS_PB_CLASH_WEIGHT` (default **0.0**) is an optional **search-time** intermolecular clash *penalty* inside `vcfunction.cpp`. It uses `posebusters_vdw_radius()` in `LIB/soft_wall.h` (RDKit-like table: N 1.60, O 1.55, Cl 1.80 Å). NativePoseQC intermolecular checks now use **the same table** via `vdw_radius()` in `ChecksGeometry.cpp`. Same *idea* (0.75 × summed vdW), **still not a PoseBusters pass** — official `bust` uses RDKit radii plus covalent radii for inorganic cofactors.

Turning the clash penalty on during search does not make `pb_pass` true. Leaving it off does not make `pb_pass` false.

---

## What you must not write

| Forbidden sentence | Why |
|--------------------|-----|
| “PoseBusters passed” after only NativePoseQC | Column names are shared; the physics is not. |
| “Docking success = 72%” from RMSD ≤ 2 Å alone | S1 is not S2. Modern claims need `success_pb`. |
| “We use PoseBusters, so waters don’t matter” | Water checks are in the pinned 27. |
| “NativePoseQC is equivalent to PoseBusters 0.6.5” | Stereo, UFF energy, sanitization, Tversky volume, cofactors, and waters are not equivalent. |
| “PoseBust computes ΔS / binding free energy” | It does not. No partition function, no tENCoM, no GPU entropy dispatch lives in this directory. |
| “`CF.pb_clash` is the PoseBusters gate” | Search penalty ≠ post-election validator. |

If `pb_backend` is `native_pose_qc`, `native_pose_qc_fallback`, `bust_cli_missing`, `skipped_*`, or `error`, you do not have a claim-ready PoseBusters result.

---

## How to run it

Install official PoseBusters (needed for any sentence that says “PoseBusters”):

```bash
python3 -m venv .venv-posebusters
.venv-posebusters/bin/pip install 'posebusters==0.6.5'
export FLEXAIDDS_POSEBUSTERS_BIN="$PWD/.venv-posebusters/bin/bust"
# Claim campaigns: refuse silent native fallback
export FLEXAIDDS_POSEBUSTERS_REQUIRE_CLI=1
```

Environment:

| Variable | Effect |
|----------|--------|
| `FLEXAIDDS_POSEBUST=0` | Skip (DatasetRunner still forces NativePoseQC as a mandatory floor unless you also disable that) |
| `FLEXAIDDS_POSEBUST_BACKEND=bust` | Official CLI (default) |
| `FLEXAIDDS_POSEBUST_BACKEND=native` | NativePoseQC only — diagnostic |
| `FLEXAIDDS_POSEBUSTERS_BIN` | Absolute path to `bust` |
| `FLEXAIDDS_POSEBUSTERS_REQUIRE_CLI=1` | Missing `bust` → fail closed, no native fallback |
| `FLEXAIDDS_INCHI_BIN` | Optional IUPAC `inchi-1` for NativePoseQC `inchi_convertible` |

C++ entry points (GoogleTest coverage in `tests/test_posebust.cpp`):

```cpp
// Filesystem path used by DatasetRunner
auto report = flexaids::posebust::evaluate_paths(
    complex_pdb, receptor_pdb, crystal_sdf, opt);

// Elected BindingMode pose — the claim gate
auto out = flexaids::posebust::validate_elected_pose(
    elected_pose_pdb, receptor_pdb, crystal_sdf, opt);
out.finalize_success_pb(success_rmsd);  // success_pb := rmsd ∧ pb_pass
```

Rebuild tests after touching this directory:

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_posebust -j "$(nproc)"
ctest --test-dir build -R PoseBustTests --output-on-failure
```

Parity vs installed `bust` (skips if the CLI is absent):

```bash
python3 tests/test_posebust_upstream_parity.py
```

---

## Code map

| File | Role |
|------|------|
| `Types.h` | Molecule graph, check keys, `success_pb` algebra (native diagnostic ≠ claim gate) |
| `PdbCoords.h` | Shared finite PDB XYZ decoder (same bytes as DatasetRunner RMSD) |
| `Loaders.cpp` | SDF V2000, FlexAID CONECT ligand extract, fail-closed topology transfer |
| `ChecksChemistry.cpp` | Load / sanitization / InChI / connectivity / formula / stereo / strain |
| `ChecksGeometry.cpp` | Bonds, angles, internal clash, ring and double-bond flatness |
| `ChecksProtein.cpp` | Ligand–protein distance, voxel volume overlap, skipped cofactor/water keys |
| `Engine.cpp` | Native orchestration, JSON sidecar, `validate_elected_pose` |
| `BustCli.cpp` | Official `bust` argv + 27-column schema pin + receipts |

---

## Limitations chemists should budget for

These remain after the 2026-08-18 honesty fixes.

1. **Native cofactor/water keys are skipped, not scored**, on apo crops. A pose that hits a crystallographic water or a heme iron can still be NativePoseQC-green (`all_passed()` ignores the skipped rows) and `bust`-red. Missing `bust` no longer copies that native green onto `pb_pass`.
2. **Volume overlap is a different estimator.** Voxels of unscaled spheres vs RDKit Tversky with 0.8-scaled vdW. The 7.5% threshold is the same *number* with different *volumes*.
3. **Official `bust` vdW/energy/stereo are still not NativePoseQC.** Aligning NativePoseQC with `soft_wall.h` does not make it PoseBusters 0.6.5.
4. **`Suite::Mol` is ligand-only.** Dock/Redock still emit the 27-key dock list when a crystal pointer is set (parity lock).
5. **SDF writer is V2000** (≤999 atoms). Fine for drug-like ligands; not a biological assembly dumper.
6. **No GPU, no ΔS, no Eigen** in this directory. Hardware-dispatch or “entropy integration” language does not describe this code.

The 2026-08-18 source audit that produced this introduction is
[`docs/audit/2026-08-18_posebust_science_and_code_audit.md`](../../docs/audit/2026-08-18_posebust_science_and_code_audit.md).
