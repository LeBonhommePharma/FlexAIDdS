# FlexAIDdS Benchmark Standard — Three Reproducible Tiers

**Status:** Normative. This document defines the *only* sanctioned ways to report a
FlexAIDdS docking-power number. Any benchmark result quoted in a paper, thesis,
issue, or commit message **must** name the tier (`TIER-1`, `TIER-2`, or `TIER-3`)
it was produced under. A bare percentage without a tier label is not a result.

**Audience:** A thesis reviewer (or independent group) who needs to know exactly
what protocol produced a number and which published method that number is
comparable to — with no methodological footnotes required.

---

## 0. Why this document exists

FlexAIDdS version numbers (v1…v103) historically mixed **three incompatible
docking protocols** under one "Astex Diverse 85" banner. The same dataset name
was used for an easy protocol that yields 94–96 % and a hard protocol that yields
22–26 %. Comparing across versions, or against the literature, was therefore
meaningless unless you reconstructed the exact env vars and input JSON of each run.

The three protocols differ along **three axes only**:

| Axis | TIER-1 | TIER-2 | TIER-3 |
|------|--------|--------|--------|
| Receptor | cognate, ligand removed (`<PDB>_apo.pdb`, same PDB ID as ligand) | **non-cognate** (different PDB's structure) | **non-cognate** |
| Binding site | oracle crystal-centroid sphere | oracle crystal-centroid sphere (from the **target** complex) | **SURFNET** (no crystal info) |
| Native pose seeding | **forbidden** | forbidden | forbidden |

Everything else — the GA engine, scoring (VCT), clustering, restart count,
softcore wall, RMSD metric, success threshold — is **held identical across all
three tiers**. The only knobs that change between tiers are the receptor source,
the site source, and (always off) native seeding. This is what makes the three
numbers a clean difficulty ladder rather than three unrelated experiments.

> **Native seeding is banned in all three tiers.** Past headline numbers
> (e.g. self-dock 94–96 %) were produced with `FLEXAIDDS_NATIVE_SEED_FRAC=0.90`,
> which floods the GA population with the crystal pose and yields **seed-echo**
> poses (RMSD ≈ 0.0 Å that are copies of the input, not search results). Roughly
> half of those "successes" are echoes. No literature physics method seeds the
> answer, so seeding makes the number incomparable. The standard fixes
> `FLEXAIDDS_NATIVE_SEED_FRAC=0.0` everywhere.

---

## 1. Common invariants (apply to every tier)

| Property | Value |
|----------|-------|
| Dataset | Astex Diverse 85 (Hartshorn et al. 2007, *J. Med. Chem.* 50, 726) |
| N targets | 85 |
| Success criterion | **RMSD < 2.0 Å** (`rmsd_hungarian` field), the universal Astex convention |
| RMSD metric | Hungarian symmetry-corrected, heavy atoms only, ligand vs deposited crystal pose |
| Reference pose | `<PDB>_ligand.sdf` (deposited crystal ligand) — used **only** for RMSD scoring, never as GA input |
| Ligand torsions | flexible (rotatable bonds perceived from SDF bond orders) |
| Headline metric (**BCR**) | **B**inding-mode **C**orrect **R**ate = (# targets with `rmsd_hungarian < 2.0`) / 85 |

The "< 2 Å = success" rule is the published convention used by rDock, GOLD, Vina,
FlexAID 2015, GNINA, SurfDock, and Uni-Mol Docking. Using it verbatim is what
makes the FlexAIDdS BCR directly stackable against those papers.

### Shared engine env block

These are identical for all three tiers. They configure the search engine and
**must not be varied** when reporting a tiered number (varying them produces a
non-standard run that may not be quoted under a tier label):

```bash
# ── FlexAIDdS shared engine configuration (TIER-invariant) ──
export FLEXAIDDS_RESTARTS=5             # independent GA restarts, pooled
export FLEXAIDDS_PARALLEL_RESTARTS=0    # restarts run serially within a job
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1  # per-dihedral evaluation scaling
export FLEXAIDDS_CONSENSUS_SCORER=1     # consensus pose selection
export FLEXAIDDS_N_ELITE=1              # GA-internal elitism (1 elite snapshot)
export FLEXAIDDS_BUDGET_SCALE=1         # GA generation budget multiplier
export FLEXAIDDS_SOFTCORE_WAL=1         # softcore (capped) wall term
export FLEXAIDDS_SOFTCORE_FLOOR=0.5     # softcore floor
export FLEXAIDDS_T_HOT=500              # hot-replica temperature (K)
export FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1# relax receptor side-chain rotamers
export FLEXAIDDS_CHAIN_NORM=1           # chain-ID normalization (post-v102 fix)
export FLEXAIDDS_NATIVE_SEED_FRAC=0.0   # NO native pose seeding (banned)
```

### Common runner invocation

```bash
RUNNER=/path/to/benchmark_datasets          # the built C++ binary, NOT the python shim
caffeinate -i "$RUNNER" \
  --benchmark "crossdock_json:<TIER_JSON>" \  # JSON differs per tier (see below)
  --output    "<OUT>/astex_diverse" \
  --threads   5 \
  --omp-threads 1 \
  --temperature 298 \
  --job-timeout-seconds 7200 \
  --cache     "<OUT>/cache/astex_diverse" \
  --mode      <oracle|autonomous>            # differs per tier
```

### The `crossdock_json` input schema

The runner reads cross-docking pairs from a JSON file via the
`crossdock_json:<file>` benchmark prefix. Each pair object is located by its
`"receptor_id"` key and the following string fields are read
(see `LIB/benchmark_datasets.cpp:388-456`):

```json
{
  "pairs": [
    {
      "receptor_id":     "1G9V",
      "ligand_id":       "1G9V",
      "receptor_pdb":    "/abs/path/1G9V_apo.pdb",
      "ligand_sdf":      "/abs/path/1G9V_ligand.sdf",
      "oracle_site_pdb": "/abs/path/1G9V_binding_site.pdb"
    }
  ]
}
```

| JSON field | Maps to | Meaning |
|------------|---------|---------|
| `receptor_id` | `entry.pdb_id` | result row label (any object lacking this key is skipped) |
| `receptor_pdb` | `entry.receptor_path` | the structure to dock **into** |
| `ligand_sdf` | `entry.ligand_path` | the molecule to dock |
| `oracle_site_pdb` | `entry.binding_site_path` | crystal-centroid site sphere; **omit for TIER-3** |
| `rmsd_ref_sdf` *(optional)* | `entry.rmsd_reference_path` | explicit RMSD reference pose (defaults to `ligand_sdf`) |

**The single field that distinguishes self-dock from cross-dock is whether
`receptor_pdb` belongs to the same PDB as `ligand_sdf`.** If the receptor is the
apo form of the ligand's own crystal → TIER-1. If it is a different PDB's
structure → TIER-2/TIER-3.

> **Current repository state (audited 2026-06-20, commit `5721b6d`):**
> The committed Astex JSON is
> `benchmarks/datasets/benchmark_astex_native_85.json`. Every pair in it has
> `receptor_id == ligand_id` and `receptor_pdb` pointing at `<PDB>_apo.pdb`
> (the holo crystal with its ligand stripped). **It is therefore a TIER-1
> self-dock specification** (`"oracle_mode": true`).
> The `benchmarks/astex_diverse/astex_diverse.json` path referenced by recent
> launch scripts (`/tmp/launch_v103_chainnorm.sh`) is a cross-paired (TIER-2/3)
> JSON that is **generated at launch time and not committed to git** — it must be
> regenerated to reproduce v50b / v102 / v103. See §6.

---

## 2. TIER-1 — Self-Dock (cognate receptor, oracle site)

**Question answered:** *Given the right pocket and the right protein conformation,
can the GA + VCT scorer find the crystal pose?* This is a pure search/scoring test
with protein flexibility removed as a variable.

- **Receptor:** the cognate crystal structure with its own ligand removed
  (`<PDB>_apo.pdb`). `receptor_id == ligand_id`.
- **Ligand:** the deposited crystal ligand (`<PDB>_ligand.sdf`), docked flexibly.
  **No native seed.**
- **Site:** oracle crystal-centroid sphere of the same complex
  (`<PDB>_binding_site.pdb`).
- **Mode:** `oracle`.

### Launch (TIER-1)

```bash
#!/usr/bin/env bash
set -euo pipefail
REPO=/Users/lp.more/Projects/FlexAIDdS
RUNNER=/path/to/benchmark_datasets
OUT="$HOME/results/tier1_selfdock_$(date +%Y%m%d)"
ORACLE_DIR="$REPO/benchmarks/astex_diverse/astex_diverse"
mkdir -p "$OUT"

# shared engine block
export FLEXAIDDS_RESTARTS=5 FLEXAIDDS_PARALLEL_RESTARTS=0 \
       FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1 FLEXAIDDS_CONSENSUS_SCORER=1 \
       FLEXAIDDS_N_ELITE=1 FLEXAIDDS_BUDGET_SCALE=1 \
       FLEXAIDDS_SOFTCORE_WAL=1 FLEXAIDDS_SOFTCORE_FLOOR=0.5 \
       FLEXAIDDS_T_HOT=500 FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1 \
       FLEXAIDDS_CHAIN_NORM=1 FLEXAIDDS_NATIVE_SEED_FRAC=0.0

# TIER-1 deltas
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_ORACLE_SITE_DIR="$ORACLE_DIR"     # oracle site ON

caffeinate -i "$RUNNER" \
  --benchmark "crossdock_json:$REPO/benchmarks/datasets/benchmark_astex_native_85.json" \
  --output "$OUT/astex_diverse" --threads 5 --omp-threads 1 \
  --temperature 298 --job-timeout-seconds 7200 \
  --cache "$OUT/cache/astex_diverse" --mode oracle
```

### Published comparators (self-dock, RMSD < 2 Å, top-1)

| Method | BCR | Class |
|--------|-----|-------|
| rDock | 87.8 % | physics |
| GOLD | 64–72 % | physics |
| FlexAID 2015 | 67.9 % | physics |
| AutoDock Vina | 56–70 % | physics |

A TIER-1 BCR is comparable to these self-docking numbers. **Do not** compare a
TIER-1 number to any cross-docking result — it is an easier task.

---

## 3. TIER-2 — Oracle Cross-Dock (non-cognate receptor, oracle site)

**Question answered:** *Given the right pocket location but the **wrong** protein
conformation, can FlexAIDdS recover the crystal pose?* This adds induced-fit /
conformational-mismatch difficulty while still handing the method the binding-site
location. **This is the FlexAIDdS headline result** (v50b: 69/85 = 81.2 %).

- **Receptor:** a non-cognate structure — an apo or alternative-holo PDB that is
  *not* the co-crystal of this ligand. `receptor_id != ligand_id`.
- **Ligand:** the target ligand from **its own** co-crystal (`<TARGET>_ligand.sdf`).
- **Site:** oracle crystal-centroid sphere taken from the **target** complex
  (`<TARGET>_binding_site.pdb`), *not* the receptor complex. The site marks where
  the answer is in the target's frame, transferred onto the non-cognate receptor.
- **Mode:** `oracle`.

### Launch (TIER-2)

Identical to TIER-1 except the input JSON pairs each ligand with a **non-cognate
receptor**, and `oracle_site_pdb` points at the **target's** binding site:

```bash
# TIER-2 deltas (everything else as TIER-1)
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_ORACLE_SITE_DIR="$ORACLE_DIR"     # oracle site ON

caffeinate -i "$RUNNER" \
  --benchmark "crossdock_json:$REPO/benchmarks/astex_diverse/astex_crossdock_85.json" \
  --output "$OUT/astex_diverse" --threads 5 --omp-threads 1 \
  --temperature 298 --job-timeout-seconds 7200 \
  --cache "$OUT/cache/astex_diverse" --mode oracle
```

The TIER-2 JSON differs from TIER-1 only in that `receptor_pdb` is a different
PDB's structure. Example pair:

```json
{ "receptor_id": "2GBP", "ligand_id": "1L2S",
  "receptor_pdb":    "/abs/.../2GBP_apo.pdb",
  "ligand_sdf":      "/abs/.../1L2S_ligand.sdf",
  "oracle_site_pdb": "/abs/.../1L2S_binding_site.pdb" }
```

### Published comparators (oracle / pocket-conditioned cross-dock, RMSD < 2 Å)

| Method | BCR | Class |
|--------|-----|-------|
| SurfDock | 77 % | AI (deep learning) |
| Uni-Mol Docking | 69.2 % | AI |
| GNINA | 56 % | physics + CNN |
| classical physics methods | < 32 % | physics |
| **FlexAIDdS v50b** | **81.2 %** | physics |

A TIER-2 BCR is the number to quote when claiming FlexAIDdS "beats all physics
methods and matches/exceeds AI on cross-docking." **It must never be reported
without the word "oracle"**, because the method is given the pocket location.

---

## 4. TIER-3 — Blind Cross-Dock (autonomous, no crystal information)

**Question answered:** *With no crystal information at all — not even the pocket —
can FlexAIDdS dock?* This is the honest fully-autonomous number (v102: 22/85 =
25.9 %).

- **Receptor:** non-cognate apo/alternative-holo PDB. `receptor_id != ligand_id`.
- **Ligand:** the target ligand.
- **Site:** **SURFNET** binding-site detection run by the engine. **No crystal
  centroid is used at any stage.**
- **Mode:** `autonomous` (enables blinding; no seed elitism).

### How TIER-3 actually suppresses the oracle (two independent gates — both required)

The runner will use an oracle site if it can find one by **either** route, so a
blind run must close **both**:

1. `FLEXAIDDS_ORACLE_SITE_DIR` must be **unset / empty**, AND
2. the input JSON must **omit `oracle_site_pdb`** on every pair.

If either is present the run is silently a TIER-2 run. (`--mode autonomous` alone
does **not** make a run blind — mode controls pose blinding and seed elitism, not
the site source; `binding_site_path` still resolves from the JSON or the env dir
regardless of mode. See `LIB/DatasetRunner.cpp:3754-3818` and
`LIB/benchmark_datasets.cpp:435`.)

> **Guard rail:** the built-in `astex` benchmark path hard-aborts on a 0-oracle
> run (`LIB/DatasetRunner.cpp:3814`, "Set FLEXAIDDS_ORACLE_SITE_DIR. Aborting").
> The `crossdock_json` path does **not** abort, so TIER-3 must be driven through
> a `crossdock_json` file with `oracle_site_pdb` stripped — that is the only
> supported way to obtain a blind FlexAIDdS run.

### Launch (TIER-3)

```bash
# TIER-3 deltas (everything else as the shared engine block)
unset FLEXAIDDS_ORACLE_SITE_DIR        # site source OFF (gate 1)
# JSON has NO oracle_site_pdb fields    # site source OFF (gate 2)

caffeinate -i "$RUNNER" \
  --benchmark "crossdock_json:$REPO/benchmarks/astex_diverse/astex_crossdock_85_blind.json" \
  --output "$OUT/astex_diverse" --threads 5 --omp-threads 1 \
  --temperature 298 --job-timeout-seconds 7200 \
  --cache "$OUT/cache/astex_diverse" --mode autonomous
```

Blind JSON pair (note absence of `oracle_site_pdb`):

```json
{ "receptor_id": "2GBP", "ligand_id": "1L2S",
  "receptor_pdb": "/abs/.../2GBP_apo.pdb",
  "ligand_sdf":   "/abs/.../1L2S_ligand.sdf" }
```

### Published comparators (blind cross-dock, RMSD < 2 Å)

| Method | BCR | Class |
|--------|-----|-------|
| classical physics methods | 20–35 % | physics |
| AI methods | *no established blind cross-dock baseline yet* | — |
| **FlexAIDdS v102** | **25.9 %** | physics |

A TIER-3 BCR is the only number that may be described as "fully autonomous" or
"blind." It is the appropriate headline for a prospective / real-world claim.

---

## 5. Reading the output — where the BCR lives

The runner writes results in two granularities. **Compute the BCR yourself from
the per-target RMSD column** — do not trust a printed "success rate" until you
have confirmed it equals the RMSD-derived count.

### Per-target file: `<out_dir>/<PDB>/result.csv`

One row per target. Header (`LIB/DatasetRunner.cpp:6805`):

```
pdb_id,best_score,rmsd_to_crystal,rmsd_hungarian,predicted_dG,predicted_dH,
predicted_TdS,shannon_entropy,search_entropy_proxy,num_poses,wall_time_s,
success,cf_native,best_cluster_rmsd,best_cluster_idx,seed_echo,pose_source,...
```

### Aggregate file: `<output>/astex_crossdock_85_results.csv`

One row per target across the whole run. Header (`LIB/DatasetRunner.cpp:7115`):

```
pdb_id,best_score,rmsd_to_crystal,rmsd_hungarian,predicted_dG,predicted_dH,
predicted_TdS,shannon_entropy,search_entropy_proxy,num_poses,wall_time_s,
success,cf_native,best_cluster_rmsd,best_cluster_idx,seed_echo,pose_source,...
```

### Authoritative fields

| Field | Use |
|-------|-----|
| **`rmsd_hungarian`** | **The BCR field.** Symmetry-corrected RMSD of the selected pose vs crystal. A target counts as correct iff `rmsd_hungarian < 2.0`. |
| `success` | Convenience flag. In the current engine `success = (docking_completed && rmsd_hungarian >= 0 && rmsd_hungarian < 2.0)` (`LIB/DatasetRunner.cpp:6738`). It *should* equal the `rmsd_hungarian < 2.0` count — verify it does and report the RMSD-derived count regardless. (Historically `success` meant only "docking ran" and overstated accuracy.) |
| `rmsd_to_crystal` | Serial-order (non-symmetry-corrected) RMSD. Diagnostic only — **not** the success metric. |
| `seed_echo` | Must be **0** for every row in a standards-compliant run. A `1` means the pose is a copy of the input ligand (only possible if native seeding leaked in) and invalidates the run. |
| `best_score` | Selected-pose CF score. Do **not** use as a correctness metric — it is decoupled from the RMSD-scored pose. |

### Canonical BCR computation

```python
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))   # *_results.csv
n = len(rows)
correct = sum(1 for r in rows if 0.0 <= float(r["rmsd_hungarian"]) < 2.0)
echoes  = sum(1 for r in rows if r.get("seed_echo") == "1")
assert echoes == 0, f"INVALID: {echoes} seed-echo poses present"
print(f"BCR = {correct}/{n} = {100*correct/n:.1f}%")
```

The aggregate `<dataset>_summary.csv` also reports `success_rate`, `mean_rmsd`,
and `median_rmsd`; treat `success_rate` as a cross-check on the computation
above, not as the source of truth.

---

## 6. Reproducibility checklist

A run is reproducible only if **all** of the following are recorded alongside the
BCR. Capture them into a `provenance.json` next to the output (the runner already
emits an oracle/config block — see `LIB/DatasetRunner.cpp:4983-5010`,
`5563`). A reviewer with this checklist can rebuild and re-run independently.

- [ ] **Tier label** — `TIER-1` / `TIER-2` / `TIER-3`.
- [ ] **Git commit SHA** of the source tree that built the binary:
      `git rev-parse HEAD` (this document was written at `5721b6d`).
- [ ] **Working tree clean** — `git status --porcelain` empty, *or* the diff
      attached. (Memory: forbidden-file edits have silently confounded past
      A/B runs — an uncommitted tree is not reproducible.)
- [ ] **Binary SHA256** of both executables:
      `shasum -a 256 benchmark_datasets FlexAIDdS`.
- [ ] **Toolchain + platform** — compiler & version, OS, arch, CMake flags
      (`Release`/`LTO`, SIMD ISA, OpenMP/Metal/CUDA on/off). Binary SHA256
      differs across platforms; this records *why*.
- [ ] **Dataset JSON hash** — `shasum -a 256 <TIER_JSON>` of the exact
      `crossdock_json` file used (TIER-1 self / TIER-2 oracle / TIER-3 blind).
      For TIER-2/3 the cross-paired JSON is generated, not committed — archive
      the generated file, do not just reference its path.
- [ ] **Full env-var snapshot** — `env | grep ^FLEXAIDDS_ | sort` captured into
      the output dir. Must show `FLEXAIDDS_NATIVE_SEED_FRAC=0.0`; for TIER-3 must
      show `FLEXAIDDS_ORACLE_SITE_DIR` **absent**.
- [ ] **Runner invocation** — the exact `--mode`, `--threads`, `--omp-threads`,
      `--temperature`, `--job-timeout-seconds` line.
- [ ] **Determinism note** — the GA is non-deterministic across runs unless
      `FLEXAIDDS_SEED` + `ga.seed` are fixed and `--omp-threads 1`. ~6 Astex
      targets sit within 0.5 Å of the 2.0 Å cutoff, so a single run has run-to-run
      BCR noise of a few targets. Report **N≥3 repeats** (mean ± range) for any
      claimed delta smaller than ~4 targets, or fix the seed and state it.
- [ ] **Output artifacts** — archive `*_results.csv`, `*_summary.csv`, and the
      per-target `result.csv` files; the BCR must be recomputable from them
      via §5.

### One-liner provenance capture

```bash
{
  echo "tier: TIER-2"
  echo "commit: $(git rev-parse HEAD)"
  echo "tree_clean: $(git status --porcelain | wc -l | tr -d ' ') changes"
  echo "binary_sha256: $(shasum -a 256 "$RUNNER" | awk '{print $1}')"
  echo "engine_sha256: $(shasum -a 256 "$(dirname "$RUNNER")/FlexAIDdS" | awk '{print $1}')"
  echo "dataset_sha256: $(shasum -a 256 "$TIER_JSON" | awk '{print $1}')"
  echo "--- env ---"; env | grep '^FLEXAIDDS_' | sort
} > "$OUT/provenance.txt"
```

---

## 7. Quick decision guide

> *"Which tier is my run, and what may I compare it to?"*

```
Is receptor the ligand's own crystal (apo)?
├─ YES → is a crystal-centroid site supplied?
│        └─ YES → TIER-1 (self-dock).   Compare to rDock/GOLD/FlexAID2015/Vina.
└─ NO (non-cognate receptor)
         ├─ crystal-centroid site supplied → TIER-2 (oracle cross-dock).
         │                                    Compare to SurfDock/Uni-Mol/GNINA.
         └─ SURFNET site, no crystal info  → TIER-3 (blind cross-dock).
                                              Compare to physics 20–35%; no AI baseline.
```

If `seed_echo` is non-zero anywhere, the run belongs to **no tier** — it is
invalid and may not be quoted.

## 7. Grand canonical / competitive validation (P4+ extension)

This section defines sanctioned reporting for concentration-dependent
competitive binding using the GrandPartitionFunction (Ξ).

**Scope**: multi-ligand per-receptor campaigns where observables are
p_bind (occupancy probability per ligand at given c), mean_occupancy,
log_Xi, apparent + intrinsic selectivity, derived ΔΔG from Ki ratios.

**Not a replacement for docking-power BCR**. Grand validation augments
thermodynamic accuracy on top of pose correctness (PoseBusters + RMSD ≤2.0 Å
still required for any "successful" binding mode feeding the ensemble Z).

### Data
- Use `benchmarks/datasets/competition_example.yaml` (or derived real sets).
- Pair with per-ligand ensemble log_Z (from StatMechEngine on GA poses).
- Synthetic exact fixtures: `benchmarks/grand_synthetic/*.json` (analytical
  ground truth, Z + c supplied, expected p/Ξ/sel known to 1e-12).
- Literature: Ki ratios → ΔΔG = RT ln(Ki_weak/Ki_tight) at 298 K for
  selectivity checks (see yaml for conversion notes).

### Metrics (must name "grand canonical")
- p_bind agreement (RMSE or max abs err vs analytical/literature at stated concs)
- selectivity error (log units; intrinsic vs apparent distinguished)
- occupancy curve match (vary c, check p(empty) and total occ)
- ΔΔG recovery from selectivity (within experimental variance)
- Repro: full manifest + concentrations + engine SHA + seed in every report

### Harness
- `python3 scripts/grand_calibrate.py --synthetic ...` (pure-Py reference impl)
- Extended `scripts/validate_benchmark_results.py --manifest competition_example.yaml`
- Future: load_results() + grand post-processing in Python package.

### Reproducibility (in addition to §6 checklist)
- Record exact per-ligand concentrations (M) and T.
- Record provenance of every log_Z (run id, binding mode indices, whether
  CCBM / vib corrections active).
- HW parity: identical grand quantities (within tol) across scalar / Metal /
  CUDA builds when fed identical input log_Z (see GPF_IMPLEMENTATION_PLAN.md
  P4 HW notes).
- Two runs with same seed + manifest + concs → bitwise-identical grand summary
  or documented fp tolerance.

### Example reporting sentence
"Grand canonical validation on competition_example tier-1 synthetic fixtures
(n=3 cases) + literature Ki-ratio dual-ligand set: p_bind RMSE < 0.01,
intrinsic selectivity error < 0.05 log units; HW parity (Metal vs scalar)
held to 1e-9 rel on log_Xi. All using RngSeed=1234, T=298 K."

**Current status (P4)**: data + harness prep complete; integration with
docking outputs pending P1–P3. See GPF_IMPLEMENTATION_PLAN.md and
docs/GrandPartitionFunction_Report.md for observables and theory.

