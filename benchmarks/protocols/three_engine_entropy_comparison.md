# Protocol: Three-engine entropy comparison (queue-ready)

**Status:** Refined for computing-queue submission. Normative for any campaign that claims  
“FlexAID 2015 vs FlexAID (entropy) vs FlexAIDdS” on the same matrix.

**Audience:** Queue operators, thesis methods, agents launching array jobs.

**Related:** `benchmarks/BENCHMARK_STANDARD.md` (tiers, no-seed rule), `AGENTS.md` (scientific guardrails),
**`benchmarks/protocols/admission_metrics_contract.md`** (normative S1/S2/S3 + claim admission; enforced by `scripts/aggregate_claim_metrics.py`).

---

## 0. Scientific question

Under a **matched cognate-pocket, no native-seed** redocking protocol, and an **identical energy matrix**, how do docking power and ranking change across:

| Arm ID | Engine | Entropy ranking |
|--------|--------|-----------------|
| **A** | FlexAID **2015-era** (JCIM paper lineage; CF/VCT, no BindingMode free-energy ranking) | Off (era default) |
| **B0** | FlexAID **current master** | **Off** (`TEMPER 0` → CF clustering forced) |
| **B** | FlexAID **current master** | **On** (`TEMPER 21` + FO clustering / engine ACF; **not** DatasetRunner Softβ S1 rescoring of CF ensembles) |
| **C0** | FlexAIDdS (pinned SHA) | **CF / consensus election only** (DatasetRunner Softβ S1 default OFF; `FLEXAIDDS_SOFTBETA_ELECTION=0`) |
| **C** | FlexAIDdS (same SHA as C0) | Full stack available; **default claim path** = elected pose under standard DatasetRunner (document exact election); thermo ledger reported separately |

**Primary contrast for entropy *ranking*:** B0→B and C0→C on **shared ensembles** where possible.  
**Primary contrast for engine generation:** A vs B0 vs C0 (CF election, same matrix, same site).

**Out of scope for this queue campaign:**

- Seeded / native-pose-inheritance runs (including any `native_pose_seeded=1` restore).
- Cross-dock or blind SURFNET (TIER-2/3) until TIER-1 completes.
- Claiming CF as experimental ΔG.

---

## 1. Hard invariants (fail preflight if violated)

### 1.1 Energy matrix — single shared pin

**All five arms use the exact same interaction matrix file.**

| Field | Value |
|-------|--------|
| Canonical file name | `MC_st0r5.2_6.dat` (or lab-standard alias if identical bytes) |
| **Pinned MD5** | **Record at submit time** with `md5 -q <matrix>`; must match on every node and every data dir |
| Known local copies (verify, do not assume) | e.g. `FlexAIDdS/WRK/MC_st0r5.2_6.dat` → `9dc93717…` (baseline-validated JCIM matrix — the MD5 the engine loads on disk); a stray `72d7c739…` packing-sweetened fork may exist in some trees, do **not** use it — **copy the 9dc9 file into every arm’s `DEPSPA` / data dir** |

**Preflight (every job):**

```bash
test "$(md5 -q "$MATRIX_PATH")" = "$MATRIX_MD5" || exit 90
```

No matrix-swap ablation. No “default M6 vs MC_*” drift. Atom-type defs (`AMINO*.def`, `NUCLEOTIDES*.def`) ship from the **same data bundle** as the matrix.

### 1.2 Protocol tier

**TIER-1 cognate redock, no seed** (aligned with `BENCHMARK_STANDARD.md`):

| Axis | Setting |
|------|---------|
| Dataset | Astex Diverse, **N = 85** |
| Receptor | `<PDB>_apo.pdb` (cognate, ligand removed) |
| Ligand | `<PDB>_ligand.sdf` flexible; **pose-blinded** before search |
| Site | Cognate pocket spheres / binding-site PDB (GetCleft or frozen lab set) — **known pocket only** |
| Native pose seed | **Forbidden** (`seed_fraction=0`, `pose_seed_enabled=false`, no `_INI` elitism) |
| Crystal pose | **RMSD / PoseBusters reference only** |
| Restarts | **5** independent, results pooled for election |
| Temperature (thermo arms) | C0/C: **298 K**; FlexAID B: **TEMPER 21** (locked) |
| Softcore / VCT extras | Document in provenance; prefer identical where both engines support the knob |

### 1.3 Search budget (matched effort)

| Parameter | Target | Notes |
|-----------|--------|--------|
| Population (**base**) | 1000 | Base chromosomes; **this is the modulated axis** for DoF (see below) |
| Generations | **2000 fixed** (claim freeze 2026-07-15; prior draft used 6000) | Do **not** inflate generations for flexible ligands |
| Restarts | 5 | Same RNG scheme: `SEED_BASE + stable_hash(pdb_id, restart_i)` |
| Wall clock / job | ≤ 3 h per target×arm (queue default); fail and retry once | |

#### 1.3.1 DoF budget modulation — **population, not generations** (normative)

Hard ligands need more search **diversity**, not longer trajectories. DatasetRunner implements this in `LIB/DatasetRunner.cpp` (Lever 1 + high-DoF budget):

| Knob | Value for three-engine / TIER-1 claim | Effect |
|------|--------------------------------------|--------|
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | **`1`** (default) | `pop_eff = pop_base × max(1, n_flex_bonds/4)`; **`n_gen` stays at base** |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | `0` | **Legacy — forbidden for new claim runs:** scales **generations** by `ceil(n_genes/4)` |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | `-1` / `off` | Fixed pop+gen — **oracle-ceiling / restore only**, not three-engine claim |
| `FLEXAIDDS_BUDGET_SCALE` | `1` (default ON) | Extra **population** multiplier when `n_genes ≥ 14` (`× max(1, n_genes/7)`); gens still fixed |

**Anti-pattern (Codex / agents):** treating CLI `--ga-population 1000 --ga-generations 6000` as a fully fixed budget and setting `EVAL_SCALE_DIHEDRAL=-1` “to match.” That **disables** the intended chromosome (pop) modulation. Receipts may still show `pop: 1000, gen: 6000` while logs show effective pop via `[EVAL-BUDGET] … pop=…`.

**Correct reading of logs:**

```text
[EVAL-SCALE]  … FIXED pop=1000 n_gen=2000     ← mode -1: NO dihedral pop-scale (wrong for claim path)
[EVAL-SCALE]  … pop_base=1000 pop_effective=N n_gen=2000  ← mode 1: correct (gens fixed, pop grows with DoF)
[EVAL-BUDGET] … budget_scale=S n_gen=2000 pop=P           ← final pop after high-DoF multiplier; n_gen must stay 2000
```

If FlexAID-2015 cannot sustain large `pop_eff`, freeze a **budget ladder** in provenance and match **total CF evaluations** where possible — still prefer **pop** ladder over **gen** ladder:

```text
budget_class: full | half | quarter
evals_nominal: pop_effective * gen_fixed * restarts
evals_actual: logged
modulation_axis: population   # never generations for claim runs
```

Do not mix budget classes or modulation axes in the headline table.

### 1.4 Success metrics (all reported; abstract uses S1)

| ID | Definition | Role |
|----|------------|------|
| **S1** | Top-1 elected pose RMSD ≤ 2.0 Å (Hungarian, heavy atoms) | **Primary / queue KPI** |
| **S2** | S1 ∧ PoseBusters pass on elected pose | Modern secondary |
| **S3** | Any emitted cluster pose RMSD ≤ 2.0 Å (BCR / sampling ceiling) | Diagnostic only |
| **S4** | Rank of native-containing mode (0 = best) | Ranking / entropy story |
| **seed_echo** | Elected path is `_INI` or protocol seed flag | Must be **0** for claim rows |

**Queue must not** summarize S3 as “success rate.”

---

## 2. Binary pins

| Arm | Source | Pin method |
|-----|--------|------------|
| A | FlexAID JCIM-2015 lineage | Git tag/commit **or** frozen binary SHA256 used for paper reproduction |
| B / B0 | `/path/to/FlexAID` **master** | `git rev-parse HEAD` + binary SHA256 |
| C / C0 | FlexAIDdS repo | `git rev-parse HEAD` + `build/FlexAIDdS` or `build_lto/FlexAIDdS` SHA256 |

**Build once, stage to read-only `$QUEUE_ROOT/bin/{A,B,C}/`.** Jobs never rebuild.

---

## 3. Arm-specific configuration

### 3.1 Shared inputs per target

Frozen under `$QUEUE_ROOT/inputs/astex_diverse/<PDB>/`:

```text
<PDB>_apo.pdb
<PDB>_ligand.sdf
<PDB>_site.pdb          # cleft / binding-site spheres (same for all arms)
```

Single generator script builds this tree once; all arms symlink it.

### 3.2 Arm A — FlexAID 2015

- Classic `.inp` / LOCCLF (or LOCCEN) from **same** site spheres.
- No BindingMode free-energy ranking.
- Election: paper-era default (lowest CF / cluster protocol as implemented).
- Output: cluster PDBs + scores → normalized via adapter.

### 3.3 Arm B0 — FlexAID master, entropy off

- Same inputs as B.
- **`TEMPER 0`** in docking input → engine forces **CF** clustering (`read_input.c`: temperature ≤ 0 overrides CLUSTA to CF).
- Election: CF rank-0 (and document frequency gate if any).

### 3.4 Arm B — FlexAID master, entropy on

- **`TEMPER 21`** (locked claim freeze 2026-07-15; operator-optimized). Clustering: FO when T>0 (generator default).
- BindingMode / colony entropy path active.
- Election: free-energy / ACF ranking as master implements; dump `dG,dH,TdS` when present.

### 3.5 Arm C0 — FlexAIDdS CF election

- DatasetRunner **no-seed** path (ORACLE_CEILING or DEFINED_CLEFT_REDOCK with seed off — **not** native flood).
- Known site env / cleft sphere file.
- Election: CF / consensus **without** thermo override for the claim column.
- PoseBusters + optional tENCoM for S2 / claim_ready columns (claim_ready not required for primary KPI).

### 3.6 Arm C — FlexAIDdS production claim path

- Same binary and search budget as C0.
- Prefer **re-rank frozen C0 ensembles** with thermo stack when possible (true entropy isolation).  
  If full re-dock is used instead, document as “search+rank coupled.”

### 3.7 Ideal entropy isolation (preferred when I/O allows)

```text
Search once per (engine_family, pdb) → freeze pose ensemble
  ├─ rank CF     → B0 / C0 metrics
  └─ rank F/S    → B  / C  metrics
```

FlexAID A has no F path → CF only.

---

## 4. Queue layout

```text
$QUEUE_ROOT/
  protocol.md                 # copy of this file + submit SHA
  provenance.json             # matrix MD5, binary SHAs, seed base, dataset path
  bin/{A,B,C}/                # staged executables + data dirs (shared matrix)
  data/                       # MC_st0r5.2_6.dat + AMINO*.def + … (one tree)
  inputs/astex_diverse/<PDB>/
  work/<ARM>/<PDB>/           # per-job scratch
  results/<ARM>/<PDB>/result.csv
  results/normalized/<ARM>.csv
  results/joined/astex85_paired.parquet  # or csv
  logs/<ARM>/<PDB>.{out,err}
  aggregate/summary.json
```

### 4.1 Job array design

| Dimension | Spec |
|-----------|------|
| Array index | `0 … 85*N_ARMS - 1` **or** nested: arm × pdb |
| Mapping | `pdb = ASTEX[i % 85]`, `arm = ARMS[i // 85]` |
| Threads | 1 job × `OMP_NUM_THREADS=4–6` (match machine; avoid oversubscription) |
| Memory | Request ≥ 4 GB / job (5 restarts serial preferred); if parallel restarts, ≥ 12 GB |
| Time | 3 h default; 6 h for high-DoF ligands if needed |
| Idempotent | Skip if `results/<ARM>/<PDB>/result.csv` exists and `seed_echo=0` and RMSD parsed |

**Serial restarts inside a job** (`PARALLEL_RESTARTS=0`) preferred on shared nodes to reduce RAM spikes.

### 4.2 Pilot before full array

| Stage | N | Arms | Gate to proceed |
|-------|---|------|-----------------|
| **P0** Preflight | 0 | — | Matrix MD5 match; binaries run `-h`/smoke; site files exist for 85 |
| **P1** Smoke | 2 PDB (1 rigid, 1 flexible) | A,B0,B,C0 | Completes; seed flags 0; RMSD finite |
| **P2** Pilot | 8 hard panel | all | S1 rates logged; no dual matrix; wall times within 2× median |
| **P3** Full | 85 | all | Aggregate + stats |

**Do not submit P3 until P2 gate passes.**

---

## 5. Normalized result schema

One row per `(arm, pdb_id)`:

```text
arm, engine_sha, matrix_md5, pdb_id,
rmsd_top1, rmsd_bcr, success_s1, success_s2, success_s3,
rank_native_mode, n_poses, n_modes,
score_top1, H, TS, F,          # NA if unavailable
pb_pass, tencom_status,        # NA on FlexAID if not run
seed_echo, native_pose_seeded, protocol_claim_eligible,
wall_s, restarts_finished, evals_actual,
budget_class
```

**Admission rule for claim table:**  
`native_pose_seeded == 0 AND seed_echo == 0 AND matrix_md5 == PIN`.

**Enforced aggregation** (claim filters + separate S1/S2/S3; S3 never primary):

```bash
python3 scripts/aggregate_claim_metrics.py <results/<ARM>> [--json out.json]
# C0 full85 after source ~/.flexaidds_env:
python3 scripts/aggregate_claim_metrics.py --c0-full85
# FAIL: python3 scripts/aggregate_claim_metrics.py <dir> --headline s3   # needs --diagnostic-only
```

Default matrix pin: `9dc93717dfed0698006d88dd6a9627bc` (or `RUN_RECEIPT.json` / `provenance.json`).  
Full contract: `benchmarks/protocols/admission_metrics_contract.md`.

---

## 6. Statistics (run after P3)

| Analysis | Method |
|----------|--------|
| Rate ± CI | Wilson or bootstrap over 85 targets |
| A vs B0, B0 vs C0, A vs C0 | McNemar on S1 (and S2) |
| B0 vs B, C0 vs C | McNemar on shared-ensemble re-rank if available; else paired re-dock |
| Election gap | mean(S3 − S1), list targets with S3=1,S1=0 |
| Rescue rate | P(S1_thermo=1 | S1_CF=0) with bootstrap CI |
| Stratify | n_rotors, MW, pocket volume (optional) |

**Headline table columns:** arm, S1%, S2%, S3%, gap, median native rank, mean wall_s.

---

## 7. Hypotheses (pre-registered)

1. **H1 (generation, CF-only):** S1(C0) ≥ S1(B0) ≥ S1(A) under identical matrix and site.  
2. **H2 (entropy ranking):** On frozen ensembles, thermo re-rank raises S1 vs CF (rescue CI excludes 0).  
3. **H3 (no seed magic):** S1 for all arms is **not** comparable to seeded ceilings (~90%+ BCR under inheritance).  
4. **H4 (PB):** S2 < S1 for modern arms; report both.

---

## 8. Preflight checklist (copy into queue submit script)

```bash
# 1) Matrix
test -f "$DATA/MC_st0r5.2_6.dat"
test "$(md5 -q "$DATA/MC_st0r5.2_6.dat")" = "$MATRIX_MD5"

# 2) Binaries staged
for a in A B C; do test -x "$QUEUE_ROOT/bin/$a/"*FlexAID* || test -x "$QUEUE_ROOT/bin/$a/"*FlexAIDdS*; done

# 3) Inputs complete
test "$(find "$QUEUE_ROOT/inputs/astex_diverse" -name '*_apo.pdb' | wc -l)" -eq 85

# 4) No seed env pollution
test "${FLEXAIDDS_NATIVE_SEED_FRAC:-0}" = "0"
unset FLEXAIDDS_FORCE_SEED 2>/dev/null || true

# 5) PoseBusters for S2 (FlexAIDdS / post)
test -x "$FLEXAIDDS_POSEBUSTERS_BIN"

# 6) Provenance written
test -f "$QUEUE_ROOT/provenance.json"
```

`provenance.json` minimum keys:

```json
{
  "protocol": "three_engine_entropy_comparison",
  "protocol_version": "1.0",
  "matrix_path": "...",
  "matrix_md5": "...",
  "arms": {
    "A":  {"git": "...", "binary_sha256": "..."},
    "B":  {"git": "...", "binary_sha256": "..."},
    "C":  {"git": "...", "binary_sha256": "..."}
  },
  "dataset": "astex_diverse_85",
  "tier": "TIER-1_cognate_nativeseed_forbidden",
  "seed_base": 20260714,
  "pop": 1000,
  "gen": 6000,
  "restarts": 5,
  "temperature_K": 298
}
```

---

## 9. Resource estimate (planning)

| Item | Estimate |
|------|----------|
| Jobs | 85 × 5 arms = **425** (or 85 × 3 if B0/B and C0/C are re-rank only → **255** search jobs) |
| Core-hours | ~0.5–2 h × 425 ≈ **200–850 core-h** (machine-dependent) |
| Disk | ~0.5–2 GB / target / arm with logs → reserve **200–500 GB** for full trees |
| Preferred | Re-rank design: **255** heavy docks + cheap CPU re-rank |

---

## 10. Explicit exclusions

| Artifact / practice | Status |
|---------------------|--------|
| `oracle_ceiling_restore_v43proto_r3` seeded campaign | **Not** arm C; diagnostic only |
| Different MD5 matrices per engine | **Forbidden** |
| Native seed fraction > 0 | **Forbidden** |
| Reporting BCR as abstract success | **Forbidden** |
| Comparing TIER-1 to cross-dock literature without label | **Forbidden** |

---

## 11. Submit sequence (operator)

1. Freeze matrix MD5 → stage identical `data/` for A/B/C.  
2. Stage binaries + write `provenance.json`.  
3. Build `inputs/astex_diverse` once.  
4. Run **P1** (2 targets) → **P2** (8 targets).  
5. Human review: seed flags, matrix hash, wall times.  
6. Submit **P3** array (85 × arms).  
7. Aggregate → McNemar → `aggregate/summary.json` + methods table.

---

## 12. One-paragraph methods blurb (for queue ticket / paper)

We compare FlexAID (2015-era CF docking), FlexAID (current master with and without BindingMode entropy ranking via TEMPER 0/298), and FlexAIDdS under a matched **Astex Diverse (N=85) cognate-pocket redocking protocol with native pose seeding disabled**. All arms use the **same interaction matrix** (single MD5-pinned `MC_st0r5.2_6.dat` and type definitions), five independent GA restarts, and Hungarian RMSD ≤ 2 Å of the elected pose as the primary success criterion, with PoseBusters as a secondary filter. Entropy effects are isolated by re-ranking frozen ensembles (CF vs free-energy) where possible; engine-generation effects are read from CF-only arms.

---

## 13. Open decisions (resolve before P3, not during)

| Decision | Default if unset |
|----------|------------------|
| Exact FlexAID-2015 commit/binary | **Block submit** until pinned |
| Matrix MD5 choice among local copies | Use the baseline-validated JCIM matrix (`9dc93717…`, the MD5 the engine loads on disk); ignore any stray `72d7c739…` packing-sweetened fork; **copy that file** to all arms |
| FlexAIDdS mode name (`oracle-ceiling` no-seed vs `defined-cleft-redock`) | Prefer **defined-cleft-redock** or no-seed oracle with pose-blind — same site file either way |
| Parallel vs serial restarts | Serial restarts in-job |
| Run B/C as re-rank only? | **Yes** if ensembles portable |

---

**Document version:** 1.0 (refined for queue)  
**Changelog:** Matrix identity elevated to hard invariant; seeded campaigns excluded; TEMPER 0 defined as B0; job array + preflight + budget classes; S1 primary KPI.
