# Audit: `e18cb6b8c` — DP vs FO clustering small-sim pilot (8 random Astex)

| Field | Value |
|-------|--------|
| **Short SHA** | `e18cb6b8c` |
| **Full SHA** | `e18cb6b8cf85d6bf5851892128d595d6df3f37ef` |
| **Subject** | Add: DP vs FO clustering small-sim pilot (8 random Astex) |
| **Author / date** | LP \<lp@thebonhomme.com\> · 2026-07-15 03:33:17 -0400 |
| **Parent** | `4277b416f7ab49aa402b3d9b85a6911394f9569e` (C0 claim clean reboot) |
| **Immediate follow-up** | `c93b866a4` — Fix: DPFO aggregator RMSD sentinels and claim CSV columns |
| **Diffstat** | 2 files, +507 / −0 |
| **Audit scope** | Pilot scripts only (no engine source in this commit). Isolation from C0 claim OUT; budget pop=200 gen=50 scientific validity; random seed / repro; FO vs DP fair compare. |
| **Audit mode** | READ-only inspection of commit + live pilot OUT (when present). **No source edits.** |
| **Auditor** | Grok Build subagent (26h-swarm) |
| **Verdict** | **CONDITIONAL PASS as plumbing pilot** · **FAIL as scientific FO-vs-DP ranking** · **aggregator at this SHA is incorrect for DatasetRunner CSVs** (fixed next commit) |

---

## 1. Summary

This commit adds two orchestration scripts for an isolated **Density Peak (DP) vs FastOPTICS (FO)** clustering pilot:

1. `scripts/run_DPFO_pilot8_small.sh` — launches 8 Astex pairs × two arms (FO, DP) at **pop=200, gen=50**, `EVAL_SCALE_DIHEDRAL=-1`, crystal seed OFF.
2. `scripts/aggregate_dpfo_pilot.py` — aggregates S1 / BCR / election_gap / **P(S1\|BCR)** and writes `DPFO_SUMMARY.{json,md}`.

**What the commit message claims:** isolated FO/DP arms for success-rate and predictive power **without touching clean C0 claim OUT**.

**What the audit finds:**

| Axis | Rating | One-line finding |
|------|--------|------------------|
| **C0 claim isolation** | **PASS** (strong) | Dedicated OUT namespace; hard refuse if path contains `C0_full85_claim`; yields to live FlexAIDdS by default. |
| **Budget scientific validity** | **FAIL for ranking** | ~0.2% of C0 claim base evals; expected near-zero S1/BCR; useful only as plumbing / packaging smoke. |
| **Random sample + seed repro** | **PARTIAL PASS** | Manifest sample seed `20260715` is fixed; GA seeds match FO↔DP per target; crystal seed OFF is real. Sample generator not in-repo. |
| **FO vs DP fair compare** | **FAIL for metrics** | Matched seeds/budget/mode, but FO packaging/election failed on live run (RMSD sentinels −1, empty `elected_pose_path`) while DP elected 8/8; triple MinPts FO still present on the binary used. |
| **Aggregator correctness (this SHA)** | **FAIL** | Missing `rmsd_to_crystal`; accepts RMSD `< 0` as success; BCR falls back to S1. Fixed in `c93b866a4`. |

**Live pilot campaign (disk, post-commit run):**  
`$FLEXAIDDS_RESULTS/campaigns/DPFO_pilot8_small_g50_p200_20260715`  
Targets: `1JD0, 1L7F, 1N46, 1OF1, 1OPK, 1OWE, 1SJ0, 2BYS` · both arms **0/8 S1 and 0/8 BCR**.

---

## 2. Commit inventory

```
scripts/aggregate_dpfo_pilot.py   | 276 +
scripts/run_DPFO_pilot8_small.sh  | 231 +
2 files changed, 507 insertions(+)
```

No C++ / CMake / DatasetRunner / claim launcher changes. Pure pilot harness.

---

## 3. Pilot launcher — `scripts/run_DPFO_pilot8_small.sh`

### 3.1 Protocol knobs (as coded)

| Knob | Value | Notes |
|------|-------|--------|
| Population | `DPFO_POP` default **200** | CLI `--ga-population` |
| Generations | `DPFO_GEN` default **50** | CLI `--ga-generations` |
| Temperature | **298** | |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | **−1** | Fixed pop+gen (no DoF pop scale). Documented in engine as oracle-ceiling / fixed-budget mode — **not** claim mode. |
| `FLEXAIDDS_BUDGET_SCALE` | **1** (on) | High-DoF (`n_genes ≥ 14`) can still inflate **population** even under EVAL_SCALE=−1. Pilot comment “fixed small budget” is only partially true. |
| Restarts | `FLEXAIDDS_RESTARTS=2` | Claim C0 uses **5**. |
| Crystal / native seed | `NATIVE_SEED_FRAC=0`, `SEED_ELITISM=0`, `unset FLEXAIDDS_FORCE_SEED` | “seed OFF” = **no crystal pose injection**, not “GA unseeded”. |
| Mode | `defined-cleft-redock` | Same family as C0 claim. |
| Clustering | `--clustering FO` or `DP`; DP also sets `FLEXAIDDS_USE_DP=1` | FO clears `USE_DP`. |
| Threads | `--threads 1 --omp-threads 1` | Good for determinism. |
| Job timeout | 1800 s / target | |
| Parallel restarts | 0 | Serial restarts. |
| Nice | 19 | Background relative to claim. |
| `--force` | always | Overwrites pilot OUT targets. |

### 3.2 Paths and defaults

| Resource | Default |
|----------|---------|
| Manifest | `$Q/inputs/astex_dpfo_pilot8_random.json` |
| Runner / binary | `$Q/bin/C/benchmark_datasets`, `$Q/bin/C/FlexAIDdS` |
| Base OUT | `$R/campaigns/DPFO_pilot8_small_g50_p200_20260715` |
| Arm OUT | `$BASE_OUT/FO`, `$BASE_OUT/DP` |
| Logs / lock | `$Q/logs/DPFO_pilot8_small.{log,lock,pid}` |

`$Q` / `$R` come from `use_icloud_benchmark_storage.sh` or env (`FLEXAIDDS_QUEUE_ROOT`, `FLEXAIDDS_RESULTS`) — defaults under **iCloud CloudDocs**, not local-first `$FLEXAIDDS_LOCAL_ROOT`.

### 3.3 Provenance block

Writes `$BASE_OUT/PROVENANCE.json` with:

- `manifest_sha256`, `pdb_ids`, `seed_sample: 20260715`
- protocol snapshot (pop/gen/T/EVAL_SCALE/arms/metrics)
- `runner_sha256`, `binary_sha256`
- note: does not touch C0 claim OUT; yields to live FlexAIDdS unless `--nowait`

Good for audit pins. Does **not** pin DatasetRunner `git_commit` itself (runner does via `RUN_RECEIPT.json`).

---

## 4. Isolation from C0 claim OUT

### 4.1 What is protected

| Mechanism | Assessment |
|-----------|------------|
| Distinct campaign id | `DPFO_pilot8_small_g50_p200_${STAMP}` vs `C0_full85_claim_g2000_popmod_20260715` |
| Hard refuse | `BASE_OUT` contains `C0_full85_claim` → exit 91 |
| Default OUT path | Under `…/campaigns/DPFO_…` — string-disjoint from claim |
| Process yield | Default waits while `pgrep` sees `/bin/C/FlexAIDdS` or `bin/C/FlexAID`; `--nowait` opts out |
| Lock | Own lock/pid under queue logs; refuse if pilot already live (exit 92) |
| Manifest | Separate JSON (8 pairs), not full85 claim manifest |

**Live disk check:** claim OUT and pilot OUT are sibling directories under `…/results/campaigns/` — no path overlap. Pilot does not write into claim tree.

### 4.2 Residual operational risks (not claim-file corruption)

1. **Shared binary/queue tree** (`$Q/bin/C/*`, `$Q/data`) — pilot and claim can race on machine resources; pilot does not rewrite claim OUT but **can starve claim** with `--nowait`.
2. **Direct iCloud OUT** — violates local-first / thin-iCloud policy used by `run_C0_claim_clean.sh`. Risk is FileProvider hang / ops pain, not silent claim overwrite.
3. **`Path.rglob("result.csv")` in aggregator** under CloudDocs — conflicts with `AGENTS.md` “never tree-walk CloudDocs” rule; can hang ops monitors.
4. **`pgrep -f '/bin/C/FlexAIDdS|…'`** is coarse: any matching binary (smoke, other campaigns) blocks or is treated as “claim”.
5. **Shared log directory** on queue root — name-prefixed, low collision risk.

**Isolation verdict: PASS for “must not touch C0 claim OUT.”** Intention and implementation both succeed on that narrow contract.

---

## 5. Budget: pop=200 gen=50 — scientific validity

### 5.1 Arithmetic (base; no DoF scale)

| Campaign | Pop | Gen | Restarts | Base evals / target | Notes |
|----------|-----|-----|----------|---------------------|--------|
| **This pilot** | 200 | 50 | 2 | **20 000** | EVAL_SCALE=−1 (fixed) |
| **C0 claim clean** | 1000 | 2000 | 5 | **10 000 000** (+ pop×DoF when EVAL_SCALE=1) | Claim path |

**Ratio (base, no DoF):** pilot ≈ **0.2%** of C0 claim eval budget per target.  
With claim pop×DoF (e.g. `n_flex_bonds=12` → pop_eff=3000), claim is larger still (~30M evals).

### 5.2 Is the budget “too small”?

**Yes, for any scientific ranking of FO vs DP on Astex success rates.**

- Search space for flexible ligands is not meaningfully covered at 200×50.
- Live pilot: **0/8 S1 and 0/8 BCR on both arms** — consistent with under-sampling, not with a clustering winner.
- BCR (oracle ceiling over emitted cluster heads) never hit ≤2 Å either → **sampling ceiling not reached**; predictive power **P(S1\|BCR)** is undefined/null when `n_bcr=0` (aggregator correctly returns `null`).
- Comment in launcher labels this a “small-sim” pilot; live `ANALYSIS_DP_vs_FO_indepth.md` also states **plumbing only**. That framing is correct.

### 5.3 EVAL_SCALE=−1 vs claim contract

Per `LIB/DatasetRunner.cpp` / `AGENTS.md`:

- Claim: **modulate population, keep generations fixed** (`EVAL_SCALE=1`).
- Mode `−1`: **fixed pop+gen** — “oracle-ceiling restore only” / freeze budget.
- Pilot intentionally freezes budget for a cheap A/B; that is fine for plumbing, but results **must not** be cited as claim-relevant clustering evidence.

### 5.4 BUDGET_SCALE still ON

`FLEXAIDDS_BUDGET_SCALE=1` remains enabled. High-DoF targets can still get extra population. Arms remain matched to each other, but “fixed 200×50” is not absolute for every ligand.

**Budget verdict: scientifically too small for FO/DP ranking; acceptable for isolation + packaging smoke only.**

---

## 6. Random sample and reproducibility

### 6.1 Target sample (“8 random Astex”)

Manifest (not in git; lives under queue inputs):

```text
name: astex_dpfo_pilot8_random_20260715
seed: 20260715
pairs (n=8): 1JD0, 1L7F, 1N46, 1OF1, 1OPK, 1OWE, 1SJ0, 2BYS
```

- Sample seed is **recorded** in manifest + pilot `PROVENANCE.json` (`seed_sample: 20260715`).
- **No generator script** in this commit or repo search hit that rebuilds the 8-subset from Astex 85 with that seed. Repro of *which* 8 PDBs requires the pinned JSON (sha256 in PROVENANCE: `b7108a2d53e269e72dc2364497070aa12e4b63bc4253a7e082fdab9e900ba3d7`).
- Paths inside manifest are **absolute CloudDocs** paths — machine-local, not portable.

### 6.2 “seed OFF” vs GA RNG seed

| Layer | Behavior |
|-------|----------|
| Crystal / native seed | OFF (`seed_fraction=0`, `pose_seed_enabled=false`, elitism 0) |
| GA `ga.seed` in `dock_config.json` | **Deterministic** via `deterministic_ga_seed(pdb_id, restart, seed_base)` with default `seed_base=0` |
| FO vs DP same target | **Identical seeds** observed on disk (e.g. 1JD0 root `374870301`, r1 `489846628` for both arms) |

So FO and DP re-run independent full GA trajectories with **matched RNG streams** (intended fair clustering A/B), not “shared ensemble, recluster only.” Under single-thread + same binary pin, this is the right design for a clustering-only contrast **if packaging succeeds on both arms**.

### 6.3 Repro gaps

1. Manifest not versioned in the git commit (external queue artifact).
2. Binary/runner pins are SHA256 of whatever is currently in `$Q/bin/C` — not built from this commit’s tree (this commit has no engine change).
3. Live `RUN_RECEIPT.json` `git_commit` differs FO vs DP (`4277b416f` vs `c93b866a4`) because HEAD moved between arm launches; **binary/runner SHA256 matched** — good, but process provenance is sloppy.
4. Metal GPU pairwise RMSD in FO clustering logs (`Metal GPU pairwise RMSD precomputed` ×3 for triple MinPts) can introduce **non-CPU determinism** across runs even with fixed seeds.
5. Aggregator at this SHA is non-repro for metrics (see §8); post-fix `c93b866a4` required for correct summary.

**Seed / repro verdict: GA seed matching FO↔DP is solid; crystal seed OFF is solid; sample selection is pinned only by external JSON; full bit-exact re-run not guaranteed (Metal, OUT on CloudDocs, `--force`).**

---

## 7. FO vs DP fair compare

### 7.1 What is fair (matched)

| Factor | FO | DP | Match? |
|--------|----|----|--------|
| Targets (8) | same manifest | same | ✓ |
| pop / gen / T | 200 / 50 / 298 | same | ✓ |
| EVAL_SCALE | −1 | −1 | ✓ |
| Restarts | 2 | 2 | ✓ |
| Mode | defined-cleft-redock | same | ✓ |
| Crystal seed | OFF | OFF | ✓ |
| GA seeds | deterministic | same values | ✓ |
| Threads | 1 | 1 | ✓ |
| Binary / matrix pin | same SHA256 | same | ✓ |
| Clustering CLI | `--clustering FO` | `--clustering DP` + `USE_DP=1` | ✓ intent |

### 7.2 What is **not** fair (live packaging / protocol)

| Issue | Severity | Evidence |
|-------|----------|----------|
| **FO election/packaging null** | **Critical** | FO: all 8 targets `rmsd_to_crystal=-1`, `best_cluster_rmsd=-1`, empty `elected_pose_path`. DP: 8/8 elected with finite RMSDs (all S1/BCR still >2 Å). S1/BCR rates on FO measure **packaging failure**, not clustering quality. |
| **Triple MinPts FO ladder** | **High** | Live FO emits dual-suffix `_{5,7,10}_*.pdb` and three Metal RMSD precomputes; stderr: `0 binding modes after clustering (minPoints=5/7/10)`. Protocol docs (`docs/implementation/3dsig_red_pair_protocol.md`) forbid multi-scale FO for science arms. Binary used predates or does not enforce single literature MinPts (`6ec671a92` landed later in the swarm window). |
| **Independent full GA, not recluster-only** | Medium (design choice) | Clustering contrast is confounded by any non-determinism in GA/Metal even with matched seeds. Ideal fair test: freeze FO ensemble, recluster with DP only. |
| **FO dual-suffix vs DP single-suffix emission** | Medium | DatasetRunner has dual-suffix enumeration (`DatasetRunner.cpp` ~1072+), but FO still failed to elect on this pilot — either no parseable heads, empty BindingPopulation, or election path bail-out. DP single-suffix path is the well-exercised packaging path. |
| **Metrics at low BCR base rate** | High for inference | With BCR=0 on both arms, **P(S1\|BCR)** cannot distinguish FO from DP. |
| **n=8** | High for inference | Binomial noise dominates any rate difference; no CI / McNemar / paired test in aggregator. |

### 7.3 Fair-compare conclusion

**Protocol knobs for a fair A/B are mostly correct; outcome metrics from the live pilot are not a fair FO-vs-DP ranking.**  
DP “wins” packaging completeness; FO loses the election pipeline. Neither arm samples natives under this budget. Do **not** promote pilot S1/BCR rates as clustering science.

---

## 8. Aggregator — `scripts/aggregate_dpfo_pilot.py` (state at `e18cb6b8c`)

### 8.1 Intended metrics (good design)

- **S1:** elected RMSD ≤ 2.0 Å (and not seed echo when column present).
- **BCR:** best cluster RMSD ≤ 2.0 Å.
- **election_gap:** BCR ∧ ¬S1.
- **predictive_power:** P(S1 \| BCR) at arm level.
- Paired FO vs DP table + markdown/JSON summary.

### 8.2 Bugs at this SHA (confirmed; fixed in `c93b866a4`)

| Bug | Effect |
|-----|--------|
| S1 keys omit **`rmsd_to_crystal`** (DatasetRunner primary column) | S1 often always false/None on claim-schema CSVs. |
| `_f` accepts **any finite float**, including **−1** sentinels | `best_cluster_rmsd=-1` ⇒ BCR would be **true** (−1 ≤ 2). Inflates BCR / corrupts P(S1\|BCR). |
| If BCR column missing, **`rmsd_bcr = rmsd_s1`** | Forces P(S1\|BCR) → 1 whenever S1 hits; predictive power becomes tautology. |
| Prefers generic `"rmsd"` before elected-specific names | Schema-fragile. |
| No use of **`success_rmsd`** flag when present | Diverges from DatasetRunner claim contract. |
| `arm_dir.rglob("result.csv")` | CloudDocs tree walk — hang risk. |
| “Keep best S1 rmsd” merge across duplicate rows | OK intent; dangerous if mixed with summary CSVs that lack per-target BCR. |

Live `DPFO_SUMMARY.*` (FO BCR false despite −1 sentinels) indicates aggregation was re-run **after** `c93b866a4` (negative RMSD rejected). **Do not trust summaries produced by the pure `e18cb6b8c` aggregator against DatasetRunner CSVs.**

### 8.3 Missing science hygiene

- No PoseBusters / S2 in headline (campaign claims “success-rate”; S2 is mandatory for benchmark success elsewhere).
- No statistical uncertainty.
- No check that `elected_pose_path` is non-empty before counting a “failure.”
- `best_score` / CF never reported (follow-up analysis notes `_f` also drops legitimate **negative CF** when requiring `v >= 0`).

---

## 9. Relation to C0 claim and storage policy

| Item | Pilot (`e18cb6b8c`) | C0 claim clean |
|------|---------------------|----------------|
| OUT | `DPFO_pilot8_…` under iCloud results | Local-first + thin iCloud mirror |
| pop/gen | 200 / 50 fixed | 1000 / 2000 + pop×DoF |
| EVAL_SCALE | −1 | 1 |
| Restarts | 2 | 5 |
| Seed crystal | OFF | OFF |
| Purpose | DP vs FO plumbing | Publishable Astex claim path |

Pilot correctly **avoids contaminating claim OUT**. It does **not** follow local-first storage; ops risk is hang, not claim pollution.

---

## 10. Findings severity table

| ID | Severity | Finding |
|----|----------|---------|
| F1 | **Critical** (science use) | Budget too small; 0/8 S1 & BCR both arms — not a FO/DP ranking. |
| F2 | **Critical** (metrics at this SHA) | Aggregator mishandles DatasetRunner RMSD columns/sentinels; fixed next commit. |
| F3 | **High** | Live FO arm packaging null vs DP elected — unfair success-rate compare. |
| F4 | **High** | FO multi-MinPts ladder on binary used; contradicts single literature MinPts protocol. |
| F5 | **Medium** | Manifest/sample not in git; absolute CloudDocs paths. |
| F6 | **Medium** | Pilot OUT on CloudDocs + `rglob`; anti-hang policy violation. |
| F7 | **Medium** | `BUDGET_SCALE` still ON under “fixed” budget. |
| F8 | **Low** | Shared queue binaries; yield depends on coarse `pgrep`. |
| F9 | **Positive** | Strong C0 OUT isolation + refuse guard + yield-to-claim default. |
| F10 | **Positive** | Matched FO/DP GA seeds, mode, pop/gen, binary pin on live run. |

---

## 11. Recommendations (for future commits — not applied here)

1. Treat `e18cb6b8c` pilot results as **plumbing only**; never cite rates vs C0 or as FO>DP/DP>FO science.
2. Always run aggregator at **`c93b866a4` or later** for this campaign schema.
3. Before any FO-vs-DP metric claim: verify FO dual-suffix election + non-sentinel RMSDs; enforce **single MinPts** FO (`6ec671a92` lineage).
4. Prefer **recluster-only** design for pure clustering fairness (one GA ensemble → FO and DP heads).
5. Scale budget toward claim-like pop/gen (or intermediate pilot e.g. pop≥1000, gen≥500, restarts≥3) before interpreting S1/BCR or P(S1\|BCR).
6. Version the 8-target manifest (or generator + seed) in-repo under `benchmarks/` without secrets.
7. Local-first pilot OUT + pin-cache for aggregation; ban `rglob` on CloudDocs (use known `*/result.csv` layout or `icloud_safe_io`).
8. Aggregator: fail-closed if `elected_pose_path` empty; report packaging_fail vs true S1 miss separately.
9. Optionally set `FLEXAIDDS_BUDGET_SCALE=0` for true iso-budget small sim.

---

## 12. Verdict

**`e18cb6b8c` achieves its engineering goal:** an isolated DP/FO pilot harness that **does not write into C0 claim OUT**, with explicit yield-to-claim behavior and provenance pins.

**It does not achieve a scientifically valid FO vs DP comparison** at pop=200 gen=50 on 8 targets: budget is ~two orders of magnitude below claim, live FO packaging failed while DP did not, and the aggregator as committed cannot be trusted on DatasetRunner `result.csv` without `c93b866a4`.

| Use | Allowed? |
|-----|----------|
| Plumbing / isolation / ops dry-run of dual-arm runner | **Yes** |
| Packaging smoke (FO dual-suffix election gaps) | **Yes** (ironically valuable) |
| Claim-adjacent success rates or clustering superiority | **No** |
| P(S1\|BCR) predictive-power science | **No** (null base rate + packaging confound) |

**Overall: CONDITIONAL PASS (isolation harness) / FAIL (scientific FO vs DP pilot as stated in purpose metrics).**

---

## 13. Evidence anchors

| Artifact | Path / ref |
|----------|------------|
| Commit | `e18cb6b8cf85d6bf5851892128d595d6df3f37ef` |
| Files | `scripts/run_DPFO_pilot8_small.sh`, `scripts/aggregate_dpfo_pilot.py` |
| Follow-up fix | `c93b866a4aed4bfcacf9376238ef9b9dda076d92` |
| Live pilot OUT | `…/results/campaigns/DPFO_pilot8_small_g50_p200_20260715/` |
| Manifest | `…/queues/three_engine_entropy_q1/inputs/astex_dpfo_pilot8_random.json` (seed 20260715) |
| Claim OUT (must stay untouched) | `…/campaigns/C0_full85_claim_g2000_popmod_20260715` |
| Engine seed/budget | `LIB/DatasetRunner.cpp` (`deterministic_ga_seed`, EVAL_SCALE modes) |
| FO dual-suffix enum | `LIB/DatasetRunner.cpp` ~1072 |
| FO MinPts protocol | `docs/implementation/3dsig_red_pair_protocol.md` §2.1 |

---

*End of audit report for swarm entry #11 (`e18cb6b8c`).*
