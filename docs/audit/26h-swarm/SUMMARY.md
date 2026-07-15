# 26h Deep Code Audit — Swarm Synthesis

**Scope:** all `91` commits on `main` from `git log --since="26 hours ago" --oneline`  
**Method:** 91 parallel background audit agents (one per commit) → per-commit reports in `docs/audit/26h-swarm/<short>.md`  
**Focus:** scientific correctness, reproducibility, ranking integrity, CF vs thermo language, ops safety  

---

## Executive verdict

**Do not treat current `main` tip as claim-ready for 3Dsig red bars or clean C0 Shannon election without fixes.**

The last 26h shipped real science progress (soft-β G̃ election intent, single literature FO MinPts, DoF=pop not gen, S1/S2/S3 contract, RUN_RECEIPT, shell harden, apo strip, A/B pilot harness) **and** several high-risk product defaults / incomplete packaging paths. Unmerged branches still hold required pieces (`SoftBetaFreeEnergy.h`, FO dual-suffix complete packaging).

### Scorecard

| Metric | Count |
|--------|------:|
| Commits audited | 91 |
| CRITICAL | 6 |
| HIGH | 16 |
| MEDIUM | 19 |
| LOW/INFO/? | 50 |
| Verdict BLOCK / REJECT | 32 |
| Verdict MERGE_WITH_FIX | 38 |
| Verdict MERGE_OK / ACCEPT | 19 |
| Ranking-impact YES/SOFT | 7 |

### By category

| Category | n |
|----------|--:|
| merge | 32 |
| science | 26 |
| build | 14 |
| other | 11 |
| ops | 5 |
| docs | 3 |

---

## P0 — Science blockers (must fix before claim tables)

### 1. Soft-β S1 election default is incomplete identity (`c82e6fc24`) — **CRITICAL**

- Algebra \( \tilde G = \tilde H - T\tilde S \equiv E_{\min}-T\ln Z \) is **correct** (soft-β \(\beta=1/T\), not \(k_B T\)).
- **soft_T hardcodes 298** when `FLEXAIDDS_ELECTION_SOFT_T` unset — **never reads dock TEMPER**. Arm B uses **TEMPER 21** → DatasetRunner and engine disagree.
- Default **ON** is a ranking contract change without Astex revalidation.
- **FO BindingMode never writes `.mcf`** → only `cluster.cpp` does → FO path collapses to \(S̃=0\), \(G̃=\mathrm{CF}\) (entropy election theater).
- FO dual-suffix enumeration only on election path; BCR/oracle elsewhere incomplete on `main`.
- Docs claim `LIB/SoftBetaFreeEnergy.h` shared identity — **file absent on main** (lives on unmerged `fix/softbeta-ranking-identity`).
- Legacy ZH rollback τ defaults 298 not historical 0.592.

**Fix:** wire dock T into election; emit `.mcf` on FO path; merge SoftBeta header + FO dual-suffix packaging; Astex pilot before default-ON claims; or keep default OFF until validated.

### 2. FO dual-suffix packaging incomplete on main

- Engine emits `prefix_minPts_rank.pdb`.
- Complete shared `enumerate_emitted_cluster_heads()` is on **unmerged** `fix/fo-dual-suffix-packaging`.
- Without it, BCR/S3/election miss FO heads → null/sentinel RMSDs (seen on DPFO pilot).

### 3. 3Dsig red-pair metric pipeline not closed (`0e39f3a0b`)

- Primary deck metric is **S_top10** + 10k bootstrap.
- Pilot parser emits `rmsd_top1`/`rmsd_bcr` only → bootstrap fail-open / not true top-10.
- Threshold **&lt;2.0** vs **≤2.0** drift across docs/code.
- Serial launcher had local-first vs iCloud policy thrash; bash 3.2 `extra[@]` launch bug (**BLOCK** on `5a24ebbd2`).

### 4. Protocol freeze drift (TEMPER / gen / budget)

| Axis | Early (`4e87c0b3c`) | Later freeze | Risk |
|------|---------------------|--------------|------|
| Arm B TEMPER | **298** | **21** | soft-β sharpness ~14× |
| Generations | **6000** | **2000** | budget comparability |
| Restarts | 5 | 3Dsig wants **10** | deck fidelity |
| DoF | fixed 1000×6000 | **pop×DoF, gen fixed** | AGENTS contract |

`68063cc9d` correctly fixed **docs** that claimed iso-budget; code already grew total evals with pop.

### 5. Ranking-relevant co-landed science (`7e79352e3` / `dea70ea88` / `26cb99276` / `83f7a8584`)

- Cluster ACF: global mix → **local log-sum-exp** (emission order change when T>0).
- No-seed cognate stack + DirectLigandIC (fairness improved; torsions still cognate).
- Oracle seed fix restores native flood when file set; engine defaults still seed-friendly if orchestration fails.
- Clash promotion default OFF (good); flag ON can move emitted RMSD without reordering modes.

---

## P1 — Ops / reproducibility

| Issue | Commits | Note |
|-------|---------|------|
| iCloud force then reverse to local-first | `f2fce514f` → `b1ae633b6` → `ab0850a6d` → `5a24ebbd2` | Final doctrine correct; receipts still lie; dual doctrine residual |
| `icloud_safe_io` timeout isolation broken | `f75cdfc50` | ProcessPoolExecutor does not kill hung FileProvider workers; pickling broken |
| Ops monitor S1 inflate | `a9cb06e64` | Hungarian RMSD alone can force S1 true |
| Fleet placeholder catastrophe | `69aa0fab6` / `2fc7189d8` / `8eaf043ae` | Gutted CMake + DatasetRunner.h; recovered later (`292a0ce6a`, `d842e3247`, `033eeb889`) — **never pin claim binaries to 8eaf043ae** |
| Metal gate unrunnable | `ffa7499cd` / `8c42517bd` | No `self-hosted-m3` runner; nm smoke can LTO-strip symbols |
| Absolute `/Users/lp.more/...` paths | `9971dff7e`, docs | AGENTS hygiene violation |

---

## P2 — Good landings (keep)

| Area | Commits | Why |
|------|---------|-----|
| Admission S1/S2/S3 + fail-closed seeds | `668cc3095` → `9dbbd9fa9` | Correct claim gates after harden |
| Shell/exec harden | `646824df1` → `4dabec565` | Provenance injection resistance |
| ProtocolConfig + RUN_RECEIPT | `205ed0887` → `c957d32b2` → `04735c31b` | Typed env + receipts (dual-truth seed_elitism residual) |
| Apo strip gate | `65afedcb2` → `fe0b961e4` | 0/85 residual ligand (live re-run) |
| CF naming clarity | `033eeb889` | + restored DatasetRunner.h |
| Single FO MinPts literature | `6ec671a92` | Right production policy (with heuristic caveats) |
| Homebrew Metal link via flexaid_core | `3e059594b` | Real link fix |
| Formula sha256 v2.0.3 | `00a7b6eae` / `342d6650d` | Digest matches live tarball |
| Dependabot pins | `155ebbb7b` etc. | Low risk |
| Ligand-centered sites | `89e4979a9` | 20/20 crystal fidelity for bad GetCleft targets |

---

## Ranking-impact commits (YES/SOFT)

- `6ec671a92` (?): Fix: run FastOPTICS once with literature MinPts (drop triple ladder)
- `c4509428d` (?): Add: FlexAID A/B0/B pilot harness for three-engine campaign
- `8eaf043ae` (?): Merge pull request #257 from LeBonhommePharma/feature/bonhomme-fleet-dataset-runner-v1
- `83f7a8584` (?): Add: PB extract hygiene + optional BindingMode clash promotion
- `dea70ea88` (?): Add: no-seed cognate docking stack + DirectLigandIC geometry

---

## CRITICAL / HIGH severity list

### CRITICAL
- `5a24ebbd2` — Fix: FlexAID --legacy for A/B/B0 pilot; local-first OUT + deferred iCloud sync → BLOCK
- `c82e6fc24` — Fix: DatasetRunner elects by 3Dsig Shannon free energy G̃=H̃−T·S̃ → MERGE_WITH_FIX
- `8eaf043ae` — Merge pull request #257 from LeBonhommePharma/feature/bonhomme-fleet-dataset-runner-v1 → BLOCK
- `711b83cc9` — Merge branch 'main' into feature/bonhomme-fleet-dataset-runner-v1 → BLOCK
- `2fc7189d8` — fix(build): Update CMakeLists.txt to include FleetRunner, metal_microbench, tests, and new targets for Bonhomme Fleet DatasetRunner → BLOCK
- `69aa0fab6` — feat: Integrate Bonhomme Fleet as first-class DatasetRunner backend + full testing & monitoring + pose viz tools → BLOCK

### HIGH
- `0e39f3a0b` — Add: 3Dsig red-pair protocol, archived bars, serial A→B0→B launcher → MERGE_WITH_FIX
- `f75cdfc50` — Fix: Production CloudDocs anti-hang I/O for ops and agents → MERGE_WITH_FIX
- `b1ae633b6` — Fix: Claim live I/O on local disk to stop iCloud fileprovider hangs → MERGE_WITH_FIX
- `6ec671a92` — Fix: run FastOPTICS once with literature MinPts (drop triple ladder) → MERGE_WITH_FIX
- `a9cb06e64` — Add: unified benchmark ops + finished-run monitor automation → BLOCK
- `f2fce514f` — Docs/Ops: force all production benchmark results onto iCloud Drive → BLOCK
- `c4509428d` — Add: FlexAID A/B0/B pilot harness for three-engine campaign → MERGE_WITH_FIX
- `8c42517bd` — Fix: Harden self-hosted Metal gate + exact release language → MERGE_WITH_FIX
- `e2b799495` — Merge pull request #261 from LeBonhommePharma/fix/audit-macos-ci-gate → MERGE_WITH_FIX
- `5a3b95430` — Fix: Route Homebrew --with-metal through HEAD until flexaid_core fix ships → MERGE_WITH_FIX
- `292a0ce6a` — Merge: origin/master into main — align main with master tip → MERGE_WITH_FIX
- `964bec0a2` — Merge pull request #258 from LeBonhommePharma/feature/bonhomme-fleet-dataset-runner → MERGE_WITH_FIX
- `9971dff7e` — Add: Astex repro launch/score scripts and handoff notes → MERGE_WITH_FIX
- `dea70ea88` — Add: no-seed cognate docking stack + DirectLigandIC geometry → MERGE_WITH_FIX
- `4e87c0b3c` — Add: queue-ready three-engine entropy comparison protocol → MERGE_WITH_FIX
- `26cb99276` — Fix: oracle-ceiling native pose seed when reference_ligand.file is set → MERGE_WITH_FIX

---

## Unmerged branches still required for science identity

| Branch | Purpose |
|--------|---------|
| `fix/softbeta-ranking-identity` | Shared `SoftBetaFreeEnergy.h` + BindingMode local G̃ ≡ ACF |
| `fix/fo-dual-suffix-packaging` | Complete FO head enumeration for election **and** BCR |
| `feat/3dsig-live-red-pair-bars` | Live barplots from new docks |
| `feat/3dsig-full85-prepare` | Full85 prepare-only launcher |

---

## Recommended action order

1. **Hotfix launch scripts** — bash 3.2 `extra[@]`, sync parser, dual-launch flock (`5a24ebbd2`).
2. **Wire soft_T ← dock TEMPER** + FO `.mcf` emission + dual-suffix packaging merge.
3. **Default-OFF or feature-flag** Shannon S1 until Astex pilot validates vs CF election.
4. **S_top10 emission** from pilot arms (mode_rmsd_0..9) before 3Dsig bootstrap claims.
5. **Pin binaries/SHA** + matrix MD5 `72d7…` in every RUN_RECEIPT; forbid claim pins to fleet-stub SHAs.
6. **Fix icloud_safe_io** kill-on-timeout or ban CloudDocs reads entirely for ops.
7. **Clean dual storage doctrine** — one AGENTS path only (local-first already written; kill force-iCloud scripts).
8. **Register self-hosted-m3** or stop claiming Metal CI validation.
9. **Hygiene** — absolute paths out of benchmarks scripts; residual CF→pKd affinity language.

---

## Process note on the swarm

- Agents were instructed report-only; some still created `docs/audit-*` branches and pushed — **hygiene debt**; consolidate reports under `docs/audit/26h-swarm/` only.
- Some reports are full deep audits; early/lost ones were re-materialized from agent completion digests (still evidence-backed from that agent’s file reads).

---

## Per-commit reports

See [INDEX.md](INDEX.md) for the full table. Each commit: `docs/audit/26h-swarm/<short>.md`.
