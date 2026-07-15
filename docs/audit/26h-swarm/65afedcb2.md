# Deep code audit — `65afedcb2`

| Field | Value |
|-------|-------|
| **Short** | `65afedcb2` |
| **Full SHA** | `65afedcb2488ff76ae83671b4dadf8c07731aa1b` |
| **Subject** | Add: Astex apo ligand-strip validation script and report CSV |
| **Author / dates** | LP · AuthorDate 2026-07-15 00:34:06 −0400 · CommitDate 00:42:02 −0400 |
| **Parents** | `342d6650d` (no prior `validate_astex_apo_strip.py`) |
| **Successor (strict gate)** | `fe0b961e49fd1f51784064aaae0b439aba4a6af7` (`fe0b961e4`) — ~3 minutes later |
| **Scope of this audit** | Early apo-strip introduction only; compare to strict-gate successor. **No source edits.** |
| **Audit date** | 2026-07-15 |
| **Overall verdict** | **ACCEPTABLE v0 — science outcome correct; operational gate incomplete** |
| **Risk rating** | **Medium** (destructive `--fix` uncapped; mission contract incomplete) |
| **Science risk** | **Low** on committed Astex 85 snapshot (0 residual cognate ligand) |

---

## 1. Executive summary

This commit introduces the **first** Astex Diverse apo residual-ligand validator:

1. `scripts/validate_astex_apo_strip.py` (471 lines, new)
2. `benchmarks/datasets/astex_apo_strip_report.csv` (85 targets + header)

It correctly answers the science question raised by prior audit work: **most `*_apo.pdb` files are byte-identical to the deposit `*.pdb` (83/85), yet that is not automatically a fail** when the cognate ligand already lives only in CIF/SDF and is absent from the deposit PDB. Residue-name matching + tight SDF coordinate matching (0.35 Å) report **`status=ok` for all 85**, with **0 ligand residue atoms** and **0 coordinate hits**.

About **three minutes later**, `fe0b961e4` hardens the same tool into the production-style **strict gate** (summary JSON, mission column aliases, MOL2 fallback, dry-run, pilot-capped write, `pass`/`pass_fail` vocabulary, CANONICAL.md documentation). That successor is the operational contract; **this commit is the prototype**.

**Bottom line:** ship-quality *science finding*; pre-production *tooling gate*. Do not treat `65afedcb2` alone as CI-ready or as a safe bulk-rewrite path.

---

## 2. Change inventory

| Path | Δ | Role |
|------|---|------|
| `scripts/validate_astex_apo_strip.py` | +471 | Validator + optional in-place strip |
| `benchmarks/datasets/astex_apo_strip_report.csv` | +86 | Frozen per-target report (85 Astex IDs) |

No C++/CMake/engine changes. No tests. No docs in this commit (`CANONICAL.md` lands in `fe0b961e4`).

---

## 3. Architecture of the early validator

### 3.1 Inputs (canonical tree)

```
benchmarks/astex_diverse/astex_diverse/<PDB>/
  <PDB>.pdb           # deposit
  <PDB>_apo.pdb       # docking receptor
  <PDB>_ligand.sdf    # cognate ligand (title = expected residue code)
```

Discovery: directories that contain `*_apo.pdb` **or** deposit `*.pdb`.

### 3.2 Detection channels

| Channel | Mechanism | Fail? |
|---------|-----------|-------|
| Residue name | SDF line-1 title vs ATOM/HETATM `resName` in apo | **fail** if any match (non-peptide) |
| Coordinates | Each apo atom within **0.35 Å** of any SDF atom | **fail** if hits ≥ `max(3, n_lig//2)` |
| Weak coords | 0 < hits < strong threshold | **warn** |
| Peptide titles | Title ∈ standard AA set → **skip** resname match; rely on coords | **warn** only if no coords available |
| Byte identity | `filecmp.cmp(apo, deposit, shallow=False)` + dual SHA-256 columns | **Informational only** — identity alone is **not** fail |
| Nonstd HETATM inventory | HETATM resnames ∉ `STD_RESIDUES`, excluding cognate title | Inventory only (does not drive status) |

### 3.3 Status machine (early)

```
missing apo or ligand SDF     → warn
ligand_atoms_in_apo > 0       → fail
strong coord match            → fail  (+ note strong_coord_match)
weak coord match              → warn
peptide_like && no lig coords → warn  (+ peptide_needs_manual_review)
else                          → ok
  (+ identity_ok_if_ligand_absent_from_pdb when identical)
```

Exit codes: `0` clean, `1` any `fail` (unless `--allow-fail`), `2` path/usage errors. **Default already fails closed on residual ligand** — there is no separate `--strict` flag yet, but behavior is “strict-ish.”

### 3.4 Optional rewrite (`--fix`)

- Strips every ATOM/HETATM whose residue name equals the SDF title.
- Writes `.bak` once if missing.
- Skips peptide-like titles unless `--force-peptide`.
- Re-analyzes after strip.
- **No dry-run, no pilot cap, no plan/preview phase.**

---

## 4. Committed report results (science)

Recomputed from the CSV at this commit and re-validated live against the current canonical tree with the early logic:

| Metric | Value |
|--------|-------|
| Targets | **85** |
| `status=ok` | **85** |
| `status=warn` / `fail` | **0 / 0** |
| Apo byte-identical to deposit | **83/85** |
| Ligand residue atoms remaining | **0/85** |
| Coord match atoms > 0 | **0/85** |
| Missing apo / ligand SDF | **0** |
| Exit code (default) | **0** |

### 4.1 Non-identical apo vs deposit (not fails)

| PDB | Ligand title | Size note | Why still `ok` |
|-----|--------------|-----------|----------------|
| **1TW6** | `ALA` (peptide-like) | apo 124 501 vs dep 126 688 | Deposit has a few more polymer atoms; peptide resname matching disabled; **0** coord hits → no residual ALA ligand evidence |
| **2BYS** | `LOB` | apo 270 220 vs dep **1 351 246** | Large historical chain trim / multi-chain prep; LOB absent by name and coords |

### 4.2 Peptide-like titles

| PDB | Title | Identical? | Status |
|-----|-------|------------|--------|
| 1TW6 | ALA | No | ok (+ `peptide_like_ligand_title`) |
| 1X8X | TYR | Yes | ok (+ peptide + identity notes) |

Only these two titles collide with `PEPTIDE_LIKE` / `STD_RESIDUES` amino acids. Correctly **not** mass-matching protein TYR/ALA.

### 4.3 Nonstd HETATM noise (not cognate)

| PDB | nonstd | Notes |
|-----|--------|-------|
| 1N1M | `HG:4` | Mercury; not cognate A3M |
| 1YV3 | `VO4:1` | Vanadate; not cognate BIT |

Inventory is useful; status correctly ignores them for the cognate-ligand gate.

### 4.4 Interpretation (correct)

For ~98% of Astex Diverse, “apo” is **not** a stripped rewrite of a holo PDB — it is the **deposit PDB already lacking cognate HETATM**, with the ligand extracted to SDF (often from CIF). Byte-identity + pass is therefore expected and must not trigger bulk rewrites. The early script documents this in notes (`identity_ok_if_ligand_absent_from_pdb`) and in the module docstring. **Science interpretation: sound.**

---

## 5. Comparison to `fe0b961e4` strict gate

`fe0b961e4` rewrites the same script (~471 → ~788 lines) and expands artifacts. Both commits agree on the **science outcome** for the live tree: **fail=0/85**, **83/85 identical**, non-identical **1TW6 + 2BYS**.

### 5.1 Feature matrix

| Capability | `65afedcb2` (early) | `fe0b961e4` (strict gate) |
|------------|---------------------|---------------------------|
| Residue-name residual detect | Yes | Yes (same core) |
| SDF coord residual detect (0.35 Å) | Yes | Yes (returns hit atoms too) |
| Peptide-safe resname skip | Yes | Yes |
| Identity ≠ automatic fail | Yes | Yes (documented in CANONICAL.md) |
| Status vocabulary | `ok` / `warn` / `fail` | `pass` / `warn` / `fail` (+ legacy `ok` counted as pass) |
| Binary `pass_fail` column | **No** | **Yes** (`warn` → pass for binary gate) |
| Mission alias `identical_sha_to_deposit` | **No** (only `apo_identical_to_deposit`) | **Yes** (both) |
| `ligand_hetatm_names` | **No** | **Yes** |
| `ligand_source` (sdf/mol2/none) | **No** (SDF only) | **Yes** |
| MOL2 fallback | **No** | **Yes** (`resolve_ligand`) |
| Summary JSON | **No** | **Yes** (`astex_apo_strip_summary.json`) |
| Explicit `--strict` | **No** (default exit 1 on fail) | **Yes** (documented mission gate; same fail exit) |
| `--fix-dry-run` | **No** | **Yes** |
| Destructive write | `--fix` unrestricted | `--write` / `--fix` **pilot-capped ≤3** unless `--all-safe` |
| `plan_strip` / `fix_planned` | **No** | **Yes** |
| Safe `relative_to` helper | **No** (can raise) | **Yes** (`_rel`) |
| CANONICAL.md documentation | **No** | **Yes** |
| Ranked `worst_offenders` | **No** | **Yes** (in JSON) |
| Unit / CI tests | **No** | **No** (still none in successor) |

### 5.2 Gate semantics

| Aspect | Early | Strict |
|--------|-------|--------|
| Residual ligand → exit 1 | **Yes (default)** | **Yes (default + `--strict`)** |
| `--allow-fail` override | Yes | Yes |
| Machine-readable summary for ops/CI | CSV only | CSV + **JSON** |
| Safe rewrite workflow | Direct `--fix` | dry-run → pilot write → optional `--all-safe` |
| Status rename `ok`→`pass` | — | Breaking for naive CSV greps on `status=ok` |

**Important nuance:** Early already exits nonzero on `fail`. The successor’s `--strict` is mostly an **explicit mission flag** and documentation surface, not a behavioral invention. The real upgrades are **contract fields, dry-run, write caps, MOL2, and JSON**.

### 5.3 Detection logic delta (core)

Residue / peptide / strong-coord thresholds are **essentially unchanged**. Successor improvements around detection:

- MOL2 path if SDF missing.
- Coord hits also populate `ligand_hetatm_names` when resname channel is empty.
- `pass_fail` collapses warn→pass for binary dashboards.

No evidence the successor changes the Astex 85 pass/fail outcome relative to early.

---

## 6. Findings (ordered by severity)

### F1 — HIGH (operational): Uncapped destructive `--fix`

**Where:** `main()` fix loop; `strip_ligand_residue()`.

**Issue:** Early `--fix` rewrites **every** `status=fail` target in one invocation (all chains, all instances of the residue name), with only a single `.bak` guard and peptide skip. There is no:

- dry-run,
- pilot limit,
- confirmation that only cognate instances are removed,
- regeneration hook for `astex_diverse_sha256.csv` / manifest.

**Impact:** If the validator ever false-fails (wrong SDF title, homonym HET code, mis-labeled peptide forced with `--force-peptide`), a single command can mutate the **canonical docking receptors** at scale.

**Mitigation in successor:** `--fix-dry-run`, `--write` pilot cap 3, `--all-safe`, `plan_strip`.

**Audit recommendation:** Treat early `--fix` as **unsafe for production ops**; prefer successor workflow only.

---

### F2 — MEDIUM: Mission / CI contract incomplete

**Missing vs stated science-priority gate (fulfilled in `fe0b961e4`):**

- `identical_sha_to_deposit`, `ligand_hetatm_names`, `pass_fail`
- `astex_apo_strip_summary.json`
- explicit `--strict`
- CANONICAL.md operator docs

**Impact:** Consumers grepping for mission field names or relying on JSON ops monitors cannot use this commit alone. Status token `ok` later becomes `pass` → brittle automation.

---

### F3 — MEDIUM: SDF-only ligand resolution

**Where:** `ligand = d / f"{pdb_id}_ligand.sdf"`; `parse_sdf_coords` V2000-only (`lines[3][0:3]`).

**Issues:**

1. No MOL2 fallback if SDF absent or corrupt.
2. V3000 SDF would silently yield **0 coords** → false `ok` if resname also mismatches.
3. Only first molecule; multi-mol SDFs ignored.

**Astex 85 reality:** All sampled titles are 3-letter codes; counts lines are V2000; all 85 ligands exist as SDF in the committed report. **Latent**, not active on this set.

**Successor:** MOL2 fallback; still V2000-oriented for SDF.

---

### F4 — MEDIUM: Cofactor-blind nonstd inventory

**Where:** `STD_RESIDUES` includes `HEM`, `NAD`, `ATP`, `GOL`, `EDO`, `PEG`, etc.

**Issue:** Nonstd inventory **suppresses** common cofactors/solvents. Residual **non-cognate** cofactors do not fail the gate (by design for apo docking). Residual **cognate** ATP/HEM would still be caught **if** the SDF title matches that code.

**Risk:** If a dataset entry’s cognate ligand is a cofactor-like code but the SDF title is wrong/empty, residual cofactor atoms are invisible to both resname and nonstd channels; only coord match remains.

**Astex 85:** No cognate titles in the cofactor subset of `STD_RESIDUES` (only peptide ALA/TYR collide with the AA subset).

---

### F5 — LOW–MEDIUM: `relative_to(REPO_ROOT)` can crash external trees

**Where:** `analyze_target` path fields.

**Issue:** `str(apo.relative_to(REPO_ROOT))` raises `ValueError` if `--tree` points outside the repo, aborting the whole run mid-target.

**Successor:** `_rel()` swallows `ValueError`.

---

### F6 — LOW: Coordinate matcher complexity / false-positive profile

**Where:** `count_coord_matches` — nested loop O(|PDB| × |lig|).

- Fine for Astex 85; not general.
- Tolerance **0.35 Å** is intentionally tight (near-identity), so protein atoms near the empty pocket almost never false-hit. **Good scientific choice.**
- Strong threshold `max(3, n_lig//2)` requires substantial overlap before fail — resists single-atom noise.

---

### F7 — LOW: Strip semantics are name-global

`strip_ligand_residue` removes **all** residues with the cognate code (every chain/instance). Correct for multi-copy cognate ligands; dangerous for rare homonym HET codes used both as ligand and as covalent modification elsewhere. No occupancy/altLoc awareness; no CONECT cleanup.

---

### F8 — LOW: No automated tests; CSV is a frozen snapshot

- Commit adds no pytest/ctest.
- Committed CSV can drift from tree until re-run.
- Running the early script defaults to **overwriting** the report path (observed during this audit; restored to HEAD). Operators should use a temp `--report` for experimental runs.

---

### F9 — INFO: Dead code / polish

`write_csv` assigns `fieldnames` from `asdict(rows[0])` then immediately overwrites from dataclass field order. Harmless.

---

### F10 — INFO: REPO_ROOT coupled to `__file__`

Correct when script lives at `scripts/validate_astex_apo_strip.py`. Copying the file outside the repo breaks default paths (reproduced with `/tmp` copy during audit). Expected for repo-root-relative tools.

---

## 7. What the early commit got right

1. **Correct scientific framing** of apo ≡ deposit for CIF-extracted ligands (83/85).
2. **Peptide-safe** resname handling (1TW6 ALA, 1X8X TYR) — avoids shredding the protein with a naive strip.
3. **Dual evidence**: resname + tight coords; identity is diagnostic only.
4. **Per-target SHA-256** of apo and deposit for audit trails / checksum cross-checks.
5. **Fail-closed default** on residual cognate ligand (`exit 1` without needing a special flag).
6. **Conservative default mode** is report-only; no bulk rewrite unless `--fix`.
7. **Frozen CSV evidence** that the full 85 was actually scanned at commit time.
8. Clean stdlib-only Python 3 — no new deps, Apache-safe.

---

## 8. Relationship to dataset / docking pipeline

Canonical docking inputs use `*_apo.pdb` (see `generate_flexaid_inp.py`, CANONICAL.md after successor). Residual cognate ligand in apo would:

- bias Voronoi CF / contact scoring (self-docking cheats),
- distort cavity/site perception,
- invalidate self-docking RMSD claims.

This validator is therefore a **data-quality gate for the three-engine Astex campaign**, not cosmetic hygiene. Early commit establishes the measurement; strict gate makes it operable.

**Self-docking vs cross-docking:** script does not reason about docking mode; it only checks cognate residual atoms in the receptor file used as apo input. Semantics remain “is the cognate ligand gone from the receptor structure?”

---

## 9. Threat model / misuse scenarios

| Scenario | Early behavior | Risk |
|----------|----------------|------|
| Report-only on Astex 85 | exit 0, CSV ok×85 | Safe |
| Residual ligand appears later | exit 1 | Correct fail-closed |
| Operator runs bare `--fix` after a false fail | Mass rewrite of apos | **High** |
| `--force-peptide` on 1TW6/1X8X | Could delete protein ALA/TYR if status were fail | High if forced; currently status ok so no strip |
| External `--tree` | Possible crash on `relative_to` | Medium UX |
| V3000 / missing SDF | Silent under-detect | Medium latent |

---

## 10. Verification performed this audit

| Check | Result |
|-------|--------|
| `git show 65afedcb2 --stat` | 2 files, +557 lines |
| Full script review (`/tmp` export of blob) | Complete |
| CSV status distribution | ok=85, fail=0, warn=0 |
| Identity split | 83 True / 2 False (1TW6, 2BYS) |
| Live re-run of early logic on current tree | ok=85, exit 0 |
| Diff vs `fe0b961e4` script | Documented in §5 |
| `fe0b961e4` summary JSON | fail=0/85, same non-identical pair |
| Non-identical size/het diffs | 1TW6 small polymer delta; 2BYS large chain trim |
| Ligand titles ∩ STD_RESIDUES | Only ALA, TYR (peptide path) |
| No engine/CMake touch | Confirmed |

**Note:** A live re-run of the early script briefly overwrote the working-tree CSV with the early schema; content was restored to `HEAD` (strict schema, 27 columns) before finishing this report.

---

## 11. Verdict matrix

| Dimension | Score | Notes |
|-----------|-------|-------|
| Scientific correctness (Astex 85) | **Strong** | 0 residual cognate; identity interpreted correctly |
| Detection design | **Good** | Peptide-safe; tight coords; dual channel |
| Operational safety of `--fix` | **Weak** | Uncapped rewrite; fixed in successor |
| Mission / CI contract | **Incomplete** | JSON, aliases, `--strict`, docs → `fe0b961e4` |
| Test coverage | **Absent** | Both early and strict lack automated tests |
| Docs | **Absent here** | Successor adds CANONICAL.md section |
| License / hygiene | **Clean** | New script + CSV only; no secrets/paths |

### Final verdict

**`65afedcb2` is a correct and valuable first measurement of Astex apo ligand residuals.**  
It should be read as **v0 of the apo-strip gate**: science findings trustworthy; rewrite path and CI contract **not** yet production-grade.

**Prefer `fe0b961e4` (and any later refinements) as the operational reference.** When citing dataset cleanliness for campaign claims, quote:

- fail = **0/85** residual cognate ligand  
- identical apo/deposit = **83/85** (non-identical: **1TW6**, **2BYS**, not fails)  
- validator lineage: **65afedcb2 → fe0b961e4**

### Residual recommendations (for future work; not applied here)

1. Keep dry-run + pilot write as the only rewrite path (already in successor).
2. Add a small pytest using synthetic mini PDBs (residual resname, peptide title, strong coords, missing SDF).
3. Optionally wire `python3 scripts/validate_astex_apo_strip.py --strict` into CI when the dataset tree is present (or as a scheduled dataset job).
4. When rewriting apos, regenerate `astex_diverse_sha256.csv` / manifest in the same change.

---

## 12. File / symbol index (early commit)

| Symbol | Purpose |
|--------|---------|
| `STD_RESIDUES` | Polymer / solvent / common cofactor keep-list for nonstd inventory |
| `PEPTIDE_LIKE` | AA codes that disable resname matching / default strip |
| `TargetReport` | CSV row dataclass |
| `sdf_title` / `parse_sdf_coords` | Ligand identity + V2000 coords |
| `parse_pdb_atoms` | ATOM/HETATM with xyz + res/chain |
| `count_coord_matches` | 0.35 Å residual geometry |
| `analyze_target` | Status machine per PDB ID |
| `strip_ligand_residue` | Destructive name-based strip + `.bak` |
| `discover_targets` | Tree walk |
| `write_csv` / `main` | CLI + report + exit policy |

---

## 13. One-line summary for INDEX / swarm rollup

> **65afedcb2** — First Astex apo residual-ligand validator + CSV (85/85 ok, 83/85 apo≡deposit); science sound; **unsafe uncapped `--fix`** and incomplete mission gate vs successor **fe0b961e4**.
