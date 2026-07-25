---
name: flexaidds-dataset-runner
description: >
  Production DatasetRunner + classic FlexAID three-engine campaign skill for FlexAIDdS.
  Use for Astex/HAP2/CASF DatasetRunner runs, entropy collapse analysis, and A/B0/B
  red-pair prep. Defers scientific policy to repo AGENTS.md and .grok/skills/flexaidds/SKILL.md.
  Natural triggers: "run DatasetRunner", "benchmark on Astex", "entropy collapse",
  "three-engine red-pair", "FlexAIDdS dataset runner".
user_invocable: true
---

# FlexAIDdS DatasetRunner Skill (in-repo launcher)

**Source of truth:** `AGENTS.md`  
**Primary science skill:** `.grok/skills/flexaidds/SKILL.md`  
**Benchmark ops:** `.agents/skills/flexaidds-benchmarking/SKILL.md`  
**Methodology:** `METHODOLOGY.md` (cite §N; do not restate numbers)

This skill is a **launcher + packaging pointer**. If anything here conflicts with the repo flexaidds skill or `AGENTS.md`, **those win**.

```bash
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null)}"
cd "$FLEXAIDDS_ROOT"
```

Never hardcode `/Users/<username>/...` in committed commands. Use `FLEXAIDDS_ROOT`, `FLEXAIDDS_LOCAL_ROOT`, `FLEXAIDDS_ICLOUD`, `FLEXAIDDS_QUEUE_ROOT`.

---

## Paths

| Path | When |
|------|------|
| **DatasetRunner / C0** | `python3 .grok/skills/flexaidds/scripts/dataset_runner.py` / `python -m flexaidds.dataset_runner` |
| **Classic red-pair A/B0/B** | `scripts/generate_flexaid_inp.py` + `scripts/run_flexaid_arm_pilot8.sh` |
| **Live Astex entropy ops** | `.agents/skills/flexaidds-benchmarking/SKILL.md` |

Do **not** mix claim rates across paths without labeling engine and flags.

---

## Non-negotiable science (summary)

1. **CF is a scoring proxy**, not experimental ΔG.
2. **Softβ S1** defaults **OFF** (`FLEXAIDDS_SOFTBETA_ELECTION=0`). Softβ cannot fix BCR=0.
3. **Arm B FO@TEMPER21 ≠ Softβ S1.**
4. Production `ga.inp`: `SHARESCL 10`, `SHAREPEK 5`, `SHAREALF 4` (never 0.20).
5. Matrix pin: `MC_st0r5.2_6.dat` MD5 `72d7c7396702331d96ff12d18f831796`.
6. Ligand integrity + native CF oracle fail closed before ranking claims.
7. **No dual-launch** of heavy GA. Serial A → B0 → B.
8. **Local-first** live OUT; iCloud thin mirror after success.
9. Modern claim success = **RMSD ≤ 2.0 Å and PoseBusters** on the same pose + on-disk `result.csv`.
10. Refuse success language without receipts (see flexaidds *Deception-proof claim contract*).

Full contract: `.grok/skills/flexaidds/SKILL.md`.

---

## DatasetRunner campaigns

```bash
cd "$FLEXAIDDS_ROOT"
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
export FLEXAIDDS_REQUIRE_BUILD=1   # recommended for claim sessions
export FLEXAIDDS_SOFTBETA_ELECTION=0

# Dry-run first
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
  --dataset astex_diverse --tier 1 --dry-run --resume --package

# Real campaign only after dry-run + pin check (local results dir)
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
  --dataset astex_diverse --tier 2 --resume --package \
  --results-dir "${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/dataset_runner"
```

Confirm **self-docking vs cross-docking** before launch (`astex_diverse` / native vs `astex_nonnative` / non_native).

---

## Classic three-engine red-pair (A / B0 / B)

```bash
cd "$FLEXAIDDS_ROOT"
source scripts/use_local_first_benchmark_storage.sh 2>/dev/null || true
bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40
# Serial only — never dual-launch
bash scripts/run_3dsig_red_pair_serial.sh --only A   # then B0, then B
bash scripts/run_benchmark_ops_monitor.sh
```

Pin binary SHA256 + matrix MD5 in `RUN_RECEIPT.json`. Primary metric: **S_top10**.

---

## Resume / claim discipline

- Skip targets with valid `result.csv` unless `--force`.
- New protocol → new OUT namespace; do not mix old pilot receipts with new bins.
- Softβ remains OFF on resume unless the receipt says ON.
- Fail-closed admission metrics for incomplete `claim_ready` (see `benchmarks/protocols/admission_metrics_contract.md`).
- Cite `METHODOLOGY.md` for parity/determinism/Astex-85/ctest.

---

## Validation before claiming “skill ready”

```bash
cd "$FLEXAIDDS_ROOT"
python3 scripts/check_repo_hygiene.py
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

## Any target / any ligand

```bash
python3 .grok/skills/flexaidds/scripts/dock_any.py --receptor r.pdb --ligand l.mol2
python3 .grok/skills/flexaidds/scripts/validate_dataset_semantics.py
```
