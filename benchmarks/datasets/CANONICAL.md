# Canonical benchmark datasets (Astex and related)

This note resolves **duplicate Astex directory trees** and documents which
path agents, launchers, and `benchmark_datasets` should treat as authoritative.

**Do not delete large local trees without a replacement path.** Prefer this
document + manifests + deprecation notes until a dedicated data-dedup PR
removes tracked duplicates safely.

---

## Astex Diverse (85) — CANONICAL

| Item | Path |
|------|------|
| **Canonical structure tree** | `benchmarks/astex_diverse/astex_diverse/<PDB>/` |
| Dataset YAML (runner config) | `benchmarks/datasets/astex_diverse.yaml` |
| PDB ID list (parity with C++) | `LIB/DatasetRunner.cpp` → `astex_diverse_codes()` |
| Tier-1 harness | `benchmarks/astex_diverse/` (`manifest.yaml`, `run.sh`, `download.sh`) |
| Checksum CSV (key files) | `benchmarks/datasets/astex_diverse_sha256.csv` |
| Manifest summary (JSON) | `benchmarks/datasets/astex_diverse_manifest.json` |
| Apo ligand-strip report (CSV) | `benchmarks/datasets/astex_apo_strip_report.csv` |
| Apo ligand-strip summary (JSON) | `benchmarks/datasets/astex_apo_strip_summary.json` |
| Reproducibility prose | `REPRODUCIBILITY.md` § Dataset |

### Per-target files (canonical)

Under `benchmarks/astex_diverse/astex_diverse/<PDB>/`:

| File | Role |
|------|------|
| `<PDB>_apo.pdb` | Receptor with cognate ligand removed (docking input) |
| `<PDB>_ligand.sdf` | Crystal ligand (RMSD reference) |
| `<PDB>_binding_site.pdb` | Oracle site / LOCCLF sphere |
| `<PDB>.pdb` / `<PDB>.cif` | RCSB deposit snapshot |
| optional `*_ligand_centered_site.pdb` | Alternate site definition |

**All production launchers and REPRODUCIBILITY.md point at this tree**
(e.g. `ORACLE_DIR=…/benchmarks/astex_diverse/astex_diverse`).

### Checksums

```bash
# Verify a single key file:
shasum -a 256 benchmarks/astex_diverse/astex_diverse/1HNN/1HNN_apo.pdb
# Compare against:
#   benchmarks/datasets/astex_diverse_sha256.csv

# Full-tree rehash / verify:
python3 scripts/generate_astex_manifest.py
python3 scripts/generate_astex_manifest.py --check
```

The committed CSV covers **all 85 PDB IDs** and the key file suffixes above
(~445 hashed files). Runtime cache used by `DatasetRunner` defaults to
`~/.flexaidds/benchmarks/` and is **not** the canonical source tree; oracle
sites may be redirected with `FLEXAIDDS_ORACLE_SITE_DIR` pointing at the
canonical tree.

### Apo strip validation (ligand residual check)

Science data-quality gate: confirm cognate ligands are **absent** from
`*_apo.pdb` (residue-name + coordinate match vs `*_ligand.sdf` / optional
MOL2). Byte-identity of apo to deposit is reported but is **not** a fail when
the ligand already lives only in CIF/SDF.

```bash
# Report only (writes CSV + summary JSON):
python3 scripts/validate_astex_apo_strip.py --strict

# Describe strip operations without writing:
python3 scripts/validate_astex_apo_strip.py --fix-dry-run

# Real strip (≤3 pilot targets unless --all-safe; writes .bak):
python3 scripts/validate_astex_apo_strip.py --write --targets 1SQ5,1HP0,1G9V
```

**Findings (canonical 85, 2026-07-15):**

| Metric | Value |
|--------|-------|
| Targets | 85 |
| `status=fail` (residual cognate ligand) | **0** |
| `status=pass` | 85 |
| Apo byte-identical to deposit PDB | **83/85** |
| Ligand residue atoms remaining in apo | **0/85** |
| Strip dry-run planned | 0 (nothing to strip) |

**Non-identical apo vs deposit (not fails):**

| PDB | Ligand title | Notes |
|-----|--------------|-------|
| `1TW6` | `ALA` (peptide-like) | Apo differs from deposit; peptide ligand handled by coord match only; no residual coords in apo |
| `2BYS` | `LOB` | Apo differs from deposit (historical chain trim / prep); no residual LOB by name or coords |

**Interpretation:** For most of the set, “apo” means the deposit **PDB**
already lacks the cognate HETATM (ligand extracted into SDF from CIF). That is
why 83/85 files are SHA-identical to deposit **and** still pass the residual-
ligand gate. Do **not** bulk-rewrite apo files without re-running this
validator and regenerating checksums. Do **not** delete data trees as part of
this check.

---

## Duplicate / secondary trees (do not treat as source of truth)

| Path | Size / status (typical) | Role | Action |
|------|-------------------------|------|--------|
| `benchmarks/astex_diverse/data/astex_diverse/` | ~85 dirs, **content differs** from canonical (apo/PDB hashes diverge; ligands often match; many targets lack `_binding_site.pdb`) | Historical second prep / nested copy, still **git-tracked** | **Deprecated.** See `benchmarks/astex_diverse/data/DEPRECATION.md`. Do not use for new benchmarks. |
| `benchmarks/astex_diverse/data/{1r1h,1sq5,1t46,2c69,2hb1}.pdb` | Tier-1 loose PDBs | Legacy tier-1 download layout | Prefer full canonical dirs; leave in place for old scripts. |
| `benchmarks/astex_diverse/structures/` | 85 dirs of **symlinks** → `../astex_diverse/<PDB>/…` | Thin alias layer for tools expecting lowercase protein names | Keep; not a second structure store. |
| `benchmarks/astex_repro/` | Often **multi-GB local results** (`full/`, `full_v*`, logs); mostly **gitignored**; only launch/docs scripts tracked | Campaign workspace | Local-only results. Not dataset source. See `.gitignore`. |
| `benchmarks/astex_entropy/` | gitignored workspace | Entropy campaign scratch | Local-only. |
| `results/**/astex_*`, `.virgin_scratch/**` | Local / gitignored | Run outputs | Not input data. |
| `python/flexaidds/dataset_runner/datasets/astex_diverse.yaml` | Package copy of YAML | Must stay in sync with `benchmarks/datasets/astex_diverse.yaml` | Config mirror, not structure data. |
| `tests/benchmarks/astex_diverse/` | Tiny fixtures | Unit/integration stubs | Not the 85-set. |

### Content note on `data/astex_diverse` vs canonical

Spot checks (audit): ligand SDFs often **byte-identical**, while
`<PDB>.pdb` / `<PDB>_apo.pdb` **differ in size and SHA256**. Canonical includes
oracle `_binding_site.pdb` for all 85; the nested `data/` copy frequently does
not. Always prefer canonical.

---

## Astex Non-Native

| Item | Path |
|------|------|
| Harness / docs | `benchmarks/astex_nonnative/` (preferred spelling) |
| Cross-dock pair JSON (legacy name) | `benchmarks/astex_non_native/` (`pairs*.json`) |
| YAML | `benchmarks/datasets/astex_nonnative.yaml` |

Structures for non-native runs are derived from the **canonical diverse tree**
plus pair lists; there is no second full 85-structure payload under
`astex_nonnative/`.

---

## iCloud / machine-local results (not in git)

Benchmark **outputs** and fleet queues often live under
`$FLEXAIDDS_ICLOUD` / `$FLEXAIDDS_RESULTS` (see `scripts/use_icloud_benchmark_storage.sh`
and agent handoff docs). Those paths are **results**, never the canonical
structure source. Do not commit iCloud-mirrored result trees into the repo.

---

## Deprecation policy

1. New code and docs must reference `benchmarks/astex_diverse/astex_diverse/`.
2. Do not expand `benchmarks/astex_diverse/data/astex_diverse/` with new preps.
3. Removing the nested `data/astex_diverse` tree from git is a **follow-up PR**
   (requires confirming no CI/script hard-codes that path, then one deletion
   commit). This audit intentionally **documents without mass deletion**.
4. Regenerating checksums after intentional structure edits: update
   `astex_diverse_sha256.csv` + `astex_diverse_manifest.json` in the same PR.
