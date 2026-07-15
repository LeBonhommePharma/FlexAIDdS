# DEPRECATED — nested Astex copy

**This directory is not the canonical Astex Diverse structure tree.**

| | |
|--|--|
| **Canonical** | `benchmarks/astex_diverse/astex_diverse/` |
| **Docs** | `benchmarks/datasets/CANONICAL.md` |
| **Checksums** | `benchmarks/datasets/astex_diverse_sha256.csv` |

## Why this still exists

Historical nested prep (`data/astex_diverse/<PDB>/…`) and a few tier-1 loose
PDBs (`1sq5.pdb`, …). Content **diverges** from the canonical tree for
receptors/apo files; ligands often match. Many entries lack oracle
`_binding_site.pdb` files required by modern LOCCLF / oracle-mode runs.

## Rules

- Do **not** point new launchers, agents, or docs at this path.
- Do **not** grow this tree with new structure preps.
- Safe removal is a separate, explicit PR after script/CI grep is clean.

See `benchmarks/datasets/CANONICAL.md` for the full audit map (including
`structures/` symlinks and gitignored `astex_repro` campaign residue).
