# Scoring-locked baseline provenance receipt

Tracks the reconstructed shared dock configs and a **legacy production-runtime**
CF receipt for the three scoring-locked targets. This is a provenance repair and
a legacy-surface reconstruction — NOT a canonical baseline, NOT a historical
root-cause claim.

## Receptor identity — protocol tier (Codex re-review at 5af6d5f8)

The receptors that reproduce the legacy campaign numbers are the **deprecated
historical second-prep** copies, which retain crystallographic waters (+metals):

```
snapshot committed here:  ops/gates/configs/legacy_runtime_receptors/{PDB}_apo.pdb
true tracked origin:      benchmarks/astex_diverse/data/astex_diverse/{PDB}/{PDB}_apo.pdb   (DEPRECATED)
byte-identical to:        ~/.flexaidds/benchmarks/astex_diverse/{PDB}/{PDB}_apo.pdb   (runtime cache)
```

Per `benchmarks/datasets/CANONICAL.md` and `benchmarks/astex_diverse/README.md`,
the **repository-canonical** receptor tree is
`benchmarks/astex_diverse/astex_diverse/{PDB}/` (waters stripped) — a *different
protocol* that gives different CF. The water retention is why the numbers move:

| PDB | legacy prod-runtime (deprecated prep, waters) | repo-canonical (waters stripped) | prior workorder |
|---|---|---|---|
| 1OQ5 | -34.229803 | -34.229803 (coincidentally equal) | -34.230 |
| 1SQ5 | -73.413683 | -73.790961 | -73.414 |
| 1YGC | -0.871396  | +7.303774  | -0.871 |

The legacy numbers match the prior campaign workorder — i.e. the campaign scored
against the deprecated water-retaining prep. That is the surface this receipt
pins. It is a **separate protocol tier** from the repo-canonical Astex prep, and
must not be called "the canonical baseline."

## Scope — what this establishes and what it does NOT

- **Establishes:** given the exact `{binary, probe_cf, config, receptor, pose}`
  below, native `cf_total` on the legacy prod-runtime prep is reproducible and
  deterministic (1OQ5 3x, 1SQ5 3x, 1YGC 2x — bit-identical, tolerance 0.0).
- **Regenerability (narrowed):** the *inputs* are now tracked (configs + receptor
  snapshots + pose hashes). The *numeric baseline* is replayable only while a
  binary matching the recorded SHA is retained — with `source_commit: null` and
  no committed binary, it is NOT source-regenerable from tracked state. KNOWN GAP.
- **Does NOT establish:** a canonical baseline (wrong prep), nor the cause of the
  historical frozen->current drift (older binary vs deleted config, unrecoverable).

## Environment (this box, current state)

| Artifact | SHA256 | Note |
|---|---|---|
| `build/FlexAIDdS` | `968377d1fbd59948896c2199886d940afb25097a5b1f405862153b304bfd7c14` | built Jul 28 01:56; source commit NOT captured (KNOWN GAP) |
| `build/probe_cf`  | `90b3ebf25a0ad2ba12a5fcfa15aa34a8b83b4aeccc59d5cfc9b077063eddcaa7` | built Jul 28 01:55; source commit NOT captured |

Working tree at receipt time is `main` `11ce273c`, but the binary predates it —
do not assume it was built from `11ce273c`.

All three configs byte-identical (shared template = tracked `1J3J_dock_config.json`):
`2f75f024952806990ae683eb35afe2e29ee5def08c2177b305ef7ebe0f39b713`

## Per-target — exact paths, hashes, expected CF, repetitions

Command (per target), paths repo-relative:
```
build/probe_cf \
  --receptor ops/gates/configs/legacy_runtime_receptors/<PDB>_apo.pdb \
  --pose     benchmarks/astex_diverse/astex_diverse/<PDB>/<PDB>_ligand.sdf \
  --ligand   benchmarks/astex_diverse/astex_diverse/<PDB>/<PDB>_ligand.sdf \
  --config   ops/gates/configs/<PDB>_dock_config.json
```

| PDB | legacy-prep receptor SHA256 | pose/ligand SHA256 | expected cf_total | reps |
|---|---|---|---|---|
| 1OQ5 | `8452399c2354bf0fb380b9e04f42791e7ca2f2398751cb3bac407bdd982eaeb6` | `c064f28cfa5d0d289e952e6bdb55b547a9bc3c6b240e5942a7f6315fb51f0481` | `-34.229803` | 3 |
| 1SQ5 | `9a8d26f19cd00e3410ebc7f3a8baeadea311008afc3d016843b957673cd41312` | `6e074aa266d0fad08b829807170aab7943fad99d76f5b5321f5f7759f4c0186d` | `-73.413683` | 3 |
| 1YGC | `cae9d69093cb684a8ea60417cfaf4553e349fd063bd47e3ff5c1fa69d427a721` | `6f92986cbc23a42fe9cf50d1878b850eb30198d93dd4a335e203f58707b00b05` | `-0.871396`  | 2 |

Tolerance: exact (0.0) — `probe_cf` emits six decimals, bit-reproducible on identical inputs.
