# Legacy production-runtime receptor snapshots — NOT canonical

These three `_apo.pdb` files are immutable snapshots of the **deprecated**
historical second-prep receptor tree:

    benchmarks/astex_diverse/data/astex_diverse/{PDB}/{PDB}_apo.pdb   (deprecated)

which is byte-identical (hash-for-hash) to the runtime cache
`~/.flexaidds/benchmarks/astex_diverse/{PDB}/{PDB}_apo.pdb` that the legacy
campaign actually scored against. This prep **retains crystallographic waters**
(and metals where present).

They are NOT the repository-canonical Astex receptors. Per
`benchmarks/datasets/CANONICAL.md` and `benchmarks/astex_diverse/README.md`, the
canonical apo tree is `benchmarks/astex_diverse/astex_diverse/{PDB}/` (waters
stripped) — which produces materially different CF (1YGC: +7.30 canonical vs
-0.871 legacy). These snapshots exist only to pin the exact structures that
produced the legacy production-runtime baseline recorded in
`../SCORING_LOCKED_BASELINE_RECEIPT.{md,json}`. Use the canonical tree for new work.
