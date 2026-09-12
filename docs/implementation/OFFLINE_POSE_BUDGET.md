# Explicit budgets for saved pose scoring

The frozen campaign `offline_elector.build_pool(prefixes)` always enumerated at
its default of 50, even though `enumerate_heads` already accepted a budget.
`scripts/offline_elector.py` now forwards `build_pool(prefixes, budget=250)` to
every restart. Omitting the argument retains 50. The accepted range is 1–5000;
the ambient `FLEXAIDDS_MAX_RESULTS` variable never changes an offline audit.

The successor preserves the frozen port's CF/DP rank order, followed by FO heads
ordered by `(min_pts, rank, path)`, with one shared budget per restart. Seed
`_INI` files and published pose copies remain excluded. Parsing, election,
sentinel and free-energy formulas retain their historical behavior. The module
docstring identifies the original campaign assumptions; it is not a universal
implementation of every engine protocol.

## Score a saved cell

The new adapter takes the existing trusted RMSD provider explicitly, so the
caller can retain its ligand matching, no-fit symmetry RMSD, sentinel rejection
and coordinate deduplication. It does not modify a completed campaign's driver.
Run from the repository root with paths supplied by the campaign owner:

```bash
python3 scripts/score_emitted_pools.py \
  --prefix /path/to/cell/TARGET \
  --prefix /path/to/cell/r1/TARGET \
  --prefix /path/to/cell/r2/TARGET \
  --reference /path/to/TARGET_ligand.sdf \
  --scoring-code /path/to/frozen/code/score_electors.py \
  --pose-budget 250 \
  --out /path/to/new/cell_budget250.json
```

The provider must expose `Reference` and `pool_tiers`. The original campaign's
provider requires NumPy, SciPy and spyRMSD; the enumerator itself uses only the
standard library. `--scoring-code` loads Python code and must name a trusted file.
The adapter passes CF as the ordering key for the top-N diagnostic (default 10).
The oracle considers the entire retained, deduplicated, scorable pool.

The JSON records the requested budget, ordered prefixes, raw/parsed/scored pool
counts, parse failures and empty prefixes. It fingerprints enumerated PDBs,
their optional `.mcf` sidecars, the reference, provider, adapter and enumerator.
Successful scoring exits 0. Missing/empty prefixes, an empty scored pool, or
missing/nonfinite CF produce diagnostic JSON and exit 2. Invalid arguments,
provider imports and filesystem errors fail without a success receipt. Existing
output files are never overwritten.

## Interpretation limits

- Increasing the reader budget cannot recover poses never written by the
  engine. The engine's result budget also affects cluster construction; a fresh
  raised-budget arm is not a stored default-budget arm plus extra output.
- `cf_truncated` retains the frozen predicate: the exact next filename
  `<prefix>_<budget>.pdb` exists. A gap at that boundary can hide higher ranks.
  `fo_truncated` means discovered FO heads were omitted. Neither flag proves an
  upstream cap is absent or non-binding.
- `n_not_in_scored_pool` includes all downstream provider omissions, including
  intended filtering/deduplication and possible matching/RMSD errors. It does
  not classify those reasons or certify that every emitted pose was valid.
- A non-vacuity proof uses identical stored inputs at both reader budgets and
  requires `n_pool_scored` to increase. This proves reader coverage, not an
  engine write budget, scientific benefit, or docking success.

## Verification

```bash
python3 -m pytest tests/test_offline_pose_budget.py -q
```

The suite covers explicit budget propagation through every restart, unchanged
default-50 behavior despite environment settings, shared CF/FO limits, later-rank
eligibility, exclusions, malformed inputs, missing restart accounting and CLI
receipts. The implementation handoff additionally carries a real stored-pose
50-versus-250 proof and equality checks against the frozen default-50 scorer.
