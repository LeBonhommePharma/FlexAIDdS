# Reproducibility Policy

This document defines what must exist before a benchmark or scientific performance claim is treated as **repository-reproducible**.

## Levels of claim maturity

### 1. Replayable from repository artifacts

A claim reaches this level only when all of the following are present:

- dataset provenance or acquisition script
- checksums or immutable identifiers
- preprocessing steps
- exact command lines
- fixed seeds where applicable
- expected outputs and metric calculation scripts
- recorded git SHA and environment details

## Seed provenance

Repository-reproducible stochastic runs must record the run-level seed. The
canonical environment variable is:

```bash
export FLEXAID_SEED=42
```

Core fallback RNG paths that previously seeded directly from
`std::random_device` should route through `LIB/RngSeed.h`. When
`FLEXAID_SEED` is set, each call site derives a deterministic stream seed from
that run seed. When it is unset, the helper falls back to `std::random_device`
for exploratory runs.

This seed is necessary provenance, not a complete determinism guarantee for
parallel algorithms. Any benchmark intended to be replayable must also record
thread counts, backend selection, compiler flags, input order, and whether
parallel scheduling can change draw order.

### 2. Preliminary

A claim is preliminary if it appears in documentation but is not yet backed by a replayable bundle in the repository.

Preliminary claims must be clearly labeled as such.

### 3. External / published

A claim may also be labeled external or published when it is validated in a peer-reviewed publication or equivalent external artifact.

## Required benchmark bundle layout

Every replayable benchmark should include a bundle under `benchmarks/` with at least:

- `README.md`
- `manifest.yaml`
- `download.sh` or equivalent acquisition instructions
- `run.sh` or equivalent execution script
- `expected/` outputs or metric snapshots
- `environment.txt` or machine metadata template

## Minimum manifest fields

Each benchmark manifest should describe:

- benchmark name
- dataset source
- dataset version or immutable identifier
- preprocessing steps
- executable and config used
- seed(s)
- output artifact paths
- metric definitions
- known limitations

## What must not happen

- no benchmark table should imply full reproducibility if the corresponding bundle is missing
- no metric should be called final if the replay scripts do not exist
- no mixed reporting of exploratory and release-grade numbers without explicit labels

## Immediate repository direction

The first reproducibility target for Core 1.0 should be a small smoke benchmark bundle that can run in CI, followed by larger manually triggered bundles for full dataset evaluation.
