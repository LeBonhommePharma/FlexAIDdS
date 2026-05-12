# Known Limitations

This file documents limitations that matter for installation trust, scientific interpretation, and release planning.

## Product-scope limitations

- The repository contains more code than the supported 1.0 surface. Not every visible feature is part of the supported release contract.
- Experimental surfaces may compile or exist in-tree without being considered support-guaranteed.

## Reproducibility limitations

- Benchmark tables and scientific claims should not be assumed replayable unless a corresponding bundle exists under `benchmarks/`.
- Some documented results may still be preliminary or manuscript-bound rather than repository-reproducible.

## Backend limitations

- Backend breadth does not imply release-grade parity.
- Scalar CPU, OpenMP, and core CPU workflows are the primary supported base.
- Backend-specific accelerators may remain optional or experimental until covered by release validation.

## Documentation limitations

- Historical README and site language may reflect the full research vision rather than the narrow supported Core 1.0 surface.
- Public-facing documentation should be interpreted together with `PRODUCT.md`, `docs/SUPPORT_MATRIX.md`, and `docs/EXPERIMENTAL_CAPABILITIES.md`.

## Security limitations

- Legacy C and C++ parsing surfaces require continued hardening and regression testing.
- Presence of an audit or policy document does not itself imply closure; fixes and automation are the real closure criteria.

### Current Security Status (May 2026)

The March 2026 buffer-overflow audit remains useful as a historical threat
inventory, but it is not a current closure report. Several high-risk patterns
have since been bounded in source, while legacy C/C++ config parsing still
requires ongoing hardening and regression tests.

**Current live status:**
- legacy GIST config path parsing now uses bounded copy helpers instead of
  unbounded `%s` reads
- buffer-safety regression coverage exists in `tests/test_buffer_safety.cpp`
- remaining unsafe string patterns should be tracked as source findings, not
  inferred from the stale March audit totals

**Recommended mitigation until the legacy parser is fully retired:**
- do not process untrusted configuration files
- use signed or repository-controlled configuration bundles for benchmarks
- keep config files small and reviewable
- run parser changes under ASAN/UBSAN where available

See [`SECURITY_HARDENING_ROADMAP.md`](SECURITY_HARDENING_ROADMAP.md) and
[`SECURITY_AUDIT_BUFFER_OVERFLOW.md`](../SECURITY_AUDIT_BUFFER_OVERFLOW.md) for
historical context. Treat those documents as audit inputs, not proof of current
security closure.

## Scientific interpretation limitations

- Entropy-aware docking can be useful without implying universal superiority on every target class.
- Claims should be evaluated against the specific benchmark bundle, dataset scope, preprocessing rules, and metric definition used.
