# Support Boundary

**Last actualized**: 2026-07-12  

The repository is a research platform. **Core 1.0** is a deliberate subset.

## Supported (Core 1.0)

- `FlexAIDdS` / `FlexAID` / `tENCoM` CLI workflows documented for core docking and vibrational entropy
- `flexaidds` Python package for results I/O, thermodynamics analysis, and documented APIs
- JSON configuration paths that map to those CLIs
- Benchmark **bundles** under `benchmarks/` with reproducibility artifacts
- Installation, support matrix, known limitations, and security policy docs

Details: [`docs/VALIDATED_CAPABILITIES.md`](../VALIDATED_CAPABILITIES.md).

## Experimental (not support-guaranteed)

- Swift packages and Apple device integration
- TypeScript / PWA dashboards
- Fleet / iCloud distributed execution layers
- NATURaL co-translational workflows (code + tests exist; product support experimental)
- GPU backends not covered by the support matrix for a given release

Details: [`docs/EXPERIMENTAL_CAPABILITIES.md`](../EXPERIMENTAL_CAPABILITIES.md).

## Interpreting docs and marketing language

README vision language may describe the full research platform. For release trust, combine:

1. Support Matrix
2. Validated Capabilities
3. Known Limitations
4. Reproducibility policy for any numeric benchmark claim
