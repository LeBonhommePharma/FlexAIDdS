# Testing Guide

**Last actualized**: 2026-07-12

This guide is the operator-facing companion to [`TEST_COVERAGE_ANALYSIS.md`](TEST_COVERAGE_ANALYSIS.md).  
Authoritative workflow rules: repository root [`AGENTS.md`](../AGENTS.md).

## Goals

1. Prove behavioral contracts for Core 1.0 (CLI engines + `flexaidds` package).
2. Protect thermodynamic identities (partition function, F/H/S consistency, audit schema).
3. Keep parsers and ranking diagnostics from regressing silently.
4. Never conflate CF/contact-function **scoring proxy** success with thermodynamic ledger correctness.

## Quick start

### C++ unit tests

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
ctest --test-dir build --output-on-failure
```

Filter examples:

```bash
ctest --test-dir build -R 'StatMech|ThermoLedger|BindingMode' --output-on-failure
ctest --test-dir build -R BufferSafety --output-on-failure
```

### Python package tests

```bash
# Recommended interpreter with scientific stack
export PYTHONPATH=python
python3.11 -m pytest python/tests/ -q --tb=line
```

Focused slices:

```bash
python3.11 -m pytest python/tests/test_thermodynamics.py python/tests/test_thermo_schema.py -q
python3.11 -m pytest python/tests/test_dataset_adapters.py python/tests/test_truncate_chain.py -q
python3.11 -m pytest python/tests/test_ranking_bias.py python/tests/test_fallback_types.py -q
```

### Skill / agent packaging

```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
python3 scripts/check_repo_hygiene.py
```

### Optional C++ bindings smoke

```bash
cmake -B build -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
# then pytest with @requires_core tests enabled
```

## Markers and optional deps

| Situation | Behavior |
|-----------|----------|
| `_core` extension missing | `@requires_core` tests skip; pure-Python StatMech/ENCoM fallbacks still tested |
| numpy / pyyaml missing | Some dataset_runner / adapter tests fail import — install deps or use env with numpy |
| PyMOL missing | Plugin code is not part of default pytest; treat as optional integration |
| GPU / Metal | Dispatch unit tests run on CPU; device kernels need hardware-gated jobs |

## What “passing” means scientifically

| Claim | Required evidence |
|-------|-------------------|
| Thermodynamic identities hold | `test_statmech`, `test_thermo_ledger`, Python thermo / schema tests |
| Audit package schema valid | `test_thermo_schema` (+ runner packaging dry-runs) |
| CF ranking bias diagnostics | `test_ranking_bias` (does **not** reimplement C++ CF) |
| Pose/mode bookkeeping | `test_binding_mode_*` |
| Config buffer safety | `test_buffer_safety` |
| Benchmark orchestration | DatasetRunner dry-run tests + skill validator |

**Not** proven by unit tests alone: published docking success rates, PoseBusters pass rates, ITC affinity accuracy. Those require benchmark bundles under `benchmarks/` and the benchmarking skill contract.

## Adding tests (checklist)

1. Prefer pure functions and small fixtures over full docking runs in CI.
2. Register every new `tests/test_*.cpp` in `CMakeLists.txt`.
3. Mirror new `python/flexaidds/<module>.py` with `python/tests/test_<module>.py` when practical.
4. Use precise names: CF score vs free energy; never assert “ΔG experimental” from uncalibrated CF.
5. Run the **affected** suite before commit; run broader suites before push (`AGENTS.md`).

## CI map

| Workflow | Role |
|----------|------|
| `ci.yml` | Multi-platform C++ ctest + pure Python slices + bindings smoke |
| `coverage.yml` | lcov on `LIB/*` |
| `sanitizers.yml` / `tsan.yml` | Memory / race instrumentation |
| `license-scan.yml` | Apache-2.0 / GPL isolation guardrails |

## Related

- Coverage inventory: [`TEST_COVERAGE_ANALYSIS.md`](TEST_COVERAGE_ANALYSIS.md)
- Thermodynamics semantics: [`thermodynamics.md`](thermodynamics.md)
- Support boundary: [`SUPPORT_MATRIX.md`](SUPPORT_MATRIX.md), [`VALIDATED_CAPABILITIES.md`](VALIDATED_CAPABILITIES.md)
- Benchmark claims: [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), [`BENCHMARKS.md`](BENCHMARKS.md)
