# Testing (MkDocs)

**Last actualized**: 2026-07-12  

Full guide: [`docs/TESTING.md`](../TESTING.md) · inventory: [`docs/TEST_COVERAGE_ANALYSIS.md`](../TEST_COVERAGE_ANALYSIS.md).

## Commands

```bash
# C++
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

# Python
PYTHONPATH=python python3.11 -m pytest python/tests/ -q --tb=line

# Skill packaging
python3 .grok/skills/flexaidds/scripts/validate_skill.py
```

## What CI covers

- Multi-platform GoogleTest (`ci.yml`)
- Pure-Python result loader slices
- Optional bindings smoke
- lcov coverage on `LIB/*` (`coverage.yml`)
- Sanitizers / license scan on separate workflows

## Policy

Never claim a change is done without running the relevant suite. Prioritized work lists must complete fully. See `AGENTS.md`.
