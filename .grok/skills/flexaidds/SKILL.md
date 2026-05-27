---
name: flexaidds
description: Expert developer assistant for the FlexAIDdS entropy-driven molecular docking engine. Deep knowledge of the C++26 core (LIB/), Python package (flexaidds), CMake build system, testing discipline (ctest + pytest), architecture, and strict "verify then commit" workflow rules. Use when working on any part of the FlexAIDdS codebase, running benchmarks, modifying LIB/ modules, Python bindings, or asking about docking/statmech/tENCoM. Slash command: /flexaidds. Automatically activates on mentions of FlexAIDdS, FlexAID, docking engine, or related modules.
---

# FlexAIDdS Development Skill

You are an expert on FlexAIDdS (FlexAID with ΔS Entropy) — an entropy-driven molecular docking engine combining genetic algorithms with statistical mechanics thermodynamics. It targets real-world psychopharmacology and drug discovery applications.

**Primary languages**: C++26 (core engine in LIB/), Python (bindings + analysis in `python/flexaidds/`), Objective-C++ (Metal GPU), CUDA (optional).

**License**: Apache-2.0 only. Never introduce GPL/AGPL dependencies. Follow the clean-room policy in `docs/licensing/`.

## Core Workflow Rules (Non-Negotiable)

- **Verify with actual execution before claiming anything is done**. Run the build or test command and show passing output. Never say "fixed", "implemented", "complete", or "should work" without evidence from running the tools.
- **Multi-step tasks (3+ actions) start with todo_write**. Define the full list with `merge: false`. Keep exactly **one** item `in_progress` at a time. Mark items `completed` immediately when finished — never batch. Re-read the todo list before ending any turn that still has pending/in-progress work.
- **After any code change, commit and push immediately**. Use conventional prefixes (`Fix:`, `Add:`, `Update:`, `Refactor:`). Do not batch multiple changes. If a git operation hangs, kill stale processes and retry.
- **Fresh builds after CMake or source changes**. Never assume a target builds. After editing CMakeLists.txt or adding .cpp/.h files under LIB/, run a clean configure + build and confirm linking succeeds. Watch for disk space issues on long builds.
- **0 test failures before any push**. Run `ctest --output-on-failure` (C++) or `pytest` (Python) after relevant changes. Fix in the same session.

## Build & Test Commands

```bash
# Full C++ release build with tests (recommended starting point)
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure

# Python package only (no C++ needed for many tasks)
cd python
pip install -e .
pytest tests/ -q

# Python bindings smoke test
cmake -B build -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
# then run relevant pytest tests that use @requires_core
```

**Key CMake options** (see CLAUDE.md for full table):
- `BUILD_TESTING=ON` — GoogleTest suite
- `BUILD_PYTHON_BINDINGS=ON` — builds `_core` extension (statmech + encom)
- `FLEXAIDS_USE_METAL=ON` — macOS GPU acceleration
- `FLEXAIDS_USE_CUDA=ON` — CUDA GPU kernels

**Main targets**: `FlexAID` (core executable), `_core` (pybind11), `tENCoM`, `benchmark_tencom`, `benchmark_vcfbatch`.

## Architecture Overview

Core pipeline:
1. Genetic Algorithm — `LIB/gaboom.cpp`
2. Voronoi Contact Scoring — `LIB/Vcontacts.cpp`
3. Statistical Mechanics — `LIB/statmech.cpp` (partition functions, free energy, entropy, WHAM, TI)
4. Binding Mode Clustering + Thermodynamics — `LIB/BindingMode.cpp`
5. Vibrational Entropy (ENCoM) — `LIB/encom.cpp`
6. Shannon Configurational Entropy — `LIB/ShannonThermoStack/` (CPU/CUDA/Metal)
7. Cavity Detection — `LIB/CavityDetect/` (SURFNET + Metal)

Major specialized modules under `LIB/`:
- `tENCoM/` — Torsional ENCoM backbone flexibility + differential entropy
- `ShannonThermoStack/` — Hardware-accelerated histogram entropy
- `LigandRingFlex/`, `ChiralCenter/`, `NATURaL/`, `CavityDetect/`

Python package (`python/flexaidds/`):
- `models.py`, `results.py` (`load_results()`), `docking.py`, `encom.py`, `tencm.py`, `thermodynamics.py`, `io.py`
- CLI entry: `python -m flexaidds <results_dir> [--json|--csv|--top N]`
- Bindings: `python/bindings/core_bindings.cpp` (exposes StatMechEngine, ENCoMEngine, etc.)

PyMOL plugin lives in `pymol_plugin/`.

## Key Files & Navigation

See the full table in CLAUDE.md. Most important entry points:
- `LIB/flexaid.h` (central header)
- `CMakeLists.txt` + `cmake/`
- `python/flexaidds/__init__.py` (public API surface)
- `tests/` (C++) and `python/tests/` (pytest, many with `@requires_core`)
- `benchmarks/` (dataset YAMLs, runners, Astex Diverse, CASF, etc.)
- `docs/` (roadmaps, implementation notes, licensing)

## Common Development Tasks

**Adding a new C++ source file**
1. Place `.cpp`/`.h` under `LIB/`
2. Add it to the `FLEXAID_SOURCES` list in `CMakeLists.txt`
3. Write corresponding GoogleTest in `tests/`
4. Enable `BUILD_TESTING=ON`, build, and run `ctest`

**Adding a new Python module**
1. Create under `python/flexaidds/`
2. Export from `__init__.py` if public
3. Add tests in `python/tests/`
4. Add fixtures in `python/conftest.py` if needed
5. For C++ exposure: add to `python/bindings/core_bindings.cpp`

**Running the full verification loop** (use this before any commit that touches core logic):
- Build with tests
- Run `ctest --output-on-failure`
- (If Python changed) `cd python && pip install -e . && pytest tests/`
- Confirm zero failures
- Commit + push immediately

## Usage

- **Slash command**: `/flexaidds` (injects this skill)
- **Automatic**: Grok will invoke this skill when the prompt mentions FlexAIDdS, FlexAID, docking, statmech, tENCoM, BindingMode, or related development tasks in this repo.
- **With other skills**: Can be combined with `/review`, `/implement`, `/check`, etc.

Always read the latest CLAUDE.md at the start of any substantial session for the complete reference. This skill + CLAUDE.md together give you production-grade context for the entire codebase.
