---
name: flexaidds
description: Expert developer assistant for the FlexAIDdS entropy-driven molecular docking engine. Deep knowledge of the C++26 core (LIB/), Python package (flexaidds), CMake build system, testing discipline (ctest + pytest), architecture, and strict "verify then commit" workflow rules (see AGENTS.md). Use when working on any part of the FlexAIDdS codebase, running benchmarks, modifying LIB/ modules, Python bindings, or asking about docking/statmech/tENCoM. Also supports on-demand generation of Nature Reviews Drug Discovery-style cover figures (dramatic entropy/enthalpy or E–E index molecular compositions with full FlexAID∆S branding). Slash command: /flexaidds. Automatically activates on mentions of FlexAIDdS, FlexAID, docking engine, or related modules.
---

# FlexAIDdS Development Skill

You are an expert on FlexAIDdS (FlexAID with ΔS Entropy) — an entropy-driven molecular docking engine combining genetic algorithms with statistical mechanics thermodynamics. It targets real-world psychopharmacology and drug discovery applications.

**Primary languages**: C++26 (core engine in LIB/), Python (bindings + analysis in `python/flexaidds/`), Objective-C++ (Metal GPU), CUDA (optional).

**License**: Apache-2.0 only. Never introduce GPL/AGPL dependencies. Follow the clean-room policy in `docs/licensing/`.

## Core Workflow Rules (Non-Negotiable)

**Authoritative version**: See `AGENTS.md` in the repository root. The rules below are the Grok-optimized summary.

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

## Figure Generation (new first-class action)

The skill now exposes automated, publication-quality NRDD-style cover figure generation via the dedicated module `.grok/skills/flexaidds/flexaidds_skill.py` and its manifest `.grok/skills/flexaidds/flexaidds_skill_manifest.json`.

The primary action is `generate_flexaids_figure`. See the module docstring and the JSON manifest for the full contract (parameters, validation, returns, dependency-injection of the image generator, reproducibility sidecar).

Example usage from within the skill/agent:
```python
from .flexaidds_skill import generate_flexaids_figure, FigureParameters

params = FigureParameters(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style="dramatic_faces")
result = generate_flexaids_figure(
    params=params,
    # Supply your platform's image tool here (Grok image_gen, DALL·E wrapper, etc.)
    image_generator=your_platform_image_callable,
    output_dir="results/figures"
)
# result["prompt"] and (if generated) result["path"] + result["metadata_path"]
```

The implementation follows the five-point professional integration pattern:
- Entry point in this skill (flexaidds_skill.py).
- Robust handler with validation, safe overrides, logging, and dependency injection for the image generator.
- Exposed via the companion JSON manifest.
- No hard dependencies on any specific image library.
- Deterministic, reproducible output (unique filenames + full metadata JSON).

Always prefer calling the real Python helper (it reuses the canonical prompt builder from `python/flexaidds/figures` for branding, typography, scientific accuracy, and E-E index handling). The skill will never hand-craft prompts.

Always read `AGENTS.md` (source of truth for rules) + the latest `CLAUDE.md` at the start of any substantial session. This skill is the Grok-optimized companion to those two files.

## NRDD Cover & Journal-Style Figure Generation (automated, reproducible)

The /flexaidds skill now supports on-demand generation of high-end Nature Reviews Drug Discovery cover-style illustrations featuring FlexAID∆S (dramatic entropy/enthalpy personification or molecular E–E gauge compositions, precise thermodynamic call-outs for TΔS / -TΔS / I_E–E index, full branding, JetBrains Mono + thebonhomme.com typography, PLIP-inspired interaction clarity on key contacts, and emphasis on the most favourable CF-contributing interactions).

This is implemented as a first-class, professional capability in the Python package (the canonical entry point) and fully documented here so the agent can invoke it reliably.

### The 5-point integration (implemented bulletproof / idiotproof / scientifically robust)

1. **Identify your skill’s entry point**  
   The primary implementation lives in the package the skill always recommends: `python/flexaidds/figures.py` (and exported via `python/flexaidds/__init__.py`). The skill description (this file) is the "manifest" that registers the capability for /flexaidds users and agents. No separate JSON/OpenAPI is needed; the SKILL.md + package docstrings serve that role.

2. **Add a figure-generation handler** (`generate_flexaids_nrdd_cover`)  
   Located at `flexaidds.figures.generate_flexaids_nrdd_cover` (and the pure `build_nrdd_cover_prompt` + `NRDDCoverParams` dataclass).  
   - Accepts typed parameters (entropy_value, enthalpy_value, index_value for I_E–E, style="dramatic_faces"|"molecular_gauge", title, subtitle, results_dir for sourcing real ensemble values, etc.).  
   - Constructs an extremely detailed prompt that faithfully reproduces the visual language, composition, call-out boxes, dramatic lighting, central ligand, floating cubes/gauge, bottom text panels, logos, and branding of the reference covers.  
   - Injects the exact numeric values, enforces JetBrains Mono / thebonhomme.com typography for all text, PLIP-style color-coded interaction lines where molecules are shown, and prioritizes "the most favourable contacts and those contributing the most to the CF".  
   - Full validation (ranges for thermodynamic quantities, index ~[-1,1], existing results_dir, etc.).  
   - Returns a dict with the ready-to-use prompt, rich metadata (every param, timestamp, git sha, scientific guardrail note, suggested image_gen + image_edit calls), and suggested aspect ratio.  
   - Never claims experimental ΔG; always uses precise FlexAID∆S terminology ("ensemble thermodynamic ledger", "visualisation only").  
   - Reproducible sidecar-friendly (save the metadata dict as .json next to the image).

   Example usage (the professional implementation the skill exposes):

   ```python
   from flexaidds.figures import generate_flexaids_nrdd_cover, NRDDCoverParams

   res = generate_flexaids_nrdd_cover(
       entropy_value=0.93,
       enthalpy_value=1.4,
       index_value=0.92,
       style="dramatic_faces",           # or "molecular_gauge"
       title="The ΔG balance",
       subtitle="Striking the right pose in drug discovery",
       results_dir="/path/to/a/docking/run",  # optional: auto-source real ledger values
   )
   # res["prompt"]          → the full engineered string
   # res["metadata"]        → everything needed for audit/reproducibility
   # res["suggested_call"]  → "image_gen(prompt=..., aspect_ratio='16:9')"
   ```

3. **Expose the function via the skill manifest**  
   This SKILL.md (and the package `__init__`) documents the action `generate_flexaids_nrdd_cover` / "create NRDD cover" with all parameters and return contract. When a user says "/flexaidds generate a cover figure with E-E = 0.92 ..." or "make the Nature Reviews illustration for these thermodynamic values", the skill activates and the agent follows the procedure above.

4. **Manage dependencies**  
   The heavy lifting for actual pixel generation is the host environment's image generation tool (in this Grok context: the `image_gen` tool, plus `image_edit` for refinements). The Python package itself has no hard dependency on any image model or GPL code. PLIP is optional for base interaction diagrams (already integrated in the same module). All prompts are plain text and work with any compatible image model (Grok, Claude, etc.).

5. **Return the result**  
   The function returns the prompt + metadata immediately. In the agent/skill execution:
   - Call the Python helper (or run the equivalent prompt construction).
   - Invoke the available image generation tool with the prompt and correct aspect (16:9 landscape or 3:2 for cover feel).
   - (Optional) follow up with image_edit for label tightening, contrast on the I_E–E gauge/cubes, etc.
   - Save the output image next to a `cover_metadata.json` (the returned metadata) for full reproducibility.
   - Report the local path (or sync_file identifier) back to the user. The client can display or download it.

### Usage examples (agent & user)
```bash
# After a docking run (pulls real values where possible)
python -c "
from flexaidds.figures import generate_flexaids_nrdd_cover
res = generate_flexaids_nrdd_cover(results_dir='results/my_run', style='molecular_gauge')
print(res['prompt'][:300])
# Then (in Grok/agent):
# image = image_gen(prompt=res['prompt'], aspect_ratio='16:9')
# print('Generated:', image['path'])
"
```

The skill will also happily generate the figures completely standalone with user-supplied thermodynamic numbers (perfect for talks, papers, or the exact reference covers you showed).

All generation is post-processing / promotional only. It never affects scientific results, ranking, or claims.

Always verify the generated image against the metadata and the original docking ledger before use in any publication context.
