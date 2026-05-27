# AGENTS.md — FlexAIDdS AI Coding Agents

**This is the single source of truth for all AI coding agents working in this repository.**

Every other instruction file (CLAUDE.md, `.grok/skills/flexaidds/SKILL.md`, ChatGPT instructions, etc.) must derive from or explicitly defer to this document. When rules conflict, **this file takes precedence**.

All agents are expected to internalize and strictly follow the rules below on every task.

---

## Core Workflow Rules (Non-Negotiable)

These rules exist because this is a complex, performance-critical scientific codebase. Sloppy execution creates expensive bugs and broken builds.

1. **Verify with actual execution before claiming success.**  
   Never say “done”, “fixed”, “implemented”, “should work”, or “complete” until you have run the relevant build and/or test commands and shown clean, passing output.  
   - C++ changes → run the full test build + `ctest --output-on-failure`.  
   - Python changes → run the affected pytest suite.  
   - CMake or source list changes → perform a fresh configure + build.

2. **Use `todo_write` for any task with 3 or more distinct actions.**  
   Open the list with `merge: false` at the start of the task. Keep **exactly one** item in `in_progress` at all times. Mark items `completed` immediately when finished — never batch completions. Re-read the current todo list before ending any turn that still has pending or in-progress work.

3. **After any code change, commit and push immediately.**  
   Do not batch multiple logical changes. Use conventional commit prefixes (`Fix:`, `Add:`, `Update:`, `Refactor:`, etc.).  
   If a git operation appears to hang, kill stale git processes (`kill $(pgrep -f git)`) and retry. Check `git config core.fsmonitor` if problems persist.

4. **Fresh builds after touching build system or sources.**  
   Never assume a target still builds after editing `CMakeLists.txt` or adding/removing `.cpp`/`.h` files under `LIB/`. Always run a clean configure + build and confirm linking succeeds. Watch for disk space issues before long builds.

5. **Zero test failures before any push.**  
   Run the full relevant test suite (`ctest --output-on-failure` for C++, `pytest` for Python) after changes that could affect behavior. Fix all failures in the same session. Never push with red tests.

6. **Prioritized lists are completed in full.**  
   When given P0/P1/P2 (or similar), finish every item before stopping. Do not declare victory after the first few.

7. **Run the thing when asked.**  
   If the user explicitly asks you to *run* a command, benchmark, or test, do it — do not spend 20+ tool calls exploring first.

---

## Essential Commands

**C++ (recommended starting point for most work)**
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure
```

**Python package only** (many tasks do not need the full C++ build)
```bash
cd python
pip install -e .
pytest tests/ -q
```

**Python bindings smoke**
```bash
cmake -B build -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
```

See `CLAUDE.md` → “Build System” for the full table of CMake options and targets.

---

## Project Snapshot

FlexAIDdS is an entropy-driven molecular docking engine that combines genetic algorithms with statistical mechanics (partition functions, free energy, vibrational + configurational entropy).

**Primary languages**: C++26 (core in `LIB/`), Python (`python/flexaidds/` + bindings), Objective-C++ (Metal), optional CUDA.

**Strict licensing**: Apache-2.0 only. Never introduce GPL/AGPL dependencies or use GPL code as inspiration (see `docs/licensing/clean-room-policy.md` and `THIRD_PARTY_LICENSES.md`).

**Key architectural layers** (in execution order):
1. Genetic Algorithm (`LIB/gaboom.cpp`)
2. Voronoi Contact Scoring (`LIB/Vcontacts.cpp`)
3. Statistical Mechanics engine (`LIB/statmech.cpp`)
4. BindingMode clustering + thermodynamic integration (`LIB/BindingMode.cpp`)
5. Vibrational entropy (ENCoM + tENCoM)
6. Shannon configurational entropy with hardware acceleration (`ShannonThermoStack/`)
7. Cavity detection (`CavityDetect/`)

Full details live in `CLAUDE.md`.

---

## Maintaining These Instructions (Production Discipline)

- `AGENTS.md` is the **authoritative source** for workflow rules and high-level constraints.
- When you change core rules or commands, update `AGENTS.md` first.
- Then propagate the important changes into:
  - `CLAUDE.md` (Claude’s richer, more detailed view)
  - `.grok/skills/flexaidds/SKILL.md` (Grok skill)
  - Any ChatGPT / custom agent instructions
- All three files should explicitly state that `AGENTS.md` is the source of truth.
- Keep the sacred “Core Workflow Rules” section as close to identical as possible across files to prevent drift.

This three-file system (AGENTS.md + CLAUDE.md + Grok skill) is deliberately designed for minimal long-term maintenance burden while giving each tool the format it works with best.

---

## When in Doubt

1. Read `AGENTS.md` first (this file).
2. Then read the tool-specific file for the agent you are currently using.
3. When the task involves running builds or tests, actually run them and show the output.

These rules exist to protect the quality and velocity of a real scientific codebase. Follow them without exception.
