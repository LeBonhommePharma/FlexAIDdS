# ChatGPT Instructions — FlexAIDdS

Use these instructions when working with ChatGPT (Custom GPT, Projects, or agent mode) on the FlexAIDdS repository.

**Source of truth**: The full authoritative rules live in `AGENTS.md` (repo root). This is a condensed, GPT-optimized version of the most important constraints. When in doubt, ask the user to paste the latest `AGENTS.md`.

---

## Core Rules (Memorize These)

- Always verify with **actual command execution** (build + tests) before claiming anything is fixed, done, or working. Show the passing output.
- For any task with 3+ steps, explicitly use a todo list (one item in progress at a time).
- After code changes, the user must commit and push immediately. Do not batch changes.
- Fresh configure + build after any CMakeLists.txt or LIB/ source change.
- Zero test failures before any suggested push.
- Complete every item on a prioritized list (P0/P1/etc.) before stopping.
- When the user says “run the command / benchmark / test”, actually run it.

---

## Project Context (Short)

FlexAIDdS = entropy-driven molecular docking (genetic algorithm + statistical mechanics thermodynamics) for drug discovery.

Key layers:
1. GA (`gaboom.cpp`)
2. Voronoi scoring (`Vcontacts.cpp`)
3. StatMech (`statmech.cpp`)
4. BindingMode clustering + thermo
5. ENCoM / tENCoM vibrational entropy
6. ShannonThermoStack (configurational entropy, GPU accelerated)
7. Cavity detection

**Languages**: C++26 (core), Python (flexaidds package + bindings), Metal, optional CUDA.

**License**: Apache-2.0 only. No GPL anything.

---

## Essential Commands

**Most common C++ workflow**
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure
```

**Python only**
```bash
cd python && pip install -e . && pytest tests/ -q
```

Full details and CMake options are in `CLAUDE.md`.

---

## Behavior Expectations

- Be extremely strict about the verification and commit rules above.
- When the user gives you a task, start by confirming you have read `AGENTS.md` (ask them to paste it if you don't have it in context).
- Prefer running real commands over guessing.
- If something is ambiguous, ask clarifying questions rather than over-exploring the codebase.

Keep responses focused. This is a serious scientific codebase — precision and discipline matter.
