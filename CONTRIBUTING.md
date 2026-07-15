# Contributing to FlexAID∆S

Thank you for helping build an entropy-driven docking engine for real-world psychopharmacology.

This document explains how to contribute while keeping the codebase fast, scientifically honest, and legally clean.

## 1. Licensing

- All contributions are accepted under the **Apache License 2.0**.
- By opening a pull request, you agree that your changes may be redistributed under Apache-2.0.
- Do **not** submit GPL/AGPL-licensed code or anything derived from such projects. See `docs/licensing/clean-room-policy.md`.

## 2. How to Contribute

1. **Fork and branch**

   Fork the repo and create a feature branch from `master`.

2. **Open an issue (recommended)**

   For new features, dependency changes, or behavior changes, open an issue describing the proposal and design.

3. **Implement your change**

   - Follow the existing C++ style in nearby files.
   - Keep performance in mind: FlexAID∆S is a docking engine, not a toy.
   - Add or update tests where behavior changes.

4. **Run tests locally**

   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
   cmake --build build -j $(nproc)
   ctest --test-dir build --output-on-failure
   ```

   If Python bindings are present, also run:

   ```bash
   pytest python/tests -v
   ```

5. **Open a pull request**

   - Summarize what you changed and why.
   - Call out any API changes, new config options, or new dependencies.
   - Mention any licensing or build-system impact.

## 3. Coding Guidelines

- **Languages**

  - C++26 for the core engine.
  - Python for high-level workflows, analysis, and CLI tooling.

- **Formatting & static analysis (scaffolding — no mass reformat)**

  Style is historically mixed (tabs in legacy GA/scoring files; 4-space indent in newer modules such as `statmech` / `encom`). Configs exist to guide **new** code without rewriting the whole tree:

  | File | Purpose |
  |------|---------|
  | `.clang-format` | C/C++ format (LLVM-based, 4 spaces, column 100) |
  | `.clang-tidy` | Small advisory check set (nullptr, override, a few bugprone/performance) |
  | `.pre-commit-config.yaml` | Optional local hooks (clang-format on staged files, YAML, large-file guard) |

  ```bash
  # Dry-run check on the default modern-module allowlist (never rewrites):
  scripts/format_check.sh

  # Check only untracked/modified C++ under LIB/, src/, tests/, python/bindings/:
  scripts/format_check.sh --all-new

  # Format a file you own in this PR (intentionally):
  clang-format -i path/to/file.cpp

  # Optional pre-commit (local):
  pip install pre-commit && pre-commit install
  pre-commit run
  ```

  **Do not** open a PR that only runs `clang-format` over all of `LIB/`. That explodes review and mixes style debt with functional changes.

- **Performance**

  - Avoid unnecessary heap allocations, hidden O(n²) loops, or branches inside hot scoring paths.
  - Respect existing parallelization (OpenMP, parallel GA) and avoid race conditions.

- **Thread safety**

  - Any change touching GA, scoring, or ensemble statistics must document its thread-safety assumptions.

- **Error handling**

  - Return clear error codes in C++ or throw standard exceptions.
  - In Python, raise meaningful exceptions instead of silently failing.

- **Astex / benchmark structures**

  Canonical Astex Diverse tree: `benchmarks/astex_diverse/astex_diverse/`.
  See `benchmarks/datasets/CANONICAL.md` and checksums in
  `benchmarks/datasets/astex_diverse_sha256.csv`. Do not treat nested
  `benchmarks/astex_diverse/data/astex_diverse/` as source of truth.

## 4. Dependencies

- Preferred licenses: Apache-2.0, BSD, MIT, MPL-2.0, PSF.
- GPL/AGPL dependencies are not allowed.
- Any new dependency must be documented in `THIRD_PARTY_LICENSES.md` and justified in the PR.

## 5. Documentation

Documentation is part of the work, not an optional extra:

- Update README or `docs/` for new features.
- Document configuration keywords, default values, and expected ranges.
- For scientific changes, briefly state the method and key references.

## 6. Code of Conduct

We attack bad ideas and fragile methods, not people. See `CODE_OF_CONDUCT.md` for behavior expectations.

Harm reduction, anti-stigma, and scientific rigor are non-negotiable here. If your contribution helps people get better drugs, faster and safer, you're in the right place.
