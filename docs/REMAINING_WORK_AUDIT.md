# Remaining Work Audit (as of 2026-04-05)

This document captures **remaining bugs, TODOs, and roadmap work** currently visible in the repository.

## 1) Open TODOs in production C++ code

### Core pipeline integration gaps

1. **DatasetRunner docking execution is still a stub**.
   - `LIB/DatasetRunner.cpp` still carries a TODO to wire the real
     `setup_direct_input() -> gaboom()` execution path and currently only records timing/placeholder output.

2. **Parallel campaign still waits on GA re-entrancy refactor**.
   - `LIB/ParallelCampaign.cpp` has a TODO indicating model-level docking should become an OpenMP parallel loop after GA becomes re-entrant.

3. **Top-hit PDB export is placeholder-only**.
   - `LIB/ParallelCampaign.cpp` has a TODO stating actual docked pose coordinates are not yet written; current output is REMARK-only metadata.

4. **Ligand reader bridge work is pending**.
   - `LIB/ProcessLigand/ProcessLigand.cpp` has two TODOs for bridging existing FlexAID readers to `BonMol` conversion:
     - `SdfReader -> BonMol`
     - `Mol2Reader -> BonMol`

5. **Direct-input header still marks Phase-2 TODO context**.
   - `LIB/direct_input.h` explicitly states this area as a "Phase 2" TODO implementation track.

## 2) Roadmap work still marked active/planned

1. **Main implementation roadmap** (`docs/IMPLEMENTATION_ROADMAP.md`)
   - Global status still marks **Phase 4 active**.
   - Phase 4 objective includes **Voronoi hydration** as planned work.
   - Timeline section still shows pending Phase-4 milestones in April/May 2026.

2. **Windows roadmap** (`docs/WINDOWS_BUILD_ROADMAP.md`)
   - File status says only Phase 1 is implemented, and explicitly labels later phases as future:
     - OpenMP enablement strategy on Windows
     - vcpkg manifest integration
     - CMake presets
     - Windows wheel distribution
     - Windows ARM64 runner support

## 3) Known bug/risk backlog documented in security/testing reports

1. **Security hardening backlog remains substantial**.
   - The retired `SECURITY_AUDIT_BUFFER_OVERFLOW.md` (removed from publication snapshot; see git history) reported:
     - 28 distinct vulnerability patterns (25+ files)
     - 7 high-severity findings
     - 14 medium-severity findings
   - It repeatedly called out unsafe string handling patterns (`strcpy`, `strcat`, `sprintf`, unbounded `sscanf`).

2. **Coverage gap backlog exists**.
   - `docs/test-coverage-analysis.md` highlights missing tests for:
     - ring conformer/sugar pucker logic
     - chirality classification
     - PyMOL visualization fallback behavior
     - ENCoM mode-file parsing and malformed-input handling

## 4) Suggested execution order

A practical short-term order based on engineering risk:

1. **Security fixes first**: address high-severity unsafe string operations and path-length overflow risks.
2. **Finish real docking plumbing**: `DatasetRunner` + `ParallelCampaign` GA integration and real pose output.
3. **Close parser/reader bridging**: `ProcessLigand` SDF/MOL2 bridge TODOs.
4. **Increase test depth on critical chemistry features**: chirality, conformers, ENCoM file parsing.
5. **Complete remaining roadmap tracks**: Phase 4 Voronoi hydration + Windows future phases.

