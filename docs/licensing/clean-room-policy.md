# Clean-Room Licensing & GPL-Avoidance Policy

**Project:** FlexAID∆S
**Primary License:** Apache License 2.0
**Maintainer:** Louis-Philippe Morency (Le Bonhomme Pharma)

FlexAID∆S is designed to stay Apache-2.0 while stealing every useful *idea* from the literature and GPL ecosystems without importing a single GPL byte. Methods are fair game; code is not.

## 1. Core Principles

- **Ideas vs. expression.** Algorithms, equations, and screening workflows described in papers (FlexAID, NRGRank, ENCoM, Boltz-ABFE, etc.) can be reimplemented from scratch under Apache-2.0. Only the original code expression is GPL.
- **One-way compatibility.** Apache code can be embedded into GPL projects, but GPL code cannot flow back into Apache without relicensing everything as GPL. FlexAID∆S refuses that collapse.
- **No viral links.** The codebase never imports, links against, or vendors GPL/AGPL libraries. GPL tools may be *run by the user* as external programs, but they are not part of FlexAID∆S.

## 2. Relationship to Upstream Projects

### 2.1 FlexAID (Apache-2.0)

FlexAID provides the original genetic algorithm search, NATURaL scoring matrix, and side-chain flexibility framework under Apache-2.0. FlexAID∆S extends this with explicit thermodynamic ensembles, entropy analysis, and modern acceleration, all under the same permissive license.

### 2.2 NRGRank (GPL-3.0)

NRGRank is GPL-3.0 and **is not a dependency** of FlexAID∆S.

- NRGRank's cube screening / ultra-HTS approach is treated as **scientific prior art**, not as software we redistribute.
- Any similar functionality in FlexAID∆S is reimplemented from published equations and independent derivations.
- No NRGRank source files, wheels, or modules are imported, copied, vendored, or linked.

Result: NRGRank's GPL terms do **not** propagate into FlexAID∆S. Apache-2.0 remains the governing license.

## 3. Allowed vs Forbidden Licenses

### 3.1 Allowed

Direct or transitive dependencies must be under:

- Apache-2.0
- BSD-2-Clause / BSD-3-Clause
- MIT
- MPL-2.0 (file-level copyleft only)
- PSF / Python-2.0
- System runtimes with permissive redistributable terms (e.g., CUDA Toolkit, Metal SDK, OpenMP runtimes)

Current examples include RDKit (BSD-3), Eigen (MPL-2.0, header-only), PyMOL (PSF), OpenMP runtime via compilers, CUDA and Metal as user-installed toolkits.

**Pose validation:** `LIB/PoseBust` is first-party Apache-2.0 (clean-room NativePoseQC). The optional official PoseBusters `bust` CLI (BSD) may be invoked as a subprocess for claim-ready `pb_pass`; it is not vendored and its source must not be copied into this tree. A standalone sibling package is published as [PoseBust](https://github.com/LeBonhommePharma/PoseBust) under Apache-2.0.

### 3.2 Forbidden

Not allowed as dependencies for code shipped in this repo:

- GPL-2.0, GPL-3.0 (without explicit linking exceptions compatible with Apache-2.0)
- AGPL (any version)
- Any license that would force FlexAID∆S itself to become GPL/AGPL or to disclose private downstream modifications

If you are unsure, open an issue and assume "forbidden" until cleared.

## 4. Contributor Rules

When you open a PR, you agree to these clean-room constraints:

1. **No GPL/AGPL code.** Do not copy, translate, or "lightly edit" code from GPL/AGPL projects. If you looked at GPL sources implementing the same idea, say so in the PR for review.
2. **Use the literature, not their repo.** Implement from papers, preprints, or high-level docs. Don't treat someone else's GitHub as a spec.
3. **Document new dependencies.** Any new library must be added to `THIRD_PARTY_LICENSES.md` with license and usage notes.
4. **Keep boundaries loose.** If users want to chain FlexAID∆S with GPL tools, that's their business. FlexAID∆S itself will not link against them.

## 5. User Rights Under Apache-2.0

With FlexAID∆S you can:

- Use it in academic, non-profit, or commercial settings.
- Modify the code and keep your modifications private.
- Embed it in closed-source products as long as you respect Apache-2.0 notice requirements.

You are **not** forced to open-source your pipelines just because you dock drugs against receptors.

## 6. Governance

Licensing and dependency questions should be filed as GitHub issues tagged `licensing`. Maintainers may reject contributions that introduce GPL/AGPL risk or muddy the Apache-2.0 story.

Think of this policy as license thermodynamics: the goal is minimum free energy—maximum freedom—for people doing serious drug work without getting trapped in viral license states.
