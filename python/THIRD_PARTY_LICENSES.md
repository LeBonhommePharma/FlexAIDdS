# Third-Party Licenses

**FlexAID∆S** incorporates and depends on the following open-source components:

---

## Direct Dependencies

### RDKit (BSD-3-Clause)

**Project:** Open-source cheminformatics and machine learning toolkit  
**License:** BSD 3-Clause License  
**Source:** https://github.com/rdkit/rdkit  
**Use in FlexAID∆S:** Molecular I/O, SMILES parsing, molecular property calculations

```
Copyright (c) 2006-2024, RDKit Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
```

---

### Eigen (MPL 2.0) — vendored

**Project:** C++ template library for linear algebra  
**Version:** 3.4.0  
**License:** Mozilla Public License 2.0  
**Source:** https://gitlab.com/libeigen/eigen  
**Vendored at:** `LIB/vendor/eigen/` (git submodule, pinned to tag `3.4.0`)  
**Use in FlexAID∆S:** SIMD vectorization, matrix operations, eigenvalue decomposition, numerical stability

Eigen is vendored directly in the repository as a git submodule so that the build works out-of-the-box on any platform without requiring `libeigen3-dev`, `brew install eigen`, or `choco install eigen`. A fresh clone with `git clone --recurse-submodules` is all that is needed.

```
Copyright (C) 2008-2024 Eigen Contributors

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
```

**Note:** MPL 2.0 is a file-level copyleft license. Modified Eigen files must remain MPL 2.0. FlexAID∆S uses Eigen as an unmodified header-only dependency, maintaining compatibility with Apache-2.0.

---

### PyMOL (Python Software Foundation License)

**Project:** Molecular visualization system  
**License:** Python Software Foundation License (permissive)  
**Source:** https://github.com/schrodinger/pymol-open-source  
**Use in FlexAID∆S:** Optional visualization backend, PyMOL plugin integration

```
Copyright (c) 2001-2024 Schrodinger, LLC
SPDX-License-Identifier: Python-2.0

Python Software Foundation License Version 2
```

---

## Build & Runtime Dependencies

### CMake (BSD-3-Clause)

**Project:** Cross-platform build system generator  
**License:** BSD 3-Clause License  
**Source:** https://cmake.org/  
**Use in FlexAID∆S:** Build configuration, hardware detection

---

### OpenMP (Various)

**Standard:** OpenMP API Specification  
**License:** Implementation-specific (GCC: GPL with Runtime Exception, LLVM: Apache-2.0/MIT)  
**Source:** https://www.openmp.org/  
**Use in FlexAID∆S:** CPU parallelization, thread management

**Note:** OpenMP runtime libraries include linking exceptions that permit use in non-GPL software. FlexAID∆S uses OpenMP directives without modifying compiler sources.

---

### CUDA Toolkit (NVIDIA EULA - Proprietary but Free)

**Project:** NVIDIA GPU computing platform  
**License:** NVIDIA End User License Agreement (free for research and commercial use)  
**Source:** https://developer.nvidia.com/cuda-toolkit  
**Use in FlexAID∆S:** Optional GPU acceleration on NVIDIA hardware

**Note:** Not open source, but freely available. Distribution of CUDA binaries restricted by NVIDIA EULA. FlexAID∆S ships CUDA-compatible source code; users install CUDA Toolkit separately.

---

### Metal Framework (Apple Inc. - Free)

**Project:** Apple GPU programming framework  
**License:** Included with macOS SDK (free)  
**Source:** https://developer.apple.com/metal/  
**Use in FlexAID∆S:** Optional GPU acceleration on Apple Silicon (M1/M2/M3)

**Note:** Proprietary but freely available on macOS. FlexAID∆S uses public Metal APIs.

---

## Scientific Methods & Parameters

### FlexAID Core Method (Apache-2.0)

**Original Project:** NRGlab/FlexAID  
**License:** Apache License 2.0  
**Source:** https://github.com/NRGlab/FlexAID  
**Papers:**
- Gaudreault & Najmanovich (2015). *J. Chem. Inf. Model.* 55(7):1323-36.  
  DOI: [10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)

**Use in FlexAID∆S:** Genetic algorithm docking engine, NATURaL scoring function parameters (40-SYBYL atom type interaction matrix), side-chain flexibility framework

**Notes:**
- FlexAID∆S is derived from the Apache-2.0 licensed FlexAID
- Energy parameters (ε_ij values) are reused from original FlexAID
- Original docking engine architecture extended with entropy calculations
- Full Apache-2.0 compatibility maintained

---

### NRGRank Method (Scientific Inspiration Only - GPL-3.0)

**Project:** NRGlab/NRGRank  
**License:** GNU General Public License v3.0  
**Source:** https://github.com/NRGlab/NRGRank (NOT INCLUDED AS DEPENDENCY)  
**Paper:** DesCôteaux T, Mailhot O, Najmanovich RJ. *NRGRank: Coarse-grained structurally-informed ultra-massive virtual screening.* bioRxiv 2025.02.17.638675v2. doi:[10.1101/2025.02.17.638675](https://doi.org/10.1101/2025.02.17.638675)  
**Supplementary data:** Zenodo [10.5281/zenodo.16861024](https://doi.org/10.5281/zenodo.16861024) (CC-BY-4.0)

**Relationship to FlexAID∆S:**
- **Published method only:** reimplemented in `LIB/CoarseScreen.*` under Apache-2.0
- **No code dependency:** NRGlab/NRGRank source is not imported, linked, vendored, or translated
- **Clean-room:** `docs/licensing/clean_room_policy.md`

**Citation in Scientific Work:**
> "Coarse cube screening follows the published NRGRank method (DesCôteaux, Mailhot, Najmanovich, bioRxiv 2025.02.17.638675v2) and was reimplemented from the paper under Apache-2.0; the GPL NRGRank codebase is not a dependency."

---

## License Compatibility Summary

| Component | License | Compatibility with Apache-2.0 | Distribution Status |
|-----------|---------|-------------------------------|---------------------|
| **FlexAID Core** | Apache-2.0 | ✅ Same license | Incorporated |
| **RDKit** | BSD-3-Clause | ✅ Permissive | Dependency |
| **Eigen** | MPL 2.0 | ✅ File-level copyleft | Header-only dependency |
| **PyMOL** | PSF | ✅ Permissive | Optional dependency |
| **OpenMP** | Runtime exception | ✅ Linking exception | Compiler built-in |
| **CUDA** | NVIDIA EULA | ⚠️ Proprietary (free) | User installs separately |
| **Metal** | Apple | ⚠️ Proprietary (free) | macOS built-in |
| **libcurl** | curl (MIT/X-inspired) | ✅ Permissive | Optional dependency (dataset runner) |
| **NRGRank** | GPL-3.0 | ❌ **NOT A DEPENDENCY** | Cited as scientific inspiration only |

**Legend:**
- ✅ Compatible: No restrictions on Apache-2.0 distribution
- ⚠️ Proprietary but free: Not open source, but freely available
- ❌ Avoided: GPL code not included to maintain Apache-2.0 purity

---

## Compliance Notes

### For Users

**No GPL obligations:** FlexAID∆S can be used in proprietary software, modified freely, and relicensed under any terms. No source disclosure required for downstream works.

**Optional dependencies:** CUDA and Metal are optional. FlexAID∆S includes CPU-only fallback implementations.

**No viral licenses:** All included dependencies are permissive or have library exceptions. No strong copyleft propagation.

### For Contributors

**CLA Required:** All contributions must be submitted under Apache License 2.0 terms.

**No GPL code:** Do not submit patches incorporating GPL-licensed code. See [docs/licensing/clean-room-policy.md](docs/licensing/clean-room-policy.md) for detailed GPL avoidance strategy.

**Dependency review:** Proposed new dependencies must use Apache-2.0, BSD, MIT, PSF, or MPL 2.0 licenses. GPL and AGPL forbidden.

### For Redistributors

**Binary distribution:** Include this file and main LICENSE in all binary packages.

**Source distribution:** Preserve all copyright notices in source headers.

**CUDA/Metal:** If distributing binaries with GPU support, note that users must accept NVIDIA EULA (CUDA) or Apple terms (Metal) separately. Do not redistribute CUDA Toolkit itself.

**FlexAID attribution:** Retain Apache-2.0 notice for FlexAID-derived components.

---

## Updates

This file is maintained alongside FlexAID∆S releases:

- **Version 1.0 (March 2026):** Initial comprehensive third-party license documentation
- Future updates will track dependency changes and new integrations

For questions about licensing compliance, see:
- Main LICENSE file (Apache-2.0 full text)
- [docs/licensing/clean-room-policy.md](docs/licensing/clean-room-policy.md) (GPL avoidance strategy)
- GitHub issues: https://github.com/LeBonhommePharma/FlexAIDdS/issues

---

**Maintained by:** Louis-Philippe Morency, PhD (Candidate)  
**Last Updated:** March 7, 2026  
**Repository:** https://github.com/LeBonhommePharma/FlexAIDdS
