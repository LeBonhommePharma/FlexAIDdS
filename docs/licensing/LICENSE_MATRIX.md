# License Matrix

Overview of components used in FlexAID∆S and their compatibility with Apache-2.0.

| Component | License | Role | Apache-2.0 Compatible? |
|-----------|---------|------|------------------------|
| FlexAID Core | Apache-2.0 | Docking engine foundation | Yes |
| `LIB/PoseBust` (NativePoseQC + BustCli) | Apache-2.0 | In-tree pose QC + optional `bust` bridge | Yes (first-party) |
| [PoseBust](https://github.com/LeBonhommePharma/PoseBust) sibling repo | Apache-2.0 | Standalone C++26 pose validation package | Yes (separate repo) |
| PoseBusters (`bust` CLI) | BSD | Authoritative `pb_pass` / S2 claim gate | Yes (optional subprocess; not vendored) |
| RDKit | BSD-3-Clause | Cheminformatics toolkit | Yes |
| Eigen | MPL-2.0 | Linear algebra (header-only) | Yes |
| PyMOL | PSF | Visualization (optional) | Yes |
| OpenMP runtime | Various + exceptions | Parallelization | Yes (with exceptions) |
| CUDA Toolkit | NVIDIA EULA | GPU backend (optional) | Yes (not OSS) |
| Metal | Apple SDK | GPU backend (optional) | Yes (not OSS) |
| NRGRank | GPL-3.0 | **Inspiration only** | Not used as dependency |

**Claim language:** NativePoseQC is diagnostic only. Official PoseBusters pass requires the upstream BSD `bust` tool when reporting S2 / `success_pb`.

For full details and legal text, see root [`THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md), [`NOTICE`](../../NOTICE), and [`clean-room-policy.md`](clean-room-policy.md).
