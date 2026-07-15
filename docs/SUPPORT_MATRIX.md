# Support Matrix

This matrix defines the **supported** combinations for `FlexAIDdS Core 1.0`.

Anything not listed here should be treated as experimental.

## Supported product surfaces

| Surface | Status |
|:--|:--|
| C++ CLI (`FlexAIDdS`) | Supported |
| Python package (`flexaidds`) | Supported |
| JSON configuration + documented CLI workflows | Supported |
| Benchmark runner / reproducibility bundle | Supported |
| Swift packages | Experimental |
| TypeScript / PWA / dashboards | Experimental |
| Bonhomme Fleet / iCloud orchestration | Experimental |
| NATURaL workflows | Experimental |

## Supported platforms

| OS | Compiler / runtime | Status | Notes |
|:--|:--|:--|:--|
| Linux | GCC >= 10 | Supported | Core target for release validation |
| Linux | Clang >= 10 | Supported | Core target for release validation |
| macOS | Apple Clang | Supported | Core CLI and Python surfaces only |
| Windows | MSVC 2022 | Supported | Core CLI and Python surfaces only |

## Backend support tiers

| Backend | Status | Release expectation |
|:--|:--|:--|
| Scalar CPU | Supported | Must build and run |
| OpenMP | Supported | Must build and run on supported matrix where available |
| AVX2 | Supported | Performance optimization, not correctness dependency |
| AVX-512 | Experimental | Optional and architecture-specific |
| CUDA | Experimental | Not required for 1.0 support guarantees |
| Metal | Experimental | Hosted CI: shader compile smoke only. Full link/runtime gate: self-hosted `self-hosted-m3` (`metal-self-hosted.yml`). Not required for 1.0 CPU support guarantees. See `docs/CI_RELEASE_GATES.md`. |
| ROCm/HIP | Experimental | Not required for 1.0 support guarantees |

### Platform CI gates (summary)

| Platform | Blocking PR gate | What is still self-hosted |
|:--|:--|:--|
| Linux CPU | `cxx_core_build` + `ctest` | Full docking campaigns |
| macOS CPU | `macos_cpu_tests` + `ctest` (no soft-fail) | Full docking campaigns |
| Metal GPU | `macos_metal_compile_smoke` (shaders) | Link + runtime via `self-hosted-m3` |

Full table and release checklist: **`docs/CI_RELEASE_GATES.md`**.

## Python support

| Python | Status |
|:--|:--|
| 3.9 | Supported |
| 3.10 | Supported |
| 3.11 | Supported |

## What “supported” means here

A supported combination must have:

- documented installation path
- CI coverage or release validation coverage
- at least one automated test or smoke test
- no known release-blocking issue at the time of release

## What “experimental” means here

Experimental means one or more of the following:

- incomplete CI coverage
- unstable API or UX
- undocumented installation path
- incomplete reproducibility coverage
- performance path exists but is not release-gated
