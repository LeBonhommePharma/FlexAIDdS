# Installing FlexAID∆S

Two separate things ship from this repo, and conflating them is the most common
install failure:

| | what it is | how you get it |
|---|---|---|
| **`FlexAIDdS`, `tENCoM`, `FlexAID`** | the native C++ docking engine | CMake, Homebrew, Docker |
| **`flexaidds` (Python)** | analysis, `load_results`, StatMech, DatasetRunner | `pip`, conda |

The Python package does **not** contain the docking engine, and installing the
engine does not give you the Python package. Most "it installed but nothing
docks" reports are this.

---

## Support matrix

| Path | macOS arm64 | macOS x86_64 | Linux x86_64 | Windows | Ships engine | Ships Python |
|---|---|---|---|---|---|---|
| CMake from source | ✅ | ✅ | ✅ | ⚠️ partial | ✅ | optional |
| Homebrew (`Formula/flexaidds.rb`) | ✅ | ✅ | — | — | ✅ | ❌ |
| Docker (`containers/Dockerfile.locked`) | ✅ via Docker | ✅ | ✅ | ✅ via Docker | ✅ | ❌ |
| Apptainer (`containers/*.def`) | — | — | ✅ | — | ✅ | ❌ |
| `pip` | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| conda (`conda/meta.yaml`) | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| Release binaries (`release.yml`) | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ |

⚠️ **Windows**: MSVC has no `/std:c++26`, so `CMakeLists.txt` caps Windows at
C++20 and only the `_core` Python extension builds correctly. The full engine
requires GCC ≥ 14, Clang ≥ 18, or AppleClang ≥ 16 (Xcode 16). Windows is
excluded from the engine CI matrix on purpose — treat Docker or WSL as the
Windows story.

### Compiler floor

The engine is C++26 (C++23 on AppleClang, which defaults down for SDK
compatibility — override with `-DFLEXAIDS_FORCE_CXX26=ON`). CMake hard-fails
below the floor rather than producing a subtly different binary:

| toolchain | minimum |
|---|---|
| GCC | 14 |
| AppleClang | 16 (Xcode 16) |
| MSVC | 19.40 (VS 2022 17.10) — engine still unsupported |
| CMake | 3.28 |

---

## Engine: CMake from source

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/FlexAIDdS --help
```

Dependencies: CMake ≥ 3.28, Ninja or Make, Eigen 3.4+, and OpenMP
(`libomp` on macOS). The build uses LTO and Unix Makefiles by default; pass
`-G Ninja` if you prefer.

## Engine: Homebrew (lowest friction on macOS)

```bash
brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS
brew trust --formula lebonhommepharma/flexaidds/flexaidds
brew install lebonhommepharma/flexaidds/flexaidds
```

CPU + OpenMP by default. Metal is opt-in and needs a real Metal toolchain:

```bash
brew install --build-from-source --with-metal lebonhommepharma/flexaidds/flexaidds
```

## Engine: Docker (most reproducible)

Always pass the commit — see [Build provenance](#build-provenance).

```bash
docker build -f containers/Dockerfile.locked \
  --build-arg FLEXAIDS_GIT_COMMIT="$(git rev-parse --short HEAD)" \
  -t flexaidds:locked-x86_64 .

docker run --rm -v "$PWD:/work" -w /work flexaidds:locked-x86_64 --help
```

The build **fails** if the resulting image has no recoverable commit. Override
with `--build-arg ALLOW_UNKNOWN_PROVENANCE=1` only for a throwaway image.

## Python: pip

```bash
pip install flexaidds          # once published
# or, from the repo:
pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
flexaidds --help
```

`pip install` is designed never to fail. The accelerated `flexaidds._core`
extension is optional: when pybind11, Eigen, and a C++26 compiler are all
present it is compiled; otherwise the install succeeds with pure-Python
fallbacks and `flexaidds.HAS_CORE_BINDINGS` is `False`. Check which you got:

```bash
python -c "import flexaidds as fd; print(fd.__version__, fd.HAS_CORE_BINDINGS)"
```

Force the pure path (useful in CI and containers): `FLEXAIDDS_SKIP_CORE=1`.
Point at a non-standard Eigen: `EIGEN_INCLUDE_DIR=/path/to/eigen`.

## Python: conda

```bash
conda env create -f environment.yml && conda activate flexaidds
# or build the recipe:
FLEXAIDS_GIT_COMMIT_OVERRIDE=$(git rev-parse --short HEAD) conda build conda/
```

`conda/meta.yaml` derives its version from `python/flexaidds/__version__.py`;
do not hardcode it.

---

## Build provenance

**A binary you cannot identify is a binary whose results you cannot defend.**
The engine writes its commit into the `REMARK` block of every pose file:

```
REMARK FLEXAID.commit=3f2ad344 FLEXAID.dirty=0 FLEXAID.seed=…
```

If that commit is empty or wrong, a relinked engine can silently invalidate a
comparison and the pose files will not show it.

CMake resolves the stamp from four sources, in order:

| source | when | `git_dirty` |
|---|---|---|
| `-DFLEXAIDS_GIT_COMMIT_OVERRIDE=<sha>` | packager passes it (conda, Docker, Homebrew stable) | `2` (unknown) |
| git checkout | normal clone / `--HEAD` install | `0` or `1`, measured |
| `.git_archival.txt` | release tarball — populated by `export-subst` at `git archive` time | `0` (archives are immutable) |
| none of the above | — | `2`, plus a loud CMake warning |

`git_dirty` is **tri-state**: `0` clean, `1` dirty, `2` unknown. It never
reports `0` on a guess. Every build also writes
`flexaidds-build-provenance.json` into the build tree, and each packaging route
installs it next to the binary:

```bash
cat build/flexaidds-build-provenance.json
flexaidds-buildinfo          # Homebrew installs
cat /opt/flexaidds/flexaidds-build-provenance.json   # Docker image
```

### Known gap

The native engine has **no `--version` flag**. Provenance is recoverable from
the provenance JSON, from `strings $(command -v FlexAIDdS) | grep 'FLEXAID.commit'`,
or from any emitted pose file — but not from the CLI. Adding
`FlexAIDdS --version` requires a change to the C++ entry point and is tracked
separately.

Tarballs cut before `.git_archival.txt` existed (**v2.0.3 and earlier**) carry
no commit. The Homebrew formula compensates by passing the release tag as an
override, so those installs report `v2.0.3` rather than nothing.

---

## Reproducibility

Where the packaging format allows a pin, it is pinned:

* `python/pyproject.toml` — static `version`, floors on `numpy`/`scipy`.
* `conda/meta.yaml` — every `run:` dependency carries a constraint; CI fails the
  recipe if a floating dep appears.
* `containers/Dockerfile.locked` — `ARG BASE_IMAGE` for a digest pin; resolve
  with `docker buildx imagetools inspect ubuntu:22.04 --format '{{.Manifest.Digest}}'`.
* `.github/workflows/*` — actions pinned by commit SHA.

## Troubleshooting

**`pip install` succeeded but `HAS_CORE_BINDINGS` is `False`.** Expected unless
you have pybind11 + Eigen + a C++26 compiler. Not an error; the pure-Python
fallbacks are correct, just slower.

**CMake: "requires GCC ≥ 14".** The floor is real. `CXX=g++-14 cmake -S . -B build …`.

**CMake warns "Build provenance is UNKNOWN".** You are building from an
extracted tarball with no `.git_archival.txt`. Pass
`-DFLEXAIDS_GIT_COMMIT_OVERRIDE=<sha>`. Do not ship results from a build that
warns.

**`brew reinstall` fails with "couldn't find remote ref".** An old `--HEAD`
install pinned a deleted branch. See the formula caveats for the reset sequence.
