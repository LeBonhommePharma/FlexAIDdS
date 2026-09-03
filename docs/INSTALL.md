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
| `-DFLEXAIDS_GIT_COMMIT_OVERRIDE=<sha>` | packager passes it (Docker; the Homebrew stable formula once its tag carries this block) | `2` (unknown) |
| git checkout | normal clone / `--HEAD` install | `0` or `1`, measured |
| `.git_archival.txt` | release tarball — populated by `export-subst` at `git archive` time | `2` (commit known, tree state not) |
| none of the above | — | `2`, plus a loud CMake warning |

`git_dirty` is **tri-state**: `0` clean, `1` dirty, `2` unknown. It never
reports `0` on a guess. `0` is reachable only from a live git checkout whose
`git status` ran and came back empty. A release tarball is immutable, but the
directory you extract it into is not, and the tarball carries no file-hash
manifest that could detect an edit — so an archive build reports `2`, not `0`.

Every build also writes `flexaidds-build-provenance.json` into the build tree,
and each packaging route installs it next to the binary:

```bash
cat build/flexaidds-build-provenance.json
flexaidds-buildinfo          # Homebrew installs
cat /opt/flexaidds/flexaidds-build-provenance.json   # Docker image
```

### Asking the binary itself

Builds from `main` after the v2.2.0 tag (commit 087e8434 onward; no tagged
release carries it yet) answer `FlexAIDdS --version` with one `key=value` per
line:

```
name=FlexAIDdS
version=2.0.3
git_commit=3f2ad344
git_commit_form=short
git_dirty=0
git_dirty_meaning=0=clean 1=dirty 2=unknown
src_provenance=git
build_type=Release
compiler=AppleClang 21.0.0.21000333
built_utc=2026-09-03T14:02:11Z
```

The JSON and `--version` are stamped from the same CMake configure — one
version string (parsed from `python/flexaidds/__version__.py`, the single
source of truth), one timestamp, one commit/dirty/provenance triple — so they
must agree. If they do not, the file and the binary describe different builds.

Tagged tarballs v2.0.3 through v2.2.0 predate `--version`, the provenance
JSON and `.git_archival.txt` alike: their own `CMakeLists.txt` governs the
build, and it neither reads the override nor writes the JSON. For those builds
the only identity source is the `REMARK` block of an emitted pose file.

The Homebrew formula passes the release tag as an override on its stable
path. That takes effect only for a tag cut after this block landed; against
v2.0.3 through v2.2.0 the override is ignored, no JSON is produced, and the
formula's own provenance check refuses to install. The stable formula path
therefore needs a newer tag; `--HEAD` installs are unaffected.

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
