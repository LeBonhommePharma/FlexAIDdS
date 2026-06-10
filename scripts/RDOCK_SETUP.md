# rDock setup for the Astex-85 head-to-head (macOS arm64)

Verified working build recipe for [CBDD/rDock](https://github.com/CBDD/rDock)
on Apple Silicon (macOS 26 "Tahoe", Apple Clang / libc++). rDock is **not**
available as a prebuilt conda/brew package for osx-arm64 (bioconda ships
linux-64 only), so it must be built from source.

> **Licensing note**: rDock is **LGPL-3.0**. It is used here only as an
> *external, standalone benchmark binary* invoked as a separate process from
> `~/Software/rDock`. No rDock source is vendored into or linked against
> FlexAIDdS, so the project's "no GPL/AGPL dependencies" rule is not violated.
> Do **not** copy rDock code into this repo.

## 1. Dependencies (Homebrew)

```bash
brew install gsl cppunit popt      # popt is the only hard runtime dep
```

cmake is not needed — rDock uses a hand-written Makefile.

## 2. Clone

```bash
git clone --depth 1 https://github.com/CBDD/rDock.git ~/Software/rDock
cd ~/Software/rDock
```

## 3. Three source patches required on modern macOS/libc++

The upstream tree does not build as-is against Apple's current libc++. All three
are mechanical, behaviour-preserving fixes:

1. **`include/VERSION` shadows libc++ `<version>`** (case-insensitive APFS).
   `-I./include` is searched first, so libc++'s `#include <version>` picks up
   rDock's `VERSION` header → an avalanche of bogus `initializer_list` /
   `streampos` / `pos_type` errors. Rename it and fix its one reference:
   ```bash
   mv include/VERSION include/RbtVersion.h
   sed -i '' 's|#include "VERSION"|#include "RbtVersion.h"|' include/RbtResources.h
   ```

2. **`std::_Ios_Openmode` is libstdc++-only** (does not exist in libc++).
   Replace with the portable `std::ios_base::openmode` in three places:
   ```bash
   sed -i '' 's/std::_Ios_Openmode/std::ios_base::openmode/g' \
     include/RbtMOEGrid.h src/lib/RbtBaseFileSink.cxx src/lib/RbtMOEGrid.cxx
   ```

3. **Malformed cxxopts float defaults in `rbcavity`** (`"5.0f"`, `"8.0f"`).
   cxxopts can't parse the trailing `f`, so `rbcavity` aborts on *every*
   invocation. Strip the suffix:
   ```bash
   sed -i '' 's/"5\.0f"/"5.0"/; s/"8\.0f"/"8.0"/' \
     src/lib/rbcavity/rbcavity_argparser.cxx
   ```

## 4. Build

The C++17-removed functional helpers (`std::unary_function`, `bind2nd`,
`ptr_fun`) are re-enabled via libc++ escape-hatch defines:

```bash
make build -j4 CXX=clang++ \
  CXX_EXTRA_FLAGS="-D_LIBCPP_ENABLE_CXX17_REMOVED_FEATURES \
                   -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION \
                   -D_LIBCPP_ENABLE_CXX17_REMOVED_BINDERS"
```

Produces `lib/libRbt.so`, `bin/rbdock`, `bin/rbcavity` (+ rbmoegrid, rblist,
rbcalcgrid).

## 5. Environment (needed every session before running)

```bash
export RBT_ROOT=~/Software/rDock
export PATH="$RBT_ROOT/bin:$PATH"
export DYLD_LIBRARY_PATH="$RBT_ROOT/lib:$DYLD_LIBRARY_PATH"
```

Verify:
```bash
rbdock        # prints usage + "RBT_ROOT: .../rDock"
rbcavity -r <system.prm> -W -d    # carves the cavity .as file
```

## 6. Run the benchmark

```bash
bash scripts/run_rdock_astex.sh                 # all 85 complexes
python scripts/parse_rdock_results.py           # score + comparison CSV
```

Smoke-tested end-to-end on 1G9V: cavity 1843 Å³ centred on the crystal ligand,
10 poses emitted, top-1 Hungarian RMSD 3.42 Å (rDock score −36.85).
