# Baseline Build Audit

Date: 2026-05-27

Commit: `ba0f0455e7ba530ac282f3f294d71754d11b6017`

Branch at audit start: `master`

Dirty files present before thermodynamic roadmap edits:

- `python/setup.py`
- `scripts/validate_sources.py`

Environment:

- OS: Darwin `25.5.0` arm64 (`LPmore.local`)
- Compiler: Apple clang `21.0.0 (clang-2100.1.1.101)`
- CMake: `4.3.2`
- OpenMP: not found by CMake
- Eigen3: found via pkg-config

## Commands Run Before Modifications

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake -S . -B build-test -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-test --parallel
```

Initial unprivileged configure failed because the confirmed checkout is outside
the session writable root and CMake could not create `build/CMakeFiles`.
The same configure command succeeded after write access was granted.

## Baseline Status

| Command | Status | Notes |
| --- | --- | --- |
| `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` | Pass | Configure generated build files. |
| `cmake --build build --parallel` | Fail | Pre-existing compile failure in `flexaid_core`. |
| `cmake -S . -B build-test -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release` | Pass | Configure generated test build files. |
| `cmake --build build-test --parallel` | Fail | Same pre-existing compile failure in `flexaid_core`. |
| `ctest --test-dir build-test --output-on-failure` | Not run | Test binaries were not fully built. |

## Pre-existing Build Failures

The baseline full build failed before roadmap code edits. Representative errors:

```text
LIB/spfunction.cpp:196:20: error: no member named 'energy' in 'FA_Global_struct'
LIB/cffunction.cpp:56:14: error: no member named 'bondlist' in 'FA_Global_struct'
LIB/cffunction.cpp:218:6: error: no member named 'nor' in 'cf_str'
LIB/BinarySnapshot.cpp:566:32: error: member access into incomplete type 'FA_Global_struct'
```

These failures are unrelated to the thermodynamic ledger implementation and
pre-existed before modifications on `feature/thermo-ledger`.

## RAM Policy For Follow-up Verification

Full builds must use constrained parallelism after this audit, for example:

```bash
cmake --build build-test --parallel 2 --target test_thermo_ledger
ctest --test-dir build-test --output-on-failure -R ThermoLedgerTests
```

