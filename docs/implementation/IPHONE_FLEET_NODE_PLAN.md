# iPhone 18 Pro as a Bonhomme Fleet node — implementation plan

**Status:** proposal, no code landed. Every chunk below is gated; nothing here changes
pose ranking, election, or claim output on the Mac path.
**Governs:** `swift/Sources/FleetScheduler/`, `swift/Sources/FlexAIDCore/`, the five
Metal bridges under `LIB/`, `CMakeLists.txt` metallib handling, and `RUN_RECEIPT`.
**Defers to:** `AGENTS.md` (conduct), `METHODOLOGY.md` (parity / determinism gates),
`docs/implementation/GPU_FITNESS_DECISION.md` (GA fitness stays on CPU).

## 0. The one design decision

The phone is a **Metal post-process helper**, never a GA cell owner.

FlexAID's hot path is CPU Voronoi (`LIB/Vcontacts.cpp` via `ic2cf` → `vcfunction`),
and `UnifiedHardwareDispatch::best_backend(KernelType::FITNESS_EVAL)` returns the CPU
backend even under a GPU override (`LIB/UnifiedHardwareDispatch.cpp`, pinned by
`FitnessEvalDefaultsToCpuBackend` / `GpuOverrideDoesNotApplyToFitnessEval` in
`tests/test_unified_dispatch.cpp`). Sending an Astex cell to a 2-P-core phone fights
that decision and loses. What the phone *can* run is every kernel the repo already
ships on Metal:

| Kernel | Bridge | Output the Mac needs | FP32-sensitive? |
|---|---|---|---|
| Shannon histogram / Boltzmann / LSE | `LIB/ShannonThermoStack/ShannonMetalBridge.h` | bin counts, weights, log-sum-exp | **yes** (sums) — see §5 |
| Pairwise RMSD matrix | `LIB/MetalRMSDBridge.h` | `n_conf × n_conf` float | low (clustering thresholds) |
| SURFNET cavity spheres | `LIB/CavityDetect/CavityDetectMetalBridge.h` | sphere centers + radii | no (geometry) |
| tENCoM contacts + Hessian assembly | `LIB/tENCoM/tencm_metal.h` | `contacts_ij`, `k`, `r0`, `H[M×M]` | contacts no; Hessian **yes** |
| GPU FastOPTICS kNN | `LIB/gpu_fast_optics_metal.h` | neighbour lists | low |
| TurboQuant | `LIB/TurboQuantMetalBridge.h` | quantised scores | diagnostic only |

Two operating modes share one node app and one work format (`WorkChunk`,
`ChunkResult` in `swift/Sources/FleetScheduler/WorkChunk.swift`):

- **Fleet mode** — chunks travel through the encrypted iCloud drop-box the scheduler
  already uses (`FleetScheduler.submitToiCloud`, `iCloudWatcher`).
- **Tethered mode** — phone on USB-C to the Mac. Transport is Ethernet-over-USB via
  `Network.framework` with Bonjour discovery; iCloud is not touched. Power comes from
  the Mac, so `DeviceCapability.isCharging == true` and the battery multiplier is 1.0.

Vapor chamber note: Apple ships the same next-generation chamber, attached directly
to A20 Pro, in **both** iPhone 18 Pro and 18 Pro Max, with one published sustained
performance claim for both sizes. Nothing in this plan is Pro Max–only.

## 1. Constraint: hot path is CPU, 6 cores (2 P-cores) vs 11 (5 P-cores)

**Mechanism**

1. Add `WorkKind` to `WorkChunk`:
   `case gaSearch, shannonHistogram, boltzmannWeights, pairwiseRMSD, cavityDetect,
   tencomContacts, tencomHessian, fastOpticsKNN`.
   `gaSearch` keeps today's semantics; everything else is a bounded kernel with
   serialised inputs.
2. Extend `DeviceCapability` (`DeviceCapability.swift`):
   - `performanceCoreCount` from `sysctlbyname("hw.perflevel0.physicalcpu")`
     (macOS and iOS both expose it on Apple silicon).
   - `deviceClass: .mac | .iPad | .iPhone`.
   - `estimateTFLOPS`: add `iPhone18 → A19 Pro` and `iPhone19 → A20 Pro` rows. Today
     the switch stops at `case 17...` and under-weights every phone since 2025.
3. Admission rule in `FleetScheduler.splitWork` and `claimChunk`:
   `deviceClass == .iPhone || performanceCoreCount < 4` ⇒ eligible only for
   non-`gaSearch` kinds. Macs remain the only GA owners.
4. If a GA chunk is ever forced onto a phone for experiments, `FXGA` must set
   `omp_set_num_threads(performanceCoreCount)` and QoS `.userInitiated`. This is an
   experiment flag (`FLEXAIDDS_PHONE_GA_EXPERIMENT=1`), default off, never a claim
   path.

**Gate G1** — `swift test --filter FleetSchedulerTests`:
- `testPhoneNeverReceivesGASearchChunk` — a fleet of one Mac and one `iPhone19,2`
  splits a `gaSearch` job to the Mac only; post-process kinds go to both, weighted.
- `testA20ProEstimatedTFLOPSExceedsA17Pro` — regression on the model table.
- Existing `testSplitWorkProportional` and `testCriticalThermalDevicesExcluded` keep
  passing unchanged.

## 2. Constraint: sustained thermals on a long serial queue

**Mechanism**

1. **Short chunks by construction.** Every non-GA kind has a wall budget of seconds to
   a few minutes. `WorkChunk.timeoutSeconds` default drops from 3600 to 300 for
   phone-eligible kinds. There is no multi-hour cell on the phone to throttle.
2. **Live thermal telemetry, not a snapshot.** `DeviceCapability.current()` reads
   `ProcessInfo.thermalState` once at split time. The node runtime observes
   `ProcessInfo.thermalStateDidChangeNotification` and republishes capability on
   every change. The scheduler already maps nominal/fair/serious/critical to
   1.0 / 0.75 / 0.4 / 0.0; the plan makes that mapping *continuous* by re-running
   admission on the next claim, not only at job creation.
3. **Yield at `.serious`.** A node finishing a chunk in `.serious` does not claim the
   next one; it waits for `.fair` or better. At `.critical` it releases its lease so
   `sweepOrphanedChunks` reassigns it (retry path already exists, with priority
   escalation).
4. **Throughput EWMA.** `deviceThroughput` in `FleetScheduler.submitResult` is a
   last-value overwrite. Replace with an exponentially weighted mean (α = 0.3) so a
   throttled phone is downweighted smoothly instead of oscillating.
5. **Duty cycle knob.** `FLEXAIDDS_NODE_DUTY_CYCLE=0.7` caps busy fraction per minute
   on phones; default 1.0 on Macs.

**Gate G2** — on-device, plugged in, screen on at minimum brightness:
- 30-minute loop of alternating `pairwiseRMSD` (10k × 10k, 30 atoms) and
  `shannonHistogram` (1M energies) chunks.
- Log `thermalState`, chunk wall time, and `ChunkResult.thermalStateAtCompletion`.
- Accept if throughput at minute 30 ≥ 60% of minute 1, `.critical` never observed,
  and no chunk is orphaned. Record the CSV under the node's `RUN_RECEIPT`
  (`fleet_nodes[].thermal_trace`).

## 3. Constraint: 12 GB unified memory vs 18 GB

**Mechanism**

1. `WorkChunk.estimatedPeakBytes` — computed by the chunk producer from inputs:
   - RMSD: `n_conf × n_conf × 4` for the matrix plus one coordinate tile.
     `MetalRMSDBridge.mm` already tiles over row ranges to keep GPU memory bounded;
     the producer must tile so that a single tile ≤ 512 MB.
   - Shannon: `n_energies × 4 + bins × 4`.
   - tENCoM Hessian: `(N−1)² × 8`. Cap N at 8 000 residues on phones (≈ 512 MB).
   - Cavity: grid cells × sphere structs; producer reports the grid extent.
2. Hard admission cap on phones: `estimatedPeakBytes ≤ 2 GiB`, checked against
   `os_proc_available_memory()` at claim time (iOS API; returns the jetsam headroom
   for the current process). Chunks above the cap are not offered to phones at all.
3. All Metal buffers on phone use `MTLResourceStorageModeShared` (already the pattern
   in `gpu_fast_optics_metal.mm`); no private-mode copies.
4. No GA on phone means no per-thread Vcontacts/`ca_rec` clones
   (`LIB/gaboom.cpp` OpenMP path) — the main memory multiplier on the Mac never
   exists on the phone.

**Gate G3**
- XCTest: a `DeviceCapability` with `availableMemoryGB: 12` and `deviceClass: .iPhone`
  rejects a chunk with `estimatedPeakBytes = 6 GiB` and accepts one at 1.5 GiB.
- On-device: run the G2 loop under Xcode's memory gauge; peak resident < 2 GiB, zero
  `EXC_RESOURCE` / jetsam events in the device console.

## 4. Constraint: iOS, not macOS — no engine binary, background suspension

**Mechanism**

1. **iOS slice of `flexaid_core`.** Add a CMake preset:
   `-DCMAKE_SYSTEM_NAME=iOS -DCMAKE_OSX_ARCHITECTURES=arm64
   -DCMAKE_OSX_DEPLOYMENT_TARGET=17.0 -DFLEXAIDS_USE_METAL=ON
   -DFLEXAIDS_USE_OPENMP=OFF -DFLEXAIDS_BUILD_CORE=ON` producing a static
   `libflexaid_core.a` (arm64-ios). OpenMP is off because Apple's iOS toolchain ships
   no `libomp` and the node runs no GA. `swift/scripts/build-core-archive.sh` gains
   `--platform ios` and writes an `.xcframework` with macOS + iOS slices so
   `Package.swift` (already `.iOS(.v17)`) links the right one via
   `FLEXAIDDS_CORE_LIB_DIR`.
2. **Metallib resolution.** Every bridge today loads a compile-time absolute path
   (`SHANNON_METALLIB_PATH` etc., set in `CMakeLists.txt` from
   `CMAKE_CURRENT_BINARY_DIR`) then falls back to `newDefaultLibrary`. Add one shared
   helper `flexaids_metallib_url(const char* name)` used by all five bridges, in order:
   `FLEXAIDDS_METALLIB_DIR` env override → `[NSBundle mainBundle]` resource → compile-time
   path → `newDefaultLibrary`. On iOS the `.metal` sources compile into the app's
   `default.metallib` through Xcode, so the last step succeeds without any absolute
   path. This also fixes the Homebrew "metallib not installed" caveat noted in
   `docs/audit/26h-swarm/1dba43f4b.md`.
3. **Foreground compute contract.** iOS will not run unbounded background compute.
   The node app is a single-screen SwiftUI "compute" view that sets
   `UIApplication.shared.isIdleTimerDisabled = true`, dims to minimum brightness, and
   documents Guided Access (triple-click) as the `caffeinate` equivalent. Opportunistic
   background work uses `BGProcessingTaskRequest` with
   `requiresExternalPower = true`, `requiresNetworkConnectivity = true`, scoped to
   chunks whose `estimatedWallSeconds ≤ 60`. Anything longer needs foreground.
4. **Tethered transport.** `NWListener` on the phone advertising Bonjour
   `_flexaidds-node._tcp` with `includePeerToPeer = true`; the Mac coordinator browses,
   prefers a path whose interface is `.wiredEthernet` (the iPhone-over-USB NCM link),
   and falls back to Wi-Fi, then to the iCloud drop-box. Frames are length-prefixed
   `WorkChunk.encrypt(using:)` payloads (ChaChaPoly, already implemented). Key exchange
   is a one-time QR code shown on the phone, scanned by the Mac.
5. **Node identity.** `DeviceCapability.modelIdentifier()` returns `hostName` on iOS;
   switch to `sysctlbyname("hw.machine")` (e.g. `iPhone19,2`) so the TFLOPS table and
   receipts see the real model.

**Gate G4**
- Hosted CI job `ios_core_compile_smoke` (sibling of `macos_metal_compile_smoke` in
  `.github/workflows/ci.yml`): configure + build the iOS preset on `macos-latest`;
  `nm -gU libflexaid_core.a | grep -E 'metal_eval_get_capabilities|ShannonMetalBridge|metal_rmsd'`
  must hit all three.
- Simulator XCTest: loopback `NWListener`/`NWConnection` round-trips a
  `shannonHistogram` chunk and returns a decodable `ChunkResult`.
- On-device: phone on USB-C, Mac sends one RMSD chunk, `ChunkResult` arrives in
  < 2 s end-to-end; path reports `.wiredEthernet`.

## 5. Constraint: not claim-safe

**Mechanism**

1. **No fitness on phone, ever.** The node build does not link `metal_eval.mm`'s
   batch entry points into any reachable path; `WorkKind.gaSearch` is refused by
   admission (§1). `GPU_FITNESS_DECISION.md` stays authoritative.
2. **FP64 stays on the Mac.** `ShannonMetalBridge.mm` computes min/max, binning,
   `exp` and sums in `float`; `CODEBASE_REVIEW_2026-08.md` flags that an Apple run
   and a Linux run report different `log_Z`/entropy, and `ParallelDock.cpp` registers
   that `log_Z` into the grand-canonical ledger. The phone therefore returns only
   **FP32-insensitive products**:
   - Shannon: integer bin counts against **Mac-supplied bin edges** (not phone-computed
     min/max). Counts are exact; the Mac computes `H` in double.
   - Boltzmann / LSE: **not delegated**. The Mac keeps these on its FP64 CPU path.
   - RMSD: float matrix, consumed only for clustering thresholds; parity-gated.
   - tENCoM: contacts (ints + float `k`, `r0`) delegated; Hessian assembly and
     diagonalisation stay on the Mac in double.
   - Cavity spheres: geometry only, never enters CF or the ledger.
3. **Provenance on every result.** `ChunkResult` gains
   `backendProvenance: {deviceModel, osVersion, backend: "metal", precision: "fp32",
   kernel, metallibSHA256}`. `FleetAggregator` sets `ScientificProvenance` to
   `proxyOnly` for any mode whose inputs include a phone-computed value that has not
   passed the parity gate below. `swift/README.md`'s claim firewall already fails closed
   on missing provenance; this extends it to "computed on a phone".
4. **Receipts.** `RUN_RECEIPT` gains `fleet_nodes[]` (model, OS, metallib SHA, kernels,
   chunk count, thermal trace path). `resolve_build.py --check` learns to pin the iOS
   archive SHA256 alongside the engine binary.
5. **Language.** Node output is labelled "Metal helper kernels (post-process)". It is
   never "docking", never "success", never a benchmark cell. `AGENTS.md`'s
   deception-proof contract applies unchanged.

**Gate G5** — offline Astex parity, per `METHODOLOGY.md` parity section:
- For ≥ 20 receipted Astex OUT directories: Shannon bin counts from phone Metal are
  **bit-identical** to CPU counts given identical edges; RMSD matrix
  `max |Δ| < 1e-3 Å`; FOPTICS/density-peak cluster labels **identical**; elected
  pose and `result.csv` **bit-identical** to the Mac-only run.
- Any mismatch fails the gate; the offending kind stays `proxyOnly` and is dropped
  from phone eligibility.
- XCTest: a `ChunkResult` with phone provenance and no parity attestation aggregates
  to `claimValidity == .proxyOnly`.

## 6. Chunk order and dependencies

| Chunk | Scope | Gate | Needs hardware? |
|---|---|---|---|
| C0 | `WorkKind`, `DeviceCapability` P-core/model/A20 rows, admission rules, EWMA, memory cap | G1, G3 (XCTest) | Mac only |
| C1 | iOS CMake preset, `.xcframework`, `flexaids_metallib_url` in all five bridges | G4 CI smoke | Mac with Xcode |
| C2 | Node app: compute screen, thermal observer, `BGProcessingTask`, telemetry | G2 | iPhone 18 Pro |
| C3 | Tethered transport (Bonjour over USB NCM), QR key exchange, iCloud fallback | G4 simulator + device | iPhone + Mac |
| C4 | Provenance fields, receipt schema, Astex parity harness | G5 | Mac + phone |
| C5 | End-to-end: one Astex target, Mac GA + phone Shannon/RMSD, `result.csv` diff | G5 bit-identical election | Mac + phone |

C0 and C1 are independent and can land in parallel. C2 and C3 depend on C1. C4 can
start after C0. C5 requires everything.

## 7. What this buys, honestly

- The Mac keeps every OpenMP thread on Voronoi while the phone absorbs Shannon
  histograms, RMSD matrices, cavity spheres and tENCoM contact discovery. Expect a
  measurable but modest wall-time gain on post-GA phases, not a 2× campaign speedup.
- Tethered mode removes iCloud FileProvider from the hot loop entirely, which
  `AGENTS.md` § Benchmark storage already treats as a hazard.
- Nothing on the phone can change an elected pose. If G5 ever shows it can, the kind
  is pulled.

## 8. Explicitly not promised

- GA cells or claim cells on the phone.
- Unattended overnight background compute on iOS.
- Any FP64 thermodynamic reduction on Apple GPUs (Metal has no `double`).
- A Pro Max–specific thermal advantage; Apple publishes one claim for both sizes.
