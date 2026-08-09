# FlexAIDdS Swift package

Native Swift wrapper over the FlexAID entropy-driven docking engine.

## Building the tests

`Sources/FlexAIDCore` is a bridge, not an implementation. It calls the real
engine — `statmech::StatMechEngine`, `BindingMode`/`BindingPopulation`, ENCoM,
tENCoM, the Shannon stack, the GA, `read_input`, `ic2cf` — which lives in
`../LIB` and is compiled by the root `CMakeLists.txt`.

SwiftPM cannot compile sources outside the package root, so those translation
units cannot be added to a SwiftPM target. Build the genuine CMake product and
`Package.swift` will link it:

```bash
swift/scripts/build-core-archive.sh     # writes swift/.build/cxxcore/swiftlink
cd swift && swift test
```

Point `FLEXAIDDS_CORE_LIB_DIR` at a different directory containing
`libflexaid_core.a` (plus an optional `flexaid_core.link` flag file) to link an
archive produced elsewhere.

Without the archive, library targets still build:

```bash
cd swift && swift build --target FlexAIDdS
cd swift && swift build --target Intelligence
```

…but the XCTest bundle fails to link with `Undefined symbols` naming the real
engine functions. That failure is intentional. The package must never ship stub
targets, synthesized symbols, or `-undefined dynamic_lookup` to make the link
succeed: a suite that links green against fabricated symbols reports a
fabricated scientific result.

The build script refuses to write into `<repo>/build`, which is reserved for
running benchmark campaigns.

## Scientific claim firewall

`Sources/FlexAIDdS/ScientificProvenance.swift` mirrors the C++ contract in
`LIB/statmech.h`. C++ is the single source of truth: the enum spellings,
schema version, `sha256:<64 hex>` receipt syntax and claim predicates are
copied, never re-derived.

Presentation rules enforced across `Sources/Intelligence` and
`Sources/FleetScheduler`:

- kcal/mol, ΔG, ΔΔG, affinity, potency and Kd/Ki wording require a provenance
  record that authorizes the corresponding claim.
- Medicinal-chemistry, SAR, lead-optimization and druggability guidance require
  `binding_physical`.
- Missing, absent or malformed provenance fails closed to `proxy_only`.
- Retained numeric compatibility fields (for example
  `CrossPlatformSelectivityAnalysis.deltaG`) travel with an explicit claim
  marker that defaults to `proxy_only`, including when decoded from payloads
  written before the marker existed.
- Bridge records carry the C++ engine's own provenance. Swift never
  manufactures a stand-in record for arbitrary or calibrated inputs.
