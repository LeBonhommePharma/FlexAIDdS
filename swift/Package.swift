// swift-tools-version: 5.9
// Package.swift — FlexAIDdS Swift Package
//
// Native Swift wrapper for the FlexAID entropy-driven molecular docking engine.
// Provides: statistical mechanics, vibrational entropy, fleet scheduling,
// HealthKit, MusicKit, Apple Intelligence integrations.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import PackageDescription
import Foundation

let repositoryRoot = URL(fileURLWithPath: #filePath)
    .deletingLastPathComponent()
    .deletingLastPathComponent()
let coreHeaderPath = repositoryRoot.appendingPathComponent("LIB").path
let tENCoMHeaderPath = repositoryRoot.appendingPathComponent("LIB/tENCoM").path

// ─── real C++ engine linkage ────────────────────────────────────────────────
//
// Sources/FlexAIDCore/*.mm is a bridge, not an implementation: it calls
// statmech::StatMechEngine, BindingMode/BindingPopulation, ENCoM, tENCoM, the
// Shannon stack, GA, read_input and ic2cf, all of which live in <repo>/LIB and
// are built by the root CMakeLists (149 translation units).
//
// SwiftPM cannot compile sources outside the package root, so those units
// cannot be added to a target here. Instead we link the genuine CMake build
// product when it exists. Build it with:
//
//     swift/scripts/build-core-archive.sh
//
// If the archive is absent, the test bundle fails to link with an honest
// "Undefined symbols" error naming the real engine functions. That failure is
// the correct outcome: stub targets, synthesized symbols and
// `-undefined dynamic_lookup` are deliberately NOT used, because a test suite
// that links green against fabricated symbols reports a fabricated result.
let coreArchiveDirectory: String? = {
    let candidate = ProcessInfo.processInfo
        .environment["FLEXAIDDS_CORE_LIB_DIR"]
        .flatMap { $0.isEmpty ? nil : $0 }
        ?? repositoryRoot
            .appendingPathComponent("swift/.build/cxxcore/swiftlink").path
    let archive = candidate + "/libflexaid_core.a"
    return FileManager.default.fileExists(atPath: archive) ? candidate : nil
}()

// Runtime dependencies of the archive itself (OpenMP, ...) as recorded by the
// build script. Read from disk so no machine-specific toolchain path is baked
// into this manifest.
let coreArchiveLinkFlags: [String] = {
    guard let directory = coreArchiveDirectory,
          let contents = try? String(
            contentsOfFile: directory + "/flexaid_core.link", encoding: .utf8)
    else { return [] }
    return contents
        .split(whereSeparator: \.isNewline)
        .map { $0.trimmingCharacters(in: .whitespaces) }
        .filter { !$0.isEmpty }
}()

let coreLinkerSettings: [LinkerSetting] = {
    var settings: [LinkerSetting] = [.linkedLibrary("c++")]
    if let directory = coreArchiveDirectory {
        settings.append(.unsafeFlags(
            ["-L\(directory)", "-lflexaid_core"] + coreArchiveLinkFlags))
    }
    return settings
}()

let package = Package(
    name: "FlexAIDdS",
    platforms: [
        .macOS(.v14),
        .iOS(.v17),
    ],
    products: [
        .library(name: "FlexAIDdS", targets: ["FlexAIDdS"]),
        .library(name: "FleetScheduler", targets: ["FleetScheduler"]),
        .library(name: "HealthIntegration", targets: ["HealthIntegration"]),
        .library(name: "MediaIntegration", targets: ["MediaIntegration"]),
        .library(name: "Intelligence", targets: ["Intelligence"]),
    ],
    targets: [
        // Layer 1: C/Obj-C++ bridge to the C++ core.
        // The engine itself is linked from the CMake archive resolved above
        // (see `coreArchiveDirectory`); this target only holds the bridge.
        .target(
            name: "FlexAIDCore",
            path: "Sources/FlexAIDCore",
            publicHeadersPath: "include",
            cxxSettings: [
                .define("FLEXAIDS_SWIFT_BRIDGE"),
                // SwiftPM forbids headerSearchPath entries outside the package
                // root. Resolve the sibling core from this manifest instead.
                .unsafeFlags([
                    "-std=c++20",
                    "-I", coreHeaderPath,
                    "-I", tENCoMHeaderPath,
                ]),
            ],
            linkerSettings: coreLinkerSettings
        ),

        // Layer 2: Swift module — actors, models, Swift-native API
        .target(
            name: "FlexAIDdS",
            dependencies: ["FlexAIDCore"],
            path: "Sources/FlexAIDdS"
        ),

        // Layer 3: Feature modules
        .target(
            name: "FleetScheduler",
            dependencies: ["FlexAIDdS", "Intelligence"],
            path: "Sources/FleetScheduler"
        ),
        .target(
            name: "HealthIntegration",
            dependencies: ["FlexAIDdS"],
            path: "Sources/HealthIntegration"
        ),
        .target(
            name: "MediaIntegration",
            dependencies: ["FlexAIDdS", "HealthIntegration"],
            path: "Sources/MediaIntegration"
        ),
        .target(
            name: "Intelligence",
            dependencies: ["FlexAIDdS", "HealthIntegration"],
            path: "Sources/Intelligence"
        ),

        // Tests
        .testTarget(
            name: "FlexAIDdSTests",
            dependencies: ["FlexAIDdS", "Intelligence", "FleetScheduler"],
            path: "Tests/FlexAIDdSTests"
        ),
    ]
)
