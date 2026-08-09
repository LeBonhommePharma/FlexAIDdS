// DeviceCapability.swift — Device capability detection for fleet scheduling
//
// Detects hardware capabilities, thermal state, battery level, and compute capacity.
// Battery-aware scheduling prevents work allocation to low-power unplugged devices.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
#if canImport(IOKit)
import IOKit.ps
#endif
#if canImport(UIKit)
import UIKit
#endif

/// Hardware capabilities of a device in the Bonhomme fleet.
public struct DeviceCapability: Sendable, Codable, Hashable {
    /// Unique device identifier (hashed for privacy)
    public let deviceID: String

    /// Device model identifier (e.g., "MacBookPro18,3")
    public let model: String

    /// Estimated TFLOPS for the device GPU
    public let estimatedTFLOPS: Double

    /// Whether the device has a discrete/integrated GPU
    public let hasGPU: Bool

    /// Available memory in GB
    public let availableMemoryGB: Double

    /// Current thermal state
    public let thermalState: ThermalState

    /// Battery level (0.0 - 1.0), nil for desktop Macs
    public let batteryLevel: Double?

    /// Whether the device is connected to power
    public let isCharging: Bool

    /// Compute share weight (0.0-1.0) based on capability, thermal, and battery
    public let computeWeight: Double

    /// Timestamp of capability snapshot
    public let snapshotAt: Date

    /// Thermal state mirroring ProcessInfo.ThermalState for Codable support.
    public enum ThermalState: String, Sendable, Codable, Hashable {
        case nominal
        case fair
        case serious
        case critical
    }

    public init(
        deviceID: String,
        model: String,
        estimatedTFLOPS: Double,
        hasGPU: Bool,
        availableMemoryGB: Double,
        thermalState: ThermalState,
        batteryLevel: Double? = nil,
        isCharging: Bool = true,
        computeWeight: Double,
        snapshotAt: Date = Date()
    ) {
        self.deviceID = deviceID
        self.model = model
        self.estimatedTFLOPS = estimatedTFLOPS
        self.hasGPU = hasGPU
        self.availableMemoryGB = availableMemoryGB
        self.thermalState = thermalState
        self.batteryLevel = batteryLevel
        self.isCharging = isCharging
        self.computeWeight = computeWeight
        self.snapshotAt = snapshotAt
    }

    /// Detect capabilities of the current device.
    public static func current() -> DeviceCapability {
        let processInfo = ProcessInfo.processInfo

        let thermalState: ThermalState = {
            switch processInfo.thermalState {
            case .nominal: return .nominal
            case .fair: return .fair
            case .serious: return .serious
            case .critical: return .critical
            @unknown default: return .nominal
            }
        }()

        let memoryGB = Double(processInfo.physicalMemory) / (1024 * 1024 * 1024)

        // Estimate TFLOPS based on available processors and Apple Silicon generation
        let cores = processInfo.processorCount
        let modelID = Self.modelIdentifier()
        let estimatedTFLOPS = Self.estimateTFLOPS(cores: cores, model: modelID)

        // Battery detection
        let (battery, charging) = Self.detectBattery()

        // Compute weight incorporating thermal + battery state
        let thermalMultiplier: Double = {
            switch thermalState {
            case .nominal: return 1.0
            case .fair: return 0.75
            case .serious: return 0.4
            case .critical: return 0.0
            }
        }()

        let batteryMultiplier: Double = {
            guard let level = battery else { return 1.0 }  // Desktop Mac
            if charging { return 1.0 }  // Plugged in: full weight
            if level < 0.10 { return 0.0 }  // Critical battery: exclude
            if level < 0.20 { return 0.25 }  // Low battery: minimal work
            if level < 0.50 { return 0.6 }   // Medium battery: reduced work
            return 1.0
        }()

        let weight = min(1.0, estimatedTFLOPS / 10.0) * thermalMultiplier * batteryMultiplier

        return DeviceCapability(
            deviceID: Self.hashedDeviceID(),
            model: modelID,
            estimatedTFLOPS: estimatedTFLOPS,
            hasGPU: true,  // All Apple Silicon has GPU
            availableMemoryGB: memoryGB,
            thermalState: thermalState,
            batteryLevel: battery,
            isCharging: charging,
            computeWeight: weight,
            snapshotAt: Date()
        )
    }

    /// Whether this device is safe to schedule work on.
    public var isAvailable: Bool {
        thermalState != .critical && computeWeight > 0
    }

    /// Human-readable status summary.
    public var statusSummary: String {
        var parts: [String] = ["\(model)"]
        parts.append("\(String(format: "%.1f", estimatedTFLOPS)) TFLOPS")
        parts.append("thermal: \(thermalState.rawValue)")
        if let level = batteryLevel {
            parts.append("battery: \(Int(level * 100))%\(isCharging ? " (charging)" : "")")
        }
        parts.append("weight: \(String(format: "%.0f", computeWeight * 100))%")
        return parts.joined(separator: ", ")
    }

    // MARK: - Private Helpers

    /// Estimate TFLOPS with Apple Silicon generation awareness.
    ///
    /// `sysctlbyname("hw.model")` returns identifiers like "Mac14,7" (M2 MacBook Pro),
    /// "Mac15,3" (M3 MacBook Pro), "MacBookAir10,1" (M1 MacBook Air), etc.
    /// We map known model number ranges to Apple Silicon generations.
    /// For iOS, `hw.machine` returns e.g. "iPhone15,4", "iPad14,1".
    private static func estimateTFLOPS(cores: Int, model: String) -> Double {
        // Extract major number from identifiers like "Mac15,3", "MacBookPro18,3", "iPhone15,4"
        let pattern = #"(\d+),\d+"#
        guard let regex = try? NSRegularExpression(pattern: pattern),
              let match = regex.firstMatch(in: model, range: NSRange(model.startIndex..., in: model)),
              let range = Range(match.range(at: 1), in: model),
              let majorNum = Int(model[range]) else {
            return Double(cores) * 0.1
        }

        let lower = model.lowercased()

        // Mac-prefixed identifiers (Mac14 = M2, Mac15 = M3, Mac16 = M4)
        if lower.hasPrefix("mac") && !lower.hasPrefix("macbook") {
            switch majorNum {
            case 16...: return Double(cores) * 0.25  // M4 series
            case 15:    return Double(cores) * 0.22  // M3 series
            case 14:    return Double(cores) * 0.19  // M2 series
            case 13:    return Double(cores) * 0.16  // M1 series
            default:    return Double(cores) * 0.1
            }
        }

        // MacBookPro / MacBookAir identifiers (e.g., MacBookPro18,3 = M1 Pro)
        if lower.hasPrefix("macbook") {
            switch majorNum {
            case 16...: return Double(cores) * 0.25  // M4 series
            case 15:    return Double(cores) * 0.22  // M3 series
            case 13...14: return Double(cores) * 0.19  // M2 series
            case 10...12: return Double(cores) * 0.16  // M1 series
            default:      return Double(cores) * 0.1
            }
        }

        // iPhone identifiers (iPhone15 = A16, iPhone16 = A17 Pro, iPhone17 = A18)
        if lower.hasPrefix("iphone") {
            switch majorNum {
            case 17...: return Double(cores) * 0.14  // A18+
            case 16:    return Double(cores) * 0.13  // A17 Pro
            case 15:    return Double(cores) * 0.12  // A16
            default:    return Double(cores) * 0.10
            }
        }

        // iPad identifiers
        if lower.hasPrefix("ipad") {
            switch majorNum {
            case 16...: return Double(cores) * 0.22  // M3+ iPads
            case 14...15: return Double(cores) * 0.19  // M2 iPads
            case 13:    return Double(cores) * 0.16  // M1 iPads
            default:    return Double(cores) * 0.12
            }
        }

        // Conservative fallback for unknown models
        return Double(cores) * 0.1
    }

    /// Detect battery level and charging status across macOS and iOS.
    private static func detectBattery() -> (level: Double?, isCharging: Bool) {
        #if os(iOS) || os(watchOS)
        UIDevice.current.isBatteryMonitoringEnabled = true
        let level = Double(UIDevice.current.batteryLevel)
        let charging = UIDevice.current.batteryState == .charging || UIDevice.current.batteryState == .full
        return (level >= 0 ? level : nil, charging)
        #elseif os(macOS)
        return detectMacOSBattery()
        #else
        return (nil, false)
        #endif
    }

    #if os(macOS)
    private static func detectMacOSBattery() -> (level: Double?, isCharging: Bool) {
        guard let snapshot = IOPSCopyPowerSourcesInfo()?.takeRetainedValue(),
              let sources = IOPSCopyPowerSourcesList(snapshot)?.takeRetainedValue() as? [Any],
              let firstSource = sources.first,
              let desc = IOPSGetPowerSourceDescription(snapshot, firstSource as CFTypeRef)?.takeUnretainedValue() as? [String: Any] else {
            return (nil, false)  // Desktop Mac without battery
        }

        let currentCap = desc[kIOPSCurrentCapacityKey] as? Int ?? 0
        let maxCap = desc[kIOPSMaxCapacityKey] as? Int ?? 100
        let isCharging = (desc[kIOPSIsChargingKey] as? Bool) ?? false
        let level = maxCap > 0 ? Double(currentCap) / Double(maxCap) : nil

        return (level, isCharging)
    }
    #endif

    private static func hashedDeviceID() -> String {
        let hostName = ProcessInfo.processInfo.hostName
        return String(hostName.hashValue, radix: 16)
    }

    private static func modelIdentifier() -> String {
        #if os(macOS)
        var size: size_t = 0
        sysctlbyname("hw.model", nil, &size, nil, 0)
        var model = [CChar](repeating: 0, count: size)
        sysctlbyname("hw.model", &model, &size, nil, 0)
        return String(cString: model)
        #else
        return ProcessInfo.processInfo.hostName
        #endif
    }
}
