// FleetExplainer.swift — Natural-language fleet scheduling explanations
//
// Explains fleet scheduling decisions, device allocation, bottlenecks,
// and actionable fixes using Apple FoundationModels.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Fleet Context

/// Device summary for LLM interpretation.
public struct DeviceSummary: Sendable, Codable {
    public let model: String
    public let tflops: Double
    public let thermalState: String   // "nominal", "fair", "serious", "critical"
    public let batteryPercent: Int?    // nil for desktop Macs
    public let isCharging: Bool
    public let computeWeight: Double   // 0-1 share of total compute
    public let chunksAssigned: Int
    public let chunksCompleted: Int
    public let chunksFailed: Int

    public init(model: String, tflops: Double, thermalState: String,
                batteryPercent: Int?, isCharging: Bool, computeWeight: Double,
                chunksAssigned: Int, chunksCompleted: Int, chunksFailed: Int) {
        self.model = model
        self.tflops = tflops
        self.thermalState = thermalState
        self.batteryPercent = batteryPercent
        self.isCharging = isCharging
        self.computeWeight = computeWeight
        self.chunksAssigned = chunksAssigned
        self.chunksCompleted = chunksCompleted
        self.chunksFailed = chunksFailed
    }
}

/// Full fleet status for LLM explanation.
public struct FleetStatusContext: Sendable, Codable {
    public let devices: [DeviceSummary]
    public let totalChunks: Int
    public let completedChunks: Int
    public let failedChunks: Int
    public let orphanedChunks: Int
    public let etaSeconds: Double?
    public let totalTFLOPS: Double
    public let refereeRecommendation: String?  // from RefereeRecommendation enum

    public init(devices: [DeviceSummary], totalChunks: Int, completedChunks: Int,
                failedChunks: Int, orphanedChunks: Int, etaSeconds: Double?,
                totalTFLOPS: Double, refereeRecommendation: String?) {
        self.devices = devices
        self.totalChunks = totalChunks
        self.completedChunks = completedChunks
        self.failedChunks = failedChunks
        self.orphanedChunks = orphanedChunks
        self.etaSeconds = etaSeconds
        self.totalTFLOPS = totalTFLOPS
        self.refereeRecommendation = refereeRecommendation
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Fleet scheduling explanation.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct FleetExplanation: Sendable, Codable {
    /// Why devices were allocated the way they are
    public var allocationRationale: String
    /// What's slowing things down
    public var bottleneckAnalysis: String
    /// Concrete actions the user can take
    public var actionItems: [String]
    /// Human-friendly ETA estimate
    public var estimatedCompletion: String
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor FleetExplainerActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a fleet scheduling advisor for distributed molecular docking. \
        All device metrics and job progress are pre-computed — DO NOT calculate. \
        Explain scheduling decisions in plain English, identify bottlenecks, \
        and suggest concrete actions (e.g., "plug in iPad", "close background apps"). \
        Be practical and user-friendly.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Explain the current fleet scheduling state.
    public func explain(context: FleetStatusContext) async throws -> FleetExplanation {
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            // Truncate to top 3 devices by compute weight
            let topDevices = Array(context.devices.sorted { $0.computeWeight > $1.computeWeight }.prefix(3))
            let truncated = FleetStatusContext(
                devices: topDevices, totalChunks: context.totalChunks,
                completedChunks: context.completedChunks, failedChunks: context.failedChunks,
                orphanedChunks: context.orphanedChunks, etaSeconds: context.etaSeconds,
                totalTFLOPS: context.totalTFLOPS, refereeRecommendation: context.refereeRecommendation
            )
            prompt = buildPrompt(context: truncated)
        }
        return try await session.respond(to: prompt, generating: FleetExplanation.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: FleetStatusContext) -> String {
        var p = "Explain this fleet status. Produce a FleetExplanation.\n"
        p += "Progress: \(context.completedChunks)/\(context.totalChunks) chunks"
        if context.failedChunks > 0 { p += " (\(context.failedChunks) failed)" }
        if context.orphanedChunks > 0 { p += " (\(context.orphanedChunks) orphaned)" }
        p += "\nFleet: \(context.devices.count) devices, \(String(format: "%.1f", context.totalTFLOPS)) TFLOPS total"
        if let eta = context.etaSeconds {
            let mins = Int(eta / 60)
            p += "\nETA: ~\(mins) minutes"
        }

        p += "\n\nDevices:"
        for d in context.devices.prefix(5) {
            p += "\n  \(d.model): \(String(format: "%.1f", d.tflops))TF, thermal=\(d.thermalState)"
            if let batt = d.batteryPercent { p += ", battery=\(batt)%\(d.isCharging ? "+" : "")" }
            p += ", weight=\(String(format: "%.0f", d.computeWeight * 100))%"
            p += ", chunks=\(d.chunksCompleted)/\(d.chunksAssigned)"
            if d.chunksFailed > 0 { p += " (\(d.chunksFailed) failed)" }
        }

        if let rec = context.refereeRecommendation {
            p += "\nReferee recommendation: \(rec)"
        }

        // Flags
        let throttled = context.devices.filter { $0.thermalState == "serious" || $0.thermalState == "critical" }
        if !throttled.isEmpty {
            p += "\nFLAG: \(throttled.count) device(s) thermally throttled"
        }
        let lowBattery = context.devices.filter { ($0.batteryPercent ?? 100) < 20 && !$0.isCharging }
        if !lowBattery.isEmpty {
            p += "\nFLAG: \(lowBattery.count) device(s) low battery and unplugged"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent fleet explanation.
public struct CrossPlatformFleetExplanation: Sendable, Codable {
    public let allocationRationale: String
    public let bottleneckAnalysis: String
    public let actionItems: [String]
    public let estimatedCompletion: String

    public init(allocationRationale: String, bottleneckAnalysis: String,
                actionItems: [String], estimatedCompletion: String) {
        self.allocationRationale = allocationRationale
        self.bottleneckAnalysis = bottleneckAnalysis
        self.actionItems = actionItems
        self.estimatedCompletion = estimatedCompletion
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic fleet explainer for non-Apple platforms.
public struct RuleBasedFleetExplainer: Sendable {

    public init() {}

    /// Explain fleet scheduling using template logic.
    public func explain(context: FleetStatusContext) -> CrossPlatformFleetExplanation {
        guard !context.devices.isEmpty else {
            return CrossPlatformFleetExplanation(
                allocationRationale: "No devices connected to the fleet.",
                bottleneckAnalysis: "No devices available for compute.",
                actionItems: ["Connect at least one device to begin distributed docking."],
                estimatedCompletion: "Cannot estimate — no devices available."
            )
        }
        // Allocation rationale
        let sorted = context.devices.sorted { $0.computeWeight > $1.computeWeight }
        var rationale = "Work distributed by compute weight: "
        rationale += sorted.prefix(3).map { "\($0.model) (\(String(format: "%.0f", $0.computeWeight * 100))%)" }.joined(separator: ", ")
        if sorted.count > 3 { rationale += ", and \(sorted.count - 3) more" }
        rationale += "."

        // Bottleneck analysis
        var bottleneck = ""
        let throttled = context.devices.filter { $0.thermalState == "serious" || $0.thermalState == "critical" }
        let lowBattery = context.devices.filter { ($0.batteryPercent ?? 100) < 20 && !$0.isCharging }
        let failingDevices = context.devices.filter { $0.chunksFailed > 0 }

        if !throttled.isEmpty {
            bottleneck += "\(throttled.map(\.model).joined(separator: ", ")) thermally throttled — reduced compute capacity. "
        }
        if !lowBattery.isEmpty {
            bottleneck += "\(lowBattery.map(\.model).joined(separator: ", ")) at low battery without power. "
        }
        if !failingDevices.isEmpty {
            bottleneck += "\(failingDevices.map(\.model).joined(separator: ", ")) experiencing chunk failures. "
        }
        if context.orphanedChunks > 0 {
            bottleneck += "\(context.orphanedChunks) orphaned chunks being reassigned. "
        }
        if bottleneck.isEmpty {
            bottleneck = "No bottlenecks detected. Fleet operating normally."
        }

        // Action items
        var actions: [String] = []
        for d in lowBattery { actions.append("Plug in \(d.model) to maintain compute capacity.") }
        for d in throttled { actions.append("Allow \(d.model) to cool down or reduce workload.") }
        if context.failedChunks > 2 { actions.append("Investigate chunk failures — may indicate memory pressure.") }
        if actions.isEmpty { actions.append("No action needed — fleet is healthy.") }

        // ETA
        let eta: String
        if let seconds = context.etaSeconds {
            if seconds < 60 { eta = "Less than 1 minute remaining." }
            else if seconds < 3600 { eta = "~\(Int(seconds / 60)) minutes remaining." }
            else { eta = "~\(String(format: "%.1f", seconds / 3600)) hours remaining." }
        } else {
            let remaining = context.totalChunks - context.completedChunks
            eta = "\(remaining) chunks remaining."
        }

        return CrossPlatformFleetExplanation(
            allocationRationale: rationale,
            bottleneckAnalysis: bottleneck,
            actionItems: actions,
            estimatedCompletion: eta
        )
    }
}
