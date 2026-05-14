// ShannonMetalBridge.h — C++ interface to Apple Metal Shannon GPU kernels
//
// Provides a C++ header-only interface. The implementation is in
// ShannonMetalBridge.mm (Objective-C++) and compiled only on APPLE targets.
// On non-Apple platforms this header is included but the implementation is
// never compiled; callers guard behind #ifdef FLEXAIDS_HAS_METAL_SHANNON.
//
// Features:
//   – Persistent device/pipeline/queue caching (no per-call init overhead)
//   – Shannon histogram kernel
//   – Boltzmann weight batch kernel
//   – Parallel sum reduction for partition function
//   – Log-sum-exp kernel for numerical stability
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#pragma once

#include <vector>
#include <string>

namespace ShannonMetalBridge {

/// Compute Shannon entropy of `energies` using a Metal GPU histogram kernel.
/// Falls back to CPU if no Metal device is available.
/// Returns entropy in nats (natural logarithm, base-e).
/// Use H_bits = H_nats / ln(2) only at reporting/convergence boundaries.
double compute_shannon_entropy_metal(const std::vector<double>& energies,
                                     int num_bins = 20);

/// Compute Boltzmann weights on Metal GPU.
/// Returns weights[i] = exp(-beta * (E[i] - E_min)).
/// E_min is computed internally.
/// sum_w is set to the sum of all weights (for partition function).
std::vector<double> compute_boltzmann_weights_metal(
    const std::vector<double>& energies,
    double beta,
    double& sum_w,
    double& E_min);

/// Compute log-sum-exp on Metal GPU.
/// Returns log(sum(exp(values))) with numerical stability.
double log_sum_exp_metal(const std::vector<double>& values);

/// Check if Metal GPU is available and cached pipeline is ready.
bool is_metal_available();

/// Get a diagnostic string describing the Metal device.
std::string metal_device_info();

} // namespace ShannonMetalBridge
