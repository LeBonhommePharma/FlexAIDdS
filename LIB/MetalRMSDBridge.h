// MetalRMSDBridge.h — C++ interface to Metal GPU pairwise RMSD kernel
//
// Provides batch pairwise RMSD computation on Apple Silicon GPU.
// Falls back to CPU when Metal is unavailable.
// The implementation is in MetalRMSDBridge.mm (Objective-C++) and
// compiled only on APPLE targets with FLEXAIDS_USE_METAL defined.
//
// Apache-2.0 (C) 2026 Le Bonhomme Pharma
#pragma once

#include <vector>
#include <cstddef>

namespace metal_rmsd {

/// Compute pairwise RMSD matrix on Metal GPU.
///
/// @param coords  Flat coordinate buffer [n_conf * 3 * n_atoms].
///                Layout: conf[i] starts at coords[i * 3 * n_atoms],
///                with atoms laid out as x0,y0,z0, x1,y1,z1, ...
/// @param n_conf  Number of conformations (poses).
/// @param n_atoms Number of atoms per conformation.
/// @param dist_matrix  Output vector of size n_conf * n_conf, row-major.
///                dist_matrix[i * n_conf + j] = RMSD(i, j).
///                Diagonal is 0.0, matrix is symmetric.
///
/// @return true if Metal GPU was used, false if CPU fallback was used.
bool compute_pairwise_rmsd_metal(const float* coords,
                                  int n_conf,
                                  int n_atoms,
                                  std::vector<float>& dist_matrix);

/// Check if Metal GPU RMSD is available and the pipeline is ready.
bool is_metal_rmsd_available();

/// Get diagnostic string describing the Metal device used for RMSD.
const char* metal_rmsd_device_info();

} // namespace metal_rmsd
