// gpu_fast_optics_metal.h — Metal GPU kNN for FastOPTICS
//
// When FLEXAIDS_USE_METAL is defined, provides a GPU-accelerated k-nearest-
// neighbour search using Apple's Metal compute framework.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab
#pragma once

#ifdef FLEXAIDS_USE_METAL

#include <vector>
#include <utility>

struct chromosome;

// Metal-accelerated kNN for FastOPTICS neighbor discovery.
// Uploads Cartesian point data to GPU, launches Metal compute pipeline,
// downloads per-point neighbour lists.
//
// points:     vector of (chromosome*, coordinates) pairs
// k:          number of nearest neighbours (typically minPts)
// nDim:       dimensionality of each point
// out_neighbors: [N][variable] neighbour indices (resized and filled)
// out_distances: [N][k] neighbour distances (resized and filled)
void metal_foptics_knn(
    const std::vector<std::pair<chromosome*, std::vector<float>>>& points,
    int k, int nDim,
    std::vector<std::vector<int>>& out_neighbors,
    std::vector<std::vector<float>>& out_distances);

#endif // FLEXAIDS_USE_METAL
