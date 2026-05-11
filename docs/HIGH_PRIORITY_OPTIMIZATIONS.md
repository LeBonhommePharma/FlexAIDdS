# High-Priority Optimizations & Enhancements for FlexAID∆S

**Branch:** `feature/high-priority-optimizations`
**Date:** May 2026
**Author:** Grok (on behalf of LeBonhommePharma research workflow)

This PR introduces concrete code diffs, refactored snippets, and implementation guidance for the three **High Priority** items from the deep code audit.

## Summary of Changes

1. **Numerical Stability (Log-Sum-Exp)**: Critical correctness fix for partition function and Shannon entropy calculations.
2. **Memory Layout (Expand SoA)**: Performance improvement in hot paths (entropy reduction, GA fitness).
3. **Python API Enhancement**: High-level `dock_and_analyze()` for better usability and structured entropy output.

All changes are designed to be low-risk, high-impact, and align with modern C++26 and scientific computing best practices.

---

## 1. Numerical Stability: Robust Log-Sum-Exp for Partition Functions & Entropy

**Problem:** Direct `exp(energy)` in ensemble calculations can under/overflow when energies span >20 kcal/mol (common in docking).

**Solution:** Add a reusable `log_sum_exp` utility and update entropy/partition function code to use it. Also provide `softmax` helper for probabilities.

### Refactored Snippet (New File)

```cpp
// LIB/utils/LogSumExp.h (new file - proposed)
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace flexaid::utils {

/**
 * @brief Numerically stable log-sum-exp.
 * log(sum(exp(x_i))) = max(x) + log(sum(exp(x_i - max(x))))
 */
inline double log_sum_exp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    double max_val = *std::max_element(values.begin(), values.end());
    double sum = 0.0;
    for (double v : values) {
        sum += std::exp(v - max_val);
    }
    return max_val + std::log(sum);
}

/**
 * @brief Stable softmax / normalized probabilities from energies.
 * p_i = exp(-beta * E_i) / Z
 */
inline std::vector<double> stable_softmax(const std::vector<double>& energies, double beta = 1.0) {
    if (energies.empty()) return {};
    std::vector<double> shifted(energies.size());
    double max_e = *std::max_element(energies.begin(), energies.end());
    for (size_t i = 0; i < energies.size(); ++i) {
        shifted[i] = -beta * (energies[i] - max_e);
    }
    double logZ = log_sum_exp(shifted);
    std::vector<double> probs(energies.size());
    for (size_t i = 0; i < energies.size(); ++i) {
        probs[i] = std::exp(shifted[i] - logZ);
    }
    return probs;
}

} // namespace flexaid::utils
```

### Usage Example in ShannonThermoStack or GrandPartitionFunction

```cpp
// Before (risky):
// double Z = 0; for(auto e : energies) Z += exp(-beta * e);

// After (stable):
#include "utils/LogSumExp.h"
auto probs = flexaid::utils::stable_softmax(energies, beta);
double entropy = 0.0;
for (double p : probs) if (p > 0) entropy -= p * std::log(p);
```

**Expected Impact:** Eliminates NaN/Inf in entropy calculations for challenging ligands. Critical for reproducibility of thermodynamic results.

**Files to modify:** `LIB/ShannonThermoStack/ShannonThermoStack.cpp`, `LIB/GrandPartitionFunction.cpp`, add new header.

---

## 2. Memory Layout: Expand Structure-of-Arrays (SoA) Patterns

**Problem:** AoS (Array of Structs) for Pose data causes poor cache locality in hot loops (entropy reduction, GA fitness evaluation, Voronoi scoring).

**Solution:** Expand existing SoA infrastructure (`AtomSoA.h`, `VoronoiCFBatch_SoA.h`) to Pose/ensemble data and key reduction kernels.

### Refactored Snippet (Proposed Extension)

```cpp
// Example extension in LIB/ (or new PoseSoA.h)
struct PoseSoA {
    std::vector<double> x, y, z;           // coordinates
    std::vector<double> energies;
    std::vector<double> boltzmann_weights;
    // Add more as needed (e.g., partial charges, radii)

    size_t size() const { return x.size(); }

    void reserve(size_t n) {
        x.reserve(n); y.reserve(n); z.reserve(n);
        energies.reserve(n); boltzmann_weights.reserve(n);
    }

    void push_back(const Pose& pose) { /* ... */ }
};

// In entropy reduction kernel:
double compute_shannon_entropy_soa(const PoseSoA& ensemble, double T) {
    auto probs = utils::stable_softmax(ensemble.energies, 1.0 / (kB * T));
    double S = 0.0;
    for (double p : probs) if (p > 1e-300) S -= p * std::log(p);
    return S;
}
```

**Recommended Hot Paths to Convert First:**
- Entropy calculation in `ShannonThermoStack`
- GA population fitness evaluation in `GAContext`
- Batch Voronoi / energy grid evaluation

**Expected Impact:** 1.5–3× speedup on CPU for large ensembles due to better cache utilization and SIMD friendliness.

**Existing files to enhance:** `LIB/AtomSoA.h`, `LIB/VoronoiCFBatch_SoA.h`, `LIB/GAContext.cpp`, `LIB/ShannonThermoStack/ShannonThermoStack.cpp`

---

## 3. Python API Enhancement: High-Level `dock_and_analyze()`

**Problem:** Current API requires multiple calls to get full thermodynamic breakdown. Users want one-call structured output (especially for Jupyter/PyMOL workflows).

**Solution:** Add a convenience function in the Python bindings that returns a rich dataclass or dict with ΔG, ΔH, -TΔS breakdown, top poses, and entropy components.

### Refactored / New Python Snippet (in python/flexaidds/)

```python
# python/flexaidds/api.py (proposed addition)
from dataclasses import dataclass
from typing import List, Optional, Dict

import flexaidds._core as core  # the pybind11 module

@dataclass
class DockingResult:
    free_energy: float
    enthalpy: float
    entropy: float          # -TΔS
    configurational_entropy: float
    vibrational_entropy: float
    hydration_entropy: float
    top_poses: List[Dict]   # list of pose dicts with coordinates, energy, etc.
    metadata: Dict

def dock_and_analyze(
    receptor: str,
    ligand: str,
    compute_entropy: bool = True,
    n_poses: int = 100,
    temperature: float = 298.15,
    **kwargs
) -> DockingResult:
    """
    High-level convenience function.
    Returns structured thermodynamic results + top poses.
    """
    raw = core.dock(receptor=receptor, ligand=ligand, n_poses=n_poses, **kwargs)
    
    if compute_entropy:
        engine = core.StatMechEngine(temperature=temperature)
        thermo = engine.compute(raw.poses)  # assumes engine exposes breakdown
        
        return DockingResult(
            free_energy=thermo.free_energy,
            enthalpy=thermo.enthalpy,
            entropy=thermo.entropy,
            configurational_entropy=thermo.configurational,
            vibrational_entropy=thermo.vibrational,
            hydration_entropy=thermo.hydration,
            top_poses=raw.top_poses,
            metadata={"n_poses": n_poses, "temperature": temperature}
        )
    else:
        # fallback without full thermo
        ...
```

**Usage (Jupyter/PyMOL friendly):**

```python
import flexaidds as fd
result = fd.dock_and_analyze("receptor.pdb", "ligand.sdf")
print(result.free_energy)
print(result.configurational_entropy)
# Visualize in PyMOL
fd.visualize_top_poses(result)
```

**Expected Impact:** Dramatically improves usability and adoption, especially for experimental collaborators and rapid iteration in psychopharmacology workflows.

**Files to modify:** `python/flexaidds/__init__.py` or new `api.py`, update pybind11 bindings if needed for richer return types.

---

## Implementation Roadmap & Testing

- [ ] Add `LIB/utils/LogSumExp.h` and integrate into entropy kernels.
- [ ] Expand SoA in 2–3 hot paths + benchmark (use existing BenchmarkRunner).
- [ ] Implement `dock_and_analyze` in Python layer + add to docs/examples.
- [ ] Add unit tests for log_sum_exp (edge cases: empty, large range, all equal).
- [ ] Update README and `docs/` with usage examples.
- [ ] Run full benchmark suite on Astex + internal psychopharm set.

**Backward Compatibility:** All changes are additive or internal optimizations. Existing API remains unchanged.

**Performance Validation:** After implementation, re-run GA + entropy benchmarks and update numbers in manuscripts.

---

This PR serves as the living document and starting point for the optimizations. Code diffs above can be applied incrementally. Full implementation PRs can be split if preferred.

**Next Steps after Merge:**
- Iterate on the SoA conversion with profiling data.
- Expose entropy component breakdown in the main C++ API as well.

Thank you for the opportunity to contribute to FlexAID∆S. Entropy matters — let's make it fast and robust.
