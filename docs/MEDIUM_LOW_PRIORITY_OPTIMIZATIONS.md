# Medium & Low Priority Optimizations for FlexAID∆S

**Branch:** `feature/medium-low-priorities`
**Date:** May 2026
**Follow-up to:** PR #193 (High-Priority Optimizations)

This PR covers the remaining **Medium** and **Low** priority items from the deep code audit of FlexAID∆S.

## Summary

**Medium Priority (Performance + Accuracy)**
1. Genetic Algorithm improvements (diversity-preserving replacement + entropy-collapse early stopping)
2. Entropy estimation enhancements (optional robust KDE mode with adaptive bandwidth)
3. Hardware acceleration usability (expose `DeviceContext` + decision heuristics)

**Low Priority (Maintainability + Robustness)**
4. Modern C++ idioms (`std::expected`, `std::mdspan`, concepts)
5. Enhanced testing (property-based + golden regression tests for entropy)
6. C++ API exposure of full entropy component breakdown

All items are designed to be incremental, low-risk, and build directly on the high-priority foundation from PR #193.

---

## Medium Priority Items

### 1. Genetic Algorithm: Diversity-Preserving Replacement + Entropy-Collapse Early Stopping

**Current State:** The GA already has good entropy-collapse monitoring for diversity. We can make it more explicit and robust.

**Proposed Changes:**

```cpp
// In GAContext.cpp / GAContext.h

enum class ReplacementStrategy {
    Elitist,
    Crowding,           // diversity-preserving
    EntropyBased        // new: use entropy collapse signal
};

void replace_population(
    std::vector<Pose>& population,
    const std::vector<Pose>& offspring,
    ReplacementStrategy strategy = ReplacementStrategy::EntropyBased
) {
    if (strategy == ReplacementStrategy::EntropyBased) {
        // Calculate current entropy of population
        double current_entropy = compute_population_entropy(population);
        if (current_entropy < diversity_threshold) {
            // Apply crowding or restricted tournament selection
            apply_crowding_replacement(population, offspring);
            return;
        }
    }
    // fallback to elitist or standard replacement
    apply_elitist_replacement(population, offspring);
}

// Early stopping based on entropy stability
bool should_stop_early(const std::vector<double>& entropy_history) {
    if (entropy_history.size() < 5) return false;
    double recent_variance = calculate_variance(
        entropy_history.end() - 5, entropy_history.end());
    return recent_variance < entropy_stability_threshold;
}
```

**Expected Impact:** Better exploration/exploitation balance, fewer wasted generations on converged (low-diversity) populations.

**Files:** `LIB/GAContext.cpp`, `LIB/GAContext.h`

---

### 2. Entropy Estimation: Optional Robust KDE Mode

**Current:** Histogram-based or simple KDE for pose probability estimation.

**Enhancement:**

```cpp
// In ShannonThermoStack or ReferenceEntropy

enum class EntropyEstimationMode {
    Histogram,      // fast, default
    RobustKDE,      // more accurate for small ensembles
    AdaptiveKDE     // automatic bandwidth selection
};

double compute_configurational_entropy(
    const std::vector<Pose>& poses,
    EntropyEstimationMode mode = EntropyEstimationMode::Histogram
) {
    if (mode == EntropyEstimationMode::RobustKDE) {
        // Use Gaussian KDE with Silverman's rule or cross-validation bandwidth
        return kde_shannon_entropy(poses);
    }
    // histogram path...
}
```

**Expected Impact:** Higher accuracy on smaller pose ensembles or when pose distribution is multimodal.

**Files:** `LIB/ShannonThermoStack/`, `LIB/ReferenceEntropy.cpp`

---

### 3. Hardware: Expose DeviceContext + Decision Heuristics

**Goal:** Make it easy for users to control and query hardware backend.

```cpp
// Proposed new API in HardwareDispatch.h

struct DeviceContext {
    enum class Type { CPU, CUDA, Metal, Auto };
    Type type;
    int device_id = 0;           // for multi-GPU
    bool enable_profiling = false;
};

class UnifiedHardwareDispatch {
public:
    static DeviceContext select_best_device(const std::string& workload_type);
    static void set_global_context(const DeviceContext& ctx);
    static DeviceContext get_current_context();
};

// Usage
DeviceContext ctx = UnifiedHardwareDispatch::select_best_device("entropy_reduction");
ctx.enable_profiling = true;
UnifiedHardwareDispatch::set_global_context(ctx);
```

**Expected Impact:** Better control on heterogeneous systems and easier performance tuning.

**Files:** `LIB/HardwareDispatch.h`, `LIB/UnifiedHardwareDispatch.cpp`

---

## Low Priority Items

### 4. Modern C++ Idioms

- Replace raw error handling with `std::expected<double, ErrorCode>` in energy/entropy functions.
- Use `std::mdspan` for multi-dimensional grids (GIST, energy maps).
- Add C++20 concepts for scoring function interfaces:

```cpp
 template<typename T>
concept ScoringFunction = requires(T f, const Pose& p) {
    { f.score(p) } -> std::convertible_to<double>;
};
```

**Impact:** Cleaner, safer, more maintainable code.

### 5. Enhanced Testing

- Add property-based tests for `log_sum_exp` and entropy functions (using rapidcheck or similar).
- Create "golden" regression tests with small systems that have known ITC-comparable entropy values.
- Fuzz GA with random seeds + edge-case ligands.

### 6. Full Entropy Component Breakdown in C++ API

Expose in `StatMechEngine` / `BindingMode`:

```cpp
struct EntropyBreakdown {
    double configurational;
    double vibrational;
    double hydration;
    double total;
};

EntropyBreakdown get_entropy_breakdown() const;
```

Already partially available in Python via PR #193 — now bring parity to C++.

---

## Implementation Roadmap

- [ ] Implement LogSumExp-style utilities first (already proposed in PR #193)
- [ ] Add `DeviceContext` and selection heuristics
- [ ] Introduce `ReplacementStrategy::EntropyBased` in GA
- [ ] Add optional KDE mode behind a flag
- [ ] Modernize selected hot functions with `std::expected` and `std::mdspan`
- [ ] Expand test suite with property-based + golden tests
- [ ] Expose `EntropyBreakdown` struct in main C++ API

**Priority Order Recommendation:**
Medium items first (biggest user-visible gains), then Low items for long-term code health.

---

This PR serves as the planning and starting-point document for the remaining optimizations. All snippets are ready for incremental application.

**Related:** See PR #193 for the High-Priority foundation (numerical stability, SoA, Python API).

Thank you — continuing to make FlexAID∆S the most thermodynamically rigorous and usable open-source docking engine available.
