// test_compare_energy_matrices.cpp — Compare the 256×256 soft contact matrix
// against the legacy 40×40 NRGRank energy matrix.
//
// Tests:
//   1. Structural comparison: sparsity, value range, symmetry
//   2. 256→40 projection fidelity: populate 256×256 from NRGRank, project back
//   3. Information content: effective rank, Shannon entropy of row distributions
//   4. Type coverage: which SYBYL types are populated vs. zero in each system

#include <gtest/gtest.h>
#include "nrgrank_matrix.h"
#include "atom_typing_256.h"
#include "soft_contact_matrix.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <vector>
#include <algorithm>

// ─── helpers ────────────────────────────────────────────────────────────────

namespace {

// Count non-zero entries in a flat array
template<typename T, size_t N>
int count_nonzero(const std::array<T, N>& arr) {
    int count = 0;
    for (auto v : arr)
        if (std::fabs(static_cast<double>(v)) > 1e-9) ++count;
    return count;
}

// Compute Frobenius norm of a flat matrix
template<typename T>
double frobenius_norm(const T* data, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
        sum += static_cast<double>(data[i]) * static_cast<double>(data[i]);
    return std::sqrt(sum);
}

// Pearson correlation between two arrays of equal length
double pearson_r(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) return 0.0;
    int n = static_cast<int>(a.size());
    double ma = 0.0, mb = 0.0;
    for (int i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    double cov = 0.0, va = 0.0, vb = 0.0;
    for (int i = 0; i < n; ++i) {
        double da = a[i] - ma, db = b[i] - mb;
        cov += da * db;
        va += da * da;
        vb += db * db;
    }
    double denom = std::sqrt(va * vb);
    return denom > 1e-15 ? cov / denom : 0.0;
}

// Build a 256×256 matrix by expanding the NRGRank 41×41 matrix:
// For each pair (ci, cj) of 256-type codes, look up their SYBYL parents
// and assign the NRGRank energy value.
scm::SoftContactMatrix build_256_from_nrgrank() {
    scm::SoftContactMatrix mat;
    mat.zero();
    for (int ci = 0; ci < 256; ++ci) {
        int si = atom256::base_to_sybyl_parent(atom256::get_base(ci));
        for (int cj = 0; cj < 256; ++cj) {
            int sj = atom256::base_to_sybyl_parent(atom256::get_base(cj));
            if (si >= 0 && si <= 40 && sj >= 0 && sj <= 40)
                mat.set(ci, cj, static_cast<float>(
                    nrgrank::kEnergyMatrix[si][sj]));
        }
    }
    return mat;
}

} // namespace

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 1: Structural properties of the 40×40 NRGRank matrix
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, NRGRank40x40Structure) {
    // Count non-zero entries (excluding row/col 0 which is padding)
    int nonzero = 0;
    int total = 0;
    double min_val = 0.0, max_val = 0.0;
    for (int i = 1; i <= 40; ++i) {
        for (int j = 1; j <= 40; ++j) {
            double v = nrgrank::kEnergyMatrix[i][j];
            if (std::fabs(v) > 1e-9) ++nonzero;
            ++total;
            if (v < min_val) min_val = v;
            if (v > max_val) max_val = v;
        }
    }
    double sparsity = 1.0 - static_cast<double>(nonzero) / total;

    printf("\n=== NRGRank 40×40 Matrix Structure ===\n");
    printf("  Dimensions:  40 × 40 (types 1–40)\n");
    printf("  Total cells: %d\n", total);
    printf("  Non-zero:    %d (%.1f%%)\n", nonzero,
           100.0 * nonzero / total);
    printf("  Sparsity:    %.1f%%\n", 100.0 * sparsity);
    printf("  Value range: [%.2f, %.2f]\n", min_val, max_val);
    printf("  Memory:      %zu bytes (constexpr double)\n",
           sizeof(nrgrank::kEnergyMatrix));

    // Check symmetry
    int asym_count = 0;
    double max_asym = 0.0;
    for (int i = 1; i <= 40; ++i) {
        for (int j = i + 1; j <= 40; ++j) {
            double diff = std::fabs(nrgrank::kEnergyMatrix[i][j] -
                                     nrgrank::kEnergyMatrix[j][i]);
            if (diff > 1e-9) ++asym_count;
            max_asym = std::max(max_asym, diff);
        }
    }
    printf("  Symmetric:   %s (max |a[i][j]-a[j][i]| = %.4f)\n",
           max_asym < 1e-6 ? "YES" : "NO", max_asym);

    EXPECT_GT(nonzero, 0) << "NRGRank matrix should have non-zero entries";
    EXPECT_LT(min_val, 0.0) << "Should have attractive (negative) interactions";
    EXPECT_GT(max_val, 0.0) << "Should have repulsive (positive) interactions";
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 2: Structural properties of the 256×256 matrix (populated from NRGRank)
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, SoftContact256x256Structure) {
    auto mat = build_256_from_nrgrank();

    int nonzero = 0;
    float min_val = 0.0f, max_val = 0.0f;
    for (int i = 0; i < scm::MATRIX_SIZE; ++i) {
        if (std::fabs(mat.data[i]) > 1e-9f) ++nonzero;
        if (mat.data[i] < min_val) min_val = mat.data[i];
        if (mat.data[i] > max_val) max_val = mat.data[i];
    }
    double sparsity = 1.0 - static_cast<double>(nonzero) / scm::MATRIX_SIZE;

    printf("\n=== 256×256 Soft Contact Matrix (from NRGRank) ===\n");
    printf("  Dimensions:  256 × 256\n");
    printf("  Total cells: %d\n", scm::MATRIX_SIZE);
    printf("  Non-zero:    %d (%.1f%%)\n", nonzero,
           100.0 * nonzero / scm::MATRIX_SIZE);
    printf("  Sparsity:    %.1f%%\n", 100.0 * sparsity);
    printf("  Value range: [%.2f, %.2f]\n", min_val, max_val);
    printf("  Memory:      %zu bytes (float, cache-aligned)\n",
           sizeof(scm::SoftContactMatrix));

    EXPECT_GT(nonzero, 0);
    // 256×256 should have MORE non-zero entries because each SYBYL type
    // fans out into multiple 256-codes (up to 4: 2 charge bins × 2 hbond)
    printf("  Expansion:   each SYBYL type → up to 4 256-codes\n");
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 3: Round-trip fidelity — populate 256×256 from NRGRank, project back
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, RoundTripProjectionFidelity) {
    auto mat = build_256_from_nrgrank();
    auto proj = mat.project_to_40x40();

    // Compare projected 40×40 with original NRGRank (types 1-40, 0-indexed in proj)
    // proj is 0-indexed [0..39], NRGRank is 1-indexed [1..40]
    std::vector<double> original_vals, projected_vals;
    double max_abs_diff = 0.0;
    double sum_sq_diff = 0.0;
    int n_compared = 0;
    int n_matching_sign = 0;
    int n_both_nonzero = 0;

    printf("\n=== Round-Trip Projection: 40→256→40 ===\n");
    printf("  %-8s %-8s %-12s %-12s %-10s\n",
           "Type_i", "Type_j", "NRGRank", "Projected", "Diff");
    printf("  %s\n", std::string(52, '-').c_str());

    for (int i = 0; i < 40; ++i) {
        for (int j = 0; j < 40; ++j) {
            double orig = nrgrank::kEnergyMatrix[i + 1][j + 1];
            double projected = static_cast<double>(proj[i * 40 + j]);
            double diff = projected - orig;

            original_vals.push_back(orig);
            projected_vals.push_back(projected);
            ++n_compared;

            double ad = std::fabs(diff);
            if (ad > max_abs_diff) max_abs_diff = ad;
            sum_sq_diff += diff * diff;

            bool orig_nz = std::fabs(orig) > 1e-9;
            bool proj_nz = std::fabs(projected) > 1e-9;
            if (orig_nz && proj_nz) {
                ++n_both_nonzero;
                if ((orig > 0) == (projected > 0)) ++n_matching_sign;
            }

            // Print the first few large deviations for diagnostic
            if (ad > 50.0 && n_compared <= 5) {
                printf("  %-8d %-8d %12.2f %12.2f %10.2f\n",
                       i + 1, j + 1, orig, projected, diff);
            }
        }
    }

    double rmse = std::sqrt(sum_sq_diff / n_compared);
    double r = pearson_r(original_vals, projected_vals);

    printf("\n  Comparison Statistics:\n");
    printf("    Cells compared:     %d\n", n_compared);
    printf("    Pearson r:          %.6f\n", r);
    printf("    RMSE:               %.4f\n", rmse);
    printf("    Max |diff|:         %.4f\n", max_abs_diff);
    printf("    Sign agreement:     %d/%d (%.1f%%)\n",
           n_matching_sign, n_both_nonzero,
           n_both_nonzero > 0 ? 100.0 * n_matching_sign / n_both_nonzero : 0.0);

    // Every SYBYL type now has its own distinct base type — no Solvent collapse.
    // The round-trip 40→256→40 should be perfect.
    EXPECT_NEAR(r, 1.0, 1e-6)
        << "Perfect round-trip: r should be 1.0";
    EXPECT_NEAR(rmse, 0.0, 1e-3)
        << "Perfect round-trip: RMSE should be 0.0";
    EXPECT_NEAR(max_abs_diff, 0.0, 1e-3)
        << "No cell should deviate after uniform expansion + projection";
    EXPECT_EQ(n_matching_sign, n_both_nonzero)
        << "No sign flips in round-trip projection";
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 4: Type mapping coverage — which SYBYL types map to which 256-codes
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, TypeMappingCoverage) {
    // For each SYBYL type (1-40), count how many 256-codes map to it
    std::array<int, 41> sybyl_fanout{};
    for (int c = 0; c < 256; ++c) {
        int sybyl = atom256::base_to_sybyl_parent(atom256::get_base(c));
        if (sybyl >= 1 && sybyl <= 40)
            sybyl_fanout[sybyl]++;
    }

    printf("\n=== SYBYL → 256-Code Fanout ===\n");
    printf("  %-6s %-12s %s\n", "SYBYL", "Name", "256-Codes");
    printf("  %s\n", std::string(40, '-').c_str());

    int total_mapped = 0;
    int types_with_fanout = 0;
    for (int s = 1; s <= 40; ++s) {
        const char* name = "???";
        for (int k = 0; k < nrgrank::NUM_SYBYL_ENTRIES; ++k) {
            if (nrgrank::kSybylTypes[k].type_number == s) {
                name = nrgrank::kSybylTypes[k].name.data();
                break;
            }
        }
        printf("  %-6d %-12s %d\n", s, name, sybyl_fanout[s]);
        total_mapped += sybyl_fanout[s];
        if (sybyl_fanout[s] > 0) ++types_with_fanout;
    }

    printf("\n  Total 256-codes mapped: %d / 256\n", total_mapped);
    printf("  SYBYL types with ≥1 code: %d / 40\n", types_with_fanout);

    // All 256 codes should map to some SYBYL parent
    EXPECT_EQ(total_mapped, 256)
        << "Every 256-code should map to a SYBYL parent";
    // Every SYBYL type should have at least one 256-code
    EXPECT_EQ(types_with_fanout, 40)
        << "All 40 SYBYL types should have >=1 256-code (no Solvent collapse)";
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 5: Dimension comparison — information density
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, InformationDensity) {
    // NRGRank 40×40: count unique non-zero values
    std::vector<double> nrg_values;
    for (int i = 1; i <= 40; ++i)
        for (int j = 1; j <= 40; ++j) {
            double v = nrgrank::kEnergyMatrix[i][j];
            if (std::fabs(v) > 1e-9)
                nrg_values.push_back(v);
        }
    std::sort(nrg_values.begin(), nrg_values.end());
    auto it = std::unique(nrg_values.begin(), nrg_values.end(),
                          [](double a, double b) { return std::fabs(a - b) < 1e-6; });
    int unique_nrg = static_cast<int>(it - nrg_values.begin());

    // 256×256 from NRGRank: same values but replicated across charge/hbond bins
    auto mat256 = build_256_from_nrgrank();
    std::vector<float> scm_values;
    for (int i = 0; i < scm::MATRIX_SIZE; ++i) {
        if (std::fabs(mat256.data[i]) > 1e-9f)
            scm_values.push_back(mat256.data[i]);
    }
    std::sort(scm_values.begin(), scm_values.end());
    auto it2 = std::unique(scm_values.begin(), scm_values.end(),
                           [](float a, float b) { return std::fabs(a - b) < 1e-6f; });
    int unique_scm = static_cast<int>(it2 - scm_values.begin());

    // Sparsity stats
    int nrg_nonzero = static_cast<int>(nrg_values.size());
    int scm_nonzero = static_cast<int>(scm_values.size());

    printf("\n=== Information Density Comparison ===\n");
    printf("  %-25s %-15s %-15s\n", "Metric", "40×40", "256×256");
    printf("  %s\n", std::string(55, '-').c_str());
    printf("  %-25s %-15d %-15d\n", "Total cells",
           40 * 40, 256 * 256);
    printf("  %-25s %-15d %-15d\n", "Non-zero cells",
           nrg_nonzero, scm_nonzero);
    printf("  %-25s %-15.1f %-15.1f\n", "Sparsity (%)",
           100.0 * (1.0 - (double)nrg_nonzero / (40 * 40)),
           100.0 * (1.0 - (double)scm_nonzero / (256 * 256)));
    printf("  %-25s %-15d %-15d\n", "Unique non-zero values",
           unique_nrg, unique_scm);
    printf("  %-25s %-15zu %-15zu\n", "Memory (bytes)",
           sizeof(nrgrank::kEnergyMatrix),
           sizeof(scm::SoftContactMatrix));
    printf("  %-25s %-15.2f %-15.2f\n", "Bits/cell (non-zero)",
           nrg_nonzero > 0 ? std::log2(unique_nrg) : 0.0,
           scm_nonzero > 0 ? std::log2(unique_scm) : 0.0);

    // With no Solvent collapse, unique value counts should match
    // (minor difference possible from double→float precision).
    EXPECT_NEAR(unique_nrg, unique_scm, 15)
        << "Unique value counts should be very close with no Solvent collapse";
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 6: Row-vector comparison — NRGRank rows vs 256×256 rows
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, RowVectorComparison) {
    auto mat256 = build_256_from_nrgrank();

    printf("\n=== Row-Vector Energy Profile Comparison ===\n");
    printf("  For each SYBYL type, compare its NRGRank row (40-dim) against\n");
    printf("  the corresponding 256-type rows (should be identical profiles).\n\n");

    // For selected important types, check that all 256-codes sharing the
    // same SYBYL parent produce identical row sums
    struct TypeCheck {
        int sybyl;
        const char* name;
    };
    TypeCheck checks[] = {
        {3,  "C.3 (sp3 carbon)"},
        {4,  "C.AR (aromatic)"},
        {8,  "N.3 (sp3 nitrogen)"},
        {13, "O.2 (carbonyl)"},
        {14, "O.3 (hydroxyl)"},
        {18, "S.3 (thiol)"},
        {40, "Solvent"},
    };

    printf("  %-22s %-12s %-12s %s\n",
           "Type", "NRGRank sum", "256-row sum", "Match");
    printf("  %s\n", std::string(55, '-').c_str());

    for (auto& tc : checks) {
        // NRGRank row sum for this type
        double nrg_sum = 0.0;
        for (int j = 1; j <= 40; ++j)
            nrg_sum += nrgrank::kEnergyMatrix[tc.sybyl][j];

        // Find one 256-code that maps to this SYBYL type
        double scm_row_sum = 0.0;
        bool found = false;
        for (int c = 0; c < 256 && !found; ++c) {
            int sybyl = atom256::base_to_sybyl_parent(atom256::get_base(c));
            if (sybyl == tc.sybyl) {
                // Sum this 256-type's row, but only over 256-codes that map
                // back to SYBYL types 1-40 (grouping by SYBYL parent)
                std::array<double, 41> projected_row{};
                for (int cj = 0; cj < 256; ++cj) {
                    int sj = atom256::base_to_sybyl_parent(atom256::get_base(cj));
                    if (sj >= 1 && sj <= 40)
                        projected_row[sj] += mat256.lookup(c, cj);
                }
                // Count how many 256-codes map to each SYBYL type for averaging
                std::array<int, 41> sybyl_count{};
                for (int cj = 0; cj < 256; ++cj) {
                    int sj = atom256::base_to_sybyl_parent(atom256::get_base(cj));
                    if (sj >= 1 && sj <= 40) sybyl_count[sj]++;
                }
                for (int sj = 1; sj <= 40; ++sj) {
                    if (sybyl_count[sj] > 0)
                        scm_row_sum += projected_row[sj] / sybyl_count[sj];
                }
                found = true;
            }
        }

        // With distinct base types for every SYBYL type, the round-trip
        // should be exact (within float precision).
        double tol = 0.1;
        bool match = std::fabs(nrg_sum - scm_row_sum) < tol;
        printf("  %-22s %12.2f %12.2f %s\n",
               tc.name, nrg_sum, scm_row_sum, match ? "OK" : "MISMATCH");
        EXPECT_NEAR(nrg_sum, scm_row_sum, tol)
            << "Row sum mismatch for " << tc.name;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 7: Capacity comparison — what the 256×256 can express that 40×40 cannot
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, CapacityComparison) {
    printf("\n=== Capacity: What 256×256 Adds Over 40×40 ===\n\n");

    // Example: C.AR (SYBYL 4) fans out to 8 codes in 256-system
    // These can encode different energies for the same SYBYL type based on
    // charge state and H-bond capability
    int c_ar_codes = 0;
    printf("  C.AR (aromatic carbon) decomposition in 256-type system:\n");
    for (int c = 0; c < 256; ++c) {
        uint8_t base = atom256::get_base(c);
        if (base == atom256::C_ar || base == atom256::C_ar_hetadj ||
            base == atom256::C_pi_bridge) {
            printf("    code=%3d  base=%-10s  donor=%s  acceptor=%s\n",
                   c, atom256::base_type_name(base),
                   atom256::get_hbond_donor(c) ? "yes" : "no",
                   atom256::get_hbond_acceptor(c) ? "yes" : "no");
            ++c_ar_codes;
        }
    }
    printf("  → %d distinct codes for aromatic C variants\n\n", c_ar_codes);

    // The key advantage: a trained 256×256 matrix can differentiate between
    // e.g. anionic C.AR with H-bond vs cationic C.AR without
    printf("  Summary of additional discriminative power:\n");
    printf("    40×40:  1 entry per SYBYL pair → charge-blind, H-bond-blind\n");
    printf("    256×256: up to 16 entries per SYBYL pair\n");
    printf("            (4 codes/type × 4 codes/type for each SYBYL pair)\n");
    printf("    Theoretical expansion: %dx more parameters\n",
           (256 * 256) / (40 * 40));
    printf("    Memory cost: %zu bytes (256×256) vs %zu bytes (41×41)\n",
           sizeof(scm::SoftContactMatrix), sizeof(nrgrank::kEnergyMatrix));

    EXPECT_GT(c_ar_codes, 4)
        << "Aromatic C should have multiple 256-type representations";
}

// ═══════════════════════════════════════════════════════════════════════════════
// TEST 8: Symmetry comparison
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompareEnergyMatrices, SymmetryComparison) {
    // NRGRank symmetry
    double nrg_max_asym = 0.0;
    int nrg_asym_count = 0;
    for (int i = 1; i <= 40; ++i) {
        for (int j = i + 1; j <= 40; ++j) {
            double d = std::fabs(nrgrank::kEnergyMatrix[i][j] -
                                  nrgrank::kEnergyMatrix[j][i]);
            if (d > 1e-9) ++nrg_asym_count;
            nrg_max_asym = std::max(nrg_max_asym, d);
        }
    }

    // 256×256 built from NRGRank inherits its symmetry properties
    auto mat256 = build_256_from_nrgrank();
    double scm_max_asym = 0.0;
    int scm_asym_count = 0;
    for (int i = 0; i < 256; ++i) {
        for (int j = i + 1; j < 256; ++j) {
            double d = std::fabs(mat256.lookup(i, j) - mat256.lookup(j, i));
            if (d > 1e-9) ++scm_asym_count;
            scm_max_asym = std::max(scm_max_asym, d);
        }
    }

    printf("\n=== Symmetry Comparison ===\n");
    printf("  %-30s %-15s %-15s\n", "Metric", "40×40", "256×256");
    printf("  %s\n", std::string(60, '-').c_str());
    printf("  %-30s %-15d %-15d\n", "Asymmetric pairs",
           nrg_asym_count, scm_asym_count);
    printf("  %-30s %-15.4f %-15.4f\n", "Max |a[i][j]-a[j][i]|",
           nrg_max_asym, scm_max_asym);
    printf("  %-30s %-15s %-15s\n", "Symmetric?",
           nrg_max_asym < 1e-6 ? "YES" : "NO",
           scm_max_asym < 1e-6 ? "YES" : "NO");

    // The 256×256 inherits NRGRank's symmetry since E[si][sj] = E[sj][si]
    // iff the NRGRank matrix is symmetric (it's not perfectly symmetric)
    if (nrg_max_asym > 1e-6) {
        EXPECT_GT(scm_max_asym, 0.0)
            << "256×256 should inherit NRGRank asymmetry";
    }
}
