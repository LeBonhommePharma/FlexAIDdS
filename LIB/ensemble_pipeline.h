#pragma once
// =============================================================================
// ensemble_pipeline.h — 4-layer FlexAIDdS ensemble generation contract
//
// [1] Frame chart consistency   Cartesian ⇄ genes identity (CI gate)
// [2] Pocket support Ω_cleft    ligandable top-K · cleft confinement · spheres
// [3] Soft-β CF sampling        SMFREE β=1/T · niche diversity · multi-start
// [4] Classic entropy election  ACF / BindingMode H−T·S (see classic_entropy_ranking.md)
//
// Layers 1–2 define the geometry of the binding integral; layer 3 samples the
// CF measure on that support; layer 4 elects the basin. Never mix physical-kB
// Boltzmann weights into gene search, never elect on a different objective
// than you sample, and never sample outside the site support.
//
// Pure helpers live here so unit tests do not need a full GA binary.
// Apache-2.0 © 2026 Le Bonhomme Pharma
// =============================================================================

#include "flexaid.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>

namespace ensemble {

// ── Layer 1: frame chart ────────────────────────────────────────────────────

/// Soft warn threshold (Å) for native-seed IC→Cartesian round-trip.
inline constexpr double kFrameChartWarnRmsdA = 1.0;
/// Strict CI / product gate (Å). Restore identity-class chart (v34_ctrl ~0.00).
inline constexpr double kFrameChartStrictRmsdA = 0.1;

/// True when FLEXAIDDS_FRAME_CHART_STRICT is 1/true/yes.
inline bool frame_chart_strict_enabled() noexcept
{
	const char* e = std::getenv("FLEXAIDDS_FRAME_CHART_STRICT");
	if (!e || !e[0]) return false;
	return e[0] == '1' || e[0] == 't' || e[0] == 'T' || e[0] == 'y' || e[0] == 'Y';
}

/// Classify a measured native-seed RMSD.
enum class FrameChartStatus { Ok, Warn, Fail };

inline FrameChartStatus classify_frame_chart_rmsd(double rmsd_A,
                                                  bool strict) noexcept
{
	if (!(rmsd_A >= 0.0) || !std::isfinite(rmsd_A))
		return FrameChartStatus::Fail;
	// Strict CI: identity-class chart (0.1 Å). Soft path still warns at 1.0 Å.
	if (strict && rmsd_A > kFrameChartStrictRmsdA)
		return FrameChartStatus::Fail;
	if (rmsd_A > kFrameChartWarnRmsdA)
		return FrameChartStatus::Warn;
	return FrameChartStatus::Ok;
}

/// Gene-limit chart invertibility residual: |IC − genetoic(ictogene(IC))|.
/// Returns max residual over samples; identity-class charts stay ≤ one bin.
template <typename GeneToIcFn, typename IcToGeneFn>
inline double gene_chart_max_residual(double min_ic, double max_ic, double del,
                                      int n_samples,
                                      GeneToIcFn genetoic_fn,
                                      IcToGeneFn ictogene_fn)
{
	if (n_samples < 2 || del <= 0.0 || max_ic <= min_ic)
		return std::numeric_limits<double>::infinity();
	double worst = 0.0;
	for (int i = 0; i < n_samples; ++i) {
		const double t = static_cast<double>(i) / static_cast<double>(n_samples - 1);
		const double ic = min_ic + t * (max_ic - min_ic);
		const auto g = ictogene_fn(ic);
		const double ic2 = genetoic_fn(g);
		worst = std::max(worst, std::abs(ic2 - ic));
	}
	return worst;
}

// ── Layer 2: pocket support Ω_cleft ─────────────────────────────────────────

struct CleftCentroid {
	double cx = 0.0, cy = 0.0, cz = 0.0;
	/// Max |center − centroid| + radius over spheres (Å).
	double extent_A = 0.0;
	int n_spheres = 0;
};

/// Centroid and radial extent of a GetCleft / CleftDetector sphere list.
inline bool cleft_centroid_extent(const sphere* spheres, CleftCentroid* out) noexcept
{
	if (!out) return false;
	out->cx = out->cy = out->cz = 0.0;
	out->extent_A = 0.0;
	out->n_spheres = 0;
	for (const sphere* s = spheres; s; s = s->prev) {
		out->cx += s->center[0];
		out->cy += s->center[1];
		out->cz += s->center[2];
		++out->n_spheres;
	}
	if (out->n_spheres <= 0) return false;
	const double inv = 1.0 / static_cast<double>(out->n_spheres);
	out->cx *= inv;
	out->cy *= inv;
	out->cz *= inv;
	for (const sphere* s = spheres; s; s = s->prev) {
		const double dx = s->center[0] - out->cx;
		const double dy = s->center[1] - out->cy;
		const double dz = s->center[2] - out->cz;
		const double reach = std::sqrt(dx * dx + dy * dy + dz * dz)
		                     + static_cast<double>(s->radius);
		if (reach > out->extent_A) out->extent_A = reach;
	}
	return true;
}

/// Cheap ligandability proxy for a SURFNET sphere cluster (no ML).
/// volume ~ n · ⟨r⟩³, enclosure ~ ⟨r⟩ / (½·bbox diagonal), product is s.
inline double ligandable_score(int n_spheres,
                               double mean_radius,
                               double bbox_diag) noexcept
{
	if (n_spheres <= 0 || mean_radius <= 0.0) return 0.0;
	const double vol = static_cast<double>(n_spheres)
	                   * mean_radius * mean_radius * mean_radius;
	const double half_diag = 0.5 * std::max(bbox_diag, mean_radius);
	const double enclosure = mean_radius / half_diag;  // compact → higher
	return vol * enclosure;
}

/// Keep at most top_k cluster labels by descending ligandable score.
/// clusters: vector of (label, score). Returns labels to keep.
inline std::vector<int> select_top_k_clefts(
	std::vector<std::pair<int, double>> clusters,
	int top_k)
{
	if (top_k <= 0) return {};
	std::stable_sort(clusters.begin(), clusters.end(),
		[](const auto& a, const auto& b) { return a.second > b.second; });
	std::vector<int> keep;
	keep.reserve(static_cast<size_t>(std::min(top_k, static_cast<int>(clusters.size()))));
	for (const auto& kv : clusters) {
		if (static_cast<int>(keep.size()) >= top_k) break;
		keep.push_back(kv.first);
	}
	return keep;
}

/// Sphere radius accepted by pocket grid builders.
inline bool valid_sphere_radius(float radius) noexcept
{
	return std::isfinite(radius) && radius > 0.3f && radius < 50.0f;
}

// ── Layer 3: soft-β SMFREE sampling ─────────────────────────────────────────

/// Soft selection β = 1/T (CF units). Physical 1/(kB·T) is forbidden for search.
inline bool soft_selection_beta(double temperature_K, double* beta_out) noexcept
{
	if (!beta_out || !(temperature_K > 0.0) || !std::isfinite(temperature_K))
		return false;
	*beta_out = 1.0 / temperature_K;
	return true;
}

/// SMFREE product path requires T>0; otherwise fitness collapses to rank-only.
inline bool smfree_requires_temperature(double temperature_K) noexcept
{
	return temperature_K > 0.0;
}

// ── Layer 4: classic entropy election (mirrors cluster.cpp gate) ────────────

/// Return emission index for rank-0 under classic vs force_cf policy.
inline int elect_rank0_index(const std::vector<double>& acf,
                             const std::vector<double>& cf,
                             bool force_cf_rank_emission,
                             unsigned temperature) noexcept
{
	const int n = static_cast<int>(acf.size());
	if (n == 0 || static_cast<int>(cf.size()) != n) return -1;
	std::vector<int> order(static_cast<size_t>(n));
	for (int i = 0; i < n; ++i) order[static_cast<size_t>(i)] = i;
	if (temperature > 0) {
		std::stable_sort(order.begin(), order.end(),
			[&](int a, int b) { return acf[static_cast<size_t>(a)] < acf[static_cast<size_t>(b)]; });
	}
	const bool classic = (temperature > 0) && !force_cf_rank_emission;
	if (!classic) {
		std::stable_sort(order.begin(), order.end(),
			[&](int a, int b) { return cf[static_cast<size_t>(a)] < cf[static_cast<size_t>(b)]; });
	}
	return order[0];
}

}  // namespace ensemble
