// get_yval.cpp — energy-matrix interpolation + opt-in LUT
// SPDX-License-Identifier: Apache-2.0

#include "get_yval.h"

#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

double get_yval_scan(struct energy_matrix* em, double relative_area)
{
	if(!em || !em->energy_values) return 0.0;
	if(em->weight) return (double)em->energy_values->y;
	const int n = em->flat_n;
	if(n == 0) return 0.0;
	const float ra = (float)relative_area;
	const float* fx = em->flat_x;
	const float* fy = em->flat_y;
	if(!fx || !fy) return 0.0;
	if(ra < fx[0]) return 0.0;
	if(ra >= fx[n-1]) return (double)fy[n-1];
	if(n < 2 || !em->flat_slope) return 0.0;
	// Binary search for the same i as `while (i < n-2 && ra >= fx[i+1]) ++i`.
	// NaN fails both early-outs (every comparison is false) and the historical
	// linear scan left i = 0; keep that.
	int i = 0;
	if (ra == ra) {
		const float* ub = std::upper_bound(fx, fx + n, ra);
		i = static_cast<int>(ub - fx) - 1;
		if (i < 0) i = 0;
		if (i > n - 2) i = n - 2;
	}
	return (double)(fy[i] + em->flat_slope[i] * (ra - fx[i]));
}

namespace {

struct Lut {
	int bins = 0;
	std::vector<double> y;
};

using LutMap = std::unordered_map<const energy_matrix*, std::shared_ptr<const Lut>>;

// Copy-on-write map of per-matrix LUTs. After a key is published, readers
// atomic-load the snapshot pointer and hash-find with no mutex (ShannonMetalBridge
// once-then-immutable pattern). g_publish_mu is only for first-insert of a
// new energy_matrix* (fixed per dock). Snapshots are immortal so a reader
// that loaded an older pointer cannot see it freed. libc++ on this SDK has
// no std::atomic<std::shared_ptr<T>> (requires trivially copyable T).
std::atomic<const LutMap*> g_luts{nullptr};
std::mutex g_publish_mu;
std::vector<std::unique_ptr<const LutMap>> g_keep_alive;

std::shared_ptr<const Lut> lut_for(energy_matrix* em, int bins)
{
	if (const LutMap* map = g_luts.load(std::memory_order_acquire)) {
		auto it = map->find(em);
		if (it != map->end() && it->second && it->second->bins == bins)
			return it->second;
	}

	auto built = std::make_shared<Lut>();
	built->bins = bins;
	built->y.resize(static_cast<size_t>(bins));
	const double denom = static_cast<double>(bins - 1);
	for (int b = 0; b < bins; ++b) {
		const double ra = static_cast<double>(b) / denom;
		built->y[static_cast<size_t>(b)] = get_yval_scan(em, ra);
	}

	std::lock_guard<std::mutex> lock(g_publish_mu);
	if (const LutMap* map = g_luts.load(std::memory_order_relaxed)) {
		auto it = map->find(em);
		if (it != map->end() && it->second && it->second->bins == bins)
			return it->second;
		auto next = std::make_unique<LutMap>(*map);
		(*next)[em] = built;
		const LutMap* raw = next.get();
		g_keep_alive.push_back(std::move(next));
		g_luts.store(raw, std::memory_order_release);
		return built;
	}
	auto next = std::make_unique<LutMap>();
	(*next)[em] = built;
	const LutMap* raw = next.get();
	g_keep_alive.push_back(std::move(next));
	g_luts.store(raw, std::memory_order_release);
	return built;
}

// Lock-free const reads after lut_for() returns. Copy y0/y1 out so a later
// COW republish cannot invalidate the caller's references.
bool lut_samples(energy_matrix* em, int bins, int i, double* y0, double* y1)
{
	auto lut = lut_for(em, bins);
	if (!lut || lut->y.empty()) return false;
	const auto& y = lut->y;
	if (i >= bins - 1) {
		*y0 = y.back();
		*y1 = y.back();
		return true;
	}
	if (i < 0) i = 0;
	*y0 = y[static_cast<size_t>(i)];
	*y1 = y[static_cast<size_t>(i + 1)];
	return true;
}

}  // namespace

double get_yval(struct energy_matrix* em, double relative_area, bool use_lut)
{
	if (!use_lut)
		return get_yval_scan(em, relative_area);
	if (!em) return 0.0;
	if (relative_area < 0.0 || relative_area > 1.0)
		return get_yval_scan(em, relative_area);
	const int bins = flexaids::get_yval_lut_bins_cached();
	const double t = relative_area * static_cast<double>(bins - 1);
	int i = static_cast<int>(t);
	if (i < 0) i = 0;
	double y0 = 0.0, y1 = 0.0;
	if (!lut_samples(em, bins, i, &y0, &y1))
		return get_yval_scan(em, relative_area);
	if (i >= bins - 1) return y0;
	const double f = t - static_cast<double>(i);
	return y0 * (1.0 - f) + y1 * f;
}

double get_yval(struct energy_matrix* em, double relative_area)
{
	return get_yval(em, relative_area, flexaids::get_yval_lut_enabled_cached());
}
