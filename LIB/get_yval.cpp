// get_yval.cpp — energy-matrix interpolation + opt-in LUT
// SPDX-License-Identifier: Apache-2.0

#include "get_yval.h"

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
	int i = 0;
	while(i < n-2 && ra >= fx[i+1]) ++i;
	return (double)(fy[i] + em->flat_slope[i] * (ra - fx[i]));
}

namespace {

struct Lut {
	int bins = 0;
	std::vector<double> y;
};

std::mutex g_mu;
std::unordered_map<const energy_matrix*, Lut> g_luts;

// Copy the two lerp samples (or the last bin) while the mutex is held so a
// concurrent first-insert of another energy_matrix* cannot rehash g_luts
// under a live reference (OpenMP vcfunction).
bool lut_samples(energy_matrix* em, int bins, int i, double* y0, double* y1)
{
	std::lock_guard<std::mutex> lock(g_mu);
	Lut& lut = g_luts[em];
	if (lut.bins != bins || (int)lut.y.size() != bins) {
		lut.bins = bins;
		lut.y.resize(static_cast<size_t>(bins));
		const double denom = static_cast<double>(bins - 1);
		for (int b = 0; b < bins; ++b) {
			const double ra = static_cast<double>(b) / denom;
			lut.y[static_cast<size_t>(b)] = get_yval_scan(em, ra);
		}
	}
	if (lut.y.empty()) return false;
	if (i >= bins - 1) {
		*y0 = lut.y.back();
		*y1 = lut.y.back();
		return true;
	}
	if (i < 0) i = 0;
	*y0 = lut.y[static_cast<size_t>(i)];
	*y1 = lut.y[static_cast<size_t>(i + 1)];
	return true;
}

}  // namespace

double get_yval(struct energy_matrix* em, double relative_area)
{
	if (!flexaids::get_yval_lut_enabled())
		return get_yval_scan(em, relative_area);
	if (!em) return 0.0;
	if (relative_area < 0.0 || relative_area > 1.0)
		return get_yval_scan(em, relative_area);
	const int bins = flexaids::get_yval_lut_bins();
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
