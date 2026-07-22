#ifndef FLEXAIDDS_CLUSTER_REP_MODE_H
#define FLEXAIDDS_CLUSTER_REP_MODE_H

// ─────────────────────────────────────────────────────────────────────────────
// Within-cluster rank-0 representative election mode (P2).
//
// Single primary env gate, DEFAULT-OFF, so every clustering backend (CF/leader,
// DensityPeak, FastOPTICS) moves together and the unset default reproduces
// pre-medoid HEAD behavior bit-for-bit:
//
//   FLEXAIDDS_CLUSTER_REP
//     unset | "lowcf"  → lowest-CF head / Representative   (DEFAULT, bit-identical)
//     "medoid"         → pure UNWEIGHTED geometric medoid  (CF-independent)
//     "bmedoid"        → Boltzmann-CF-weighted medoid       (HEAD variant; ablation only)
//     "center"         → density-peak / OPTICS center       (where the backend has one)
//
// Rationale: within-target Spearman(CF, RMSD) ≈ 0, so CF is orthogonal to pose
// correctness. The pure geometric medoid selects the pose most central to a
// cluster's geometry — CF plays no role. "bmedoid" re-injects CF via Boltzmann
// weights and is retained only for ablation. "lowcf" is the historical default.
//
// Back-compat: the legacy FLEXAIDDS_MEDOID_REFINE var (which HEAD 3e674479c
// defaulted ON — a constraint violation this work order corrects) is honored
// ONLY when explicitly set non-zero, aliasing to "bmedoid" with a one-time
// deprecation notice. When both vars are unset the mode is LOWCF (no medoid).
// FLEXAIDDS_CLUSTER_REP always takes precedence over the legacy var.
// ─────────────────────────────────────────────────────────────────────────────

#include <cstdlib>
#include <cstdio>
#include <strings.h>   // strcasecmp

namespace flexaids {

enum class ClusterRepMode { LOWCF, MEDOID, BMEDOID, CENTER };

inline ClusterRepMode cluster_rep_mode()
{
	const char* v = std::getenv("FLEXAIDDS_CLUSTER_REP");
	if (v && *v) {
		if (!strcasecmp(v, "lowcf"))   return ClusterRepMode::LOWCF;
		if (!strcasecmp(v, "medoid"))  return ClusterRepMode::MEDOID;
		if (!strcasecmp(v, "bmedoid")) return ClusterRepMode::BMEDOID;
		if (!strcasecmp(v, "center"))  return ClusterRepMode::CENTER;
		fprintf(stderr,
		        "WARNING: FLEXAIDDS_CLUSTER_REP='%s' unrecognized "
		        "(expected lowcf|medoid|bmedoid|center); using default 'lowcf'.\n",
		        v);
		return ClusterRepMode::LOWCF;
	}

	// Legacy alias — only when explicitly enabled. Default-ON behavior is removed.
	const char* legacy = std::getenv("FLEXAIDDS_MEDOID_REFINE");
	if (legacy && std::atoi(legacy) != 0) {
		static bool warned = false;
		if (!warned) {
			fprintf(stderr,
			        "WARNING: FLEXAIDDS_MEDOID_REFINE is deprecated; treating as "
			        "FLEXAIDDS_CLUSTER_REP=bmedoid (Boltzmann-weighted medoid, "
			        "ablation only). Prefer FLEXAIDDS_CLUSTER_REP.\n");
			warned = true;
		}
		return ClusterRepMode::BMEDOID;
	}

	return ClusterRepMode::LOWCF;
}

inline const char* cluster_rep_mode_name(ClusterRepMode m)
{
	switch (m) {
		case ClusterRepMode::LOWCF:   return "lowcf";
		case ClusterRepMode::MEDOID:  return "medoid";
		case ClusterRepMode::BMEDOID: return "bmedoid";
		case ClusterRepMode::CENTER:  return "center";
	}
	return "lowcf";
}

} // namespace flexaids

#endif // FLEXAIDDS_CLUSTER_REP_MODE_H
