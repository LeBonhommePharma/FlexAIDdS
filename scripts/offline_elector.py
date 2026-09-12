#!/usr/bin/env python3
"""
offline_elector.py -- faithful offline reimplementation of the two electors.

Historical port from the frozen elector_oos_20260911 campaign. Source locations
and arm constants below document that extraction, not the current checkout or
every protocol. This successor makes the per-restart budget explicit.

TRANSCRIBED FROM SOURCE (campaign's FlexAIDdS main, LIB/):

  default elector : DatasetRunner.cpp:1452 select_pose_freq_gated_pooled()
  min-CF elector  : DatasetRunner.cpp:1629 (FLEXAIDDS_POOLED_ELECTION=mincf
                    short-circuit, committed 6dfe81ab)
  enumerator      : DatasetRunner.cpp:163  enumerate_emitted_cluster_heads()
  soft-beta G     : SoftBetaFreeEnergy.h:68 free_energy()
  prefix order    : DatasetRunner.cpp:7759 / 7796-7807  (r0, r1, r2)

ELECTION IS POST-HOC OVER FILES ON DISK. Both electors read only
<prefix>_<rank>.pdb ("REMARK CF=", "Frequency:") and the <prefix>_<rank>.mcf
sidecar. No engine state is consulted. pooled_election.h:91-93 states this
directly: "the engine never reads this variable -- it is read by the HARNESS
only. The gate re-elects from poses already on disk."

ARM-SPECIFIC CONSTANTS, read from the arms' own RUN_RECEIPT.json:
  election_shannon_free_energy = true   -> Soft-beta objective ON
  election_soft_T              = 0.0    -> falls through to dock T = 300.0
  election_include_singletons  = false  -> but use_shannon_G forces it TRUE
                                           (DatasetRunner.cpp:1485-1486)
  seed_elitism (override)      = 0      -> DEFINED_CLEFT_REDOCK forces seeds
                                           OFF (DatasetRunner.cpp:8062-8066)
  cluster_spread_max           = 0.0    -> spread guard OFF
  freqsel                      = false  -> macro-cluster reranking OFF
  cf_window_selector           = false
"""
import math
import os

# Successor to the frozen campaign port; election/scoring math is unchanged.
# The caller owns the budget. Ambient engine environment never changes an audit.

POSE_BUDGET = 50          # benchmark_pose_limit(), FLEXAIDDS_MAX_RESULTS unset
POSE_BUDGET_MAX = 5000
SENTINEL_CF = 1e3         # engine flags sentinel poses with CF >= 1e3


def validate_pose_budget(budget):
    """Reject invalid offline settings instead of silently scoring a smaller pool."""
    if type(budget) is not int or not 1 <= budget <= POSE_BUDGET_MAX:
        raise ValueError("pose budget must be an integer in [1, 5000]")
    return budget


# ---------------------------------------------------------------- enumeration
def enumerate_heads(prefix, budget=POSE_BUDGET):
    """Port of enumerate_emitted_cluster_heads(). Returns [(path, rank, min_pts)].

    CF/DP single-suffix first (rank order), then FastOPTICS dual-suffix ordered
    by (min_pts, rank, path), the whole thing bounded by one shared budget.
    """
    budget = validate_pose_budget(budget)
    out, seen = [], set()

    def try_add(path, rank, min_pts):
        if not path or path in seen or not os.path.exists(path):
            return False
        seen.add(path)
        out.append((path, rank, min_pts))
        return True

    cf_found = 0
    for pi in range(budget):
        if try_add("%s_%d.pdb" % (prefix, pi), pi, -1):
            cf_found += 1
    cf_truncated = os.path.exists("%s_%d.pdb" % (prefix, budget))

    fo = []
    d = os.path.dirname(prefix) or "."
    base = os.path.basename(prefix)
    pm = base + "_"
    if os.path.isdir(d):
        for fname in os.listdir(d):
            full = os.path.join(d, fname)
            if not os.path.isfile(full):
                continue
            if len(fname) < len(pm) + 6 or not fname.startswith(pm):
                continue
            if not fname.endswith(".pdb"):
                continue
            mid = fname[len(pm):-4]
            us = mid.find("_")
            if us in (-1, 0) or us + 1 >= len(mid):
                continue
            if mid.find("_", us + 1) != -1:
                continue
            if not all(mid[k].isdigit() for k in range(len(mid)) if k != us):
                continue
            try:
                min_pts, rank = int(mid[:us]), int(mid[us + 1:])
            except ValueError:
                continue
            if min_pts < 0 or rank < 0:
                continue
            fo.append((full, rank, min_pts))
    fo.sort(key=lambda h: (h[2], h[1], h[0]))
    fo_kept = 0
    for h in fo:
        if len(out) >= budget:
            break
        if try_add(*h):
            fo_kept += 1

    return out, {"budget": budget, "cf_found": cf_found, "fo_found": len(fo), "fo_kept": fo_kept,
                 "enumerated": len(out), "cf_truncated": bool(cf_truncated),
                 "fo_truncated": fo_kept < len(fo)}


# ------------------------------------------------------------------- parsing
def parse_pose(path):
    """Port of the parse_pose lambda (DatasetRunner.cpp:1542).

    C++ scans the whole file: FIRST "REMARK CF=" wins, LAST "Frequency:" wins.
    Measured on these arms the REMARK block is lines 1-65 with exactly one of
    each, so a header read is exact -- but we VERIFY we cleared the REMARK block
    (saw a coordinate record) and fall back to a full scan if not.
    Returns dict or None when no finite CF could be read.
    """
    cf, freq, have_cf, saw_atom = math.inf, 1, False, False
    with open(path, "r", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if not have_cf and "REMARK CF=" in line:
                p2 = line.find("CF=")
                if p2 != -1:
                    try:
                        cf = float(line[p2 + 3:].split()[0])
                        have_cf = True
                    except (ValueError, IndexError):
                        pass
            elif "Frequency:" in line:
                p2 = line.find("Frequency:")
                try:
                    freq = int(line[p2 + 10:].split()[0])
                except (ValueError, IndexError):
                    pass
            if line.startswith(("ATOM  ", "HETATM")):
                saw_atom = True
            if i >= 200 and saw_atom and have_cf:
                break
    if not have_cf or not math.isfinite(cf):
        return None

    member_cfs = []
    mcf = path[:-4] + ".mcf" if path.endswith(".pdb") else path + ".mcf"
    if os.path.exists(mcf):
        with open(mcf, "r", errors="ignore") as fh:
            for ml in fh:
                try:
                    v = float(ml)
                except ValueError:
                    continue
                if math.isfinite(v):
                    member_cfs.append(v)
    return {"path": path, "cf": cf, "freq": freq, "member_cfs": member_cfs}


def build_pool(prefixes, budget=POSE_BUDGET):
    """The `poses` vector of select_pose_freq_gated_pooled(), in engine order."""
    budget = validate_pose_budget(budget)
    poses, enum_stats = [], []
    for ri, prefix in enumerate(prefixes):
        heads, st = enumerate_heads(prefix, budget=budget)
        st["restart"] = ri
        st["prefix"] = os.fspath(prefix)
        st["parsed"] = 0
        st["parse_errors"] = []
        enum_stats.append(st)
        for path, rank, min_pts in heads:
            p = parse_pose(path)
            if p is None:
                st["parse_errors"].append({"path": path, "reason": "missing_or_nonfinite_cf"})
                continue
            st["parsed"] += 1
            p.update(restart=ri, rank=rank, min_pts=min_pts, is_seed=False)
            poses.append(p)
    return poses, enum_stats


# ------------------------------------------------------------------- scoring
def soft_free_energy(pose, soft_T):
    """Port of SoftBetaFreeEnergy.h free_energy(): G = Emin - T*ln Z."""
    energies = [e for e in pose["member_cfs"] if math.isfinite(e)]
    if not energies:
        energies = [pose["cf"]]          # singleton: S=0, G=CF
    finite = [e for e in energies if math.isfinite(e)]
    if not finite:
        return math.inf
    T = soft_T if soft_T > 1e-12 else 1e-12
    Emin = min(finite)
    Z = sum(math.exp(-(e - Emin) / T) for e in finite)
    if not (Z > 0.0) or not math.isfinite(Z):
        return Emin
    return Emin - T * math.log(Z)


# ------------------------------------------------------------------ electors
def elect_default(poses, soft_T=300.0, use_shannon_G=True):
    """Port of select_pose_freq_gated_pooled() under THESE arms' config.

    include_singletons is TRUE because use_shannon_G is true, so the freq>1 gate
    admits everything; seed elitism is forced OFF; spread guard and freqsel are
    OFF. What remains is: drop degenerate CF, score, take the minimum.
    """
    if not poses:
        return None, {}
    scored = [p for p in poses if abs(p["cf"]) > 1e-9]
    pool = scored if scored else list(poses)

    include_singletons = bool(use_shannon_G)
    chosen = [p for p in pool if p["freq"] > 1 or include_singletons]
    if not chosen:
        chosen = pool

    scored_chosen = []
    for p in chosen:
        s = soft_free_energy(p, soft_T) if use_shannon_G else (
            p["cf"] if math.isfinite(p["cf"]) else math.inf)
        if not math.isfinite(s):
            continue
        scored_chosen.append((s, p))
    if not scored_chosen:
        return None, {}
    scored_chosen.sort(key=lambda t: t[0])          # stable; ties keep enum order
    best_s, best_p = scored_chosen[0]
    n_tied = sum(1 for s, _ in scored_chosen if s == best_s)
    return best_p, {"score": best_s, "n_candidates": len(scored_chosen),
                    "n_dropped_degenerate": len(poses) - len(pool),
                    "n_tied_at_best": n_tied}


def elect_mincf(poses):
    """Port of the FLEXAIDDS_POOLED_ELECTION=mincf short-circuit.

    NOTE: it sits ABOVE the degenerate-CF drop and the frequency gate, so it
    sees every finite-CF head. Sentinels (CF >= 1e3) need no exclusion -- an
    argmin cannot pick one while a scored head is negative. Ties resolve to the
    first candidate in enumeration order.
    """
    best = None
    for p in poses:
        if not math.isfinite(p["cf"]):
            continue
        if best is None or p["cf"] < best["cf"]:
            best = p
    if best is None:
        return None, {}
    n_tied = sum(1 for p in poses
                 if math.isfinite(p["cf"]) and p["cf"] == best["cf"])
    return best, {"score": best["cf"], "n_candidates":
                  sum(1 for p in poses if math.isfinite(p["cf"])),
                  "n_tied_at_best": n_tied}


def cell_prefixes(cell_dir, pdb_id, n_restarts=3):
    """Port of the all_prefixes construction (DatasetRunner.cpp:7759/7796)."""
    out = []
    for ri in range(n_restarts):
        d = cell_dir if ri == 0 else os.path.join(cell_dir, "r%d" % ri)
        out.append(os.path.join(d, pdb_id))
    return out
