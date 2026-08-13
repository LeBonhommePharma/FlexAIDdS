#!/usr/bin/env python3
"""
FlexAIDdS CANONICAL SCORER — the single source of truth for every agent.

WHY THIS EXISTS
  Five agents computing their own RMSD produce five incomparable numbers. This project
  has already been burned three ways:
    * result.csv "rmsd_to_crystal"  = ORDERED RMSD (over-strict; 6.84 A where the correct
      symmetry-corrected value is 0.98 A on 1TZ8)
    * result.csv "rmsd_hungarian"   = element-blocked assignment (OVER-PERMISSIVE; inflated
      the measured pool ceiling from the true 48.8% to 57.8%)
    * result.csv "success"          = "docking ran", NOT "docking succeeded"
  The only correct metric is spyRMSD graph automorphism using the crystal SDF bond block,
  no superposition. That is what this script computes, and nothing else.

USAGE
  python3 score_canonical.py --run  /path/to/<campaign_dir>     # score a campaign/pilot
  python3 score_canonical.py --frozen <frozen_benchmark.csv>    # score the frozen poses
  optional: --denominator 84 (default)  --cache <dir>  --json out.json

OUTPUT
  A canonical block. Quote these numbers verbatim; do not recompute them your own way.
"""
import argparse, collections, csv, glob, json, os, re, sys

def _cache_default():
    env = os.environ.get("FLEXAIDDS_CACHE_V2")
    if env:
        return env
    results = os.environ.get("FLEXAIDDS_RESULTS")
    if results:
        return os.path.join(results, "cache_v2", "astex_diverse")
    return ""

CACHE_DEFAULT = _cache_default()
DENOM_DEFAULT = 84   # 2HR7 excluded: its "ligand" is PEG (CCD P33); no cognate ligand exists.

def die(msg):
    sys.stderr.write(f"[SCORE-FAIL] {msg}\n"); sys.exit(2)

def sdf_block(path):
    """Crystal ligand: heavy-atom coords, elements, and the bond block (for automorphism)."""
    L = open(path, errors="ignore").read().split("\n")
    na, nb = int(L[3][0:3]), int(L[3][3:6])
    xyz, el, keep = [], [], []
    for i in range(na):
        f = L[4+i]
        e = f[31:34].strip()
        if e.upper() == "H":
            keep.append(None); continue
        keep.append(len(xyz))
        xyz.append([float(f[0:10]), float(f[10:20]), float(f[20:30])]); el.append(e)
    bonds = []
    for j in range(nb):
        f = L[4+na+j]
        a, b = int(f[0:3])-1, int(f[3:6])-1
        if keep[a] is not None and keep[b] is not None:
            bonds.append((keep[a], keep[b], int(f[6:9])))
    return xyz, el, bonds

def pose_xyz(path):
    """Docked ligand = the TRAILING HETATM group (pose PDBs contain the full complex)."""
    het = [l for l in open(path, errors="ignore") if l.startswith("HETATM")]
    if not het: return None
    tag = (het[-1][17:20], het[-1][21:22])
    grp = [l for l in het if (l[17:20], l[21:22]) == tag]
    out = []
    for l in grp:
        e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if e == "H": continue
        out.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out or None

def rmsd_fn():
    try:
        import numpy as np
        from spyrmsd import rmsd as sr
        from spyrmsd.molecule import Molecule
    except ImportError:
        die("spyrmsd not installed. This script refuses to fall back to a weaker metric.")
    import numpy as np
    def f(P, Xc, elc, bonds):
        if len(P) != len(Xc): return None
        A = np.zeros((len(Xc), len(Xc)), dtype=int)
        for a, b, _ in bonds: A[a, b] = A[b, a] = 1
        Z = np.array([_z(e) for e in elc])
        try:
            return float(sr.symmrmsd(np.array(P), np.array(Xc), Z, Z, A, A, minimize=False))
        except Exception:
            return None
    return f

_PT = {"H":1,"C":6,"N":7,"O":8,"F":9,"P":15,"S":16,"CL":17,"BR":35,"I":53,
       "NA":11,"MG":12,"K":19,"CA":20,"MN":25,"FE":26,"CO":27,"NI":28,"CU":29,"ZN":30,"SE":34,"B":5,"SI":14}
def _z(e): return _PT.get(e.strip().upper(), 6)

def score_run(run_dir, cache, denom):
    rmsd = rmsd_fn()
    RUN = os.path.join(run_dir, "run") if os.path.isdir(os.path.join(run_dir, "run")) else run_dir
    tdirs = sorted(d for d in os.listdir(RUN) if os.path.isdir(os.path.join(RUN, d)))
    if not tdirs: die(f"no target directories under {RUN}")
    per = {}
    for t in tdirs:
        lig = os.path.join(cache, t, f"{t}_ligand.sdf")
        if not os.path.exists(lig): continue
        Xc, elc, bonds = sdf_block(lig)
        pool = []
        for p in glob.glob(os.path.join(RUN, t, "**", "*.pdb"), recursive=True):
            if p.endswith("_INI.pdb"): continue
            P = pose_xyz(p)
            if not P: continue
            cf = None
            m = re.search(r"REMARK CF=\s*(-?[\d.]+)", open(p, errors="ignore").read())
            if m: cf = float(m.group(1))
            if cf is None or cf > 1e3: continue          # sentinel CF -> not a real pose
            r = rmsd(P, Xc, elc, bonds)
            if r is not None: pool.append((cf, r, p))
        if pool:
            per[t] = {"n_poses": len(pool),
                      "min_cf_rmsd": min(pool, key=lambda x: x[0])[1],
                      "ceiling":     min(x[1] for x in pool)}
    return per

def score_frozen(csv_path, denom):
    byt = collections.defaultdict(list)
    for r in csv.DictReader(open(csv_path)):
        try: cf = float(r["CF_total"])
        except (TypeError, ValueError, KeyError): continue
        if cf > 1e3: continue
        byt[r["target"]].append((cf, float(r["rmsd_spyrmsd"])))
    return {t: {"n_poses": len(v),
                "min_cf_rmsd": min(v, key=lambda x: x[0])[1],
                "ceiling": min(x[1] for x in v)} for t, v in byt.items()}

def report(per, denom, label):
    n_mc = sum(1 for v in per.values() if v["min_cf_rmsd"] < 2.0)
    n_ce = sum(1 for v in per.values() if v["ceiling"] < 2.0)
    gap  = sorted(t for t, v in per.items() if v["ceiling"] < 2.0 <= v["min_cf_rmsd"])
    print("=" * 72)
    print(f"CANONICAL SCORE — {label}")
    print(f"  metric      : spyRMSD graph automorphism, no superposition")
    print(f"  denominator : N={denom} (2HR7 excluded — ligand is PEG, no cognate ligand)")
    print(f"  targets with scorable poses: {len(per)}   (missing count as FAILURES)")
    print("-" * 72)
    print(f"  min-CF election   {n_mc}/{denom} = {100*n_mc/denom:5.1f}%")
    print(f"  pool ceiling      {n_ce}/{denom} = {100*n_ce/denom:5.1f}%")
    print(f"  selection gap     {len(gap)} targets (sub-2A pose exists, min-CF misses it)")
    print("-" * 72)
    print("  REFERENCE (parent campaign astex84_dG_20260809_141245):")
    print("    as-run T=300 election  15/84 = 17.9%")
    print("    min-CF election        26/84 = 31.0%")
    print("    pool ceiling           41/84 = 48.8%")
    print("    PUBLISHED BAR                  45.2%   (top-10 bar 66.7%)")
    print("-" * 72)
    print("  DO NOT SUM GAINS FROM SEPARATE FIXES. min-CF over all poses already contains")
    print("  void recovery and cluster-representative effects. Re-elect once, report once.")
    print("=" * 72)
    return {"n_min_cf": n_mc, "n_ceiling": n_ce, "denominator": denom,
            "pct_min_cf": round(100*n_mc/denom, 1), "pct_ceiling": round(100*n_ce/denom, 1),
            "selection_gap_targets": gap, "per_target": per}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run"); ap.add_argument("--frozen")
    ap.add_argument("--cache", default=CACHE_DEFAULT)
    ap.add_argument("--denominator", type=int, default=DENOM_DEFAULT)
    ap.add_argument("--json")
    a = ap.parse_args()
    if not a.run and not a.frozen: die("give --run <dir> or --frozen <csv>")
    if a.run:
        if not os.path.isdir(a.run): die(f"not a directory: {a.run}")
        if not a.cache:
            die("give --cache <astex_diverse dir> or set FLEXAIDDS_CACHE_V2 / FLEXAIDDS_RESULTS")
        per, label = score_run(a.run, a.cache, a.denominator), os.path.basename(a.run.rstrip("/"))
    else:
        if not os.path.exists(a.frozen): die(f"no such file: {a.frozen}")
        per, label = score_frozen(a.frozen, a.denominator), os.path.basename(a.frozen)
    out = report(per, a.denominator, label)
    if a.json: json.dump(out, open(a.json, "w"), indent=1)
