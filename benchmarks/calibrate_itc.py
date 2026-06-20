#!/usr/bin/env python3
"""
calibrate_itc.py — Calibrate FlexAIDdS thermodynamic terms against ITC data.

Fits two physically-interpretable scaling factors so that FlexAIDdS's
enthalpy and entropy terms become quantitatively comparable to ITC
(Isothermal Titration Calorimetry) measurements:

    dH_pred  = alpha * (T_eff * H_vct_raw)   -> fit alpha against dH_ITC
    TdS_pred = beta  * TdS_vib               -> fit beta  against TdS_ITC
    dG_pred  = dH_pred - TdS_pred            -> validated against dG_ITC (no extra fit)

This is a 2-parameter fit (alpha, beta) that preserves the physical meaning of
each term:  alpha rescales the contact-enthalpy temperature (T_eff), beta
rescales the tENCoM vibrational-entropy magnitude.

Because FlexAIDdS already emits  H_vct = T_eff * H_vct_raw  in its output,
fitting alpha on H_vct is equivalent to rescaling T_eff:

    T_eff_calibrated = alpha * T_eff_current
    tds_vib_scale    = beta

Inputs
------
1. An ITC reference CSV with columns:
       pdb_id, dH_kcal_mol, TdS_kcal_mol, dG_kcal_mol
   (TdS_kcal_mol is the entropy *contribution* T*dS, with dG = dH - TdS.)

2. A directory of FlexAIDdS results. Per pdb_id the script looks, in order, for:
       a. an aggregate  *_results.csv  written with FLEXAIDDS_THERMO_CSV=1
          (columns g_bind,h_vct,h_vct_raw,n_heavy,tds_shannon,tds_vib)
       b. a per-target  <pdb_id>/result.csv
       c. a per-target  <pdb_id>/stdout.log  with a  [THERMO] ...  line
       d. a per-target  <pdb_id>.json / <pdb_id>/result.json

No specific complexes are hard-coded; the script calibrates against whatever
overlap exists between the ITC CSV and the available FlexAIDdS results.

Usage
-----
    python benchmarks/calibrate_itc.py \
        --itc-csv  itc_reference.csv \
        --results  ~/.../results/v89_scorpio_20260617 \
        [--t-eff-current 1.0] \
        [--itc-units kcal|kJ] \
        [--out calibration_report.json]
"""

import argparse
import csv
import glob
import json
import os
import re
import sys

# Reuse the native-format parsers from fetch_itc_data (same directory).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import fetch_itc_data as fid
except ImportError:
    fid = None

KJ_PER_KCAL = 4.184

# ── [THERMO] stdout line (matches DatasetRunner / SCORPIO correlation script) ──
THERMO_RE = re.compile(
    r"\[THERMO\]\s+"
    r"G_bind=(-?\d+(?:\.\d+)?)\s+"
    r"H_vct=(-?\d+(?:\.\d+)?)\s+"
    r"H_vct_raw=(-?\d+(?:\.\d+)?)\s+"
    r"n_heavy=(\d+)\s+"
    r"TdS_shannon=(-?\d+(?:\.\d+)?)\s+"
    r"TdS_vib=(-?\d+(?:\.\d+)?)"
)


# ─────────────────────────────────────────────────────────────────────────────
# ITC reference loading
# ─────────────────────────────────────────────────────────────────────────────
def load_itc_csv(path, units="kcal"):
    """Return {PDBID_UPPER: {dH, TdS, dG}} in kcal/mol."""
    scale = 1.0 if units == "kcal" else 1.0 / KJ_PER_KCAL
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        # tolerate a couple of column-name aliases
        cols = {c.lower().strip(): c for c in reader.fieldnames or []}

        def pick(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_id = pick("pdb_id", "pdbid", "pdb", "receptor_id")
        c_dh = pick("dh_kcal_mol", "dh", "delta_h", "dh_kcalmol")
        c_ds = pick("tds_kcal_mol", "tds", "tdelta_s", "tds_kcalmol")
        c_dg = pick("dg_kcal_mol", "dg", "delta_g", "dg_kcalmol")
        if c_id is None or c_dh is None:
            sys.exit(f"ERROR: {path} must have at least 'pdb_id' and 'dH_kcal_mol' columns; "
                     f"found {reader.fieldnames}")

        def fval(row, col):
            if col is None:
                return None
            v = (row.get(col) or "").strip()
            if v == "" or v.upper() in ("NA", "NAN", "NULL"):
                return None
            try:
                return float(v) * scale
            except ValueError:
                return None

        for row in reader:
            pid = (row.get(c_id) or "").strip().upper()
            if not pid:
                continue
            out[pid] = {
                "dH":  fval(row, c_dh),
                "TdS": fval(row, c_ds),
                "dG":  fval(row, c_dg),
                "source": "csv",
            }
    return out


def _group_rows(rows):
    """Average duplicate measurements per PDB id; track contributing source(s)."""
    by_pdb = {}
    for r in rows:
        pid = (r.get("pdb_id") or "").strip().upper()
        if not pid:
            continue
        by_pdb.setdefault(pid, []).append(r)

    def avg(g, key):
        vals = [g_[key] for g_ in g if g_.get(key) is not None]
        return sum(vals) / len(vals) if vals else None

    out = {}
    for pid, g in by_pdb.items():
        srcs = sorted({g_.get("source", "?") for g_ in g})
        out[pid] = {
            "dH":  avg(g, "dH_kcal_mol"),
            "TdS": avg(g, "TdS_kcal_mol"),
            "dG":  avg(g, "dG_kcal_mol"),
            "source": srcs[0] if len(srcs) == 1 else "mixed",
            "n_meas": len(g),
        }
    return out


def load_itc_reference(path, source="unified", units="kcal"):
    """Load ITC reference -> {PDBID_UPPER: {dH,TdS,dG,source}} in kcal/mol.

    source:
      unified  -- unified CSV (pdb_id,ligand_smiles,dH_kcal_mol,TdS_kcal_mol,dG_kcal_mol,T_K,source,doi)
      scorpio  -- native scorpio_itc_raw.csv (kcal/mol)
      bindingdb-- native BindingDB ITC TSV (kJ/mol)
      csv      -- simple canonical CSV (pdb_id,dH_kcal_mol,TdS_kcal_mol,dG_kcal_mol)
    """
    if source == "csv":
        return load_itc_csv(path, units=units)

    if fid is None:
        sys.exit("ERROR: fetch_itc_data.py must sit beside calibrate_itc.py for "
                 f"--source {source}.")

    if source == "unified":
        rows = list(fid.parse_unified_csv(path))
    elif source == "scorpio":
        rows = list(fid.parse_scorpio(path))
    elif source == "bindingdb":
        rows = list(fid.parse_bindingdb(path))
    else:
        sys.exit(f"ERROR: unknown --source {source}")

    return _group_rows(rows)


# ─────────────────────────────────────────────────────────────────────────────
# FlexAIDdS result extraction
# ─────────────────────────────────────────────────────────────────────────────
def _coerce(d):
    """Normalize one record to the canonical thermo dict."""
    def g(*keys):
        for k in keys:
            if k in d and d[k] not in (None, "", "NA", "NaN"):
                try:
                    return float(d[k])
                except (TypeError, ValueError):
                    pass
        return None

    rec = {
        "G_bind":      g("g_bind", "G_bind", "thermo_G_bind"),
        "H_vct":       g("h_vct", "H_vct", "thermo_H_vct"),
        "H_vct_raw":   g("h_vct_raw", "H_vct_raw", "thermo_H_vct_raw"),
        "n_heavy":     g("n_heavy", "n_genes", "thermo_n_heavy"),
        "TdS_shannon": g("tds_shannon", "TdS_shannon", "thermo_TdS_shannon"),
        "TdS_vib":     g("tds_vib", "TdS_vib", "thermo_TdS_vib"),
    }
    return rec


def _from_thermo_log(path):
    try:
        with open(path) as f:
            for line in f:
                m = THERMO_RE.search(line)
                if m:
                    return {
                        "G_bind":      float(m.group(1)),
                        "H_vct":       float(m.group(2)),
                        "H_vct_raw":   float(m.group(3)),
                        "n_heavy":     int(m.group(4)),
                        "TdS_shannon": float(m.group(5)),
                        "TdS_vib":     float(m.group(6)),
                    }
    except OSError:
        pass
    return None


def _rows_from_csv(path):
    """Yield (pdb_id_upper, thermo_dict) for every data row that has h_vct."""
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                rec = _coerce(row)
                if rec["H_vct"] is None and rec["H_vct_raw"] is None:
                    continue
                pid = (row.get("pdb_id") or row.get("pdbid") or "").strip().upper()
                if pid:
                    yield pid, rec
    except OSError:
        return


def load_flexaidds_results(results_dir):
    """Return {PDBID_UPPER: thermo_dict}. Aggregate CSV wins; per-target fills gaps."""
    thermo = {}

    # (a) aggregate *_results.csv written with FLEXAIDDS_THERMO_CSV=1
    for agg in sorted(glob.glob(os.path.join(results_dir, "*_results.csv"))):
        for pid, rec in _rows_from_csv(agg):
            thermo.setdefault(pid, rec)

    # (b/c/d) per-target subdirectories
    for entry in sorted(os.listdir(results_dir)):
        ep = os.path.join(results_dir, entry)
        if not os.path.isdir(ep):
            continue
        pid = entry.upper()
        if pid in thermo:
            continue
        # per-target result.csv
        rcsv = os.path.join(ep, "result.csv")
        if os.path.isfile(rcsv):
            for _, rec in _rows_from_csv(rcsv):
                thermo[pid] = rec
                break
            if pid in thermo:
                continue
        # stdout.log [THERMO]
        rec = _from_thermo_log(os.path.join(ep, "stdout.log"))
        if rec:
            thermo[pid] = rec
            continue
        # result.json
        for jp in (os.path.join(ep, "result.json"), ep + ".json"):
            if os.path.isfile(jp):
                try:
                    with open(jp) as f:
                        rec = _coerce(json.load(f))
                    if rec["H_vct"] is not None or rec["H_vct_raw"] is not None:
                        thermo[pid] = rec
                        break
                except (OSError, json.JSONDecodeError):
                    pass
    return thermo


# ─────────────────────────────────────────────────────────────────────────────
# Statistics (numpy-free fallbacks so the script runs anywhere)
# ─────────────────────────────────────────────────────────────────────────────
def pearson(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / (sxx ** 0.5 * syy ** 0.5)


def rmse(pred, obs):
    n = len(pred)
    if n == 0:
        return float("nan")
    return (sum((p - o) ** 2 for p, o in zip(pred, obs)) / n) ** 0.5


def fit_scale_through_origin(x, y):
    """alpha minimizing RMSE(alpha*x, y) with no intercept: alpha = <x,y>/<x,x>."""
    sxx = sum(a * a for a in x)
    if sxx <= 0:
        return float("nan")
    return sum(a * b for a, b in zip(x, y)) / sxx


# ─────────────────────────────────────────────────────────────────────────────
# Plotting (per source: dH, TdS, dG predicted vs ITC)
# ─────────────────────────────────────────────────────────────────────────────
def make_plots(rows, alpha, beta, plot_dir):
    """Write pred-vs-ITC scatter PNGs (3 panels) for each source and 'all'."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plot] matplotlib not available — skipping plots", file=sys.stderr)
        return []

    os.makedirs(plot_dir, exist_ok=True)

    def panel(ax, pred, obs, title):
        pts = [(p, o) for p, o in zip(pred, obs) if p is not None and o is not None]
        if len(pts) >= 1:
            xs = [p for p, _ in pts]; ys = [o for _, o in pts]
            ax.scatter(xs, ys, s=28, alpha=0.75, edgecolor="k", linewidth=0.3)
            lo = min(min(xs), min(ys)); hi = max(max(xs), max(ys))
            ax.plot([lo, hi], [lo, hi], "--", color="grey", linewidth=0.8)
            r = pearson(xs, ys); e = rmse(xs, ys)
            ax.set_title(f"{title}\nr={r:+.3f} RMSE={e:.2f} n={len(pts)}", fontsize=9)
        else:
            ax.set_title(f"{title}\n(no data)", fontsize=9)
        ax.set_xlabel("predicted (kcal/mol)", fontsize=8)
        ax.set_ylabel("ITC (kcal/mol)", fontsize=8)

    groups = {"all": rows}
    for r in rows:
        groups.setdefault(r.get("source", "?"), []).append(r)

    written = []
    for src, g in groups.items():
        dHp = [alpha * r["H_vct"] if r["H_vct"] is not None else None for r in g]
        TdSp = [beta * r["TdS_vib"] if r["TdS_vib"] is not None else None for r in g]
        dGp = [(alpha * r["H_vct"] - beta * r["TdS_vib"])
               if (r["H_vct"] is not None and r["TdS_vib"] is not None) else None for r in g]
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
        panel(axes[0], dHp, [r["itc_dH"] for r in g], f"ΔH [{src}]")
        panel(axes[1], TdSp, [r["itc_TdS"] for r in g], f"TΔS [{src}]")
        panel(axes[2], dGp, [r["itc_dG"] for r in g], f"ΔG [{src}]")
        fig.tight_layout()
        out = os.path.join(plot_dir, f"itc_calibration_{src}.png")
        fig.savefig(out, dpi=130); plt.close(fig)
        written.append(out)
    return written


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--itc-csv", required=True,
                    help="ITC reference file (format depends on --source).")
    ap.add_argument("--source", choices=["unified", "scorpio", "bindingdb", "csv"],
                    default="unified",
                    help="ITC reference format: unified table (default), native "
                         "scorpio_itc_raw.csv, native BindingDB ITC TSV, or simple "
                         "canonical CSV (pdb_id,dH_kcal_mol,TdS_kcal_mol,dG_kcal_mol).")
    ap.add_argument("--results", required=True,
                    help="Directory of FlexAIDdS results (aggregate CSV and/or per-target dirs)")
    ap.add_argument("--itc-units", choices=["kcal", "kJ"], default="kcal",
                    help="Units for --source csv only (default kcal). unified/scorpio are "
                         "already kcal/mol; bindingdb is converted from kJ automatically.")
    ap.add_argument("--t-eff-current", type=float, default=1.0,
                    help="Current T_eff used in the docking runs; reported "
                         "T_eff_calibrated = alpha * t_eff_current (default 1.0).")
    ap.add_argument("--plot-dir", default=None,
                    help="If set, write per-source ΔH/TΔS/ΔG correlation PNGs here.")
    ap.add_argument("--out", default=None, help="Write JSON calibration report here.")
    args = ap.parse_args()

    itc = load_itc_reference(args.itc_csv, source=args.source, units=args.itc_units)
    thermo = load_flexaidds_results(args.results)
    print(f"ITC reference complexes : {len(itc)}  (source={args.source})")
    print(f"FlexAIDdS results w/ thermo: {len(thermo)}")

    # ── build matched rows ────────────────────────────────────────────────────
    H, dH = [], []          # enthalpy proxy (H_vct) vs ITC dH
    V, TdS = [], []         # TdS_vib vs ITC TdS
    rows = []
    for pid in sorted(set(itc) & set(thermo)):
        i, t = itc[pid], thermo[pid]
        # enthalpy proxy: prefer H_vct (= T_eff * H_vct_raw); else reconstruct
        h_vct = t["H_vct"]
        if h_vct is None and t["H_vct_raw"] is not None:
            h_vct = args.t_eff_current * t["H_vct_raw"]
        row = {"pdb_id": pid, "H_vct": h_vct, "TdS_vib": t["TdS_vib"],
               "itc_dH": i["dH"], "itc_TdS": i["TdS"], "itc_dG": i["dG"],
               "source": i.get("source", args.source)}
        rows.append(row)
        if h_vct is not None and i["dH"] is not None:
            H.append(h_vct); dH.append(i["dH"])
        if t["TdS_vib"] is not None and i["TdS"] is not None:
            V.append(t["TdS_vib"]); TdS.append(i["TdS"])

    if not rows:
        sys.exit("ERROR: no overlap between ITC CSV pdb_ids and FlexAIDdS results.")

    print(f"Matched complexes        : {len(rows)} "
          f"(dH pairs={len(H)}, TdS pairs={len(V)})\n")

    # ── fit alpha (enthalpy) and beta (entropy) ───────────────────────────────
    alpha = fit_scale_through_origin(H, dH) if len(H) >= 2 else float("nan")
    beta = fit_scale_through_origin(V, TdS) if len(V) >= 2 else float("nan")

    # ── predictions & correlations ────────────────────────────────────────────
    dH_pred = [alpha * h for h in H]
    TdS_pred = [beta * v for v in V]
    r_dH = pearson(dH_pred, dH) if len(H) >= 2 else float("nan")  # == pearson(H,dH)
    r_TdS = pearson(TdS_pred, TdS) if len(V) >= 2 else float("nan")
    rmse_dH = rmse(dH_pred, dH)
    rmse_TdS = rmse(TdS_pred, TdS)

    # dG: predict where all three terms are available, validate against ITC dG
    Gp, Go = [], []
    for row in rows:
        if (row["H_vct"] is not None and row["TdS_vib"] is not None
                and row["itc_dG"] is not None and alpha == alpha and beta == beta):
            gp = alpha * row["H_vct"] - beta * row["TdS_vib"]
            Gp.append(gp); Go.append(row["itc_dG"])
    r_dG = pearson(Gp, Go) if len(Gp) >= 2 else float("nan")
    rmse_dG = rmse(Gp, Go) if Gp else float("nan")

    # ── report ────────────────────────────────────────────────────────────────
    t_eff_cal = alpha * args.t_eff_current
    print("── Fitted parameters ──────────────────────────────────────────")
    print(f"  alpha (enthalpy scale)       = {alpha:.6f}")
    print(f"  beta  (entropy scale)        = {beta:.6f}")
    print(f"  T_eff_current                = {args.t_eff_current:.6f}")
    print(f"  T_eff_calibrated = alpha*T   = {t_eff_cal:.6f}")
    print(f"  tds_vib_scale = beta         = {beta:.6f}")
    print("── Correlations / errors (kcal/mol) ───────────────────────────")
    print(f"  dH   : Pearson r = {r_dH:+.3f}   RMSE = {rmse_dH:.3f}   (n={len(H)})")
    print(f"  TdS  : Pearson r = {r_TdS:+.3f}   RMSE = {rmse_TdS:.3f}   (n={len(V)})")
    print(f"  dG   : Pearson r = {r_dG:+.3f}   RMSE = {rmse_dG:.3f}   (n={len(Gp)})")
    # ── per-source breakdown (validation that fitted alpha/beta generalize) ────
    per_source = {}
    srcs = sorted({r["source"] for r in rows})
    if len(srcs) > 1:
        print("── Per-source correlations (fixed global alpha/beta) ──────────")
    for src in srcs:
        g = [r for r in rows if r["source"] == src]
        sH = [(alpha * r["H_vct"], r["itc_dH"]) for r in g
              if r["H_vct"] is not None and r["itc_dH"] is not None]
        sV = [(beta * r["TdS_vib"], r["itc_TdS"]) for r in g
              if r["TdS_vib"] is not None and r["itc_TdS"] is not None]
        sG = [(alpha * r["H_vct"] - beta * r["TdS_vib"], r["itc_dG"]) for r in g
              if r["H_vct"] is not None and r["TdS_vib"] is not None and r["itc_dG"] is not None]
        m = {
            "n": len(g),
            "dH":  {"pearson_r": pearson([a for a, _ in sH], [b for _, b in sH]),
                    "rmse": rmse([a for a, _ in sH], [b for _, b in sH]), "n": len(sH)},
            "TdS": {"pearson_r": pearson([a for a, _ in sV], [b for _, b in sV]),
                    "rmse": rmse([a for a, _ in sV], [b for _, b in sV]), "n": len(sV)},
            "dG":  {"pearson_r": pearson([a for a, _ in sG], [b for _, b in sG]),
                    "rmse": rmse([a for a, _ in sG], [b for _, b in sG]), "n": len(sG)},
        }
        per_source[src] = m
        if len(srcs) > 1:
            print(f"  [{src:<10}] dH r={m['dH']['pearson_r']:+.3f}(n={m['dH']['n']})  "
                  f"TdS r={m['TdS']['pearson_r']:+.3f}(n={m['TdS']['n']})  "
                  f"dG r={m['dG']['pearson_r']:+.3f}(n={m['dG']['n']})")

    print("── Suggested env-var settings for calibrated runs ─────────────")
    print(f"  export FLEXAIDDS_THERMO_CALIBRATE=1")
    print(f"  export FLEXAIDDS_ALPHA_ENTHALPY={alpha:.6f}")
    print(f"  export FLEXAIDDS_BETA_ENTROPY={beta:.6f}")

    plots = []
    if args.plot_dir:
        plots = make_plots(rows, alpha, beta, args.plot_dir)
        for p in plots:
            print(f"  plot: {p}")

    report = {
        "itc_csv": os.path.abspath(args.itc_csv),
        "results_dir": os.path.abspath(args.results),
        "itc_units": args.itc_units,
        "n_itc": len(itc),
        "n_results": len(thermo),
        "n_matched": len(rows),
        "fit": {
            "alpha_enthalpy": alpha,
            "beta_entropy": beta,
            "t_eff_current": args.t_eff_current,
            "t_eff_calibrated": t_eff_cal,
            "tds_vib_scale": beta,
        },
        "metrics": {
            "dH":  {"pearson_r": r_dH,  "rmse": rmse_dH,  "n": len(H)},
            "TdS": {"pearson_r": r_TdS, "rmse": rmse_TdS, "n": len(V)},
            "dG":  {"pearson_r": r_dG,  "rmse": rmse_dG,  "n": len(Gp)},
        },
        "per_source": per_source,
        "plots": plots,
        "source": args.source,
        "matched": rows,
    }
    if args.out:
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nWrote {args.out}")

    return report


if __name__ == "__main__":
    main()
