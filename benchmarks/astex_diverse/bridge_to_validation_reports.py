#!/usr/bin/env python3
"""
Bridge: convert the benchmark_datasets aggregate `astex_diverse_results.csv`
into per-complex `<ENTRY>_validation_report.json` files in the schema that the
skill's analyze_entropy_collapse.py consumes (it only discovers *_validation_report.json).

Faithful mapping:
  predicted_dG         -> validation.delta_g_pred_kcal
  search_entropy_proxy -> validation.entropy_shannon_collapse
                          (falls back to legacy shannon_entropy; NULL when 0 == failure sentinel)
  predicted_dH   -> validation.delta_h_kcal              (NULL when dH==0 and TdS==0 == decomposition fell back)
  predicted_TdS  -> validation.t_delta_s_kcal            ( "" )
  rmsd_to_crystal-> validation.pose_rmsd  + an {ENTRY}.log "RMSD: x" line for the analyzer's enrichment path
Failures (shannon==0 / rmsd>=999) are written with null entropy so they are
recorded but excluded from collapse statistics and the weak-collapse outlier set.
"""
import csv, json, os, sys

CSV = os.path.expanduser(sys.argv[1])
OUTROOT = os.path.dirname(CSV)  # the astex/ dir; each complex has a subdir here

n_real = n_fail = 0
with open(CSV) as f:
    rows = list(csv.DictReader(f))

for r in rows:
    entry = r["pdb_id"].strip()
    if not entry:
        continue
    def fnum(k):
        try:
            return float(r[k])
        except (TypeError, ValueError):
            return None
    dg  = fnum("predicted_dG")
    sh  = fnum("search_entropy_proxy")
    if sh is None:
        sh = fnum("shannon_entropy")
    dh  = fnum("predicted_dH")
    tds = fnum("predicted_TdS")
    rmsd = fnum("rmsd_to_crystal")
    success = (r.get("success", "0").strip() == "1")

    failure = (sh is None or sh == 0.0) or (rmsd is not None and rmsd >= 999.0)
    # Honest null mapping
    shannon_val = None if (sh is None or sh == 0.0) else sh
    if (dh == 0.0 and tds == 0.0):
        dh_val = tds_val = None       # decomposition fell back to zero -> missing
    else:
        dh_val, tds_val = dh, tds

    d = os.path.join(OUTROOT, entry)
    os.makedirs(d, exist_ok=True)
    report = {
        "dataset": "astex_diverse",
        "entry_id": entry,
        "flexaidds_binary": "build/benchmark_datasets (v7rerun ga-2000)",
        "validation": {
            "status": "thermo_parsed" if not failure else "docking_failed_or_incomplete",
            "delta_g_pred_kcal": dg if not failure else None,
            "delta_g_exp_kcal": None,
            "delta_delta_g": None,
            "entropy_shannon_collapse": shannon_val,
            "tencom_vibrational": None,
            "delta_h_kcal": dh_val if not failure else None,
            "t_delta_s_kcal": tds_val if not failure else None,
            "pose_rmsd": rmsd,
            "num_poses": r.get("num_poses"),
            "wall_time_s": r.get("wall_time_s"),
            "success": success,
            "notes": "bridged from astex_diverse_results.csv (ga-2000, 1000x2000, 11 workers)",
        },
        "provenance": {
            "source_csv": CSV,
            "ga_population": 1000, "ga_generations": 2000,
            "skill_version": "flexaidds-dataset-runner bridge",
        },
    }
    with open(os.path.join(d, f"{entry}_validation_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    # enrichment hook so the analyzer can pick up pose RMSD via its {entry}.log path
    if rmsd is not None:
        with open(os.path.join(d, f"{entry}.log"), "w") as fh:
            fh.write(f"RMSD: {rmsd}\n")
            if shannon_val is not None:
                fh.write(f"SHANNON_ENTROPY_COLLAPSE: {shannon_val}\n")

    if failure:
        n_fail += 1
    else:
        n_real += 1

print(f"wrote {n_real+n_fail} validation reports: {n_real} with entropy data, {n_fail} failures (null entropy)")
