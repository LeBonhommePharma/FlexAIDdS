"""
Offline test suite for the ITC calibration methodology.

Exercises benchmarks/fetch_itc_data.py and benchmarks/calibrate_itc.py end to
end using only synthetic fixtures — no network, no docking runs, no numpy or
matplotlib. Safe for CI.

Run:
    pytest benchmarks/tests/test_itc_calibration.py
    # or, without pytest:
    python benchmarks/tests/test_itc_calibration.py
"""

import csv
import json
import math
import os
import sys

# Make the benchmarks/ scripts importable regardless of pytest rootdir.
_HERE = os.path.dirname(os.path.abspath(__file__))
_BENCH = os.path.dirname(_HERE)
sys.path.insert(0, _BENCH)

import fetch_itc_data as fid          # noqa: E402
import calibrate_itc as cal           # noqa: E402


# ─────────────────────────────────────────────────────────────────────────────
# tiny fixture writers
# ─────────────────────────────────────────────────────────────────────────────
def _write(path, text):
    with open(path, "w", newline="") as f:
        f.write(text)
    return path


def make_scorpio_csv(path):
    return _write(path,
        "pdb_id,ligand_name,protein,dG_kcalmol,dH_kcalmol,TdS_kcalmol,temperature_K,source_url\n"
        "1abc,LigA,ProtA,-8.0,-10.0,-2.0,298.15,http://x/1\n"
        "1abc,LigA,ProtA,-8.2,-11.0,-2.8,298.15,http://x/2\n"   # dup -> averaged
        "2xyz,LigB,ProtB,-6.0,-3.0,3.0,298.15,http://x/3\n")


def make_bindingdb_tsv(path):
    # -T Delta_S0 is -(T*dS); dG = dH - TdS must hold after sign flip.
    cols = ["Ligand SMILES", "Delta_G0 (kJ/mol)", "Delta_H0 (kJ/mol)",
            "-T Delta_S0 (kJ/mol)", "pH", "Temp (C)", "Article DOI",
            "PDB ID(s) for Ligand-Target Complex"]
    rows = [
        # self-consistent: dH=-20, TdS=+10, dG=dH-TdS=-30 kcal (=-125.52 kJ)
        ["CCO", "-125.52", "-83.68", "-41.84", "7", "25.00 C", "10.1/x", "3ABC,9ZZZ"],
        ["c1ccccc1", "", "-20.92", "-8.368", "7.4", "25", "10.1/y", "4DEF"],  # dG derived
    ]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        w.writerows(rows)
    return path


def make_flexaidds_results(root):
    """A results dir with one aggregate CSV (thermo cols) + one per-target THERMO log."""
    os.makedirs(root, exist_ok=True)
    # aggregate CSV as written by DatasetRunner with FLEXAIDDS_THERMO_CSV=1
    _write(os.path.join(root, "synth_results.csv"),
        "pdb_id,best_score,rmsd_to_crystal,g_bind,h_vct,h_vct_raw,n_heavy,tds_shannon,tds_vib\n"
        "1ABC,-5,1.2,-9.0,20.0,2.0,20,1.0,4.0\n"
        "2XYZ,-4,2.0,-6.0,6.0,0.6,15,0.5,-6.0\n")
    # per-target stdout.log with a [THERMO] line (fallback extraction path)
    d = os.path.join(root, "3ABC")
    os.makedirs(d, exist_ok=True)
    _write(os.path.join(d, "stdout.log"),
        "noise\n[THERMO] G_bind=-12.0 H_vct=30.0 H_vct_raw=3.0 n_heavy=25 "
        "TdS_shannon=2.0 TdS_vib=8.0 D_vib=0.1 compensation=0.5\nmore noise\n")
    return root


# ─────────────────────────────────────────────────────────────────────────────
# fetch_itc_data unit tests
# ─────────────────────────────────────────────────────────────────────────────
def test_complete_dg_fills_each_missing_term():
    assert fid._complete_dg(-10.0, -2.0, None) == (-10.0, -2.0, -8.0)   # dG=dH-TdS
    assert fid._complete_dg(-10.0, None, -8.0) == (-10.0, -2.0, -8.0)   # TdS=dH-dG
    assert fid._complete_dg(None, -2.0, -8.0) == (-10.0, -2.0, -8.0)    # dH=dG+TdS


def test_first_pdb_picks_valid_4char_code():
    assert fid._first_pdb("3ABC,9ZZZ") == "3ABC"
    assert fid._first_pdb("") == ""
    assert fid._first_pdb("1abc") == "1ABC"


def test_scorpio_parser(tmp_path):
    rows = list(fid.parse_scorpio(make_scorpio_csv(str(tmp_path / "s.csv"))))
    assert len(rows) == 3
    r = rows[0]
    assert r["pdb_id"] == "1ABC" and r["source"] == "scorpio"
    assert abs(r["dH_kcal_mol"] - (-10.0)) < 1e-9


def test_bindingdb_parser_units_and_sign(tmp_path):
    rows = list(fid.parse_bindingdb(make_bindingdb_tsv(str(tmp_path / "b.tsv"))))
    assert len(rows) == 2
    r0 = rows[0]
    # -83.68 kJ/mol -> -20 kcal/mol
    assert abs(r0["dH_kcal_mol"] - (-20.0)) < 1e-6
    # -T*dS = -41.84 kJ -> TdS = +10 kcal/mol
    assert abs(r0["TdS_kcal_mol"] - (10.0)) < 1e-6
    # dG = dH - TdS = -20 - 10 = -30
    assert abs(r0["dG_kcal_mol"] - (-30.0)) < 1e-6
    assert r0["pdb_id"] == "3ABC"
    # row 1 had no dG -> derived from dH,TdS
    r1 = rows[1]
    assert abs(r1["dG_kcal_mol"] - (r1["dH_kcal_mol"] - r1["TdS_kcal_mol"])) < 1e-9


def test_chembl_parser_conversion_and_filter():
    page = {"activities": [
        {"standard_type": "dH", "standard_value": "-41.84",
         "standard_units": "kJ.mol-1", "canonical_smiles": "CCO"},
        {"standard_type": "dH", "standard_value": "-5.0",
         "standard_units": "kcal.mol-1"},
        {"standard_type": "Kd", "standard_value": "100", "standard_units": "nM"},
    ]}
    rows = list(fid.parse_chembl_json(page))
    assert len(rows) == 2                                   # Kd filtered out
    assert abs(rows[0]["dH_kcal_mol"] - (-10.0)) < 1e-6     # kJ -> kcal
    assert abs(rows[1]["dH_kcal_mol"] - (-5.0)) < 1e-6      # kcal passthrough
    assert all(r["source"] == "chembl" for r in rows)


def test_pdbbind_index_dg_only(tmp_path):
    idx = _write(str(tmp_path / "INDEX_general_PL_data.2020"),
        "# comment line\n"
        "1abc  2.00  2015  6.00  Kd=1uM  // comment\n")
    rows = list(fid.parse_pdbbind_index(idx))
    assert len(rows) == 1
    r = rows[0]
    assert r["pdb_id"] == "1ABC"
    assert r["dH_kcal_mol"] is None and r["TdS_kcal_mol"] is None
    # dG = -RT ln10 * pK ; pK=6 at 298 K -> ~ -8.18 kcal/mol
    assert -8.4 < r["dG_kcal_mol"] < -7.9


def test_dedup_report_flags_variance(tmp_path):
    rows = list(fid.parse_scorpio(make_scorpio_csv(str(tmp_path / "s.csv"))))
    conflicts = fid.dedup_report(rows)
    ids = {c["pdb_id"] for c in conflicts}
    assert "1ABC" in ids                    # measured twice
    c = next(c for c in conflicts if c["pdb_id"] == "1ABC")
    assert c["n"] == 2 and c["sd_dH"] > 0


# ─────────────────────────────────────────────────────────────────────────────
# calibrate_itc unit tests
# ─────────────────────────────────────────────────────────────────────────────
def test_pearson_and_rmse_known_values():
    assert abs(cal.pearson([1, 2, 3], [2, 4, 6]) - 1.0) < 1e-9
    assert abs(cal.pearson([1, 2, 3], [6, 4, 2]) + 1.0) < 1e-9
    assert abs(cal.rmse([1, 2, 3], [1, 2, 3])) < 1e-9
    assert abs(cal.rmse([0, 0], [3, 4]) - math.sqrt(12.5)) < 1e-9


def test_fit_scale_recovers_known_alpha():
    x = [1.0, 2.0, 3.0, 4.0]
    y = [2.5 * v for v in x]
    assert abs(cal.fit_scale_through_origin(x, y) - 2.5) < 1e-9


def test_group_rows_averages_duplicates():
    rows = [
        {"pdb_id": "1abc", "dH_kcal_mol": -10.0, "TdS_kcal_mol": -2.0,
         "dG_kcal_mol": -8.0, "source": "scorpio"},
        {"pdb_id": "1ABC", "dH_kcal_mol": -12.0, "TdS_kcal_mol": -2.0,
         "dG_kcal_mol": -10.0, "source": "bindingdb"},
    ]
    g = cal._group_rows(rows)
    assert set(g) == {"1ABC"}
    assert abs(g["1ABC"]["dH"] - (-11.0)) < 1e-9
    assert g["1ABC"]["source"] == "mixed" and g["1ABC"]["n_meas"] == 2


def test_thermo_log_extraction(tmp_path):
    log = _write(str(tmp_path / "stdout.log"),
        "[THERMO] G_bind=-12.0 H_vct=30.0 H_vct_raw=3.0 n_heavy=25 "
        "TdS_shannon=2.0 TdS_vib=8.0\n")
    rec = cal._from_thermo_log(log)
    assert rec["H_vct"] == 30.0 and rec["TdS_vib"] == 8.0 and rec["n_heavy"] == 25


def test_load_flexaidds_results_aggregate_and_log(tmp_path):
    root = make_flexaidds_results(str(tmp_path / "res"))
    thermo = cal.load_flexaidds_results(root)
    assert {"1ABC", "2XYZ", "3ABC"} <= set(thermo)
    assert thermo["1ABC"]["H_vct"] == 20.0            # from aggregate CSV
    assert thermo["3ABC"]["TdS_vib"] == 8.0           # from stdout.log


def test_end_to_end_calibration_recovers_scaling(tmp_path):
    """Construct ITC data that is exactly alpha*H_vct / beta*TdS_vib so the fit
    must recover alpha and beta, and dG must correlate."""
    alpha_true, beta_true = 0.5, 1.5
    thermo = {
        "1ABC": {"H_vct": 20.0, "TdS_vib": 4.0},
        "2XYZ": {"H_vct": 6.0,  "TdS_vib": -6.0},
        "3ABC": {"H_vct": 30.0, "TdS_vib": 8.0},
    }
    # write a unified ITC csv with dH = alpha*H_vct, TdS = beta*TdS_vib
    itc_csv = str(tmp_path / "itc.csv")
    with open(itc_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb_id", "ligand_smiles", "dH_kcal_mol", "TdS_kcal_mol",
                    "dG_kcal_mol", "T_K", "source", "doi"])
        for pid, t in thermo.items():
            dH = alpha_true * t["H_vct"]
            TdS = beta_true * t["TdS_vib"]
            w.writerow([pid, "", dH, TdS, dH - TdS, 298.15, "synthetic", ""])

    # write a results dir whose thermo matches `thermo`
    root = str(tmp_path / "res")
    os.makedirs(root, exist_ok=True)
    with open(os.path.join(root, "synth_results.csv"), "w", newline="") as f:
        f.write("pdb_id,g_bind,h_vct,h_vct_raw,n_heavy,tds_shannon,tds_vib\n")
        for pid, t in thermo.items():
            f.write(f"{pid},0,{t['H_vct']},0,10,0,{t['TdS_vib']}\n")

    itc = cal.load_itc_reference(itc_csv, source="unified")
    res = cal.load_flexaidds_results(root)
    H, dH, V, TdS = [], [], [], []
    for pid in set(itc) & set(res):
        H.append(res[pid]["H_vct"]); dH.append(itc[pid]["dH"])
        V.append(res[pid]["TdS_vib"]); TdS.append(itc[pid]["TdS"])
    alpha = cal.fit_scale_through_origin(H, dH)
    beta = cal.fit_scale_through_origin(V, TdS)
    assert abs(alpha - alpha_true) < 1e-6
    assert abs(beta - beta_true) < 1e-6
    # dG prediction must correlate perfectly with synthetic dG
    Gp = [alpha * h - beta * v for h, v in zip(H, V)]
    Go = [d - t for d, t in zip(dH, TdS)]
    assert abs(cal.pearson(Gp, Go) - 1.0) < 1e-6


def test_plots_optional_when_matplotlib_absent(tmp_path):
    """make_plots must degrade gracefully (return []) if matplotlib is missing,
    and must produce files if it is present."""
    rows = [{"H_vct": 20.0, "TdS_vib": 4.0, "itc_dH": 10.0, "itc_TdS": 6.0,
             "itc_dG": 4.0, "source": "synthetic"}]
    out = cal.make_plots(rows, 0.5, 1.5, str(tmp_path / "plots"))
    try:
        import matplotlib  # noqa: F401
        assert out and all(os.path.isfile(p) for p in out)
    except ImportError:
        assert out == []


# ─────────────────────────────────────────────────────────────────────────────
# no-pytest fallback runner
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import tempfile
    import types

    class _TmpPath:
        """Minimal stand-in for pytest's tmp_path (supports / and str())."""
        def __init__(self, base):
            self._base = base
        def __truediv__(self, other):
            return _TmpPath(os.path.join(self._base, str(other)))
        def __str__(self):
            return self._base

    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and isinstance(v, types.FunctionType)]
    passed = failed = 0
    for t in tests:
        needs_tmp = t.__code__.co_argcount == 1
        with tempfile.TemporaryDirectory() as d:
            try:
                t(_TmpPath(d)) if needs_tmp else t()
                print(f"  PASS  {t.__name__}")
                passed += 1
            except Exception as e:  # noqa: BLE001
                print(f"  FAIL  {t.__name__}: {type(e).__name__}: {e}")
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)
