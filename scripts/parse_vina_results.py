#!/usr/bin/env python3
# =============================================================================
# parse_vina_results.py — score AutoDock Vina Astex output & build comparison
# =============================================================================
# Parallel to parse_rdock_results.py. Reads the per-complex Vina PDBQT output
# from run_vina_astex.sh, takes the TOP-1 pose (MODEL 1 = best affinity),
# computes the element-matched Hungarian RMSD vs the crystal ligand, applies the
# same sub-2 A success criterion, and writes:
#
#   vina_astex_results.csv      Vina-only results (engine-result schema)
#   vina_vs_flexaidds.csv       side-by-side vs FlexAIDdS
#
# Shares the SDF reader + Hungarian RMSD with parse_rdock_results.py so the
# RMSD metric is byte-for-byte identical across all three engines.
# =============================================================================
import argparse
import csv
import os
import sys

# Reuse the proven SDF parser + Hungarian RMSD from the rDock parser.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_rdock_results import _read_sdf_records, hungarian_rmsd  # noqa: E402

DEFAULT_VINA_DIR = os.path.expanduser(
    "~/flexaidds_benchmark_results/vina_astex")
DEFAULT_FLEXAIDDS_CSV = os.path.expanduser(
    "~/flexaidds_benchmark_results/astex_diverse_results.csv")
DEFAULT_ASTEX_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "benchmarks", "astex_diverse", "astex_diverse")
SUCCESS_THRESHOLD = 2.0  # Angstrom

# AutoDock atom-type (PDBQT cols 78-79) -> element, for Hungarian matching.
AD_TYPE_TO_ELEM = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
    "H": "H", "HD": "H", "HS": "H", "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "MG": "Mg", "ZN": "Zn", "CA": "Ca", "FE": "Fe", "MN": "Mn",
}


def _elem_from_pdbqt(line):
    """Derive element from a PDBQT ATOM/HETATM line via its AutoDock type."""
    ad = line[77:79].strip().upper() if len(line) >= 79 else ""
    if ad in AD_TYPE_TO_ELEM:
        return AD_TYPE_TO_ELEM[ad]
    # Fallback: PDB element column / atom name first letters
    el = line[76:78].strip() if len(line) >= 78 else ""
    if el:
        return el.capitalize()
    name = line[12:16].strip()
    return (name[0:1] or "C").capitalize()


def top1_pose_pdbqt(out_pdbqt):
    """Return (atoms, affinity, num_models) for Vina MODEL 1 (best pose).

    atoms: list of (element, x, y, z) heavy atoms only.
    affinity: kcal/mol from 'REMARK VINA RESULT' of MODEL 1.
    """
    atoms = []
    affinity = None
    num_models = 0
    in_first = False
    captured = False
    with open(out_pdbqt, errors="replace") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                num_models += 1
                in_first = (num_models == 1)
                continue
            if line.startswith("ENDMDL"):
                if in_first:
                    captured = True
                in_first = False
                continue
            if in_first and line.startswith("REMARK VINA RESULT"):
                try:
                    affinity = float(line.split()[3])
                except (IndexError, ValueError):
                    pass
            elif in_first and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                elem = _elem_from_pdbqt(line)
                if elem.upper() != "H":
                    atoms.append((elem, x, y, z))
    # Single-model files may have no MODEL/ENDMDL records.
    if num_models == 0 and atoms == [] and not captured:
        with open(out_pdbqt, errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT") and affinity is None:
                    try:
                        affinity = float(line.split()[3])
                    except (IndexError, ValueError):
                        pass
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    except ValueError:
                        continue
                    elem = _elem_from_pdbqt(line)
                    if elem.upper() != "H":
                        atoms.append((elem, x, y, z))
        num_models = 1 if atoms else 0
    return (atoms or None), affinity, num_models


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vina-dir", default=DEFAULT_VINA_DIR,
                    help="Dir with <CODE>/out.pdbqt (default: %(default)s)")
    ap.add_argument("--astex-dir", default=DEFAULT_ASTEX_DIR)
    ap.add_argument("--flexaidds-csv", default=DEFAULT_FLEXAIDDS_CSV)
    ap.add_argument("--out-dir", default=None,
                    help="Where to write CSVs (default: --vina-dir)")
    ap.add_argument("--threshold", type=float, default=SUCCESS_THRESHOLD)
    args = ap.parse_args()

    out_dir = args.out_dir or args.vina_dir
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.isdir(args.vina_dir):
        print(f"ERROR: vina dir not found: {args.vina_dir}", file=sys.stderr)
        print("Run scripts/run_vina_astex.sh first.", file=sys.stderr)
        sys.exit(1)

    codes = sorted(d for d in os.listdir(args.astex_dir)
                   if os.path.isdir(os.path.join(args.astex_dir, d)))

    rows = []
    for code in codes:
        out_pdbqt = os.path.join(args.vina_dir, code, "out.pdbqt")
        crystal = os.path.join(args.astex_dir, code, f"{code}_ligand.sdf")

        affinity = rmsd = None
        num_poses = 0
        status = "ok"

        if not os.path.isfile(crystal):
            status = "no_crystal_ref"
        elif not os.path.isfile(out_pdbqt):
            status = "no_vina_output"
        else:
            try:
                ref_atoms, _ = next(_read_sdf_records(crystal))
            except StopIteration:
                ref_atoms = None
            pose_atoms, affinity, num_poses = top1_pose_pdbqt(out_pdbqt)
            if pose_atoms is None:
                status = "no_pose"
            elif not ref_atoms:
                status = "bad_crystal_ref"
            else:
                rmsd = hungarian_rmsd(ref_atoms, pose_atoms)
                if rmsd is None:
                    status = "atom_mismatch"

        success = int(rmsd is not None and rmsd < args.threshold)
        rows.append({
            "pdb_id": code,
            "vina_affinity": "" if affinity is None else f"{affinity:.4f}",
            "rmsd_to_crystal": "" if rmsd is None else f"{rmsd:.4f}",
            "num_poses": num_poses,
            "success": success,
            "status": status,
        })

    # ---- Vina-only CSV -----------------------------------------------------
    vina_csv = os.path.join(out_dir, "vina_astex_results.csv")
    fields = ["pdb_id", "vina_affinity", "rmsd_to_crystal",
              "num_poses", "success", "status"]
    with open(vina_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n_eval = sum(1 for r in rows if r["rmsd_to_crystal"] != "")
    n_succ = sum(r["success"] for r in rows)
    print(f"Vina: {n_succ}/{n_eval} sub-{args.threshold}A "
          f"(top-1/MODEL 1, Hungarian RMSD)  -> {vina_csv}")

    # ---- Side-by-side vs FlexAIDdS ----------------------------------------
    fa = {}
    if os.path.isfile(args.flexaidds_csv):
        with open(args.flexaidds_csv, newline="") as fh:
            for r in csv.DictReader(fh):
                fa[r["pdb_id"]] = r
    else:
        print(f"WARN: FlexAIDdS CSV not found ({args.flexaidds_csv}).",
              file=sys.stderr)

    cmp_csv = os.path.join(out_dir, "vina_vs_flexaidds.csv")
    cmp_fields = ["pdb_id", "flexaidds_rmsd", "flexaidds_success",
                  "vina_rmsd", "vina_success", "winner", "vina_status"]
    fa_succ_total = 0
    with open(cmp_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cmp_fields)
        w.writeheader()
        for r in rows:
            far = fa.get(r["pdb_id"], {})
            fa_rmsd = far.get("rmsd_to_crystal", "")
            try:
                fa_succ = int(float(fa_rmsd) < args.threshold) if fa_rmsd != "" else 0
            except ValueError:
                fa_succ = 0
            fa_succ_total += fa_succ
            vn_succ = r["success"]
            if vn_succ and not fa_succ:
                winner = "Vina"
            elif fa_succ and not vn_succ:
                winner = "FlexAIDdS"
            elif fa_succ and vn_succ:
                winner = "both"
            else:
                winner = "neither"
            w.writerow({
                "pdb_id": r["pdb_id"],
                "flexaidds_rmsd": fa_rmsd,
                "flexaidds_success": fa_succ,
                "vina_rmsd": r["rmsd_to_crystal"],
                "vina_success": vn_succ,
                "winner": winner,
                "vina_status": r["status"],
            })

    print(f"FlexAIDdS: {fa_succ_total}/{len(rows)} sub-{args.threshold}A "
          f"(recomputed from RMSD)")
    print(f"Side-by-side -> {cmp_csv}")


if __name__ == "__main__":
    main()
