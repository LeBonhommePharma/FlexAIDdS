#!/usr/bin/env python3
"""
diag_vh_1jd0.py — Virtual-H geometry diagnostic for a single target.

Usage:
    python3 scripts/diag_vh_1jd0.py [PDB_ID]
    python3 scripts/diag_vh_1jd0.py 1JD0

Runs FlexAIDdS on the given target (default 1JD0) with FLEXAIDS_VH_DEBUG=1
for ONE generation (n_gen=1, pop=1) to trigger the vH assignment dump and
native-pose H-bond breakdown without running a full docking.

Outputs:
  - [vHdbg] lines: per-atom vH assignments (N.am bond_cnt=0 → VHG_NONE bug)
  - [hbdbg] lines: per-pair H-bond energies for the native pose
  - Summary table: N.am counts by VHG_NONE vs VHG_AMIDE
"""

import os, sys, subprocess, re, shutil, tempfile, pathlib

# ── Config ────────────────────────────────────────────────────────────────────
BINARY   = "/tmp/FlexAIDdS_vH"           # current stamped binary
REPO     = pathlib.Path(__file__).parent.parent
DATASET  = REPO / "benchmarks" / "astex_diverse" / "astex_diverse"
PDB_ID   = sys.argv[1].upper() if len(sys.argv)>1 else "1JD0"

# Locate receptor and ligand from dataset
def find_inputs(pdb_id):
    """Return (receptor_pdb, ligand_sdf) for a given Astex Diverse target."""
    base = DATASET / pdb_id
    rec  = base / f"{pdb_id}_apo.pdb"
    lig  = base / f"{pdb_id}_ligand.sdf"
    if not rec.exists():
        rec = next(base.glob("*_apo.pdb"), None) or next(base.glob("*.pdb"), None)
    if not lig.exists():
        lig = next(base.glob("*ligand*.sdf"), None)
    return rec, lig

def run_debug(pdb_id):
    rec, lig = find_inputs(pdb_id)
    if not rec or not lig or not rec.exists() or not lig.exists():
        print(f"[diag] ERROR: inputs not found for {pdb_id}")
        print(f"  receptor: {rec}")
        print(f"  ligand:   {lig}")
        sys.exit(1)

    print(f"[diag] {pdb_id}: receptor={rec.name}  ligand={lig.name}")
    print(f"[diag] binary: {BINARY}")
    print(f"[diag] Running with FLEXAIDS_VH_DEBUG=1, n_gen=1 …")

    env = os.environ.copy()
    env["FLEXAIDS_VH_DEBUG"] = "1"

    with tempfile.TemporaryDirectory() as td:
        cmd = [
            BINARY,
            "-rec",  str(rec),
            "-lig",  str(lig),
            "-ref",  str(lig),        # use crystal as reference for RMSD
            "-n_gen", "1",
            "-pop",   "1",
            "-n_restart", "1",
            "-out",  os.path.join(td, "diag_out"),
        ]
        result = subprocess.run(
            cmd, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=120
        )

    lines = result.stdout.splitlines()

    # ── Parse vHdbg assignment lines ─────────────────────────────────────────
    vhdbg = [l for l in lines if l.startswith("[vHdbg]")]
    hbdbg = [l for l in lines if l.startswith("[hbdbg]")]

    print(f"\n{'='*70}")
    print(f"  vH ASSIGNMENT DUMP  ({len(vhdbg)} donor atoms)")
    print(f"{'='*70}")
    for l in vhdbg:
        print(l)

    print(f"\n{'='*70}")
    print(f"  H-BOND ENERGY PAIRS  ({len(hbdbg)} pairs scored non-zero)")
    print(f"{'='*70}")
    for l in hbdbg:
        print(l)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  KEY METRICS")
    print(f"{'='*70}")

    n_am_amide = sum(1 for l in vhdbg if "VHG_AMIDE" in l and "N.am" not in l.split("VHG_AMIDE")[0])
    # Use the summary lines instead
    for l in vhdbg:
        if "N.am with VHG_AMIDE" in l or "N.am with VHG_NONE" in l \
           or "Donors with" in l:
            print(l)

    # H-bond total
    total_hb = 0.0
    for l in hbdbg:
        m = re.search(r'E=(-?[\d.]+)', l)
        if m: total_hb += float(m.group(1))
    print(f"\n  Total H-bond energy (native pose, n_gen=1): {total_hb:.4f}")

    # Fallback count
    n_fallback = sum(1 for l in hbdbg if "fallback" in l)
    n_explicit = sum(1 for l in hbdbg if "explicit-H" in l)
    n_virtual  = sum(1 for l in hbdbg if "virtual-H" in l)
    print(f"  Pairs via explicit-H : {n_explicit}")
    print(f"  Pairs via virtual-H  : {n_virtual}")
    print(f"  Pairs via 0.3 fallback: {n_fallback}  (<-- if large, VHG_NONE bug confirmed)")

    if n_fallback > 0 and n_virtual == 0:
        print(f"\n  !! CONFIRMED BUG: all H-bond angle terms from 0.3 fallback.")
        print(f"     receptor N.am atoms have bond[]=0 (PDB, no topology) → VHG_NONE.")
        print(f"     Fix: populate receptor bond connectivity by distance after PDB load.")
    elif n_virtual > 0:
        print(f"\n  vH is firing for {n_virtual} pairs. Bug is elsewhere (geometry wrong?).")

    return result.returncode

if __name__ == "__main__":
    if not os.path.exists(BINARY):
        print(f"[diag] ERROR: binary not found at {BINARY}")
        print(f"  Run: cp build_lto/FlexAIDdS /tmp/FlexAIDdS_vH")
        sys.exit(1)
    rc = run_debug(PDB_ID)
    sys.exit(rc)
