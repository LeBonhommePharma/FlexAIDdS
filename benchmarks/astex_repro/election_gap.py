#!/usr/bin/env python3
"""Election-gap failure-mode analysis for the FlexAIDdS Astex-85 corrected run.

Answers Grok 4.5's #1 next question on MY OWN consistent pipeline:
  "On election-gap systems, does the ACF score of the true best-cluster pose
   systematically lose to a high-frequency / high-contact decoy?"

For each target it computes, over the pooled poses of all restarts:
  - elected  = RMSD of the rank-0 (elected) pose  -> what the method reports
  - best     = RMSD of the closest-to-native pose in the pool -> ceiling
  - gap      = elected fails (>2 A) BUT best succeeds (<2 A)  -> pure mis-ranking
RMSD is spyrmsd graph-automorphism symmetry-corrected (PoseBusters-grade),
identical instrument to score_reference.py. Pose selection is keyed on PDB
serial >= 90000 (robust to numeric/residue ligand names).

Usage: python3 election_gap.py [RUN_DIR] [CACHE_DIR]
Writes election_gap.csv (per-target) and prints the summary + ranked failure list.
"""
import sys, os, glob, csv
import numpy as np
from rdkit import Chem
from spyrmsd import rmsd as srmsd
from spyrmsd.molecule import Molecule

RUN   = sys.argv[1] if len(sys.argv)>1 else "benchmarks/astex_repro/full"
CACHE = sys.argv[2] if len(sys.argv)>2 else "benchmarks/astex_diverse/astex_diverse"
CUT   = 2.0

def crystal(pdb):
    """(atomic numbers, coords, adjacency) heavy-atom only, from crystal SDF bond block."""
    sdf = f"{CACHE}/{pdb}/{pdb}_ligand.sdf"
    m = Chem.MolFromMolFile(sdf, removeHs=True, sanitize=False)
    if m is None: return None
    mol = Molecule.from_rdkit(m)
    heavy = mol.atomicnums != 1
    return mol.atomicnums[heavy], mol.coordinates[heavy], mol.adjacency_matrix[np.ix_(heavy,heavy)]

def pose_coords(path, nat):
    """Docked-ligand coords in file order; keyed on serial>=90000. Nx3 or None."""
    xyz=[]
    for ln in open(path):
        if not ln.startswith(("HETATM","ATOM")): continue
        try: s=int(ln[6:11])
        except ValueError: continue
        if s>=90000:
            xyz.append([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    return np.array(xyz) if len(xyz)==nat else None

def elected_pose(tdir, pdb):
    """The elected rank-0 pose path: base <pdb>_0.pdb if present, else r1/<pdb>_0.pdb."""
    for c in (f"{tdir}/{pdb}_0.pdb", f"{tdir}/r1/{pdb}_0.pdb"):
        if os.path.exists(c): return c
    g=sorted(glob.glob(f"{tdir}/{pdb}_[0-9].pdb")) or sorted(glob.glob(f"{tdir}/r*/{pdb}_[0-9].pdb"))
    return g[0] if g else None

rows=[]
for tdir in sorted(glob.glob(f"{RUN}/*")):
    pdb=os.path.basename(tdir)
    if not os.path.isdir(tdir) or len(pdb)!=4: continue
    cr=crystal(pdb)
    if cr is None: continue
    anum, ccoord, adj = cr
    poses=sorted(glob.glob(f"{tdir}/{pdb}_[0-9].pdb"))+sorted(glob.glob(f"{tdir}/r*/{pdb}_[0-9].pdb"))
    if not poses: continue
    rmsds=[]
    for p in poses:
        pc=pose_coords(p, len(anum))
        if pc is None: continue
        r=srmsd.symmrmsd(ccoord, pc, anum, anum, adj, adj, minimize=False, center=False)
        rmsds.append((float(r), p))
    if not rmsds: continue
    best_r, best_p = min(rmsds, key=lambda x:x[0])   # ceiling: closest to native in pool
    ep = elected_pose(tdir, pdb)
    er = next((r for r,p in rmsds if p==ep), rmsds[0][0]) if ep else rmsds[0][0]
    rows.append(dict(pdb=pdb, n_poses=len(rmsds),
                     elected_rmsd=round(er,3), best_rmsd=round(best_r,3),
                     elected_ok=int(er<CUT), best_ok=int(best_r<CUT),
                     election_gap=int(er>=CUT and best_r<CUT)))

with open(f"{RUN}/../election_gap.csv","w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

n=len(rows)
el=sum(r["elected_ok"] for r in rows); ce=sum(r["best_ok"] for r in rows)
gap=sum(r["election_gap"] for r in rows)
print(f"targets scored: {n}")
print(f"ELECTED success (<2A): {el}/{n} = {100*el/n:.1f}%")
print(f"CEILING success (best-in-pool <2A): {ce}/{n} = {100*ce/n:.1f}%")
print(f"ELECTION GAP (pool has <2A but elected fails): {gap}/{n} = {100*gap/n:.1f}%")
print(f"  => of the {ce-el} ceiling-minus-elected targets, {gap} are pure mis-ranking\n")
print("=== election-gap targets (elected fails, native-like pose was IN the pool) ===")
for r in sorted([r for r in rows if r["election_gap"]], key=lambda x:x["best_rmsd"]):
    print(f"  {r['pdb']}  elected {r['elected_rmsd']:5.2f} A   best-in-pool {r['best_rmsd']:5.2f} A   (n={r['n_poses']})")
