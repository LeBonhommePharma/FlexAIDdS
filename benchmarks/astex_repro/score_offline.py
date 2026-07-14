#!/usr/bin/env python3
"""Offline poster-metric scorer for the Astex-85 GetCleft cognate reproduction.
Metric (user-confirmed, pooled reading): for each target, pool ALL poses across
all restarts; success = any pose has symmetry-corrected (Hungarian) RMSD < 2.0 Å
to the crystal ligand. Disclosed deviation from strict poster wording
('10 best results'), chosen as the more forgiving reading."""
import numpy as np, glob, re, os, sys, csv
from scipy.optimize import linear_sum_assignment

ROOT=os.path.dirname(os.path.abspath(__file__))
CACHE=f"{ROOT}/../astex_diverse/astex_diverse"
FULL=f"{ROOT}/full"

def norm_el(s):
    """Canonical element token: strip digits/spaces, Title-case so SDF 'CL'
    and pose 'Cl' both -> 'Cl' (halogens/two-letter elements). Prevents the
    element-blocked Hungarian from treating CL!=Cl as a count mismatch."""
    x=re.sub(r'[^A-Za-z]','',s)
    return x[:1].upper()+x[1:].lower() if x else x

def crystal(pdb):
    """Return (xyz, elements, ligname). ligname = SDF title line (the docked
    ligand's residue name in FlexAID pose output, e.g. RQ3/SOX/HUP)."""
    p=f"{CACHE}/{pdb}/{pdb}_ligand.sdf"
    if not os.path.exists(p): return None,None,None
    L=open(p).read().splitlines()
    ligname=L[0].strip()
    n=int(L[3][:3]); xyz=[]; els=[]
    for i in range(4,4+n):
        xyz.append([float(L[i][0:10]),float(L[i][10:20]),float(L[i][20:30])])
        els.append(norm_el(L[i][31:34]))
    return np.array(xyz), els, ligname

def pose(path, ligname):
    """Docked ligand atoms carry the SDF-title residue name (serial 90001+),
    emitted in the SAME atom order as the crystal SDF (verified: 0 first-letter
    mismatches across RQ3/STC/DEX/CP6). We return coordinates only and take
    element identity from the crystal by position — FlexAID's pose element
    column (77-78) is unreliable (Cl->L truncation, H->Du dummies)."""
    tag=f" {ligname} "
    xyz=[]
    for ln in open(path):
        if ln.startswith(("HETATM","ATOM")) and tag in ln:
            xyz.append([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    return np.array(xyz)

def hrmsd(cx,ce,px):
    """Symmetry-corrected (element-blocked Hungarian) RMSD. Pose atoms are in
    crystal order, so element identity `ce` applies to both; assignment is
    optimized within each element block (accounts for topological symmetry)."""
    ce=np.array(ce)
    if len(cx)!=len(px): return None
    n=len(cx); asg=np.full(n,-1)
    for el in set(ce):
        ci=np.where(ce==el)[0]
        D=np.linalg.norm(cx[ci][:,None]-px[ci][None],axis=2)
        r,c=linear_sum_assignment(D)
        for a,b in zip(r,c): asg[ci[a]]=ci[b]
    d=cx-px[asg]; return float(np.sqrt((d*d).sum(1).mean()))

def cf_of(p):
    for ln in open(p):
        if ln.startswith("REMARK CF="):
            try: return float(ln.split("=")[1])
            except: return None
    return None

rows=[]
targets=sorted([d for d in os.listdir(FULL) if os.path.isdir(f"{FULL}/{d}") and re.match(r'^[0-9A-Z]{4}$',d)])
for t in targets:
    cx,ce,ligname=crystal(t)
    if cx is None: continue
    # pool all poses across base + r* dirs
    poses=sorted(set(glob.glob(f"{FULL}/{t}/r*/{t}_[0-9].pdb")+glob.glob(f"{FULL}/{t}/{t}_[0-9].pdb")))
    best=None; best_cf=None; nvalid=0; near=None
    for p in poses:
        px=pose(p, ligname)
        if len(px)==0: continue
        r=hrmsd(cx,ce,px)
        if r is None: continue
        nvalid+=1
        if best is None or r<best: best=r; near=p
    rows.append(dict(pdb=t, n_poses=nvalid, min_rmsd=best,
                     success=(best is not None and best<2.0), near=near))

with open(f"{ROOT}/poster_metric_results.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["pdb","n_poses","min_rmsd","success","near"])
    w.writeheader(); w.writerows(rows)

done=[r for r in rows if r["min_rmsd"] is not None]
succ=sum(1 for r in done if r["success"])
print(f"targets scored: {len(done)}/{len(targets)}")
if done:
    print(f"success (<2A pooled): {succ}/{len(done)} = {100*succ/len(done):.1f}%")
    print(f"median min-RMSD: {np.median([r['min_rmsd'] for r in done]):.2f} A")
