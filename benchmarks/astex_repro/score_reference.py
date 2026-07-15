#!/usr/bin/env python3
"""Reference-grade poster-metric scorer — Astex-85 GetCleft cognate reproduction.

Metric (user-confirmed, pooled reading): for each target, pool ALL poses across
all 10 restarts; success = any pose has symmetry-corrected heavy-atom RMSD
< 2.0 Å to the crystal ligand. Disclosed as the more forgiving of the two
poster-wording readings ("<2 Å in the 10 best results").

RMSD engine: spyrmsd (Meli & Biggin 2020, J. Cheminform. 12:49) — the
graph-isomorphism symmetry correction used by PoseBusters. NO superposition
(minimize=False): poses are already in the receptor frame, so this is a true
in-place RMSD, not an aligned one.

Reproducibility guarantees:
  * Connectivity comes from the crystal SDF *bond block* (authoritative,
    parsed directly) — never from bond perception or the pose file.
  * Atomic numbers from the crystal SDF atom block.
  * Heavy-atom only (standard docking-RMSD convention; H omitted).
  * Pose coordinates read by atom position — pose atom order == crystal SDF
    order (verified: 0 first-letter element mismatches on RQ3/STC/DEX/CP6).
    The pose PDB element column (77-78) is IGNORED (unreliable: Cl->L
    truncation, H->Du dummy atoms).
  * Falls back to element-blocked Hungarian ONLY if spyrmsd raises (logged).
"""
import numpy as np, glob, re, os, csv, sys
from spyrmsd import rmsd as spyr
from scipy.optimize import linear_sum_assignment

ROOT=os.path.dirname(os.path.abspath(__file__))
CACHE=f"{ROOT}/../astex_diverse/astex_diverse"
FULL=f"{ROOT}/full"
CUTOFF=2.0

# minimal symbol -> atomic number (covers all Astex ligand elements)
Z={'H':1,'C':6,'N':7,'O':8,'F':9,'NA':11,'MG':12,'P':15,'S':16,'CL':17,
   'K':19,'CA':20,'MN':25,'FE':26,'CO':27,'NI':28,'CU':29,'ZN':30,
   'BR':35,'I':53,'B':5,'SE':34}

def parse_sdf(path):
    """Return (coords Nx3, atomic_numbers N, adjacency NxN) from an SDF/MOL
    V2000 block, using the explicit bond list. No chemistry perception."""
    L=open(path).read().splitlines()
    counts=L[3]
    na=int(counts[0:3]); nb=int(counts[3:6])
    xyz=np.zeros((na,3)); anum=np.zeros(na,dtype=int)
    for i in range(na):
        ln=L[4+i]
        xyz[i]=[float(ln[0:10]),float(ln[10:20]),float(ln[20:30])]
        sym=ln[31:34].strip().upper()
        anum[i]=Z.get(sym,0)
    adj=np.zeros((na,na),dtype=int)
    for j in range(nb):
        ln=L[4+na+j]
        a=int(ln[0:3])-1; b=int(ln[3:6])-1
        adj[a,b]=1; adj[b,a]=1
    return xyz,anum,adj

def ligname_of(path):
    return open(path).readline().strip()

def pose_coords(path, ligname, nat):
    """Coordinates of the docked ligand, in file order. Returns Nx3 or None.

    The docked ligand is keyed on PDB serial >= 90000 (cols 6-11), a FlexAID
    output invariant that is robust to (a) numeric ligand names like '675'
    whose substring collides with coordinate fields, and (b) standard-residue
    ligand names like 'ALA'/'TYR' that collide with receptor residues. The
    ligname substring is NOT used for selection."""
    xyz=[]
    for ln in open(path):
        if not ln.startswith(("HETATM","ATOM")): continue
        try: serial=int(ln[6:11])
        except ValueError: continue
        if serial>=90000:
            xyz.append([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    if len(xyz)!=nat: return None
    return np.array(xyz)

def hungarian_fallback(cx,anum,px):
    """Element-blocked Hungarian RMSD (over-permissive; fallback only)."""
    asg=np.full(len(cx),-1)
    for el in set(anum.tolist()):
        idx=np.where(anum==el)[0]
        D=np.linalg.norm(cx[idx][:,None]-px[idx][None],axis=2)
        r,c=linear_sum_assignment(D)
        for a,b in zip(r,c): asg[idx[a]]=idx[b]
    d=cx-px[asg]; return float(np.sqrt((d*d).sum(1).mean()))

rows=[]
targets=sorted(d for d in os.listdir(FULL)
               if os.path.isdir(f"{FULL}/{d}") and re.match(r'^[0-9A-Z]{4}$',d))
for t in targets:
    sdf=f"{CACHE}/{t}/{t}_ligand.sdf"
    if not os.path.exists(sdf):
        rows.append(dict(pdb=t,n_poses=0,min_rmsd=None,success=False,engine="",near="")); continue
    lig=ligname_of(sdf)
    cxyz_all,anum_all,adj_all=parse_sdf(sdf)
    heavy=anum_all>1
    cxyz=cxyz_all[heavy]; anum=anum_all[heavy]
    adj=adj_all[np.ix_(heavy,heavy)]
    nat_all=len(anum_all)

    poses=sorted(set(glob.glob(f"{FULL}/{t}/r*/{t}_[0-9].pdb")
                    +glob.glob(f"{FULL}/{t}/{t}_[0-9].pdb")))
    coord_list=[]; path_list=[]
    for p in poses:
        pc=pose_coords(p,lig,nat_all)
        if pc is None: continue
        coord_list.append(pc[heavy]); path_list.append(p)

    if not coord_list:
        rows.append(dict(pdb=t,n_poses=0,min_rmsd=None,success=False,engine="",near="")); continue

    engine="spyrmsd"
    try:
        vals=spyr.symmrmsd(cxyz, coord_list, anum, anum, adj, adj,
                           center=False, minimize=False)
        vals=list(np.atleast_1d(vals))
    except Exception as e:
        engine="hungarian"
        vals=[hungarian_fallback(cxyz,anum,pc) for pc in coord_list]

    vals=np.array(vals,dtype=float)
    bi=int(np.argmin(vals)); best=float(vals[bi])
    rows.append(dict(pdb=t,n_poses=len(coord_list),min_rmsd=best,
                     success=bool(best<CUTOFF),engine=engine,near=path_list[bi]))

with open(f"{ROOT}/poster_metric_reference.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["pdb","n_poses","min_rmsd","success","engine","near"])
    w.writeheader(); w.writerows(rows)

done=[r for r in rows if r["min_rmsd"] is not None]
succ=sum(1 for r in done if r["success"])
nh=sum(1 for r in done if r["engine"]=="hungarian")
print(f"targets scored: {len(done)}/{len(targets)}")
if done:
    print(f"success (<{CUTOFF} A, pooled): {succ}/{len(done)} = {100*succ/len(done):.1f}%")
    print(f"median min-RMSD: {np.median([r['min_rmsd'] for r in done]):.2f} A")
    if nh: print(f"WARNING: {nh} targets used Hungarian fallback (spyrmsd failed)")
