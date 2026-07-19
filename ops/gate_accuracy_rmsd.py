#!/usr/bin/env python
"""METHODOLOGY.md §3 — Astex-85 accuracy scorer (top-1 rank-0 RMSD, spyrmsd, 2.0 A).
Run in the `python` conda env (spyrmsd 0.9.0).
Usage: gate_accuracy_rmsd.py <dock_root> <cache_dir> <targets_csv>
  dock_root : dir with <tag>/<TARGET>/d_0.pdb  (rank-0 elected pose)
  cache_dir : dir with <TARGET>/<TARGET>_ligand.sdf  (reference)
Emits per-target RMSD and success (RMSD<=2.0) for each engine tag; compares tags.
"""
import sys, os, glob, numpy as np
from spyrmsd import rmsd as spyr, io as spyrio

def load_mol(path):
    m = spyrio.loadmol(path)
    return m

def sym_rmsd(ref, pose):
    ref.strip(); pose.strip()  # heavy atoms
    cref = ref.coordinates; cpose = pose.coordinates
    return spyr.symmrmsd(cref, cpose, ref.atomicnums, pose.atomicnums,
                         ref.adjacency_matrix, pose.adjacency_matrix)

def main():
    dock_root, cache, tcsv = sys.argv[1], sys.argv[2], sys.argv[3]
    targets = [t.strip() for t in tcsv.replace(',',' ').split()]
    tags = sorted([d for d in os.listdir(dock_root) if os.path.isdir(os.path.join(dock_root,d))])
    print(f"tags={tags} targets={targets}")
    res = {}
    for tag in tags:
        res[tag] = {}
        for T in targets:
            ref_sdf = os.path.join(cache, T, f"{T}_ligand.sdf")
            pose = os.path.join(dock_root, tag, T, "d_0.pdb")  # rank-0 elected
            if not (os.path.exists(ref_sdf) and os.path.exists(pose)):
                res[tag][T] = None; continue
            try:
                r = sym_rmsd(load_mol(ref_sdf), load_mol(pose))
                res[tag][T] = float(r)
            except Exception as e:
                res[tag][T] = f"ERR:{e}"
    # report
    print(f"\n{'target':8} " + " ".join(f"{t:>12}" for t in tags))
    for T in targets:
        row = " ".join((f"{res[t][T]:12.3f}" if isinstance(res[t].get(T),float) else f"{str(res[t].get(T)):>12}") for t in tags)
        print(f"{T:8} {row}")
    print("\nsuccess@2.0A:")
    for t in tags:
        s = sum(1 for T in targets if isinstance(res[t].get(T),float) and res[t][T]<=2.0)
        n = sum(1 for T in targets if isinstance(res[t].get(T),float))
        print(f"  {t}: {s}/{n}")
    # flips
    if len(tags)==2:
        a,b=tags
        print(f"\nflips (success change {a}->{b}):")
        for T in targets:
            va,vb=res[a].get(T),res[b].get(T)
            if isinstance(va,float) and isinstance(vb,float):
                sa,sb=va<=2.0,vb<=2.0
                if sa!=sb: print(f"  {T}: {a}={va:.2f}({'PASS' if sa else 'fail'}) -> {b}={vb:.2f}({'PASS' if sb else 'fail'})  <-- FLIP")
        print("  (no flips = accuracy preserved)")

if __name__=="__main__": main()
