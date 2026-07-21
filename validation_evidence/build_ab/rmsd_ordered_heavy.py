#!/usr/bin/env python3
"""Ordered heavy-atom RMSD: pose ligand (RQ3/HETATM) vs crystal SDF."""
import sys, math, re
from pathlib import Path

def heavy_coords_pdb_ligand(path, resnames=("RQ3","LIG","UNL","MOL")):
    pts = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        res = line[17:20].strip()
        # FlexAID ligand often RQ3; also accept serial >= 90000
        serial = int(line[6:11]) if line[6:11].strip().isdigit() else 0
        if res not in resnames and serial < 90000:
            continue
        el = line[76:78].strip() if len(line) >= 78 else ""
        name = line[12:16].strip()
        if (el and el.upper() == "H") or name.upper().startswith("H"):
            continue
        try:
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        pts.append((x,y,z))
    return pts

def heavy_coords_sdf(path):
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines) and not re.match(r"^\s*\d+\s+\d+", lines[i]):
        i += 1
    if i >= len(lines):
        return []
    n_atoms = int(lines[i].split()[0])
    pts = []
    for j in range(i+1, i+1+n_atoms):
        parts = lines[j].split()
        if len(parts) < 4: continue
        x,y,z = float(parts[0]), float(parts[1]), float(parts[2])
        el = parts[3]
        if el.upper() == "H": continue
        pts.append((x,y,z))
    return pts

def rmsd(a,b):
    n = min(len(a), len(b))
    if n == 0: return float("nan")
    s = sum((a[i][0]-b[i][0])**2 + (a[i][1]-b[i][1])**2 + (a[i][2]-b[i][2])**2 for i in range(n))
    return math.sqrt(s/n)

pose, xtal = sys.argv[1], sys.argv[2]
pa = heavy_coords_pdb_ligand(pose)
xa = heavy_coords_sdf(xtal)
r = rmsd(pa, xa)
print(f"rmsd_ordered_heavy={r:.4f}")
print(f"n_pose={len(pa)} n_xtal={len(xa)} n_used={min(len(pa),len(xa))}")
