#!/usr/bin/env python3
"""
make_fa_matrix_v2_science.py — Priority-1 VCT matrix corrections (science fix v2).

Applies only the lowest-risk entries from docs/VCT_MATRIX_AUDIT.md §6 Priority 1:
  [2-4]  C.2 x C.ar   -149.4 -> -65
  [40-13] O burial     +33.99 -> +90
  [40-14] O.3 burial   +43.24 -> +90
  [40-15] carboxylate  +29.56 -> +90

Preserves all other entries verbatim from the canonical MC_st0r5.2_6.dat.
"""
from __future__ import annotations

import os
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "MC_st0r5.2_6.dat"
DST = sys.argv[2] if len(sys.argv) > 2 else "MC_st0r5.2_6_v2_science.dat"

# 1-based type indices (upper triangle)
OVERRIDES: dict[tuple[int, int], float] = {
    (2, 4): -65.0,
    (40, 13): 90.0,
    (40, 14): 90.0,
    (40, 15): 90.0,
}

N = 40
order = [(ii, jj) for ii in range(1, N + 1) for jj in range(ii, N + 1)]


def fmt(v: float) -> str:
    if v == 0:
        return "0"
    return f"{v:.4g}"


changes = []
out_lines = []

with open(SRC) as f:
    raw = [ln.rstrip("\n") for ln in f if ln.strip() != ""]

if len(raw) != len(order):
    sys.exit(f"ERROR: expected {len(order)} lines, got {len(raw)} from {SRC}")

for (a, b), line in zip(order, raw):
    eq = line.index("=")
    prefix = line[: eq + 1]
    val = float(line[eq + 1 :].strip())
    key = (a, b)
    if key in OVERRIDES:
        newv = OVERRIDES[key]
        if abs(newv - val) > 1e-6:
            changes.append((a, b, val, newv))
        val = newv
    pad = max(1, 9 - len(fmt(val)))
    out_lines.append(f"{prefix}{' ' * pad}{fmt(val)}")

os.makedirs(os.path.dirname(DST) or ".", exist_ok=True)
with open(DST, "w") as f:
    f.write("\n".join(out_lines) + "\n")

clog = os.path.splitext(DST)[0] + "_CHANGES.txt"
with open(clog, "w") as f:
    f.write("MC_st0r5.2_6_v2_science.dat — Priority-1 corrections\n")
    f.write(f"source: {SRC}\n")
    f.write(f"changed: {len(changes)} of 820\n\n")
    for a, b, old, new in changes:
        f.write(f"  {a:2d}-{b:2d}  {old:10.4g} -> {new:10.4g}\n")

print(f"wrote {DST} ({len(out_lines)} lines, {len(changes)} changed)")
print(f"wrote {clog}")