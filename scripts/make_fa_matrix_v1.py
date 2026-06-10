#!/usr/bin/env python3
"""
make_fa_matrix_v1.py — derive FA_matrix_v1.dat from the canonical MC_st0r5.2_6.dat
by conservative, physicochemically-grounded per-category rescaling.

Canonical FlexAID VCT atom-type numbering (see LIB/Mol2Reader.cpp / nrgrank_matrix.h):
   1=C.1  2=C.2  3=C.3   4=C.ar  5=C.cat
   6=N.1  7=N.2  8=N.3   9=N.4  10=N.ar 11=N.am 12=N.pl3
  13=O.2 14=O.3 15=O.co2 16=O.ar
  17=S.2 18=S.3 19=S.o   20=S.o2 21=S.ar
  22=P.3 23=F   24=Cl    25=Br   26=I    27=Se
  28..38 metals   39=DUMMY   40=SOLVENT

The .dat stores the upper triangle (i<=j) of a 40x40 symmetric matrix, one
"i- j = value" line per pair, 820 lines total. read_emat.cpp ignores the i-j
label and assigns by loop order, so we preserve the original prefix verbatim and
substitute only the numeric value (4 significant figures, matching source style).

Rescale categories (each pair assigned to exactly one, first match wins):
  diagonal (i==i)                         x1.00  (structural baseline; UNCHANGED)
    exception: (4,4) C.ar*C.ar            x1.07  (pi-stacking / London dispersion)
  halogen{F,Cl,Br,I} x {C,N,O}            x0.80  (-20%; spherical shell cannot
                                                  enforce 170 deg halogen-bond
                                                  geometry -> downscale to curb
                                                  sphere-edge over-counting)
  aromatic-aromatic cross (4,10)          x1.07  (C.ar*N.ar pi-system contact)
  H-bond donor*acceptor N{9-12}xO{13-16}  x1.10  (directional H-bond; if in-sphere
                                                  it is more likely a true H-bond
                                                  than the radial baseline assumed;
                                                  applied to favourable terms only)
  hydrophobic C*C  {1,2,3,4} off-diag     x1.05  (burial entropy / dispersion;
                                                  favourable terms only)
  polar-apolar  N/O{9-16} x C{1,2,3,4}    x0.95  (desolvation penalty for burying
                                                  a polar group on apolar carbon;
                                                  favourable terms only)

Hard constraint honoured: |Delta| <= 25% of each element's absolute value; signs
never inverted; matrix stays negative-dominant. C.cat (5) is a charged carbon and
is deliberately excluded from the hydrophobic and desolvation categories.
"""
import re, sys, os

SRC = sys.argv[1] if len(sys.argv) > 1 else "build/MC_st0r5.2_6.dat"
DST = sys.argv[2] if len(sys.argv) > 2 else "build/FA_matrix_v1.dat"

HALO   = {23, 24, 25, 26}
CARB_ALL = {1, 2, 3, 4, 5}      # all carbons (for halogen partner set)
CARB_HPHO = {1, 2, 3, 4}        # neutral carbons eligible for hydrophobic / desolv
NDON   = {9, 10, 11, 12}        # H-bond-capable nitrogen
OACC   = {13, 14, 15, 16}       # oxygen (acceptor / hydroxyl donor)
CNO    = CARB_ALL | NDON | OACC # halogen partners

def factor(a, b, v):
    """Return (multiplier, category) for upper-triangle pair (a<=b) with value v."""
    if a == b:
        if a == 4:
            return 1.07, "arom_self(C.ar.C.ar)"
        return 1.00, "diagonal"
    # halogen x C/N/O  (apply to both signs: shrink magnitude 20%)
    if (a in HALO and b in CNO) or (b in HALO and a in CNO):
        return 0.80, "halogen.CNO"
    # aromatic-aromatic cross C.ar x N.ar
    if {a, b} == {4, 10}:
        return 1.07, "arom_cross(C.ar.N.ar)"
    # H-bond donor*acceptor  N x O  (favourable only)
    if (a in NDON and b in OACC) or (b in NDON and a in OACC):
        return (1.10, "hbond.N.O") if v < 0 else (1.00, "hbond.N.O(skip+)")
    # hydrophobic C x C  (favourable only)
    if a in CARB_HPHO and b in CARB_HPHO:
        return (1.05, "hydrophobic.C.C") if v < 0 else (1.00, "hydrophobic.C.C(skip+)")
    # polar-apolar  N/O x neutral C  (favourable only)
    if ((a in (NDON | OACC) and b in CARB_HPHO) or
        (b in (NDON | OACC) and a in CARB_HPHO)):
        return (0.95, "polar_apolar.NO.C") if v < 0 else (1.00, "polar_apolar(skip+)")
    return 1.00, "unchanged"

def fmt(v):
    if v == 0:
        return "0"
    s = f"{v:.4g}"
    return s

changes = []
out_lines = []
i = j = None
# enumerate pairs in the same order read_emat consumes them: i in 0..N-1, j in i..N-1
order = []
N = 40
for ii in range(1, N + 1):
    for jj in range(ii, N + 1):
        order.append((ii, jj))

with open(SRC) as f:
    raw = [ln.rstrip("\n") for ln in f if ln.strip() != ""]
assert len(raw) == len(order) == 820, f"line/order mismatch {len(raw)} {len(order)}"

max_pct = 0.0
for (a, b), line in zip(order, raw):
    eq = line.index("=")
    prefix = line[:eq + 1]
    val = float(line[eq + 1:].strip())
    mult, cat = factor(a, b, val)
    newv = val * mult
    if mult != 1.0 and val != 0.0:
        pct = (newv - val) / abs(val) * 100.0
        max_pct = max(max_pct, abs(pct))
        changes.append((a, b, val, newv, mult, cat, pct))
    out_lines.append(f"{prefix}{' ' * (9 - len(fmt(newv)))}{fmt(newv)}")

with open(DST, "w") as f:
    f.write("\n".join(out_lines) + "\n")

print(f"wrote {DST}  ({len(out_lines)} lines)")
print(f"changed elements: {len(changes)}   max |Delta| = {max_pct:.2f}%")
assert max_pct <= 25.0 + 1e-6, "VIOLATION: an element changed by more than 25%"
# category tally
from collections import Counter
tally = Counter(c[5] for c in changes)
for k, v in sorted(tally.items()):
    print(f"  {k:28s} {v}")

# emit a changelog companion
clog = os.path.join(os.path.dirname(DST), "FA_matrix_v1_CHANGES.txt")
with open(clog, "w") as f:
    f.write("FA_matrix_v1.dat — changes vs canonical MC_st0r5.2_6.dat\n")
    f.write("source md5: derived from build/MC_st0r5.2_6.dat (canonical, v23 Arm A)\n")
    f.write(f"changed elements: {len(changes)} of 820   max |Delta| = {max_pct:.2f}% (cap 25%)\n\n")
    f.write(f"{'pair':>9}  {'category':28s} {'old':>10} {'new':>10} {'mult':>6} {'d%':>7}\n")
    for a, b, val, newv, mult, cat, pct in changes:
        f.write(f"  {a:2d}-{b:2d}  {cat:28s} {val:10.4g} {newv:10.4g} {mult:6.2f} {pct:+7.1f}\n")
print(f"wrote {clog}")
