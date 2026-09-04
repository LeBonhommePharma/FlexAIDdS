#!/usr/bin/env python3
# Apache-2.0 (c) 2026 Le Bonhomme Pharma
"""Continuous Tanimoto over FlexAIDdS atom-type-pair contact-surface profiles.

WHAT A PROFILE IS
-----------------
The FlexAID contact function contracts a rich per-type-pair profile into one
scalar, ``CF.com = sum eps(type_i, type_j) * S(i, j)``.  With the engine gate
``FLEXAIDDS_CONTACT_PROFILE=1`` set, ``LIB/vcfunction.cpp`` keeps ``S(i, j)``
BEFORE that contraction and every emitted pose gets a sidecar next to it::

    <pose>.pdb  ->  <pose>.cprof.csv

Each sidecar holds the total Voronoi contact surface, in square angstroms,
accumulated per UNORDERED atom-type pair, split into

    area_inter  different-molecule contacts  -> the ligand/receptor INTERFACE,
                including ligand contacts against a flexed side chain
    area_intra  same-molecule contacts       -> ligand-internal non-bonded
                contacts AND flexed-side-chain-versus-receptor contacts (both
                endpoints of the latter are receptor residues).  A mixed
                conformer channel; do not read it as "ligand strain".
    area_total  = area_inter + area_intra, exactly

plus ``cf_pair``, the CF contribution the engine actually credited to that same
pair.  Rows with no signal at all are omitted from the file; ``# ntypes`` in the
header is enough to rebuild the dense ``ntypes*(ntypes+1)/2`` vector with zeros,
which is what this script does.  For the 40-type MC_st0r5.2_6 matrix that is a
820-dimensional vector, matching the interaction-matrix dimensionality.

WHY THIS DESCRIPTOR
-------------------
The vector is RECEPTOR-FRAME AGNOSTIC.  It names atom types and areas; it never
names a receptor coordinate, atom index, residue or rotamer.  That is exactly
the property a flexible-receptor experiment needs, because the receptor frame is
what differs between the rigid and the flexible arm, and it is what makes a
coordinate-space overlap metric unattributable: you cannot tell whether an
overlap means the pose is wrong or the receptor moved.  Two profiles can be
compared even when the two receptors are not in the same conformation.

ORACLE STATUS -- READ THIS BEFORE USING A NUMBER FROM THIS SCRIPT
-----------------------------------------------------------------
* Comparing the profiles of TWO DOCKED POSES is ORACLE-FREE.  Nothing about the
  answer is consumed.  This is legitimate inside a production run: pose-vs-pose
  Tanimoto is a self-similarity / ensemble-diversity measure.

* Comparing a docked pose against the profile of the NATIVE complex (the
  ``<prefix>_native.cprof.csv`` written by ``LIB/native_score.cpp`` under
  ``FLEXAIDDS_SCORE_NATIVE``) is an **ORACLE METRIC**.  It consumes the crystal
  answer.  It is valid for benchmark analysis and for diagnosis -- "did the pose
  reproduce the native interaction profile even though the RMSD is large?" --
  and it is **NOT usable as a production scoring, ranking or selection term**.
  Any use of the native profile inside the docking loop, inside pose election,
  or inside any quantity that feeds a reported success rate, is a leak.  This
  script refuses nothing and enforces nothing; the discipline is yours.

THE ZERO-EPSILON BLIND SPOT
---------------------------
About three quarters of the interaction-matrix entries are exactly 0.0 and
several whole atom-type rows are entirely zero, so the atoms on those rows are
invisible to CF no matter how much surface they bury.  The engine deliberately
does NOT filter accumulation on the matrix, so those pairs appear here with
``area_total > 0`` and ``cf_pair == 0.0``.  ``--blind-spot`` reports exactly that
surface: it is the one thing CF itself can never show you, because the missing
contribution is by construction zero.

SANITY CHECK BEFORE TRUSTING A PROFILE
--------------------------------------
The sidecar header carries ``# cf_pair_sum``, the sum of the per-pair CF
contributions.  The engine accumulates those with the very same value it adds to
``cfs->com``, so ``cf_pair_sum`` must reproduce the pose's ``REMARK CF.com`` to
float32 precision.  If it does not, the profile does not belong to the pose it
sits beside and every number derived from it is meaningless.  Check it once per
new run configuration::

    grep '^# cf_pair_sum' pose_0.cprof.csv
    grep '^REMARK CF.com'  pose_0.pdb

CONTINUOUS TANIMOTO
-------------------
For non-negative vectors a and b::

    T(a, b) = (a . b) / (|a|^2 + |b|^2 - a . b)

T = 1 iff a == b (for non-zero vectors); T = 0 iff the supports are disjoint.
Both the RAW vector and the SUM-NORMALISED profile (a / sum(a)) are reported,
because they mean different things and only one of them is usually the question:

    raw         magnitude AND pattern.  Two poses burying the same TYPES in the
                same proportions but one burying twice the surface score < 1.
    normalised  pattern only.  Scale-free: it answers "is this the same kind of
                interface?" regardless of how much of it there is.

A pose that half-buries the native interface has a high normalised T and a low
raw T; a pose that buries the right amount of the wrong types has the opposite.
Reporting one without the other hides which of those two happened.

USAGE
-----
  contact_profile_tanimoto.py A.cprof.csv B.cprof.csv [--column area_inter]
  contact_profile_tanimoto.py --matrix run_dir/*.cprof.csv
  contact_profile_tanimoto.py --blind-spot A.cprof.csv
  ... add --json for machine-readable output.

Pure standard library.  Python >= 3.9.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

AREA_COLUMNS = ("area_total", "area_inter", "area_intra", "cf_pair")
DEFAULT_COLUMN = "area_inter"


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
class Profile:
    """One parsed .cprof.csv sidecar."""

    def __init__(self, path: str, ntypes: int, meta: Dict[str, str],
                 rows: Dict[Tuple[int, int], Dict[str, float]]):
        self.path = path
        self.ntypes = ntypes
        self.meta = meta
        self.rows = rows

    # -- packed upper-triangular geometry ---------------------------------- #
    @property
    def npairs(self) -> int:
        return self.ntypes * (self.ntypes + 1) // 2

    def pair_index(self, ti: int, tj: int) -> int:
        """Mirror of flexaids::contact_profile::pair_index (1-based types)."""
        lo, hi = (ti - 1, tj - 1) if ti <= tj else (tj - 1, ti - 1)
        return lo * self.ntypes - (lo * (lo - 1)) // 2 + (hi - lo)

    def pair_of_index(self, k: int) -> Tuple[int, int]:
        for ti in range(1, self.ntypes + 1):
            for tj in range(ti, self.ntypes + 1):
                if self.pair_index(ti, tj) == k:
                    return (ti, tj)
        raise IndexError(k)

    def vector(self, column: str) -> List[float]:
        """Dense packed vector for `column`, zero-filled for omitted rows."""
        if column not in AREA_COLUMNS:
            raise ValueError("unknown column %r (choose from %s)"
                             % (column, ", ".join(AREA_COLUMNS)))
        v = [0.0] * self.npairs
        for (ti, tj), cols in self.rows.items():
            v[self.pair_index(ti, tj)] = cols.get(column, 0.0)
        return v

    def label(self) -> str:
        return self.meta.get("pose_file") or os.path.basename(self.path)

    def is_oracle(self) -> bool:
        return self.meta.get("source", "") == "native_crystal_pose"


def read_profile(path: str) -> Profile:
    meta: Dict[str, str] = {}
    rows: Dict[Tuple[int, int], Dict[str, float]] = {}
    header: Optional[List[str]] = None

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                if "=" in body:
                    key, _, val = body.partition("=")
                    meta[key.strip()] = val.strip()
                continue
            fields = [f.strip() for f in line.split(",")]
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise ValueError("%s: malformed row %r" % (path, line))
            rec = dict(zip(header, fields))
            ti, tj = int(rec["type_i"]), int(rec["type_j"])
            key = (ti, tj) if ti <= tj else (tj, ti)
            rows[key] = {c: float(rec[c]) for c in AREA_COLUMNS if c in rec}

    if header is None:
        raise ValueError("%s: no CSV header line found" % path)
    if "ntypes" not in meta:
        raise ValueError("%s: header is missing '# ntypes ='; cannot size the "
                         "dense vector" % path)
    ntypes = int(meta["ntypes"])
    if ntypes <= 0:
        raise ValueError("%s: ntypes = %d" % (path, ntypes))
    for (ti, tj) in rows:
        if not (1 <= ti <= ntypes and 1 <= tj <= ntypes):
            raise ValueError("%s: type pair (%d,%d) outside 1..%d"
                             % (path, ti, tj, ntypes))
    return Profile(path, ntypes, meta, rows)


# --------------------------------------------------------------------------- #
# Metric
# --------------------------------------------------------------------------- #
def tanimoto(a: Sequence[float], b: Sequence[float]) -> float:
    """Continuous (Tanimoto / Jaccard) coefficient.

    T(a, b) = (a . b) / (|a|^2 + |b|^2 - a . b)

    Returns 0.0 when both vectors are all-zero (no interface to compare), which
    is the convention that keeps a degenerate pose from scoring as a perfect
    match against another degenerate pose.
    """
    if len(a) != len(b):
        raise ValueError("dimension mismatch: %d vs %d" % (len(a), len(b)))
    dot = sum(x * y for x, y in zip(a, b))
    na = sum(x * x for x in a)
    nb = sum(y * y for y in b)
    den = na + nb - dot
    if den <= 0.0:
        return 0.0
    return dot / den


def normalise(v: Sequence[float]) -> List[float]:
    """Sum-normalised profile: pattern with magnitude divided out."""
    s = sum(v)
    if s <= 0.0:
        return [0.0] * len(v)
    return [x / s for x in v]


def compare(pa: Profile, pb: Profile, column: str) -> Dict[str, object]:
    if pa.ntypes != pb.ntypes:
        raise ValueError("ntypes mismatch: %s has %d, %s has %d — these "
                         "profiles come from different interaction matrices "
                         "and are not comparable"
                         % (pa.path, pa.ntypes, pb.path, pb.ntypes))
    va, vb = pa.vector(column), pb.vector(column)
    return {
        "a": pa.path,
        "b": pb.path,
        "a_label": pa.label(),
        "b_label": pb.label(),
        "column": column,
        "ntypes": pa.ntypes,
        "npairs": pa.npairs,
        "tanimoto_raw": tanimoto(va, vb),
        "tanimoto_normalised": tanimoto(normalise(va), normalise(vb)),
        "sum_a": sum(va),
        "sum_b": sum(vb),
        "support_a": sum(1 for x in va if x != 0.0),
        "support_b": sum(1 for x in vb if x != 0.0),
        "support_shared": sum(1 for x, y in zip(va, vb) if x != 0.0 and y != 0.0),
        "oracle": bool(pa.is_oracle() or pb.is_oracle()),
    }


# --------------------------------------------------------------------------- #
# Blind spot
# --------------------------------------------------------------------------- #
def blind_spot(p: Profile) -> Dict[str, object]:
    """Surface the contact function priced at exactly nothing.

    A pair with area_total > 0 and cf_pair == 0.0 is buried surface that made no
    contribution to CF -- either because the interaction-matrix entry for that
    type pair is exactly 0.0, or because one of the two atom-type rows is
    entirely zero. This is invisible from CF alone: the missing contribution is
    by construction zero.
    """
    dead: List[Tuple[Tuple[int, int], float]] = []
    live_area = 0.0
    dead_area = 0.0
    for key, cols in sorted(p.rows.items()):
        area = cols.get("area_total", 0.0)
        cf = cols.get("cf_pair", 0.0)
        if area <= 0.0:
            continue
        if cf == 0.0:
            dead.append((key, area))
            dead_area += area
        else:
            live_area += area
    total = live_area + dead_area
    return {
        "path": p.path,
        "label": p.label(),
        "ntypes": p.ntypes,
        "area_scored": live_area,
        "area_unscored": dead_area,
        "area_total": total,
        "fraction_unscored": (dead_area / total) if total > 0.0 else 0.0,
        "unscored_pairs": [{"type_i": k[0], "type_j": k[1], "area_total": a}
                           for k, a in sorted(dead, key=lambda t: -t[1])],
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _oracle_banner(out) -> None:
    print("# NOTE: one side of this comparison is the NATIVE complex profile.",
          file=out)
    print("#       This is an ORACLE metric: it consumes the crystal answer.",
          file=out)
    print("#       Valid for benchmark analysis and diagnosis; NOT usable as a",
          file=out)
    print("#       production scoring, ranking or selection term.", file=out)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="contact_profile_tanimoto.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("profiles", nargs="*", metavar="PROFILE.cprof.csv")
    ap.add_argument("--column", default=DEFAULT_COLUMN, choices=AREA_COLUMNS,
                    help="which column to build the vector from "
                         "(default: %s, the ligand/receptor interface)"
                         % DEFAULT_COLUMN)
    ap.add_argument("--matrix", action="store_true",
                    help="all-vs-all over every PROFILE given")
    ap.add_argument("--blind-spot", action="store_true",
                    help="report buried surface whose CF contribution is "
                         "exactly 0.0, per profile")
    ap.add_argument("--json", action="store_true",
                    help="machine-readable output")
    args = ap.parse_args(argv)

    if not args.profiles:
        ap.error("need at least one .cprof.csv")

    try:
        loaded = [read_profile(p) for p in args.profiles]
    except (OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    if args.blind_spot:
        results = [blind_spot(p) for p in loaded]
        if args.json:
            json.dump(results, sys.stdout, indent=2)
            sys.stdout.write("\n")
            return 0
        for r in results:
            print("%s" % r["label"])
            print("  scored surface   : %12.3f A^2" % r["area_scored"])
            print("  UNSCORED surface : %12.3f A^2  (%.1f%% of contact surface,"
                  " %d type pairs)"
                  % (r["area_unscored"], 100.0 * r["fraction_unscored"],
                     len(r["unscored_pairs"])))
            for row in r["unscored_pairs"][:15]:
                print("      types %2d-%2d  %10.3f A^2  cf_pair = 0.0"
                      % (row["type_i"], row["type_j"], row["area_total"]))
            if len(r["unscored_pairs"]) > 15:
                print("      ... %d more" % (len(r["unscored_pairs"]) - 15))
        return 0

    if args.matrix:
        pairs = [(i, j) for i in range(len(loaded))
                 for j in range(i + 1, len(loaded))]
    else:
        if len(loaded) != 2:
            ap.error("give exactly two profiles, or use --matrix")
        pairs = [(0, 1)]

    try:
        results = [compare(loaded[i], loaded[j], args.column) for i, j in pairs]
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    if args.json:
        json.dump(results, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if any(r["oracle"] for r in results):
        _oracle_banner(sys.stdout)
    print("# column = %s   (raw = magnitude+pattern, norm = pattern only)"
          % args.column)
    print("%-38s %-38s %9s %9s" % ("A", "B", "T_raw", "T_norm"))
    for r in results:
        print("%-38s %-38s %9.6f %9.6f"
              % (os.path.basename(str(r["a"]))[:38],
                 os.path.basename(str(r["b"]))[:38],
                 r["tanimoto_raw"], r["tanimoto_normalised"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
