"""Pure helpers for diagnosing CF ranking bias against near-native poses.

These functions operate on already-extracted CF/RMSD tables from real dock
artifacts (REMARK CF on pose PDBs). They do not reimplement the C++ scorer;
they quantify whether ranking on shipped CF values systematically disfavors
lower-RMSD poses.

Definitions
-----------
- Lower CF is better (FlexAID complementarity function convention).
- ``near_native_missed_by_top1``: some pose has RMSD <= threshold, but the
  CF-best pose does not — algorithmic ranking bias against near-natives.
- ``scoring_pathology_gap``: CF_top1 - CF_seed < -pathology_cut ⇒ elected
  decoy scores much better than crystal/seed (seed disfavored).
- ``search_never_beats_seed``: no pose has CF < CF_seed (search failure).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence


@dataclass(frozen=True)
class PoseCF:
    """One pose with CF score and optional RMSD to reference."""

    name: str
    cf: float
    rmsd: Optional[float] = None
    cf_com: Optional[float] = None
    cf_wal: Optional[float] = None


def rank_by_cf(poses: Sequence[PoseCF]) -> List[PoseCF]:
    """Return poses sorted by CF ascending (best first)."""
    return sorted(poses, key=lambda p: p.cf)


def rank_of_best_rmsd(poses: Sequence[PoseCF]) -> Optional[int]:
    """1-based CF-rank of the lowest-RMSD pose, or None if no RMSDs."""
    with_rmsd = [p for p in poses if p.rmsd is not None]
    if not with_rmsd:
        return None
    best = min(with_rmsd, key=lambda p: p.rmsd)  # type: ignore[arg-type]
    for i, p in enumerate(rank_by_cf(poses)):
        if p.name == best.name and p.cf == best.cf:
            return i + 1
    return None


def near_native_missed_by_top1(
    poses: Sequence[PoseCF],
    rmsd_cut: float = 2.0,
) -> bool:
    """True if any pose has RMSD <= cut but CF top-1 does not."""
    if not poses:
        return False
    oracle = any(p.rmsd is not None and p.rmsd <= rmsd_cut for p in poses)
    if not oracle:
        return False
    top1 = rank_by_cf(poses)[0]
    return not (top1.rmsd is not None and top1.rmsd <= rmsd_cut)


def scoring_pathology_gap(
    cf_top1: float,
    cf_seed: float,
    pathology_cut: float = 5.0,
) -> bool:
    """True if top-1 CF is much better (more negative) than seed/native.

    gap = cf_top1 - cf_seed; pathology when gap < -pathology_cut.
    """
    return (cf_top1 - cf_seed) < -pathology_cut


def search_never_beats_seed(poses: Sequence[PoseCF], cf_seed: float) -> bool:
    """True if no pose has better (lower) CF than the seed/INI."""
    if not poses:
        return True
    return all(p.cf >= cf_seed for p in poses)


def spearman_cf_rmsd(poses: Sequence[PoseCF]) -> Optional[float]:
    """Spearman rho(CF, RMSD). Positive ⇒ better CF (lower) pairs with lower RMSD."""
    pairs = [(p.cf, p.rmsd) for p in poses if p.rmsd is not None]
    n = len(pairs)
    if n < 3:
        return None
    xs = [a for a, _ in pairs]
    ys = [b for _, b in pairs]

    def _ranks(vals: Sequence[float]) -> List[float]:
        order = sorted(range(n), key=lambda i: vals[i])
        r = [0.0] * n
        for rank, i in enumerate(order):
            r[i] = float(rank + 1)
        return r

    rx, ry = _ranks(xs), _ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    denx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    deny = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if denx == 0.0 or deny == 0.0:
        return None
    return num / (denx * deny)


def wal_over_abs_com(cf_wal: Optional[float], cf_com: Optional[float]) -> Optional[float]:
    """|CF.wal| / |CF.com| — large ratio ⇒ soft-penetration / clash-contact fight."""
    if cf_wal is None or cf_com is None or cf_com == 0.0:
        return None
    return abs(cf_wal) / abs(cf_com)
