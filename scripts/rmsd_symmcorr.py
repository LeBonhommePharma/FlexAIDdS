#!/usr/bin/env python3
"""Symmetry-corrected elected-pose RMSD — the single claim-metric producer (#365).

WHY THIS FILE EXISTS
  "Success at 2.0 A" had two producers with two meanings. The engine gates
  `success_rmsd` on `rmsd_to_crystal`, which is ORDERED/SERIAL by design
  (`LIB/DatasetRunner.cpp`: the identity atom mapping). METHODOLOGY.md §claim,
  `docs/swarm/2026-08-13/score_canonical.py` and `ops/gate_accuracy_rmsd.py`
  all mean spyRMSD graph automorphism. Both definitions were live.

  This module makes the symmetry-corrected value a first-class, joinable
  quantity so the claim gates can converge on it. It does NOT redefine
  `rmsd_to_crystal`, which keeps its documented serial meaning.

WHAT IT IS NOT
  It is not a sixth RMSD implementation. METHODOLOGY.md §0 pins five and says
  "do not invent a sixth". The symmetry correction here is performed by
  `spyrmsd.rmsd.symmrmsd` under the *exact* invocation contract of §0 method 2
  (`benchmarks/astex_repro/score_reference.py`, "the strongest instrument"):

    * connectivity from the crystal SDF *bond block* — never bond perception,
      never the pose file
    * atomic numbers from the crystal SDF atom block
    * heavy atoms only
    * center=False, minimize=False  → true in-place RMSD in the receptor
      frame, never superposed
    * ligand atoms selected on PDB serial >= 90000 (FlexAID output invariant)

  Method 2 is a top-level script that executes on import, so it cannot be
  called as a library. This module is that same contract, callable, and
  `tests/test_rmsd_symmcorr.py` pins the agreement.

FAIL-CLOSED
  There is NO fallback to element-blocked Hungarian and NO fallback to the
  serial value. Element-only Hungarian is over-permissive — the repo measured
  it inflating the pool ceiling from 48.8% to 57.8% — and silently substituting
  a weaker metric is precisely the defect being removed. A row that cannot be
  scored gets an empty `rmsd_symmcorr` and a machine-readable `status`.

ATTRIBUTABILITY
  Every row carries the pose artifact's SHA-256 and whether it verified against
  `pose_sha256` from result.csv, plus the producer identity and the pinned
  spyrmsd version. A number that cannot be tied to a specific pose is not a
  claim number.

READ-ONLY
  Reads result.csv and pose artifacts; writes only to --out. Never modifies a
  result.csv or anything under the results tree.

USAGE
  python3 scripts/rmsd_symmcorr.py <campaign_dir> --out symmcorr.csv
  python3 scripts/rmsd_symmcorr.py <campaign_dir> --cache <astex_diverse_dir> --out s.csv

  Then join it into the claim table:
  python3 scripts/aggregate_claim_metrics.py <campaign_dir> --symmcorr symmcorr.csv

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import sys
from pathlib import Path
from typing import Any

METRIC_NAME = "spyrmsd_graph_automorphism_inplace"
PRODUCER = "scripts/rmsd_symmcorr.py"
LIGAND_SERIAL_MIN = 90000  # FlexAID output invariant (METHODOLOGY.md §0 method 2)

# Minimal symbol -> atomic number, identical to §0 method 2's table.
_Z = {
    "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "NA": 11, "MG": 12, "P": 15,
    "S": 16, "CL": 17, "K": 19, "CA": 20, "MN": 25, "FE": 26, "CO": 27,
    "NI": 28, "CU": 29, "ZN": 30, "BR": 35, "I": 53, "B": 5, "SE": 34,
}

SIDECAR_FIELDS = (
    "pdb_id",
    "rmsd_symmcorr",
    "status",
    "metric",
    "producer",
    "spyrmsd_version",
    "pose_sha256",
    "sha_verified",
    "n_heavy",
    "rmsd_to_crystal_serial",
    "pose_artifact",
    "crystal_sdf",
)


class SymmCorrUnavailable(RuntimeError):
    """spyrmsd is not importable. Callers must fail closed, never substitute."""


def _spyrmsd():
    try:
        import numpy  # noqa: F401
        from spyrmsd import rmsd as spyr
        import spyrmsd as _pkg
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SymmCorrUnavailable(
            "spyrmsd is not installed. This producer refuses to fall back to a "
            "weaker metric (element-blocked Hungarian is over-permissive; the "
            "serial metric is over-strict). Install spyrmsd and re-run."
        ) from exc
    return spyr, getattr(_pkg, "__version__", "unknown")


# ── inputs ───────────────────────────────────────────────────────────────────


def parse_sdf(path: str | os.PathLike[str]):
    """(coords, atomic_numbers, adjacency) from an SDF/MOL V2000 block.

    Uses the explicit bond list. No chemistry perception, no bond guessing.
    """
    import numpy as np

    lines = Path(path).read_text(errors="ignore").splitlines()
    counts = lines[3]
    n_atoms = int(counts[0:3])
    n_bonds = int(counts[3:6])
    xyz = np.zeros((n_atoms, 3))
    anum = np.zeros(n_atoms, dtype=int)
    for i in range(n_atoms):
        ln = lines[4 + i]
        xyz[i] = [float(ln[0:10]), float(ln[10:20]), float(ln[20:30])]
        anum[i] = _Z.get(ln[31:34].strip().upper(), 0)
    adj = np.zeros((n_atoms, n_atoms), dtype=int)
    for j in range(n_bonds):
        ln = lines[4 + n_atoms + j]
        a = int(ln[0:3]) - 1
        b = int(ln[3:6]) - 1
        adj[a, b] = 1
        adj[b, a] = 1
    return xyz, anum, adj


def read_pose_artifact(path: str) -> tuple[bytes, str] | tuple[None, None]:
    """Return (raw_bytes, actual_path) for a pose PDB, transparently degzipping.

    Campaign trees compress elected poses after the run; the SHA-256 in
    result.csv is over the *uncompressed* bytes, so hashing must happen after
    decompression for the receipt to verify.
    """
    if not path:
        return None, None
    if os.path.exists(path):
        return Path(path).read_bytes(), path
    if os.path.exists(path + ".gz"):
        with gzip.open(path + ".gz", "rb") as fh:
            return fh.read(), path + ".gz"
    return None, None


def pose_ligand_coords(raw: bytes):
    """Docked-ligand heavy-atom coordinates, in file order.

    Selection is on PDB serial >= 90000 (METHODOLOGY.md §0 method 2): robust to
    numeric ligand names and to ligand names that collide with standard
    residues. The pose PDB element column (77-78) is deliberately NOT used for
    selection — it is unreliable (Cl->L truncation, H->Du dummies).
    """
    import numpy as np

    out = []
    for ln in raw.decode("utf-8", "ignore").splitlines():
        if not ln.startswith(("HETATM", "ATOM")):
            continue
        try:
            serial = int(ln[6:11])
        except ValueError:
            continue
        if serial >= LIGAND_SERIAL_MIN:
            out.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    if not out:
        return None
    return np.array(out)


# ── the metric ───────────────────────────────────────────────────────────────


def symmcorr_rmsd(crystal_sdf: str, pose_raw: bytes) -> dict[str, Any]:
    """Symmetry-corrected in-place RMSD of a docked pose against a crystal SDF.

    Returns {"rmsd": float|None, "status": str, "n_heavy": int}. Fail-closed:
    on any problem `rmsd` is None and `status` names the reason. Never returns
    a Hungarian or serial value in place of the symmetry-corrected one.
    """
    import numpy as np

    spyr, _ = _spyrmsd()

    try:
        cx_all, an_all, adj_all = parse_sdf(crystal_sdf)
    except Exception:
        return {"rmsd": None, "status": "sdf_parse_error", "n_heavy": 0}

    heavy = an_all > 1
    cx = cx_all[heavy]
    an = an_all[heavy]
    adj = adj_all[np.ix_(heavy, heavy)]
    n_heavy = int(heavy.sum())
    if n_heavy == 0:
        return {"rmsd": None, "status": "sdf_no_heavy_atoms", "n_heavy": 0}

    pose = pose_ligand_coords(pose_raw)
    if pose is None:
        return {"rmsd": None, "status": "no_ligand_serial_ge_90000", "n_heavy": n_heavy}

    # The pose may or may not carry the SDF's hydrogens. Accept exactly the two
    # unambiguous cardinalities; anything else is a mismatch and fails closed.
    if len(pose) == len(an_all):
        pose = pose[heavy]
    elif len(pose) != n_heavy:
        return {
            "rmsd": None,
            "status": f"atom_count_mismatch_pose{len(pose)}_vs_heavy{n_heavy}",
            "n_heavy": n_heavy,
        }

    try:
        val = spyr.symmrmsd(cx, pose, an, an, adj, adj, center=False, minimize=False)
    except Exception as exc:
        # NO Hungarian fallback. See module docstring.
        return {
            "rmsd": None,
            "status": f"spyrmsd_error:{type(exc).__name__}",
            "n_heavy": n_heavy,
        }
    return {"rmsd": float(np.atleast_1d(val)[0]), "status": "ok", "n_heavy": n_heavy}


def serial_rmsd(crystal_sdf: str, pose_raw: bytes) -> float | None:
    """Ordered/identity-mapping RMSD — the `rmsd_to_crystal` definition.

    Provided ONLY so callers and tests can assert the invariant
    symmcorr <= serial on identical inputs. Never used as a claim metric.
    """
    import numpy as np

    try:
        cx_all, an_all, _ = parse_sdf(crystal_sdf)
    except Exception:
        return None
    heavy = an_all > 1
    cx = cx_all[heavy]
    pose = pose_ligand_coords(pose_raw)
    if pose is None:
        return None
    if len(pose) == len(an_all):
        pose = pose[heavy]
    elif len(pose) != len(cx):
        return None
    d = cx - pose
    return float(np.sqrt((d * d).sum(1).mean()))


# ── campaign walk ────────────────────────────────────────────────────────────


def default_cache() -> str:
    env = os.environ.get("FLEXAIDDS_CACHE_V2")
    if env:
        return env
    results = os.environ.get("FLEXAIDDS_RESULTS")
    if results:
        return os.path.join(results, "cache_v2", "astex_diverse")
    home_results = os.path.expanduser("~/flexaidds_results")
    if os.path.isdir(home_results):
        return os.path.join(home_results, "cache_v2", "astex_diverse")
    return ""


def find_result_csvs(campaign_dir: str) -> list[Path]:
    root = Path(campaign_dir)
    run = root / "run"
    base = run if run.is_dir() else root
    return sorted(base.glob("*/result.csv"))


def score_campaign(campaign_dir: str, cache: str) -> list[dict[str, Any]]:
    _, spy_version = _spyrmsd()
    rows: list[dict[str, Any]] = []
    for csv_path in find_result_csvs(campaign_dir):
        try:
            with csv_path.open(newline="") as fh:
                recs = [dict(r) for r in csv.DictReader(fh)]
        except Exception:
            continue
        if not recs:
            continue
        rec = recs[0]
        pdb_id = (rec.get("pdb_id") or csv_path.parent.name).strip()
        declared_sha = (rec.get("pose_sha256") or "").strip()
        serial_csv = (rec.get("rmsd_to_crystal") or "").strip()

        sdf = os.path.join(cache, pdb_id, f"{pdb_id}_ligand.sdf")
        pose_path = (rec.get("elected_pose_path") or "").strip()
        if not pose_path:
            pose_path = str(csv_path.parent / "elected_pose.pdb")

        raw, actual = read_pose_artifact(pose_path)
        out: dict[str, Any] = {
            "pdb_id": pdb_id,
            "rmsd_symmcorr": "",
            "status": "",
            "metric": METRIC_NAME,
            "producer": PRODUCER,
            "spyrmsd_version": spy_version,
            "pose_sha256": "",
            "sha_verified": "0",
            "n_heavy": "",
            "rmsd_to_crystal_serial": serial_csv,
            "pose_artifact": actual or pose_path,
            "crystal_sdf": sdf,
        }
        if raw is None:
            out["status"] = "no_pose_artifact"
            rows.append(out)
            continue
        actual_sha = hashlib.sha256(raw).hexdigest()
        out["pose_sha256"] = actual_sha
        out["sha_verified"] = "1" if (declared_sha and actual_sha == declared_sha) else "0"
        if declared_sha and actual_sha != declared_sha:
            # The artifact is not the pose the claim row is about. Fail closed.
            out["status"] = "sha_mismatch"
            rows.append(out)
            continue
        if not os.path.exists(sdf):
            out["status"] = "no_crystal_sdf"
            rows.append(out)
            continue
        res = symmcorr_rmsd(sdf, raw)
        out["status"] = res["status"]
        out["n_heavy"] = res["n_heavy"]
        if res["rmsd"] is not None:
            out["rmsd_symmcorr"] = f"{res['rmsd']:.4f}"
        rows.append(out)
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("campaign_dir", help="Campaign dir (with run/<TARGET>/result.csv)")
    ap.add_argument("--cache", default=default_cache(), help="astex_diverse cache dir")
    ap.add_argument("--out", required=True, help="Sidecar CSV to write")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.campaign_dir):
        print(f"not a directory: {args.campaign_dir}", file=sys.stderr)
        return 2
    if not args.cache or not os.path.isdir(args.cache):
        print(
            "give --cache <astex_diverse dir> or set FLEXAIDDS_CACHE_V2 / "
            "FLEXAIDDS_RESULTS",
            file=sys.stderr,
        )
        return 2
    try:
        rows = score_campaign(args.campaign_dir, args.cache)
    except SymmCorrUnavailable as exc:
        print(f"[SYMMCORR-FAIL] {exc}", file=sys.stderr)
        return 2

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(SIDECAR_FIELDS))
        w.writeheader()
        w.writerows(rows)

    ok = [r for r in rows if r["status"] == "ok"]
    scored = [float(r["rmsd_symmcorr"]) for r in ok]
    n_pass = sum(1 for v in scored if v <= 2.0)
    print(f"targets seen    : {len(rows)}")
    print(f"scored (ok)     : {len(ok)}")
    print(f"sha verified    : {sum(1 for r in rows if r['sha_verified'] == '1')}")
    print(f"<= 2.0 A        : {n_pass}/{len(ok)}")
    if len(ok) != len(rows):
        print("unscored reasons:")
        reasons: dict[str, int] = {}
        for r in rows:
            if r["status"] != "ok":
                reasons[r["status"]] = reasons.get(r["status"], 0) + 1
        for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
            print(f"  {k}: {v}")
    print(f"wrote           : {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
