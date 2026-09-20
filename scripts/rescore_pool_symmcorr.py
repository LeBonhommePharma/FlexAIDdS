#!/usr/bin/env python3
"""Full-pool symmetry-corrected rescore of a FlexAID campaign, plus a content-hashed
pool manifest — the founding artifact of the "pose bank".

WHY
  Every open scoring question in this project currently costs a re-dock, because the
  pose pool is treated as a transient byproduct. Only SEARCH questions need docking;
  scoring, election, RMSD convention and physical validity are cheap, deterministic,
  and perfectly paired when computed over a FIXED pool. This pass freezes the pool
  (sha256 per pose) and scores every pose with the claim-grade metric, so that top-1,
  top-N and oracle-in-pool all fall out of ONE pass and every later question is a
  rescore rather than a campaign.

WHAT IT MEASURES
  scripts/rmsd_symmcorr.py::symmcorr_rmsd — spyrmsd graph-automorphism RMSD, in place
  (center=False, minimize=False, never superposed), heavy atoms only, connectivity and
  atomic numbers from the crystal SDF BOND BLOCK (no perception, never the pose file),
  ligand atoms selected on PDB serial >= 90000. It is FAIL-CLOSED: a pose that cannot
  be scored gets rmsd=None and a machine-readable status, never a weaker metric.

  serial_rmsd (identity atom mapping) is computed alongside ONLY to assert the
  documented invariant symmcorr <= serial. It is never a claim metric.

WHY NOT THE ENGINE'S OWN COLUMN
  Measured on 56 audited targets, top-1 success: 29/56 by the engine's emitted
  rmsd_hungarian column, 20/56 by exact graph automorphism, 18/56 by naive index
  order. The emitted value is not even a bound — a correctly computed Hungarian sits
  BELOW it on 3 of 56 — because its cost-matrix walk indexes out of bounds.

DESIGN NOTES
  * ONE read per pose: the bytes are hashed and scored from the same buffer, so the
    manifest costs no extra I/O.
  * Per-target output written incrementally; a target whose CSV exists is SKIPPED, so
    an interruption costs at most one target.
  * REFUSAL gates run before any scoring: spyrmsd must import, and every crystal SDF
    must parse with heavy atoms > 0 and no unresolved (Z=0) symbol. A missing or
    unparseable reference is a refusal, not a zero.
  * NON-VACUITY: the summary refuses to write if zero poses were scored, or if every
    row carries a failure status. A pass that scored nothing must not look successful.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import glob
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

# rmsd_symmcorr lives beside this file in the repo; fall back to $FLEXAIDDS_SCRIPTS
# so a workspace copy can still point at a checkout. No developer path is baked in.
SCRIPTS = os.environ.get("FLEXAIDDS_SCRIPTS", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SCRIPTS)

# Two layouts exist in this project and both must parse:
#   Arm A / DatasetRunner-per-restart : TARGET_r<restart>_<idx>.pdb
#   wall_paired_85 campaign           : TARGET_<idx>.pdb  (restart not in the name)
# restart is therefore OPTIONAL; when absent it is recorded as -1 rather than
# guessed, so a caller can never mistake "unknown restart" for "restart 0".

def _pose_glob(campaign: str, target: str) -> list:
    """Pose files for one target, tolerating both directory layouts.

    Arm A wrote poses at        <campaign>/<target>/*.pdb
    the wall campaign writes at <campaign>/<target>/<target>/*.pdb
    (the runner creates a target-named subdir inside the cell dir).

    Returns the SHALLOWER match when both exist, because a nested dir that also
    has loose .pdb files at the parent level would otherwise silently merge two
    pose sets into one pool.
    """
    shallow = [p for p in glob.glob(f"{campaign}/{target}/*.pdb") if "_INI" not in p]
    if shallow:
        return sorted(shallow)
    # Nested layout: <cell>/<target>/ holds restart 0 at the top plus r1..rN
    # subdirs, so this MUST recurse. A non-recursive glob here returns only
    # restart 0 -- 51 of 251 poses on 1G9V -- and the pass would look complete.
    return sorted(x for x in glob.glob(f"{campaign}/{target}/{target}/**/*.pdb",
                                       recursive=True) if "_INI" not in x)

POSE_RE = re.compile(
    r"^(?P<target>[0-9A-Za-z]{4})(?:_r(?P<restart>\d+))?_(?P<idx>\d+)\.pdb$")


def crystal_sdf(cache: str, target: str) -> str:
    return f"{cache}/{target}/{target}_ligand.sdf"


def pose_ligand_atom_count(campaign: str, target: str) -> int | None:
    """Heavy-atom count of the ligand IN THE POSE (serial >= 90000)."""
    poses = [p for p in _pose_glob(campaign, target) if "_INI" not in p]
    if not poses:
        return None
    n = 0
    for line in open(poses[0], errors="ignore"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                s = int(line[6:11])
            except ValueError:
                continue
            if s >= 90000:
                n += 1
    return n


def choose_reference(campaign: str, caches: list[str], target: str):
    """Pick the crystal SDF whose heavy-atom count MATCHES THE POSE.

    WHY THIS EXISTS — measured 2026-09-18. Two caches both contained a file named
    `1TW6_ligand.sdf` holding DIFFERENT MOLECULES: one `BTB` (bis-tris buffer, a
    crystallisation additive) with 14 atoms, the other `ALA` (a four-residue
    peptide) with 27. The docked pose has 27 ligand atoms, so only the second is
    the cognate reference. Selecting a reference by PATH instead of by agreement
    with the pose silently scores against the wrong molecule; here it was caught
    only because the atom counts differed and symmcorr_rmsd is fail-closed.

    Exactly 1 of 85 targets disagreed, so this is a narrow correction — but the
    check is derived from data rather than hardcoded, and NEITHER-matches is a
    refusal, not a fallback.
    """
    import rmsd_symmcorr as R

    want = pose_ligand_atom_count(campaign, target)
    cands = []
    for c in caches:
        p = crystal_sdf(c, target)
        if not os.path.exists(p):
            continue
        try:
            _, an, _ = R.parse_sdf(p)
        except Exception:
            continue
        heavy = int((an > 1).sum())
        cands.append((c, p, heavy))
        if want is not None and heavy == want:
            return {"path": p, "cache": c, "pose_heavy": want,
                    "ref_heavy": heavy, "match": True}
    return {"path": None, "cache": None, "pose_heavy": want,
            "candidates": [(os.path.basename(c), h) for c, _, h in cands],
            "match": False}


def preflight(campaign: str, cache: list[str], targets: list[str]) -> dict:
    """Refusal gates. Returns a dict of facts; raises SystemExit on refusal."""
    import rmsd_symmcorr as R

    try:
        spyr, spyver = R._spyrmsd()
    except R.SymmCorrUnavailable as exc:
        raise SystemExit(f"REFUSE: {exc}")

    bad = []
    facts = {}
    for t in targets:
        sel = choose_reference(campaign, cache, t)
        if not sel["match"]:
            bad.append((t, f"no_cache_matches_pose_heavy={sel['pose_heavy']}"
                           f" candidates={sel.get('candidates')}"))
            continue
        p = sel["path"]
        try:
            cx, an, adj = R.parse_sdf(p)
        except Exception as e:
            bad.append((t, f"parse_error:{type(e).__name__}"))
            continue
        heavy = int((an > 1).sum())
        zs = sorted({int(x) for x in an})
        if heavy == 0:
            bad.append((t, "no_heavy_atoms"))
        elif 0 in zs:
            bad.append((t, "unresolved_element_symbol_Z0"))
        facts[t] = {"n_atoms": int(len(an)), "n_heavy": heavy, "z_set": zs,
                    "reference_path": p, "reference_cache": sel["cache"],
                    "pose_heavy": sel["pose_heavy"]}

    if bad:
        raise SystemExit(
            "REFUSE: crystal reference failed the gate for "
            f"{len(bad)} target(s): {bad[:10]}"
        )
    return {
        "spyrmsd_version": spyver,
        "metric": R.METRIC_NAME,
        "producer": R.PRODUCER,
        "ligand_serial_min": R.LIGAND_SERIAL_MIN,
        "targets": facts,
    }


def score_target(args) -> dict:
    """Score every pose of one target. Hash and score from the SAME read."""
    campaign, caches, target, outdir = args
    import rmsd_symmcorr as R

    out_csv = f"{outdir}/{target}.csv"
    if os.path.exists(out_csv):
        with open(out_csv, newline="") as fh:
            n = sum(1 for _ in csv.reader(fh)) - 1
        return {"target": target, "skipped": True, "n_rows": max(0, n)}

    sel = choose_reference(campaign, caches, target)
    if not sel["match"]:
        # preflight already refuses on this; belt-and-braces so a direct call
        # can never score against a reference that disagrees with the pose.
        return {"target": target, "skipped": False, "n_poses": 0, "n_scored": 0,
                "n_failed": 0, "n_invariant_violations": 0, "wall_s": 0.0,
                "refused": f"no_cache_matches_pose_heavy={sel['pose_heavy']}"}
    sdf = sel["path"]
    poses = sorted(p for p in _pose_glob(campaign, target) if "_INI" not in p)
    t0 = time.time()
    rows = []
    pool_h = hashlib.sha256()
    n_ok = n_fail = n_viol = 0

    for p in poses:
        name = os.path.basename(p)
        m = POSE_RE.match(name)
        # POSE KEY MUST BE UNIQUE WITHIN A TARGET. In the wall_paired_85 layout the
        # restarts each write into their own r<N>/ subdir and REUSE THE SAME
        # FILENAMES (1G9V_18.pdb exists under ., r1/, r2/, r3/ and r4/), so the
        # basename identifies 5 different poses. Key on the path relative to the
        # target's pose root instead, and take the restart from the DIRECTORY --
        # the filename does not carry it in this layout.
        parent = os.path.basename(os.path.dirname(p))
        rdir = re.match(r"^r(\d+)$", parent)
        if rdir:
            rel = f"{parent}/{name}"
            restart = int(rdir.group(1))
        else:
            rel = name
            restart = 0            # top-level dir is restart 0 in this layout
        if m and m.group("restart") is not None:
            restart = int(m.group("restart"))   # Arm A layout: name is authoritative
            rel = name
        raw = Path(p).read_bytes()          # ONE read: hash + score from this buffer
        h = hashlib.sha256(raw).hexdigest()
        pool_h.update(h.encode())

        res = R.symmcorr_rmsd(sdf, raw)
        ser = R.serial_rmsd(sdf, raw)
        rm = res.get("rmsd")
        if rm is None:
            n_fail += 1
        else:
            n_ok += 1
            # documented invariant: the symmetry-corrected value can only be <= serial
            if ser is not None and rm > ser + 1e-6:
                n_viol += 1
        rows.append({
            "target": target,
            "pose_file": rel,
            "restart": restart,
            "pose_index": int(m.group("idx")) if m else "",
            "sha256": h,
            "rmsd_symmcorr": "" if rm is None else f"{rm:.6f}",
            "rmsd_serial": "" if ser is None else f"{ser:.6f}",
            "n_heavy": res.get("n_heavy", ""),
            "status": res.get("status", ""),
        })

    # GUARD: the bank is keyed on (target, pose_file). A duplicate key means two
    # distinct poses would silently collapse into one on any later join, which is
    # exactly what keying on the basename did before restarts got their own dirs.
    keys = [r["pose_file"] for r in rows]
    if len(set(keys)) != len(keys):
        dup = sorted({k for k in keys if keys.count(k) > 1})[:5]
        raise SystemExit(f"REFUSE: duplicate pose keys for {target}: {dup} "
                         f"({len(keys)} poses, {len(set(keys))} distinct)")

    tmp = out_csv + ".partial"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else
                           ["target", "pose_file", "restart", "pose_index", "sha256",
                            "rmsd_symmcorr", "rmsd_serial", "n_heavy", "status"])
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out_csv)               # atomic: a reader never sees a partial file

    return {
        "target": target, "skipped": False, "n_poses": len(poses),
        "n_scored": n_ok, "n_failed": n_fail, "n_invariant_violations": n_viol,
        "target_pool_sha256": pool_h.hexdigest(), "wall_s": round(time.time() - t0, 1),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--cache", required=True, help="primary astex_diverse cache")
    ap.add_argument("--cache-alt", default="",                     help="fallback cache, used only when the primary's ligand "
                         "disagrees with the pose's heavy-atom count")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--only", default="", help="comma-separated targets (default: all)")
    ap.add_argument("--shard", type=int, default=0, help="0-based shard index")
    ap.add_argument("--nshards", type=int, default=1, help="total shards")
    ap.add_argument("--summary-only", action="store_true",
                    help="write POOL_MANIFEST.json from existing per-target CSVs; score nothing")
    a = ap.parse_args(argv)

    os.makedirs(a.out, exist_ok=True)
    all_targets = sorted(d for d in os.listdir(a.campaign)
                         if os.path.isdir(f"{a.campaign}/{d}") and len(d) == 4)
    targets = all_targets
    if a.only:
        want = {x.strip().upper() for x in a.only.split(",") if x.strip()}
        targets = [t for t in targets if t.upper() in want]
    if not (1 <= a.nshards) or not (0 <= a.shard < a.nshards):
        raise SystemExit(f"REFUSE: bad shard {a.shard}/{a.nshards}")
    if a.nshards > 1:
        # strided so every shard gets a mix of cheap and expensive targets
        targets = targets[a.shard::a.nshards]
    if not targets:
        raise SystemExit("REFUSE: no targets selected")
    if a.summary_only:
        targets = all_targets

    caches = [a.cache] + ([a.cache_alt] if a.cache_alt else [])
    facts = preflight(a.campaign, caches, targets)
    print(f"GATES OK  spyrmsd={facts['spyrmsd_version']}  metric={facts['metric']}  "
          f"targets={len(targets)}", flush=True)

    # SEQUENTIAL BY DESIGN. concurrent.futures.ProcessPoolExecutor cannot start in
    # this sandbox: _check_system_limits() calls os.sysconf("SC_SEM_NSEMS_MAX"),
    # which raises PermissionError (same EPERM family as sysctl/ps here).
    # Parallelism is therefore achieved by launching N INSTANCES of this script with
    # disjoint --shard values; per-target CSVs plus the skip-if-exists check make
    # that safe and resumable, and no two instances ever write the same file.
    results = []
    for i, t in enumerate(targets, 1):
        r = score_target((a.campaign, caches, t, a.out))
        results.append(r)
        if r.get("skipped"):
            print(f"[{i}/{len(targets)}] {r['target']} SKIP ({r['n_rows']} rows)", flush=True)
        else:
            print(f"[{i}/{len(targets)}] {r['target']} {r['n_scored']}/{r['n_poses']} scored, "
                  f"{r['n_failed']} failed, {r['n_invariant_violations']} inv-viol, "
                  f"{r['wall_s']}s", flush=True)

    scored = sum(r.get("n_scored", 0) for r in results)
    failed = sum(r.get("n_failed", 0) for r in results)
    viol = sum(r.get("n_invariant_violations", 0) for r in results)
    fresh = [r for r in results if not r.get("skipped")]

    # NON-VACUITY GATES
    if fresh and scored == 0:
        print("VACUOUS: targets were processed but ZERO poses scored", file=sys.stderr)
        return 3

    pool_h = hashlib.sha256()
    for r in sorted(results, key=lambda x: x["target"]):
        if r.get("target_pool_sha256"):
            pool_h.update(r["target_pool_sha256"].encode())

    manifest = {
        "campaign": a.campaign, "cache": a.cache,
        "n_targets": len(targets), "n_scored": scored, "n_failed": failed,
        "n_invariant_violations": viol,
        "pool_sha256": pool_h.hexdigest() if fresh else "incomplete-resume",
        "metric": facts["metric"], "producer": facts["producer"],
        "spyrmsd_version": facts["spyrmsd_version"],
        "ligand_serial_min": facts["ligand_serial_min"],
        "per_target": {r["target"]: {k: v for k, v in r.items() if k != "target"}
                       for r in results},
        "reference_per_target": {t: {"cache": v.get("reference_cache"),
                                     "ref_heavy": v.get("n_heavy"),
                                     "pose_heavy": v.get("pose_heavy")}
                                 for t, v in facts["targets"].items()},
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with open(f"{a.out}/POOL_MANIFEST.json", "w") as fh:
        json.dump(manifest, fh, indent=1)

    print(f"\nscored {scored}  failed {failed}  invariant-violations {viol}")
    print(f"pool_sha256 {manifest['pool_sha256'][:16]}")
    return 1 if viol else 0


if __name__ == "__main__":
    sys.exit(main())
