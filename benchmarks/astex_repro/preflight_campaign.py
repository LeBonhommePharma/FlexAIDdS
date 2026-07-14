#!/usr/bin/env python3
"""Campaign preflight (Chunk 3): site QC, chain catalog, binary hygiene."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path


def sdf_centroid(path: Path):
    lines = path.read_text().splitlines()
    coords = []
    for i, line in enumerate(lines):
        if "V2000" in line:
            n = int(line[:3])
            for j in range(i + 1, i + 1 + n):
                p = lines[j].split()
                coords.append(tuple(map(float, p[:3])))
            break
    if not coords:
        return None
    n = len(coords)
    return tuple(sum(c[i] for c in coords) / n for i in range(3))


def pdb_centroid(path: Path):
    xs = []
    for line in path.open(errors="ignore"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                xs.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            except ValueError:
                pass
    if not xs:
        return None
    n = len(xs)
    return tuple(sum(a[i] for a in xs) / n for i in range(3))


def chains_atom(path: Path) -> set[str]:
    return {line[21] for line in path.open(errors="ignore") if line.startswith("ATOM")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--engine-dir", type=Path, required=True)
    ap.add_argument("--codes-file", type=Path, required=True)
    ap.add_argument("--site-max-delta", type=float, default=8.0)
    args = ap.parse_args()

    errors = []
    warnings = []

    flex = args.engine_dir / "FlexAIDdS"
    runner = args.engine_dir / "benchmark_datasets"
    for b in (flex, runner):
        if not b.is_file():
            errors.append(f"missing binary {b}")
        else:
            # PGO instrumented?
            try:
                out = subprocess.check_output(["nm", str(b)], stderr=subprocess.DEVNULL, text=True)
                if "llvm_profile" in out.lower() or "__llvm_prf" in out.lower():
                    errors.append(f"PGO symbols in {b.name} — refuse Stage-1 instrumented binary")
            except (subprocess.CalledProcessError, FileNotFoundError):
                pass

    # CF-primary / CF-stagnation strings
    try:
        s = subprocess.check_output(["strings", str(runner)], text=True, errors="ignore")
        if "CF-PRIMARY" not in s and "cf_primary" not in s.lower():
            # may still be OK if only compile-time env
            warnings.append("runner strings lack CF-PRIMARY log tag (rebuild?)")
        if "pop_effective" not in s:
            warnings.append("runner may lack pop-scaling")
    except Exception as e:
        warnings.append(f"strings check failed: {e}")

    codes = [c.strip() for c in args.codes_file.read_text().splitlines() if c.strip()]
    cache = args.cache / "astex_diverse" if (args.cache / "astex_diverse").is_dir() else args.cache

    # catalog
    cat_path = Path(__file__).resolve().parent / "chain_catalog.yaml"
    catalog = {}
    if cat_path.exists():
        try:
            import yaml  # optional
            catalog = yaml.safe_load(cat_path.read_text()).get("targets", {}) or {}
        except Exception:
            # minimal parse
            cur = None
            for line in cat_path.read_text().splitlines():
                if line.strip().startswith("1G9V") or line.strip().startswith("2BYS"):
                    cur = line.strip().rstrip(":")
                if "keep_chains" in line and cur:
                    # keep_chains: ["A", "B", "C"]
                    import re
                    catalog[cur] = {"keep_chains": re.findall(r'"([A-Z])"', line)}

    for code in codes:
        d = cache / code
        if not d.is_dir():
            warnings.append(f"{code}: missing cache dir")
            continue
        apo = d / f"{code}_apo.pdb"
        lig = d / f"{code}_ligand.sdf"
        site = d / f"{code}_binding_site.pdb"
        if apo.exists() and code in catalog:
            want = set(catalog[code].get("keep_chains", []))
            have = chains_atom(apo)
            if want and not want.issubset(have):
                errors.append(f"{code}: apo missing catalog chains {want - have} (have {have})")
            if want and have - want:
                warnings.append(f"{code}: apo has extra chains {have - want} (catalog {want})")
        if lig.exists() and site.exists():
            lc, sc = sdf_centroid(lig), pdb_centroid(site)
            if lc and sc:
                delta = math.dist(lc, sc)
                if delta > args.site_max_delta:
                    # ligand-centered repair is OK if present
                    lc_site = d / f"{code}_ligand_centered_site.pdb"
                    if lc_site.exists():
                        sc2 = pdb_centroid(lc_site)
                        if sc2 and math.dist(lc, sc2) <= args.site_max_delta:
                            warnings.append(
                                f"{code}: binding_site {delta:.1f}A from ligand; "
                                f"ligand_centered_site OK"
                            )
                        else:
                            errors.append(
                                f"{code}: site {delta:.1f}A from ligand and no good repair"
                            )
                    else:
                        errors.append(f"{code}: site {delta:.1f}A from ligand (>{args.site_max_delta})")

    for w in warnings:
        print(f"WARN: {w}")
    for e in errors:
        print(f"ERROR: {e}")
    if errors:
        print(f"PREFLIGHT FAIL ({len(errors)} errors)")
        return 1
    print(f"PREFLIGHT OK ({len(codes)} codes, {len(warnings)} warnings)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
