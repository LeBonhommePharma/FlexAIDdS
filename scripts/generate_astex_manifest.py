#!/usr/bin/env python3
"""Regenerate Astex Diverse checksum manifests for the canonical tree.

Canonical root: benchmarks/astex_diverse/astex_diverse/<PDB>/

Writes:
  benchmarks/datasets/astex_diverse_sha256.csv
  benchmarks/datasets/astex_diverse_manifest.json

Usage (from repo root):
  python3 scripts/generate_astex_manifest.py
  python3 scripts/generate_astex_manifest.py --check   # verify only; exit 1 on drift
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CANONICAL = REPO / "benchmarks" / "astex_diverse" / "astex_diverse"
OUT_CSV = REPO / "benchmarks" / "datasets" / "astex_diverse_sha256.csv"
OUT_JSON = REPO / "benchmarks" / "datasets" / "astex_diverse_manifest.json"

KEY_SUFFIXES = (
    "{id}.pdb",
    "{id}.cif",
    "{id}_apo.pdb",
    "{id}_ligand.sdf",
    "{id}_binding_site.pdb",
    "{id}_ligand_centered_site.pdb",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect() -> tuple[list[str], list[dict], int]:
    if not CANONICAL.is_dir():
        raise SystemExit(f"canonical tree missing: {CANONICAL}")
    targets = sorted(p.name for p in CANONICAL.iterdir() if p.is_dir())
    entries: list[dict] = []
    n_files = 0
    for tid in targets:
        tdir = CANONICAL / tid
        files: dict[str, dict] = {}
        for tmpl in KEY_SUFFIXES:
            name = tmpl.format(id=tid)
            fp = tdir / name
            if not fp.exists():
                continue
            target = fp.resolve() if fp.is_symlink() else fp
            if not target.is_file():
                continue
            files[name] = {
                "sha256": sha256_file(target),
                "bytes": target.stat().st_size,
            }
            if fp.is_symlink():
                files[name]["symlink"] = True
            n_files += 1
        entries.append({"pdb_id": tid, "files": files})
    return targets, entries, n_files


def render_csv(entries: list[dict]) -> str:
    rows = ["pdb_id,relative_path,bytes,sha256"]
    for e in entries:
        for name, meta in sorted(e["files"].items()):
            rel = f"benchmarks/astex_diverse/astex_diverse/{e['pdb_id']}/{name}"
            rows.append(f"{e['pdb_id']},{rel},{meta['bytes']},{meta['sha256']}")
    return "\n".join(rows) + "\n"


def render_json(targets: list[str], n_files: int) -> str:
    summary = {
        "dataset": "astex_diverse",
        "canonical_root": "benchmarks/astex_diverse/astex_diverse",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        ),
        "n_targets": len(targets),
        "n_hashed_files": n_files,
        "key_suffixes": list(KEY_SUFFIXES),
        "reference": {
            "hartshorn_2007": "https://doi.org/10.1021/jm061277y",
            "yaml": "benchmarks/datasets/astex_diverse.yaml",
            "canonical_doc": "benchmarks/datasets/CANONICAL.md",
            "reproducibility": "REPRODUCIBILITY.md",
        },
        "targets": targets,
    }
    return json.dumps(summary, indent=2) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify committed CSV matches a fresh hash of the canonical tree",
    )
    args = ap.parse_args()

    targets, entries, n_files = collect()
    csv_text = render_csv(entries)
    json_text = render_json(targets, n_files)

    if args.check:
        if not OUT_CSV.is_file():
            print(f"FAIL: missing {OUT_CSV}", file=sys.stderr)
            return 1
        committed = OUT_CSV.read_text()
        if committed != csv_text:
            print(
                f"FAIL: {OUT_CSV} drift vs canonical tree "
                f"({n_files} files / {len(targets)} targets hashed)",
                file=sys.stderr,
            )
            c_lines = committed.splitlines()
            n_lines = csv_text.splitlines()
            drift = 0
            for a, b in zip(c_lines, n_lines):
                if a != b:
                    print(f"  - {a[:120]}", file=sys.stderr)
                    print(f"  + {b[:120]}", file=sys.stderr)
                    drift += 1
                    if drift >= 5:
                        break
            if len(c_lines) != len(n_lines):
                print(
                    f"  line count committed={len(c_lines)} fresh={len(n_lines)}",
                    file=sys.stderr,
                )
            return 1
        print(
            f"OK: {OUT_CSV} matches canonical tree "
            f"({n_files} files, {len(targets)} targets)"
        )
        return 0

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    OUT_CSV.write_text(csv_text)
    OUT_JSON.write_text(json_text)
    print(f"wrote {OUT_CSV.relative_to(REPO)} ({n_files} files)")
    print(f"wrote {OUT_JSON.relative_to(REPO)} ({len(targets)} targets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
