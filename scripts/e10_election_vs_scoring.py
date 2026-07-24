#!/usr/bin/env python3
"""E10 offline diagnostic: election vs scoring / sampling split on frozen heads.

No GA re-dock. Reads existing result.csv + ranked pose PDBs (REMARK CF) under
campaign target directories and reports:

  - election_gap: seed_echo=0, BCR ≤ 2.5 Å, elected RMSD > 2.0 Å
  - among heads: does any near-native (RMSD ≤ 2.5) exist with CF better than elected?
  - soft-β size bias proxy: REMARK soft_beta_G vs REMARK CF for rank-0

Optional: --probe-cf path (not required). smina not required.

Usage:
  python3 scripts/e10_election_vs_scoring.py \\
    --target-dir ~/flexaidds_results/v_autonomous_*/1G9V \\
    --out-json workorders/e10_1g9v.json

  python3 scripts/e10_election_vs_scoring.py \\
    --campaign-dir ~/flexaidds_results/v_autonomous_20260724_160919 \\
    --out-csv /tmp/e10.csv --out-md /tmp/e10.md

Exit: 0 always when analysis written; 2 on usage error.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional


RMSD_OK = 2.0
BCR_NEAR = 2.5


def _f(x: Any) -> float:
    if x is None or x == "" or str(x).upper() == "NA":
        return float("nan")
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def parse_pdb_remarks(path: Path) -> dict[str, float]:
    """Extract REMARK CF / soft_beta_G / frequency from a FlexAID pose PDB."""
    out: dict[str, float] = {}
    if not path.is_file():
        return out
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return out
    for line in text.splitlines():
        if not line.startswith("REMARK"):
            continue
        m = re.search(r"REMARK\s+CF=(-?[\d.]+)", line)
        if m and "CF." not in line and "CF=" in line:
            # First bare CF= wins as pose CF
            if "cf" not in out:
                out["cf"] = float(m.group(1))
        m = re.search(r"REMARK\s+soft_beta_G\s*=\s*(-?[\d.]+)", line)
        if m:
            out["soft_beta_G"] = float(m.group(1))
        m = re.search(r"REMARK\s+frequency\s*=\s*([\d.]+)", line)
        if m:
            out["frequency"] = float(m.group(1))
        m = re.search(r"REMARK\s+CF\.com=(-?[\d.]+)", line)
        if m:
            out["cf_com"] = float(m.group(1))
        m = re.search(r"Average CF:(-?[\d.]+)\s+Frequency:(\d+)", line)
        if m:
            out["cluster_avg_cf"] = float(m.group(1))
            out["cluster_freq"] = float(m.group(2))
    return out


def list_ranked_poses(target_dir: Path, pdb: str) -> list[tuple[int, Path]]:
    """Return (rank_index, path) for {pdb}_{i}.pdb excluding _INI."""
    found: list[tuple[int, Path]] = []
    for p in target_dir.glob(f"{pdb}_*.pdb"):
        name = p.name
        if "_INI" in name or name.endswith("_ini.pdb"):
            continue
        m = re.match(rf"{re.escape(pdb)}_(\d+)\.pdb$", name)
        if not m:
            continue
        found.append((int(m.group(1)), p))
    found.sort(key=lambda x: x[0])
    return found


def load_result_row(target_dir: Path) -> dict[str, str]:
    for cand in (target_dir / "result.csv",):
        if cand.is_file():
            with cand.open(newline="", encoding="utf-8", errors="replace") as fh:
                rows = list(csv.DictReader(fh))
                if rows:
                    return rows[0]
    # nested r*/
    for sub in sorted(target_dir.iterdir()) if target_dir.is_dir() else []:
        if sub.is_dir() and (sub / "result.csv").is_file():
            with (sub / "result.csv").open(newline="", encoding="utf-8", errors="replace") as fh:
                rows = list(csv.DictReader(fh))
                if rows:
                    return rows[0]
    return {}


@dataclass
class TargetE10:
    pdb: str
    path: str
    seed_echo: str = ""
    rmsd_elected: float = float("nan")
    bcr: float = float("nan")
    elected_cf: float = float("nan")
    elected_soft_beta_G: float = float("nan")
    elected_freq: float = float("nan")
    n_heads: int = 0
    best_head_rmsd_proxy: float = float("nan")  # min rank index among heads — use BCR from csv
    best_cf_among_heads: float = float("nan")
    n_heads_cf_better_than_elected: int = 0
    election_gap: bool = False  # BCR near-native, elected fails
    size_bias_suspect: bool = False  # |soft_beta_G - cf| large with high freq
    notes: str = ""
    heads: list[dict[str, Any]] = field(default_factory=list)


def analyze_target(target_dir: Path, max_heads: int = 40) -> TargetE10:
    pdb = target_dir.name
    row = load_result_row(target_dir)
    if row.get("pdb_id"):
        pdb = str(row["pdb_id"]).strip() or pdb

    rmsd = _f(row.get("rmsd_hungarian"))
    if not math.isfinite(rmsd):
        rmsd = _f(row.get("rmsd_to_crystal"))
    bcr = _f(row.get("best_cluster_rmsd"))
    elect_cf = _f(row.get("elected_cf"))
    if not math.isfinite(elect_cf):
        elect_cf = _f(row.get("best_score"))
    se = str(row.get("seed_echo", "")).strip()

    ranked = list_ranked_poses(target_dir, pdb)
    heads_out: list[dict[str, Any]] = []
    best_cf = float("nan")
    n_better = 0

    elected_path = target_dir / f"{pdb}_0.pdb"
    if not elected_path.is_file() and row.get("elected_pose_path"):
        ep = Path(str(row["elected_pose_path"]))
        if ep.is_file():
            elected_path = ep
    elected_rem = parse_pdb_remarks(elected_path)
    if not math.isfinite(elect_cf) and "cf" in elected_rem:
        elect_cf = elected_rem["cf"]
    soft_g = elected_rem.get("soft_beta_G", float("nan"))
    freq = elected_rem.get("frequency", elected_rem.get("cluster_freq", float("nan")))

    for rank, path in ranked[:max_heads]:
        rem = parse_pdb_remarks(path)
        cf = rem.get("cf", float("nan"))
        heads_out.append(
            {
                "rank": rank,
                "path": str(path),
                "cf": cf,
                "soft_beta_G": rem.get("soft_beta_G", float("nan")),
                "frequency": rem.get("frequency", rem.get("cluster_freq", float("nan"))),
                "cf_com": rem.get("cf_com", float("nan")),
            }
        )
        if math.isfinite(cf):
            if not math.isfinite(best_cf) or cf < best_cf:
                best_cf = cf
            if math.isfinite(elect_cf) and cf < elect_cf - 1e-6:
                n_better += 1

    election_gap = (
        se in ("0", "0.0", "")
        and math.isfinite(bcr)
        and bcr <= BCR_NEAR
        and math.isfinite(rmsd)
        and rmsd > RMSD_OK
    )
    # Size bias: soft_beta_G much lower than pose CF with huge frequency
    size_bias = False
    if math.isfinite(soft_g) and math.isfinite(elect_cf) and math.isfinite(freq):
        if freq >= 10 and soft_g < elect_cf - 50.0:
            size_bias = True

    notes = []
    if election_gap:
        notes.append(
            f"election_gap: BCR={bcr:.2f}Å near-native pool but elected rmsd={rmsd:.2f}Å"
        )
    if size_bias:
        notes.append(
            f"size_bias_suspect: soft_beta_G={soft_g:.1f} vs CF={elect_cf:.1f} freq={freq:.0f}"
        )
    if n_better > 0:
        notes.append(f"{n_better} heads have better (lower) CF than elected")

    return TargetE10(
        pdb=pdb,
        path=str(target_dir),
        seed_echo=se,
        rmsd_elected=rmsd,
        bcr=bcr,
        elected_cf=elect_cf,
        elected_soft_beta_G=soft_g,
        elected_freq=freq,
        n_heads=len(heads_out),
        best_cf_among_heads=best_cf,
        n_heads_cf_better_than_elected=n_better,
        election_gap=election_gap,
        size_bias_suspect=size_bias,
        notes="; ".join(notes),
        heads=heads_out,
    )


def discover_targets(campaign_dir: Path) -> list[Path]:
    out: list[Path] = []
    if not campaign_dir.is_dir():
        return out
    for child in sorted(campaign_dir.iterdir()):
        if not child.is_dir() or child.name.startswith("."):
            continue
        # has result or ranked poses
        if (child / "result.csv").is_file() or list(child.glob("*_0.pdb")):
            out.append(child)
    return out


def render_md(results: list[TargetE10], *, title: str) -> str:
    lines = [
        f"# {title}",
        "",
        "E10 offline diagnostic (no re-dock). Metrics from `result.csv` + REMARK CF on heads.",
        "",
        f"Targets analyzed: **{len(results)}**",
        f"Election-gap (BCR≤{BCR_NEAR}, elected>{RMSD_OK}, seed_echo=0): "
        f"**{sum(1 for r in results if r.election_gap)}**",
        f"Size-bias suspects (soft_beta_G ≪ CF with high freq): "
        f"**{sum(1 for r in results if r.size_bias_suspect)}**",
        "",
        "| PDB | rmsd_h | BCR | elect_CF | soft_β_G | freq | gap? | size_bias? | notes |",
        "|-----|--------|-----|----------|----------|------|------|------------|-------|",
    ]
    for r in results:
        lines.append(
            f"| {r.pdb} | {_fmt(r.rmsd_elected)} | {_fmt(r.bcr)} | {_fmt(r.elected_cf)} | "
            f"{_fmt(r.elected_soft_beta_G)} | {_fmt(r.elected_freq, 0)} | "
            f"{'Y' if r.election_gap else 'n'} | {'Y' if r.size_bias_suspect else 'n'} | "
            f"{r.notes[:80]} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- **election_gap**: sampling found a near-native basin (BCR) but rank-0 election failed — "
            "primary target for **E1b ACF_STRICT** / election fixes.",
            "- **size_bias_suspect**: REMARK soft_beta_G far below pose CF with large cluster frequency "
            "— consistent with ACF ≈ Emin − T ln Z multiplicity inflation.",
            "- Incomplete / UNCITABLE campaigns remain mechanism evidence only.",
            "",
        ]
    )
    return "\n".join(lines)


def _fmt(x: float, nd: int = 2) -> str:
    if not math.isfinite(x):
        return "—"
    if nd == 0:
        return f"{x:.0f}"
    return f"{x:.{nd}f}"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-dir", type=Path, action="append", default=None)
    ap.add_argument("--campaign-dir", type=Path, default=None)
    ap.add_argument("--max-heads", type=int, default=40)
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--out-csv", type=Path, default=None)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)

    dirs: list[Path] = []
    if args.target_dir:
        dirs.extend(args.target_dir)
    if args.campaign_dir:
        dirs.extend(discover_targets(args.campaign_dir))
    if not dirs:
        print("error: pass --target-dir and/or --campaign-dir", file=sys.stderr)
        return 2

    results = [analyze_target(d.expanduser(), max_heads=args.max_heads) for d in dirs]
    payload = {
        "n_targets": len(results),
        "n_election_gap": sum(1 for r in results if r.election_gap),
        "n_size_bias_suspect": sum(1 for r in results if r.size_bias_suspect),
        "targets": [asdict(r) for r in results],
    }

    if args.out_json:
        args.out_json.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.out_json.expanduser().write_text(
            json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8"
        )
        print(f"Wrote {args.out_json}")
    if args.out_csv:
        args.out_csv.expanduser().parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "pdb",
            "rmsd_elected",
            "bcr",
            "elected_cf",
            "elected_soft_beta_G",
            "elected_freq",
            "election_gap",
            "size_bias_suspect",
            "n_heads",
            "n_heads_cf_better_than_elected",
            "notes",
        ]
        with args.out_csv.expanduser().open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in results:
                w.writerow(asdict(r))
        print(f"Wrote {args.out_csv}")
    if args.out_md:
        args.out_md.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.out_md.expanduser().write_text(
            render_md(results, title="E10 election vs scoring diagnostic"),
            encoding="utf-8",
        )
        print(f"Wrote {args.out_md}")

    print(
        f"E10: n={len(results)} election_gap={payload['n_election_gap']} "
        f"size_bias_suspect={payload['n_size_bias_suspect']}"
    )
    for r in results:
        if r.election_gap or r.size_bias_suspect:
            print(f"  {r.pdb}: rmsd={r.rmsd_elected:.2f} bcr={r.bcr:.2f} {r.notes}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
