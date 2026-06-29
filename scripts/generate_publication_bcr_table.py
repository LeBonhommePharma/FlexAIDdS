#!/usr/bin/env python3
"""Generate publication-ready Astex Diverse BCR comparison tables (CSV/JSON/LaTeX).

Search-only BCR excludes successes where seed_echo=1 (crystal INI copies).
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = Path.home() / "Documents/PhD/Programs/FlexAIDdS/results"
OUT_DIR = REPO / "benchmarks/publication"
N_DENOM = 85
SUCCESS_RMSD = 2.0

FLEXAIDDS_RUNS: list[dict[str, Any]] = [
    {
        "run_id": "v109",
        "label": "v109 tier1 consensus5r",
        "results_subdir": "v109_20260626_tier1_consensus5r",
        "protocol_summary": "oracle-ceiling",
        "publishable_tier": "non-standard (oracle-ceiling)",
        "notes": "Internal record; 46/80 passes are seed_echo ini_elitism",
        "vct_r0": None,
        "consensus_scorer": "1",
        "native_seed_frac": "0.0",
        "git_commit": "15b536f8",
        "literature_comparable": "no",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "v127",
        "label": "v127 optB full85",
        "results_subdir": "v127_20260629_0139_optB_full85",
        "protocol_summary": "oracle-ceiling + Fix-B logsumexp",
        "publishable_tier": "non-standard (oracle-ceiling)",
        "notes": "Fix-B selector; native_seed_frac=0.9",
        "vct_r0": "4",
        "consensus_scorer": "1",
        "native_seed_frac": "0.90",
        "git_commit": "82ad51f4",
        "literature_comparable": "partial",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "v128",
        "label": "v128 v50b repro",
        "results_subdir": "v128_20260629_0254_v50b_repro",
        "protocol_summary": "v50b baseline reproduction",
        "publishable_tier": "non-standard (pinned old binary)",
        "notes": "Missing post-v50b bug fixes; 2BYS timeout",
        "vct_r0": None,
        "consensus_scorer": "1",
        "native_seed_frac": "0.90",
        "git_commit": "efc4f5d",
        "literature_comparable": "no",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "v130",
        "label": "v130 sulfo+expB full85",
        "results_subdir": "v130_20260629_0548_sulfo_expB_full85",
        "protocol_summary": "oracle-ceiling + sulfo remap + 1HNN expB",
        "publishable_tier": "non-standard (oracle-ceiling); best search-only",
        "notes": "Only 4 seed_echo passes; 81.2% search-only",
        "vct_r0": "4",
        "consensus_scorer": "1",
        "native_seed_frac": "0.90",
        "git_commit": "bf8cf1d2",
        "literature_comparable": "partial",
        "recommended_for_main_text": True,
    },
    {
        "run_id": "v131",
        "label": "v131 r07 nofixb holo+sulfo+expB",
        "results_subdir": "v131_20260629_0835_r07_nofixb_full85",
        "protocol_summary": "r0=7, consensus OFF, holo+sulfo+expB",
        "publishable_tier": "in progress",
        "notes": "Chase v109 record with fairer selector",
        "vct_r0": "7",
        "consensus_scorer": "0",
        "native_seed_frac": "0.90",
        "git_commit": "8864bd17",
        "literature_comparable": "partial",
        "recommended_for_main_text": False,
    },
]

LITERATURE_ROWS: list[dict[str, Any]] = [
    {
        "run_id": "gold_chemplp_all",
        "label": "GOLD 5 + ChemPLP",
        "protocol_summary": "Cognate redock; all binding sites; top-ranked; 25 repeats",
        "publishable_tier": "TIER-1 comparable (cognate)",
        "notes": "No crystal pose seeding; 6 Å site; waters excluded",
        "headline_pass": 69,
        "source": "CCDC GOLD workcase (Astex Diverse symposium)",
        "doi": "",
        "recommended_for_main_text": True,
    },
    {
        "run_id": "gold_chemplp_best",
        "label": "GOLD 5 + ChemPLP (best site)",
        "protocol_summary": "Cognate redock; best site per target selected post hoc",
        "publishable_tier": "TIER-1 comparable (optimistic site)",
        "notes": "Best-site selection inflates vs all-sites metric",
        "headline_pass": 74,
        "source": "CCDC GOLD workcase",
        "doi": "",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "rdock",
        "label": "rDock",
        "protocol_summary": "Cognate self-dock",
        "publishable_tier": "TIER-1 comparable (cognate)",
        "notes": "Verify against Ruiz-Carmona rDock paper before thesis citation",
        "headline_pass": 75,
        "headline_bcr_pct": 87.8,
        "source": "FlexAIDdS BENCHMARK_STANDARD.md comparator",
        "doi": "",
        "recommended_for_main_text": True,
    },
    {
        "run_id": "vina_hartshorn",
        "label": "AutoDock Vina",
        "protocol_summary": "Cognate redock",
        "publishable_tier": "TIER-1 comparable (cognate)",
        "notes": "Classic baseline ~58%",
        "headline_pass": 49,
        "headline_bcr_pct": 57.65,
        "source": "Hartshorn 2007 follow-up",
        "doi": "10.1021/jm061277y",
        "recommended_for_main_text": True,
    },
    {
        "run_id": "flexaid2015_full",
        "label": "FlexAID 2015",
        "protocol_summary": "Native self-dock",
        "publishable_tier": "TIER-1 lineage",
        "notes": "Predecessor engine; non-native focus in paper",
        "headline_pass": 58,
        "headline_bcr_pct": 67.9,
        "source": "BENCHMARK_STANDARD.md",
        "doi": "10.1021/acs.jcim.5b00078",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "flexaidds_v50b_t2",
        "label": "FlexAIDdS v50b (oracle cross-dock)",
        "protocol_summary": "TIER-2 oracle cross-dock",
        "publishable_tier": "TIER-2 only",
        "notes": "Non-cognate receptor; NOT cognate redock",
        "headline_pass": 69,
        "headline_bcr_pct": 81.18,
        "source": "BENCHMARK_STANDARD.md",
        "doi": "",
        "recommended_for_main_text": False,
    },
    {
        "run_id": "surfdock_t2",
        "label": "SurfDock",
        "protocol_summary": "TIER-2 oracle cross-dock",
        "publishable_tier": "TIER-2 only",
        "notes": "AI method; pocket given",
        "headline_pass": 65,
        "headline_bcr_pct": 77.0,
        "source": "BENCHMARK_STANDARD.md",
        "doi": "",
        "recommended_for_main_text": False,
    },
]

PUBLISHABILITY: dict[str, dict[str, Any]] = {
    "v109": {
        "main_text": False,
        "supplementary": True,
        "headline_claim_allowed": False,
        "search_only_claim": "40.0% — below AutoDock Vina (~58%)",
        "vs_gold_chemplp": "Headline appears +13 pp but 57.5% of successes are crystal seed copies",
        "vs_rdock": "Not comparable on headline; search-only far below rDock 87.8%",
        "must_disclose": [
            "oracle-ceiling mode (IC anchored to crystal ligand)",
            "46/80 successes are ini_elitism with seed_echo=1",
            "median success RMSD = 0.00 Å (artifact of seed election)",
            "consensus_scorer=1",
        ],
        "eli5": "Looks like 94% but half the wins are 'you already knew the answer' copies.",
    },
    "v127": {
        "main_text": False,
        "supplementary": True,
        "headline_claim_allowed": False,
        "search_only_claim": "62.4% — between Vina and GOLD",
        "vs_gold_chemplp": "Search-only below GOLD 81%",
        "vs_rdock": "Search-only below rDock 87.8%",
        "must_disclose": [
            "Fix-B logsumexp selector (regresses 1SG0/1T9B/1GPK vs v109)",
            "native_seed_frac=0.90",
            "25 seed_echo passes",
            "r0=4 (not v109 r0=7)",
        ],
        "eli5": "High score but still cheats a bit with seeded poses; selector changes hurt some targets.",
    },
    "v128": {
        "main_text": False,
        "supplementary": True,
        "headline_claim_allowed": False,
        "search_only_claim": "50.6% — near Vina",
        "vs_gold_chemplp": "Below GOLD on both metrics",
        "vs_rdock": "Well below rDock",
        "must_disclose": [
            "Pinned v50b binary — missing 6+ bug fixes",
            "2BYS timeout (data bug)",
            "Not a current-engine result",
        ],
        "eli5": "Old engine replay — floor, not ceiling.",
    },
    "v130": {
        "main_text": True,
        "supplementary": True,
        "headline_claim_allowed": "Only with oracle-ceiling qualifier",
        "search_only_claim": "81.2% — matches GOLD ChemPLP 81% (all sites)",
        "vs_gold_chemplp": "Search-only ≈ tie with GOLD ChemPLP; headline +5 pp with caveat",
        "vs_rdock": "Search-only below rDock 87.8% by ~6.6 pp",
        "must_disclose": [
            "oracle-ceiling mode",
            "Only 4/73 successes are seed_echo (5.5%)",
            "pose_source: 47 ga_cluster, 22 bcr_gate, 4 ini_elitism",
            "1TW6 still used 4-chain PDB (fixed in v131 JSON)",
            "sulfo remap + 1HNN expB site bundled",
            "r0=4, consensus_scorer=1",
        ],
        "eli5": "Best honest run: 81% found by search, almost no crystal copy wins.",
    },
    "v131": {
        "main_text": False,
        "supplementary": False,
        "status": "in_progress",
        "expected_claim": "Target ≥80/85 search-only with r0=7 + consensus off + holo pipeline",
        "must_disclose": [
            "r0=7 (v109-matched restart count)",
            "consensus_scorer=0 (no consensus pose selection)",
            "1TW6_holo.pdb + sulfo remap + 1HNN expB",
            "native_seed_frac=0.90 until TIER-1 run completes",
        ],
        "eli5": "The run designed to beat the real record fairly — wait for full 85 before citing.",
    },
    "gold_chemplp_all": {
        "main_text": True,
        "role": "Primary physics comparator (cognate)",
        "eli5": "Industry standard genetic algorithm docking, no answer seeding.",
    },
    "rdock": {
        "main_text": True,
        "role": "Strongest cited cognate comparator in repo (87.8%)",
        "eli5": "Fast physics docking — the bar to beat for real search.",
    },
}


def _pct(num: int, denom: int) -> float:
    return round(100.0 * num / denom, 2) if denom else 0.0


def analyze_run(results_dir: Path) -> dict[str, Any] | None:
    if not results_dir.is_dir():
        return None
    rows: list[dict[str, str]] = []
    for target_dir in sorted(results_dir.iterdir()):
        rc = target_dir / "result.csv"
        if not rc.is_file():
            continue
        with rc.open(newline="") as fh:
            rows.append(next(csv.DictReader(fh)))
    if not rows:
        return None

    successes = [r for r in rows if r.get("success") in ("1", "True", "true")]
    search_only = [r for r in successes if str(r.get("seed_echo", "0")) != "1"]
    seed_echo = [r for r in successes if str(r.get("seed_echo", "0")) == "1"]

    def count_pose(src: str) -> int:
        return sum(1 for r in successes if r.get("pose_source") == src)

    rmsds = [float(r["rmsd_hungarian"]) for r in successes if r.get("rmsd_hungarian")]
    rmsds.sort()

    n = len(rows)
    hp = len(successes)
    sp = len(search_only)
    return {
        "n_completed": n,
        "complete": n >= N_DENOM,
        "headline_pass": hp,
        "headline_bcr_pct": _pct(hp, N_DENOM),
        "search_only_pass": sp,
        "search_only_bcr_pct": _pct(sp, N_DENOM),
        "seed_echo_passes": len(seed_echo),
        "seed_echo_frac_of_successes_pct": _pct(len(seed_echo), hp) if hp else 0.0,
        "ini_elitism_passes": count_pose("ini_elitism"),
        "ga_cluster_passes": count_pose("ga_cluster"),
        "bcr_gate_passes": count_pose("bcr_gate"),
        "mean_rmsd_success_angstrom": round(sum(rmsds) / len(rmsds), 3) if rmsds else None,
        "median_rmsd_success_angstrom": round(rmsds[len(rmsds) // 2], 3) if rmsds else None,
    }


def build_rows(results_root: Path) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for meta in FLEXAIDDS_RUNS:
        stats = analyze_run(results_root / meta["results_subdir"])
        row = {
            "category": "FlexAIDdS",
            **{k: v for k, v in meta.items() if k != "results_subdir"},
            "n_denominator": N_DENOM,
            "n_completed": stats["n_completed"] if stats else 0,
            "complete": stats["complete"] if stats else False,
        }
        if stats:
            row.update(stats)
        else:
            row.update(
                {
                    "headline_pass": 0,
                    "headline_bcr_pct": 0.0,
                    "search_only_pass": 0,
                    "search_only_bcr_pct": 0.0,
                    "seed_echo_passes": 0,
                    "seed_echo_frac_of_successes_pct": 0.0,
                    "ini_elitism_passes": 0,
                    "ga_cluster_passes": 0,
                    "bcr_gate_passes": 0,
                    "mean_rmsd_success_angstrom": None,
                    "median_rmsd_success_angstrom": None,
                }
            )
        out.append(row)

    for lit in LITERATURE_ROWS:
        hp = lit["headline_pass"]
        pct = lit.get("headline_bcr_pct", _pct(hp, N_DENOM))
        out.append(
            {
                "category": "Literature",
                "run_id": lit["run_id"],
                "label": lit["label"],
                "protocol_summary": lit["protocol_summary"],
                "publishable_tier": lit["publishable_tier"],
                "notes": lit["notes"],
                "n_denominator": N_DENOM,
                "n_completed": N_DENOM,
                "complete": True,
                "headline_pass": hp,
                "headline_bcr_pct": pct,
                "search_only_pass": hp,
                "search_only_bcr_pct": pct,
                "seed_echo_passes": 0,
                "seed_echo_frac_of_successes_pct": 0.0,
                "ini_elitism_passes": 0,
                "ga_cluster_passes": 0,
                "bcr_gate_passes": 0,
                "mean_rmsd_success_angstrom": None,
                "median_rmsd_success_angstrom": None,
                "vct_r0": None,
                "consensus_scorer": None,
                "native_seed_frac": None,
                "git_commit": None,
                "literature_comparable": "yes",
                "recommended_for_main_text": lit["recommended_for_main_text"],
                "source": lit.get("source"),
                "doi": lit.get("doi", ""),
            }
        )
    return out


CSV_FIELDS = [
    "category",
    "run_id",
    "label",
    "protocol_summary",
    "publishable_tier",
    "notes",
    "n_denominator",
    "n_completed",
    "complete",
    "headline_pass",
    "headline_bcr_pct",
    "search_only_pass",
    "search_only_bcr_pct",
    "seed_echo_passes",
    "seed_echo_frac_of_successes_pct",
    "ini_elitism_passes",
    "ga_cluster_passes",
    "bcr_gate_passes",
    "mean_rmsd_success_angstrom",
    "median_rmsd_success_angstrom",
    "vct_r0",
    "consensus_scorer",
    "native_seed_frac",
    "git_commit",
    "literature_comparable",
    "recommended_for_main_text",
]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_json(path: Path, rows: list[dict[str, Any]]) -> None:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Astex Diverse 85 (Hartshorn 2007)",
        "success_criterion": f"rmsd_hungarian < {SUCCESS_RMSD} Angstrom (top-1)",
        "search_only_definition": "success AND seed_echo != 1",
        "disclosure_requirements": [
            "Report tier (TIER-1 cognate / TIER-2 oracle cross-dock / TIER-3 blind)",
            "Report headline BCR and search-only BCR side by side",
            "Report seed_echo pass count and pose_source breakdown",
            "State oracle-ceiling vs autonomous mode",
            "State FLEXAIDDS_NATIVE_SEED_FRAC and whether IC is crystal-anchored",
            "Do not compare oracle-ceiling headline to GOLD/Vina without search-only column",
        ],
        "rows": rows,
        "publishability_by_run": PUBLISHABILITY,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def _tex_pct(v: float) -> str:
    return f"{v:.1f}"


def write_tex(path: Path, rows: list[dict[str, Any]]) -> None:
    flex = [r for r in rows if r["category"] == "FlexAIDdS" and r["complete"]]
    flex.sort(key=lambda r: r["search_only_bcr_pct"], reverse=True)
    lit = [r for r in rows if r["category"] == "Literature" and r["run_id"] in (
        "gold_chemplp_all", "gold_chemplp_best", "rdock", "vina_hartshorn", "flexaid2015_full"
    )]

    lines = [
        "% Astex Diverse 85 — publication comparison table",
        "% Generated by scripts/generate_publication_bcr_table.py",
        f"% Metric: top-1 RMSD_Hungarian < {SUCCESS_RMSD:.1f} A",
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Astex Diverse Set ($N{=}85$) binding-mode correctness (BCR): headline vs.\ search-only (excluding \texttt{seed\_echo} crystal-input poses). Literature rows use cognate redocking protocols; FlexAIDdS rows use oracle-ceiling mode unless noted.}",
        r"\label{tab:astex-diverse-bcr}",
        r"\small",
        r"\begin{tabular}{@{}llcccccc@{}}",
        r"\toprule",
        r"\textbf{Method} & \textbf{Protocol tier} &",
        r"\multicolumn{2}{c}{\textbf{Headline BCR}} &",
        r"\multicolumn{2}{c}{\textbf{Search-only BCR}} &",
        r"\textbf{Seed-echo} & \textbf{Comparable?} \\",
        r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
        r" & & Pass & \% & Pass & \% & passes & \\",
        r"\midrule",
        r"\multicolumn{8}{l}{\textit{Published cognate redocking comparators}} \\",
    ]
    for r in lit:
        comp = "yes" if r["literature_comparable"] == "yes" else "no"
        star = ""
        lines.append(
            f"{r['label']} & {r['publishable_tier']} & {r['headline_pass']} & "
            f"{_tex_pct(r['headline_bcr_pct'])} & {r['search_only_pass']} & "
            f"{_tex_pct(r['search_only_bcr_pct'])} & 0 & {comp} \\\\"
        )
    lines += [
        r"\midrule",
        r"\multicolumn{8}{l}{\textit{FlexAIDdS oracle-ceiling runs (this work)}} \\",
    ]
    for r in flex:
        star = r"$^\star$" if r.get("recommended_for_main_text") else ""
        comp = r.get("literature_comparable", "no")
        lines.append(
            f"{r['run_id']}{star} & {r['publishable_tier'].split(';')[0]} & "
            f"{r['headline_pass']} & {_tex_pct(r['headline_bcr_pct'])} & "
            f"{r['search_only_pass']} & {_tex_pct(r['search_only_bcr_pct'])} & "
            f"{r['seed_echo_passes']} & {comp} \\\\"
        )
    # partial v131
    v131 = next((r for r in rows if r["run_id"] == "v131"), None)
    if v131 and not v131["complete"]:
        n_done = v131["n_completed"]
        lines.append(
            f"v131 (partial, $n={n_done}$) & in progress & "
            f"{v131['headline_pass']} & {_tex_pct(v131['headline_bcr_pct'])} & "
            f"{v131['search_only_pass']} & {_tex_pct(v131['search_only_bcr_pct'])} & "
            f"{v131['seed_echo_passes']} & partial \\\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}\small",
        r"\item $^\star$ Recommended for main-text comparison: lowest seed-echo contamination among complete full-85 runs.",
        r"\item \textbf{Search-only BCR}: successes with \texttt{seed\_echo}$\neq$\texttt{1} (pose found by GA/selector, not crystal INI copy).",
        r"\item FlexAIDdS oracle-ceiling uses oracle binding-site spheres and may anchor the internal-coordinate frame to the crystal ligand; this is \emph{not} identical to GOLD/rDock cognate redocking.",
        r"\item Sources: GOLD ChemPLP/rDock/Vina from comparative Astex literature; FlexAIDdS from per-target \texttt{result.csv} ($n{=}85$).",
        r"\end{tablenotes}",
        r"\end{table*}",
        "",
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Pose provenance for successful targets (FlexAIDdS full-85 runs).}",
        r"\label{tab:astex-pose-source}",
        r"\small",
        r"\begin{tabular}{@{}lrrrr@{}}",
        r"\toprule",
        r"\textbf{Run} & \textbf{ini\_elitism} & \textbf{ga\_cluster} & \textbf{bcr\_gate} & \textbf{seed\_echo} \\",
        r"\midrule",
    ]
    for r in flex:
        lines.append(
            f"{r['run_id']} & {r['ini_elitism_passes']} & {r['ga_cluster_passes']} & "
            f"{r['bcr_gate_passes']} & {r['seed_echo_passes']} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = build_rows(args.results_root)
    write_csv(args.out_dir / "astex_diverse_bcr_comparison.csv", rows)
    write_json(args.out_dir / "astex_diverse_bcr_comparison.json", rows)
    write_tex(args.out_dir / "astex_diverse_bcr_comparison.tex", rows)
    print(f"Wrote {args.out_dir}/astex_diverse_bcr_comparison.{{csv,json,tex}}")
    for r in rows:
        if r["category"] == "FlexAIDdS":
            status = "complete" if r["complete"] else f"partial {r['n_completed']}/85"
            print(
                f"  {r['run_id']:5s} {status:14s} headline {r['headline_pass']}/85 "
                f"search-only {r['search_only_pass']}/85 seed_echo {r['seed_echo_passes']}"
            )


if __name__ == "__main__":
    main()