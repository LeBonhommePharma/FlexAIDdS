"""
RE-DOCK Command-Line Interface
===============================

Human-in-the-loop CLI for managing distributed docking campaigns.

Commands
--------
  redock init <targets.json>   Initialize campaign from target list
  redock dispatch              Dispatch next generation to workers
  redock ingest <result.json>  Ingest completed chunk results
  redock status                Campaign progress summary
  redock analyze               Van't Hoff decomposition + Shannon entropy

Options
-------
  --campaign-dir PATH    Campaign working directory (default: ./campaign)
                         **On this M3 Pro machine, always point at iCloud**
                         (e.g. ~/Library/Mobile\ Documents/com~apple~CloudDocs/FlexAIDdS/re-dock-campaigns/my-campaign
                          or use paths from ~/.flexaidds_env)

  --t-min FLOAT          Minimum temperature in K (default: 298)
  --t-max FLOAT          Maximum temperature in K (default: 600)
  --n-replicas INT       Number of temperature replicas (default: 8)
  --fit-dcp              Include ΔCp in Van't Hoff fit

Invocation on this machine (M3 Pro + iCloud-only results policy)
----------------------------------------------------------------
Because the directory is named re-dock/ (hyphen), standard "python -m" does not
work directly. Use one of:

  # Recommended for Codex/local resume work (all outputs must go to iCloud)
  PYTHONPATH=. python -c '
import sys
sys.path.insert(0, "benchmarks/re-dock")
from cli import main
main()
' --campaign-dir ~/Library/Mobile\ Documents/.../my-campaign ...

Or:
  cd benchmarks/re-dock && PYTHONPATH=../.. python cli.py --campaign-dir /iCloud/path/...

All persistent outputs (checkpoints, results, visualizations) **must** land on
the 2 TB iCloud Drive per the Storage Invariant.

Le Bonhomme Pharma / Najmanovich Research Group
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .orchestrator import BenchmarkCampaign
from .thermodynamics import shannon_entropy_of_ensemble


def cmd_init(args: argparse.Namespace) -> None:
    """Initialize a new benchmark campaign."""
    campaign = BenchmarkCampaign(
        campaign_dir=args.campaign_dir,
        T_min=args.t_min,
        T_max=args.t_max,
        n_replicas=args.n_replicas,
    )
    campaign.initialize(args.targets)
    print(f"Campaign initialized: {campaign.campaign_id}")
    print(f"  Targets: {len(campaign.targets)}")
    print(f"  Replicas: {campaign.n_replicas}")
    print(f"  Temperatures: {', '.join(f'{T:.1f}K' for T in campaign.temperatures)}")
    print(f"  Checkpoint: {campaign.campaign_dir / 'checkpoint.json'}")


def cmd_dispatch(args: argparse.Namespace) -> None:
    """Dispatch the next generation of docking jobs."""
    checkpoint = Path(args.campaign_dir) / "checkpoint.json"
    if not checkpoint.exists():
        print("Error: No campaign found. Run 'redock init' first.", file=sys.stderr)
        sys.exit(1)

    campaign = BenchmarkCampaign.load_checkpoint(str(checkpoint))
    generation = campaign.current_generation + 1

    chunks = campaign.dispatch_generation(generation)
    print(f"Generation {generation}: dispatched {len(chunks)} chunks")

    # Write worker scripts to campaign dir
    scripts_dir = campaign.campaign_dir / f"gen_{generation}"
    scripts_dir.mkdir(exist_ok=True)
    for chunk in chunks:
        script_path = scripts_dir / f"{chunk.chunk_id}.py"
        with open(script_path, "w") as f:
            f.write(chunk.to_worker_script())

    print(f"  Worker scripts written to: {scripts_dir}")


def cmd_ingest(args: argparse.Namespace) -> None:
    """Ingest completed chunk results."""
    checkpoint = Path(args.campaign_dir) / "checkpoint.json"
    if not checkpoint.exists():
        print("Error: No campaign found.", file=sys.stderr)
        sys.exit(1)

    campaign = BenchmarkCampaign.load_checkpoint(str(checkpoint))

    with open(args.result) as f:
        result = json.load(f)

    chunk_id = result.get("chunk_id", "")
    campaign.process_chunk_result(chunk_id, result)

    # Run exchange round for this target
    pdb_id = result.get("pdb_id", "")
    if pdb_id:
        exchanges = campaign.run_exchange_round(pdb_id)
        accepted = sum(1 for e in exchanges if e["accepted"])
        print(f"Ingested {chunk_id} | Exchanges: {accepted}/{len(exchanges)} accepted")

    campaign.save_checkpoint()


def cmd_status(args: argparse.Namespace) -> None:
    """Show campaign status summary."""
    checkpoint = Path(args.campaign_dir) / "checkpoint.json"
    if not checkpoint.exists():
        print("No active campaign found.", file=sys.stderr)
        sys.exit(1)

    campaign = BenchmarkCampaign.load_checkpoint(str(checkpoint))

    print(f"Campaign: {campaign.campaign_id}")
    print(f"Generation: {campaign.current_generation}")
    print(f"Targets: {len(campaign.targets)}")
    print(f"Replicas: {campaign.n_replicas}")
    print(f"Temperature range: {campaign.T_min:.1f}K – {campaign.T_max:.1f}K")

    total_chunks = len(campaign.chunks)
    completed = sum(1 for c in campaign.chunks if c.status == "completed")
    print(f"Chunks: {completed}/{total_chunks} completed")

    for target in campaign.targets:
        pid = target.pdb_id
        if pid in campaign.replicas:
            n_poses = sum(len(r.poses) for r in campaign.replicas[pid])
            print(f"  {pid}: {n_poses} poses collected")
            if pid in campaign.vant_hoff_results:
                vr = campaign.vant_hoff_results[pid]
                print(f"    ΔH°={vr.delta_H_kcal:.2f} kcal/mol, "
                      f"ΔS°={vr.delta_S_cal:.1f} cal/mol·K, "
                      f"ΔG°={vr.delta_G_kcal:.2f} kcal/mol (R²={vr.r_squared:.3f})")


def cmd_analyze(args: argparse.Namespace) -> None:
    """Run Van't Hoff analysis and Shannon entropy on all targets."""
    checkpoint = Path(args.campaign_dir) / "checkpoint.json"
    if not checkpoint.exists():
        print("Error: No campaign found.", file=sys.stderr)
        sys.exit(1)

    campaign = BenchmarkCampaign.load_checkpoint(str(checkpoint))

    for target in campaign.targets:
        pid = target.pdb_id
        print(f"\n{'='*60}")
        print(f"Target: {pid} ({target.receptor_name})")
        print(f"{'='*60}")

        # Van't Hoff
        vr = campaign.run_vant_hoff(pid, fit_dCp=args.fit_dcp)
        if vr:
            print(f"  Van't Hoff (R²={vr.r_squared:.4f}):")
            print(f"    ΔH° = {vr.delta_H_kcal:.2f} kcal/mol")
            print(f"    ΔS° = {vr.delta_S_cal:.1f} cal/(mol·K)")
            print(f"    ΔG° = {vr.delta_G_kcal:.2f} kcal/mol at {vr.T_ref:.1f}K")
            if vr.delta_Cp_cal != 0:
                print(f"    ΔCp = {vr.delta_Cp_cal:.1f} cal/(mol·K)")
            if not math.isnan(target.exp_dG_kcal):
                print(f"    Exp  = {target.exp_dG_kcal:.2f} kcal/mol "
                      f"(error: {abs(vr.delta_G_kcal - target.exp_dG_kcal):.2f})")

        # Shannon entropy per replica
        if pid in campaign.replicas:
            print(f"  Shannon entropy by replica:")
            for replica in campaign.replicas[pid]:
                if replica.poses:
                    energies = [p.energy_kcal for p in replica.poses]
                    S = shannon_entropy_of_ensemble(energies, replica.temperature)
                    print(f"    T={replica.temperature:.1f}K: "
                          f"S_config={S:.4f} kcal/(mol·K), "
                          f"n_poses={len(replica.poses)}, "
                          f"exchange_rate={replica.acceptance_ratio:.2f}")

    campaign.save_checkpoint()


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog="redock",
        description="RE-DOCK: Replica Exchange Distributed Orchestrated Docking Kit",
    )
    parser.add_argument(
        "--campaign-dir", default="./campaign",
        help="Campaign working directory (default: ./campaign)",
    )
    parser.add_argument(
        "--t-min", type=float, default=298.0,
        help="Minimum temperature in K (default: 298)",
    )
    parser.add_argument(
        "--t-max", type=float, default=600.0,
        help="Maximum temperature in K (default: 600)",
    )
    parser.add_argument(
        "--n-replicas", type=int, default=8,
        help="Number of temperature replicas (default: 8)",
    )

    sub = parser.add_subparsers(dest="command", help="Available commands")

    # init
    p_init = sub.add_parser("init", help="Initialize campaign from target list")
    p_init.add_argument("targets", help="Path to targets JSON file")

    # dispatch
    sub.add_parser("dispatch", help="Dispatch next generation to workers")

    # ingest
    p_ingest = sub.add_parser("ingest", help="Ingest completed chunk results")
    p_ingest.add_argument("result", help="Path to result JSON file")

    # status
    sub.add_parser("status", help="Campaign progress summary")

    # analyze
    p_analyze = sub.add_parser("analyze", help="Van't Hoff + Shannon analysis")
    p_analyze.add_argument(
        "--fit-dcp", action="store_true",
        help="Include ΔCp in Van't Hoff fit",
    )

    return parser


def main() -> None:
    """CLI entry point."""
    import math  # noqa: F811 — needed in cmd_analyze scope

    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    commands = {
        "init": cmd_init,
        "dispatch": cmd_dispatch,
        "ingest": cmd_ingest,
        "status": cmd_status,
        "analyze": cmd_analyze,
    }

    commands[args.command](args)


if __name__ == "__main__":
    main()
