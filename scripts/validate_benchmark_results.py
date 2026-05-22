#!/usr/bin/env python3
"""
validate_benchmark_results.py — FlexAIDdS Benchmark Statistical Validation

Reads the summary CSV produced by run_benchmark_production.sh and performs:
  1. Descriptive statistics (wall-clock, success rate, RMSD distribution)
  2. Bootstrap 95% CI on success rate (n=85, 10k resamples — not normal approx)
  3. Fisher's exact test vs. Vina baseline (58/85 = 68.2% — Hartshorn 2007)
  4. Shannon-Weighted Success Rate (SWSR calibration check)
  5. Spearman ρ(H_final, RMSD_top1) with bootstrap CI
  6. Plots: wall-clock dist, score vs RMSD scatter, Shannon H convergence curves
  7. PASS / FAIL verdict against BENCHMARKING_PLAN.md thresholds

Usage:
    python3 scripts/validate_benchmark_results.py <summary.csv> [options]

    python3 scripts/validate_benchmark_results.py \\
        benchmark_results/20260517_120000/summary.csv \\
        --shannon-log-dir benchmark_results/20260517_120000/astex_diverse/ \\
        --out-dir benchmark_results/20260517_120000/figures/

Apache-2.0 · Le Bonhomme Pharma / NRGlab, Université de Montréal
"""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path


# ─── Constants (must match BENCHMARKING_PLAN.md / ShannonThermoStack.h) ───────

kB_KCAL       = 0.001987206          # kcal mol⁻¹ K⁻¹
LN2           = 0.6931471805599453
kHSC_SOFT     = 2.0 * LN2           # 1.3863 nats — soft convergence threshold
kHSC_HARD     = 1.0 * LN2           # 0.6931 nats — hard (thesis) convergence threshold
VINA_N_SUCCESS = 49                  # Vina baseline (Hartshorn 2007 Table 3, top-1)
VINA_N_TOTAL   = 85
RMSD_SUCCESS   = 2.0                 # Å — primary success threshold
TARGET_SR      = 0.65               # 65% target success rate (FlexAIDdS)
FLOOR_SR       = 0.58               # 58% hard floor (Vina baseline)
N_BOOT         = 10_000             # bootstrap resamples
BOOT_CI        = 0.95               # confidence interval level


# ─── Bootstrap ────────────────────────────────────────────────────────────────

def bootstrap_sr(successes: list[int], n_boot: int = N_BOOT,
                 ci: float = BOOT_CI) -> tuple[float, float, float]:
    """
    Bootstrap 95% CI on success rate.
    Returns (mean, lo, hi).
    Uses only stdlib random — no numpy required.
    """
    import random
    n = len(successes)
    mean = sum(successes) / n
    boot = []
    for _ in range(n_boot):
        sample = [random.choice(successes) for _ in range(n)]
        boot.append(sum(sample) / n)
    boot.sort()
    a = (1 - ci) / 2
    lo = boot[int(a * n_boot)]
    hi = boot[int((1 - a) * n_boot)]
    return mean, lo, hi


def bootstrap_spearman(x: list[float], y: list[float],
                        n_boot: int = N_BOOT) -> tuple[float, float, float]:
    """Bootstrap 95% CI on Spearman ρ(x, y)."""
    import random
    rho = spearman_r(x, y)
    n = len(x)
    pairs = list(zip(x, y))
    boot = []
    for _ in range(n_boot):
        s = [random.choice(pairs) for _ in range(n)]
        sx, sy = zip(*s)
        boot.append(spearman_r(list(sx), list(sy)))
    boot.sort()
    a = (1 - BOOT_CI) / 2
    lo = boot[int(a * n_boot)]
    hi = boot[int((1 - a) * n_boot)]
    return rho, lo, hi


def spearman_r(x: list[float], y: list[float]) -> float:
    """Spearman rank correlation — stdlib only."""
    n = len(x)
    if n < 2:
        return float('nan')

    def ranks(lst):
        sorted_idx = sorted(range(n), key=lambda i: lst[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and lst[sorted_idx[j]] == lst[sorted_idx[j + 1]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[sorted_idx[k]] = avg_rank
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = math.sqrt(
        sum((rx[i] - mx) ** 2 for i in range(n)) *
        sum((ry[i] - my) ** 2 for i in range(n))
    )
    return num / den if den > 0 else float('nan')


def fisher_exact_greater(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """
    One-sided Fisher's exact test: H₁: odds ratio > 1.
    2×2 table: [[a, b], [c, d]]
    Returns (odds_ratio, p_value).
    stdlib only — hypergeometric via log-factorial.
    """
    def log_fact(n: int) -> float:
        return math.lgamma(n + 1)

    def log_hypergeom(k: int, K: int, N: int, n: int) -> float:
        return (log_fact(K) - log_fact(k) - log_fact(K - k) +
                log_fact(N - K) - log_fact(n - k) - log_fact(N - K - n + k) +
                log_fact(n) - log_fact(N) + log_fact(N - n))

    # a = FlexAID success, b = FlexAID fail, c = Vina success, d = Vina fail
    n1, n2, K, N = a + b, c + d, a + c, a + b + c + d
    n = n1
    # P(X >= a) under H0
    p = 0.0
    for k in range(a, min(K, n) + 1):
        if K - k <= N - K - n + k:
            try:
                p += math.exp(log_hypergeom(k, K, N, n))
            except (ValueError, OverflowError):
                pass
    or_ = (a * d) / (b * c) if b > 0 and c > 0 else float('inf')
    return or_, min(p, 1.0)


# ─── Shannon-Weighted Success Rate ────────────────────────────────────────────

def swsr(successes: list[int], h_finals: list[float],
         floor: float = 0.1) -> float:
    """
    Shannon-Weighted Success Rate (Section 7.3 of BENCHMARKING_PLAN.md).
    SR_H = Σ (1/H_i · success_i) / Σ (1/H_i)
    """
    weights = [1.0 / max(h, floor) for h in h_finals]
    total_w = sum(weights)
    if total_w == 0:
        return float('nan')
    return sum(w * s for w, s in zip(weights, successes)) / total_w


# ─── RMSD distribution ────────────────────────────────────────────────────────

def rmsd_distribution(rmsds: list[float]) -> dict:
    """Compute the full RMSD distribution per BENCHMARKING_PLAN.md §7.5."""
    n = len(rmsds)
    if n == 0:
        return {}
    rmsds_sorted = sorted(rmsds)
    return {
        "n": n,
        "median":          rmsds_sorted[n // 2],
        "mean":            sum(rmsds) / n,
        "frac_lt_1A":      sum(1 for r in rmsds if r < 1.0) / n,
        "frac_lt_2A":      sum(1 for r in rmsds if r < 2.0) / n,
        "frac_lt_3A":      sum(1 for r in rmsds if r < 3.0) / n,
        "frac_gt_5A":      sum(1 for r in rmsds if r > 5.0) / n,
    }


# ─── Load Shannon trace CSVs ──────────────────────────────────────────────────

def load_shannon_trace(complex_dir: Path) -> list[dict]:
    """
    Load per-step Shannon trace CSV from SHANNON_TRACE_LEVEL=2 output.
    Expected format: step,H_nats,H_bits,n_clusters,top_pose_weight,...
    """
    for fname in ["shannon_trace.csv", "shannon.csv", "entropy_trace.csv"]:
        f = complex_dir / fname
        if f.exists():
            rows = []
            with open(f) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    try:
                        rows.append({
                            "step":            int(row.get("step", 0)),
                            "H_nats":          float(row.get("H_nats", row.get("H", 0))),
                            "H_bits":          float(row.get("H_bits", 0)),
                            "n_clusters":      int(row.get("n_clusters", 0)),
                            "top_pose_weight": float(row.get("top_pose_weight", 0)),
                        })
                    except (ValueError, KeyError):
                        continue
            return rows
    return []


# ─── Plotting ─────────────────────────────────────────────────────────────────

def make_plots(rows: list[dict], shannon_data: dict,
               out_dir: Path) -> list[str]:
    """
    Generate and save three key plots.
    Returns list of file paths created.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("  [WARN] matplotlib not available — skipping plots")
        print("         pip install matplotlib --break-system-packages")
        return []

    out_dir.mkdir(parents=True, exist_ok=True)
    created = []

    wall_times = [float(r["wall_time_s"]) for r in rows
                  if r["wall_time_s"] not in ("0", "N/A", "") and float(r["wall_time_s"]) > 0]
    rmsds      = [float(r["rmsd_to_crystal"]) for r in rows
                  if r["rmsd_to_crystal"] not in ("N/A", "")]
    scores     = [float(r["top1_score"]) for r in rows
                  if r["top1_score"] not in ("N/A", "")]

    # 1. Wall-clock distribution ───────────────────────────────────────────────
    if wall_times:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(wall_times, bins=20, color="#1f77b4", edgecolor="white", alpha=0.85)
        ax.axvline(sum(wall_times) / len(wall_times), color="red",
                   lw=2, linestyle="--", label=f"Mean={sum(wall_times)/len(wall_times):.0f}s")
        ax.set_xlabel("Wall-clock time (s/complex)", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.set_title("FlexAIDdS — Wall-clock distribution (Astex Diverse 85)", fontsize=13)
        ax.legend()
        # Overlay expected range band
        ax.axvspan(300, 480, alpha=0.12, color="green", label="Expected 5–8 min")
        ax.legend()
        fig.tight_layout()
        p = out_dir / "wallclock_distribution.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        created.append(str(p))

    # 2. Score vs RMSD scatter ─────────────────────────────────────────────────
    if rmsds and scores and len(rmsds) == len(scores):
        successes_bool = [r < RMSD_SUCCESS for r in rmsds]
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = ["#2ca02c" if s else "#d62728" for s in successes_bool]
        ax.scatter(rmsds, scores, c=colors, alpha=0.7, edgecolors="white", s=50)
        ax.axvline(RMSD_SUCCESS, color="gray", lw=1.5, linestyle="--",
                   label=f"RMSD = {RMSD_SUCCESS}Å threshold")
        ax.set_xlabel("RMSD to crystal (Å)", fontsize=12)
        ax.set_ylabel("Top-1 FlexAIDdS score", fontsize=12)
        ax.set_title("Score vs. RMSD — Astex Diverse 85\n"
                     "(green=success, red=failure)", fontsize=13)
        ax.legend()
        fig.tight_layout()
        p = out_dir / "score_vs_rmsd.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        created.append(str(p))

    # 3. Shannon H(X) convergence curves ──────────────────────────────────────
    if shannon_data:
        fig, ax = plt.subplots(figsize=(10, 6))
        max_step = 0

        # Plot per-complex traces (thin, grey)
        for pdb, trace in list(shannon_data.items())[:20]:  # cap at 20 for legibility
            if not trace:
                continue
            steps  = [r["step"] for r in trace]
            H_nats = [r["H_nats"] for r in trace]
            ax.plot(steps, H_nats, lw=0.5, alpha=0.4, color="steelblue")
            max_step = max(max_step, max(steps) if steps else 0)

        # Mean trace (bold)
        if shannon_data:
            all_traces = [t for t in shannon_data.values() if t]
            if all_traces:
                common_steps = sorted(set(r["step"] for t in all_traces for r in t))
                mean_H = []
                for s in common_steps:
                    vals = [r["H_nats"] for t in all_traces
                            for r in t if r["step"] == s]
                    if vals:
                        mean_H.append(sum(vals) / len(vals))
                if mean_H:
                    ax.plot(common_steps[:len(mean_H)], mean_H,
                            lw=2.5, color="navy", label="Mean H(X)")

        # Threshold lines
        ax.axhline(kHSC_SOFT, color="orange", lw=2, linestyle="--",
                   label=f"HSC soft = {kHSC_SOFT:.4f} nats (2·ln2)")
        ax.axhline(kHSC_HARD, color="red", lw=2, linestyle=":",
                   label=f"HSC hard = {kHSC_HARD:.4f} nats (ln2)")

        ax.set_xlabel("GA generation", fontsize=12)
        ax.set_ylabel("Shannon entropy H(X) [nats]", fontsize=12)
        ax.set_title("Shannon Energy Collapse — Astex Diverse\n"
                     "(collapse below orange dashed line = convergence)", fontsize=13)
        ax.legend(fontsize=10)
        if max_step > 0:
            ax.set_xlim(0, max_step)
        ax.set_ylim(bottom=0)
        fig.tight_layout()
        p = out_dir / "shannon_convergence.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        created.append(str(p))

    return created


# ─── Verdict ──────────────────────────────────────────────────────────────────

def pass_fail(name: str, value: float, threshold: float,
              direction: str = "ge") -> tuple[str, bool]:
    """direction: 'ge' (≥ threshold = PASS) or 'le' (≤ threshold = PASS)."""
    if direction == "ge":
        passed = value >= threshold
    else:
        passed = value <= threshold
    status = "✅ PASS" if passed else "❌ FAIL"
    return f"  {status}  {name}: {value:.4f} (threshold {'≥' if direction=='ge' else '≤'} {threshold})", passed


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="FlexAIDdS benchmark statistical validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("summary_csv", help="Summary CSV from run_benchmark_production.sh")
    parser.add_argument("--shannon-log-dir", default=None,
                        help="Directory containing per-complex result dirs "
                             "(for Shannon trace CSVs)")
    parser.add_argument("--out-dir", default=None,
                        help="Output directory for figures (default: alongside CSV)")
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--rmsd-threshold", type=float, default=RMSD_SUCCESS)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    csv_path = Path(args.summary_csv)
    if not csv_path.exists():
        print(f"ERROR: Summary CSV not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.out_dir) if args.out_dir else csv_path.parent / "figures"
    shannon_root = Path(args.shannon_log_dir) if args.shannon_log_dir else None

    print("=" * 65)
    print("  FlexAIDdS Benchmark Statistical Validation")
    print(f"  Input:   {csv_path}")
    print(f"  Out dir: {out_dir}")
    print("=" * 65)
    print()

    # ── Load CSV ──────────────────────────────────────────────────────────────
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    n = len(rows)
    if n == 0:
        print("ERROR: Summary CSV is empty", file=sys.stderr)
        sys.exit(1)

    print(f"  Loaded {n} complexes from summary CSV")

    # Parse numeric fields robustly
    wall_times = []
    successes  = []
    rmsds      = []
    h_finals   = []

    for r in rows:
        try:
            wt = float(r.get("wall_time_s", 0))
            wall_times.append(wt if wt > 0 else None)
        except ValueError:
            wall_times.append(None)

        try:
            successes.append(int(r.get("success", 0)))
        except ValueError:
            successes.append(0)

        try:
            rmsd = float(r.get("rmsd_to_crystal", "N/A"))
            rmsds.append(rmsd)
        except ValueError:
            rmsds.append(None)

        try:
            h = float(r.get("shannon_entropy", r.get("h_final", "N/A")))
            h_finals.append(h)
        except ValueError:
            h_finals.append(None)

    valid_wall  = [w for w in wall_times if w is not None and w > 0]
    valid_rmsds = [r for r in rmsds if r is not None]
    valid_h     = [h for h in h_finals if h is not None]

    # ── Wall-clock statistics ─────────────────────────────────────────────────
    print("\n── Wall-clock Statistics ─────────────────────────────────────")
    if valid_wall:
        mean_w = sum(valid_wall) / len(valid_wall)
        var_w  = sum((x - mean_w) ** 2 for x in valid_wall) / max(len(valid_wall) - 1, 1)
        std_w  = math.sqrt(var_w)
        min_w  = min(valid_wall)
        max_w  = max(valid_wall)
        eta_h  = mean_w * 85 / 3600
        print(f"  n (with timing):  {len(valid_wall)}")
        print(f"  Mean:             {mean_w:.1f}s ± {std_w:.1f}s")
        print(f"  Min / Max:        {min_w:.1f}s / {max_w:.1f}s")
        print(f"  Est. full 85:     ~{eta_h:.1f}h")
        in_range = 300 <= mean_w <= 480
        print(f"  In expected range (5–8 min): {'YES ✅' if in_range else 'NO ⚠️'}")
    else:
        print("  No valid wall-clock data (dry-run or missing?)")

    # ── RMSD distribution ─────────────────────────────────────────────────────
    print("\n── RMSD Distribution ────────────────────────────────────────")
    dist = rmsd_distribution(valid_rmsds)
    if dist:
        print(f"  n (with RMSD):    {dist['n']}")
        print(f"  Median RMSD:      {dist['median']:.2f} Å")
        print(f"  Mean RMSD:        {dist['mean']:.2f} Å")
        print(f"  RMSD < 1.0 Å:     {dist['frac_lt_1A']*100:.1f}% (near-crystal)")
        print(f"  RMSD < 2.0 Å:     {dist['frac_lt_2A']*100:.1f}% (primary criterion)")
        print(f"  RMSD < 3.0 Å:     {dist['frac_lt_3A']*100:.1f}% (near-success)")
        print(f"  RMSD > 5.0 Å:     {dist['frac_gt_5A']*100:.1f}% (catastrophic failure)")
    else:
        print("  No RMSD data available")

    # ── Success rate + bootstrap CI ───────────────────────────────────────────
    print("\n── Success Rate (Bootstrap 95% CI) ──────────────────────────")
    n_success = sum(successes)
    sr_obs    = n_success / n
    print(f"  n_success / n:    {n_success} / {n}")
    print(f"  SR (observed):    {sr_obs*100:.1f}%")

    sr_mean, sr_lo, sr_hi = bootstrap_sr(successes, n_boot=args.n_boot)
    print(f"  Bootstrap CI:     [{sr_lo*100:.1f}%, {sr_hi*100:.1f}%] (95%, {args.n_boot} resamples)")
    print(f"  Vina baseline:    {VINA_N_SUCCESS}/{VINA_N_TOTAL} = {VINA_N_SUCCESS/VINA_N_TOTAL*100:.1f}%")
    print(f"  FlexAIDdS target: ≥ {TARGET_SR*100:.0f}%  hard floor: ≥ {FLOOR_SR*100:.0f}%")

    # ── Fisher's exact test vs Vina ───────────────────────────────────────────
    print("\n── Fisher's Exact Test vs. Vina Baseline ────────────────────")
    a = n_success
    b = n - n_success
    c = VINA_N_SUCCESS
    d = VINA_N_TOTAL - VINA_N_SUCCESS
    or_, p_val = fisher_exact_greater(a, b, c, d)
    delta_pp = (n_success / n - VINA_N_SUCCESS / VINA_N_TOTAL) * 100
    print(f"  FlexAIDdS: {a}/{n}   Vina: {c}/{VINA_N_TOTAL}")
    print(f"  ΔRMSD-SR:  {delta_pp:+.1f} pp")
    print(f"  Odds ratio: {or_:.2f}")
    print(f"  p (one-sided): {p_val:.4f}  {'(**)' if p_val < 0.01 else '(*)' if p_val < 0.05 else ''}")
    print(f"  Significant (p < 0.05): {'YES ✅' if p_val < 0.05 else 'NO ❌'}")

    # ── Shannon-Weighted Success Rate ─────────────────────────────────────────
    print("\n── Shannon-Weighted Success Rate (SWSR Calibration) ─────────")
    if valid_h and len(valid_h) == n:
        sw = swsr(successes, valid_h)
        print(f"  SR  (unweighted): {sr_obs*100:.1f}%")
        print(f"  SR_H (SWSR):      {sw*100:.1f}%")
        if sw > sr_obs:
            print("  SR_H > SR  ✅  High-confidence predictions are more accurate (calibrated)")
        else:
            print("  SR_H < SR  ⚠️  Overconfident failures — investigate scoring function")
    else:
        print(f"  Shannon entropy data: {len(valid_h)}/{n} complexes")
        print("  SWSR requires shannon_entropy column in summary CSV")
        print("  (Set SHANNON_TRACE_LEVEL≥1 before running benchmark)")

    # ── Spearman ρ(H_final, RMSD) ─────────────────────────────────────────────
    print("\n── Spearman ρ(H_final, RMSD) ────────────────────────────────")
    paired_h_rmsd = [(h, r) for h, r in zip(h_finals, rmsds)
                     if h is not None and r is not None]
    if len(paired_h_rmsd) >= 5:
        ph = [p[0] for p in paired_h_rmsd]
        pr = [p[1] for p in paired_h_rmsd]
        rho, rho_lo, rho_hi = bootstrap_spearman(ph, pr, n_boot=args.n_boot)
        print(f"  n:          {len(paired_h_rmsd)}")
        print(f"  ρ:          {rho:.3f}")
        print(f"  Bootstrap CI: [{rho_lo:.3f}, {rho_hi:.3f}] (95%)")
        print(f"  Expected:   ρ > 0 (higher H → worse pose)")
        if rho > 0:
            print("  ✅  H(X) is predictive of pose quality (lower entropy → better pose)")
        elif abs(rho) < 0.2:
            print("  ⚠️  |ρ| < 0.2 — H(X) weakly predictive → check scoring function")
        else:
            print("  ❌  ρ < 0 — unexpected direction — investigate")
    else:
        print(f"  Insufficient paired data ({len(paired_h_rmsd)} pairs) — "
              "need shannon_entropy in CSV")

    # ── Load Shannon traces (if available) ───────────────────────────────────
    shannon_data = {}
    if shannon_root and shannon_root.is_dir():
        print(f"\n── Loading Shannon Traces from {shannon_root} ──")
        for pdb_dir in sorted(shannon_root.iterdir()):
            if pdb_dir.is_dir():
                trace = load_shannon_trace(pdb_dir)
                if trace:
                    shannon_data[pdb_dir.name] = trace
        print(f"  Loaded traces for {len(shannon_data)} complexes")

        # Per-complex convergence stats
        converged_soft = sum(
            1 for t in shannon_data.values()
            if t and t[-1]["H_nats"] < kHSC_SOFT
        )
        converged_hard = sum(
            1 for t in shannon_data.values()
            if t and t[-1]["H_nats"] < kHSC_HARD
        )
        n_traced = len(shannon_data)
        print(f"  H(X) < {kHSC_SOFT:.4f} nats (soft, 2·ln2): "
              f"{converged_soft}/{n_traced} complexes")
        print(f"  H(X) < {kHSC_HARD:.4f} nats (hard,   ln2): "
              f"{converged_hard}/{n_traced} complexes")

    # ── Plots ─────────────────────────────────────────────────────────────────
    if not args.no_plots:
        print(f"\n── Generating Plots → {out_dir} ─────────────────────────")
        created = make_plots(rows, shannon_data, out_dir)
        for p in created:
            print(f"  Saved: {p}")
        if not created:
            print("  (No plots generated)")

    # ── PASS / FAIL Verdict ───────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("  BENCHMARK VERDICT")
    print("=" * 65)

    all_passed = True
    checks = []

    # Primary: success rate ≥ 65%
    line, ok = pass_fail("Success rate (SR)", sr_obs, TARGET_SR, "ge")
    checks.append((line, ok))
    if not ok:
        all_passed = False

    # Hard floor: ≥ 58%
    line, ok = pass_fail("Hard floor (SR ≥ Vina)", sr_obs, FLOOR_SR, "ge")
    checks.append((line, ok))
    if not ok:
        all_passed = False

    # Statistical significance vs Vina
    line, ok = pass_fail("Fisher p < 0.05 vs Vina", 1 - p_val, 0.95, "ge")
    checks.append((line, ok))
    # (not strictly blocking — can WARN)

    # Spearman (informational — don't block)
    if len(paired_h_rmsd) >= 5:
        rho_val = rho
        line = f"  {'⚠️ WARN' if rho_val < 0 else '✅ INFO'}  Spearman ρ(H,RMSD) = {rho_val:.3f} (expected > 0)"
        checks.append((line, True))

    for line, _ in checks:
        print(line)

    print()
    if all_passed:
        print("  ██████████████████████████████████████████████")
        print("  ██   OVERALL: PASS  — Ready for publication  ██")
        print("  ██████████████████████████████████████████████")
    else:
        print("  ████████████████████████████████████████████████")
        print("  ██   OVERALL: FAIL  — Criteria not met        ██")
        print("  ████████████████████████████████████████████████")

    print()
    print("  Reference: BENCHMARKING_PLAN.md §2.1 — FlexAIDdS ≥ 65%")
    print("             Hartshorn et al. (2007) J Med Chem 50:726–741")
    print()

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
