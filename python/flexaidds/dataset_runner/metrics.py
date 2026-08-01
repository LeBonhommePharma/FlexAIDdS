"""Benchmark metric calculations for FlexAIDdS DatasetRunner.

All functions are pure-computation: they accept plain Python / NumPy inputs and
return scalar floats (or structured dicts).  They carry no FlexAIDdS-specific
imports so they can be used in isolation or in CI runners.

Metric catalogue
----------------
* entropy_rescue_rate   — Shannon Energy Collapse: primary ΔS discriminator
* enrichment_factor     — EF at any percentile
* log_auc               — Logarithmic AUC for early enrichment
* scoring_power         — Pearson r + RMSE vs experimental affinities
* docking_power         — % top-ranked poses with RMSD < threshold
* target_specificity_zscore — Z-score of binder scores vs random background
* hit_rate_top_n        — Fraction of true binders in top-N predictions
* bootstrap_ci          — 95% CI via non-parametric bootstrapping
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Data carrier
# ---------------------------------------------------------------------------


@dataclass
class PoseScore:
    """Score for a single docked pose from one target/ligand pair.

    Attributes:
        target_id:          Unique target identifier (e.g. PDB code or gene name).
        ligand_id:          Unique ligand/compound identifier (e.g. CID, SMILES key).
        pose_rank:          1-indexed rank in the docking output list.
        rmsd:               RMSD (Å) to the reference crystal structure; -1 if unknown.
        enthalpy_score:     Contact-function score *without* entropy correction (kcal/mol).
                            Lower = better binding.
        entropy_correction: TΔS contribution (kcal/mol, positive = entropy-favoured).
        total_score:        enthalpy_score − entropy_correction (final ranking score).
        is_active:          Whether the compound is a known binder / true positive.
        exp_affinity:       Experimental ΔG (kcal/mol); None if unavailable.
        structural_state:   Receptor source: ``"holo"``, ``"apo"``, or ``"af2"``.
    """

    target_id: str
    ligand_id: str
    pose_rank: int
    rmsd: float
    enthalpy_score: float
    entropy_correction: float
    total_score: float
    is_active: bool
    exp_affinity: Optional[float] = None
    structural_state: str = "holo"
    ensemble_log_Z: Optional[float] = None  # P3 grand canonical: real log_Z from ensemble (for conc-weighted Xi)


# ---------------------------------------------------------------------------
# Entropy rescue rate (Shannon Energy Collapse)
# ---------------------------------------------------------------------------


def entropy_rescue_rate(
    poses: Sequence[PoseScore],
    rmsd_threshold: float = 2.0,
    rank_threshold: int = 3,
) -> float:
    """Shannon Energy Collapse metric — primary ΔS discriminator.

    Measures the fraction of docking targets where the entropy correction
    *rescues* the crystallographic binding mode: the mode is outside the top
    ``rank_threshold`` positions under enthalpy-only scoring but moves inside
    the top ``rank_threshold`` once the entropy term is applied.

    A rescue event for target *t* requires all of the following:
    1. The crystal pose (lowest RMSD across all poses for *t*) has
       ``rmsd < rmsd_threshold``.
    2. Its enthalpy rank > ``rank_threshold`` (missed by enthalpy alone).
    3. Its entropy-corrected rank ≤ ``rank_threshold`` (rescued by ΔS).

    The denominator is only targets where condition 1 **and** 2 hold (i.e.
    enthalpy alone would have missed the pose).  Returns 0 if the denominator
    is zero.

    Args:
        poses:            All poses across all targets.
        rmsd_threshold:   RMSD cut-off for a pose to count as "correct" (Å).
        rank_threshold:   Top-N threshold for ranking success.

    Returns:
        Rescue rate in [0, 1].
    """
    from collections import defaultdict

    # Group by target
    by_target: dict[str, list[PoseScore]] = defaultdict(list)
    for p in poses:
        by_target[p.target_id].append(p)

    n_missed = 0
    n_rescued = 0

    for target_poses in by_target.values():
        # Find the crystal pose: lowest RMSD that is ≥ 0
        valid = [p for p in target_poses if p.rmsd >= 0]
        if not valid:
            continue
        crystal = min(valid, key=lambda p: p.rmsd)
        if crystal.rmsd >= rmsd_threshold:
            continue  # no reference pose within threshold

        # Rank by enthalpy (lower = better)
        sorted_enthalpy = sorted(target_poses, key=lambda p: p.enthalpy_score)
        enthalpy_rank = next(
            (i + 1 for i, p in enumerate(sorted_enthalpy) if p is crystal), None
        )
        if enthalpy_rank is None or enthalpy_rank <= rank_threshold:
            continue  # enthalpy already found it — not a rescue candidate

        n_missed += 1

        # Rank by total score (lower = better)
        sorted_total = sorted(target_poses, key=lambda p: p.total_score)
        total_rank = next(
            (i + 1 for i, p in enumerate(sorted_total) if p is crystal), None
        )
        if total_rank is not None and total_rank <= rank_threshold:
            n_rescued += 1

    return n_rescued / n_missed if n_missed > 0 else 0.0


# ---------------------------------------------------------------------------
# Virtual screening: enrichment factor
# ---------------------------------------------------------------------------


def enrichment_factor(
    scores: Sequence[float],
    labels: Sequence[bool],
    fraction: float = 0.01,
) -> float:
    """Enrichment factor at a given fraction of the ranked list.

    EF(α) = (hits_in_top_α / N_top_α) / (total_hits / N)

    Args:
        scores:   Predicted scores; lower = more favourable.
        labels:   True binary labels (True = active).
        fraction: Top fraction of the list to consider (default 1%).

    Returns:
        EF value (e.g. 10 means 10× enrichment over random).
    """
    n = len(scores)
    if n == 0:
        return 0.0
    pairs = sorted(zip(scores, labels), key=lambda x: x[0])
    n_top = max(1, int(math.ceil(n * fraction)))
    hits_total = sum(1 for lbl in labels if lbl)
    if hits_total == 0:
        return 0.0
    hits_top = sum(1 for _, lbl in pairs[:n_top] if lbl)
    random_rate = hits_total / n
    return (hits_top / n_top) / random_rate


# ---------------------------------------------------------------------------
# Logarithmic AUC
# ---------------------------------------------------------------------------


def log_auc(
    scores: Sequence[float],
    labels: Sequence[bool],
    min_fraction: float = 0.001,
    max_fraction: float = 1.0,
) -> float:
    """Logarithmic AUC for early-enrichment emphasis.

    Integrates the ROC curve on a log-scale x-axis (false-positive rate) from
    ``min_fraction`` to ``max_fraction``.  This emphasises performance at very
    low false-positive rates where prospective screening operates.

    The random baseline under this metric is 1/ln(max/min) ≈ 0.145 for the
    default [0.001, 1.0] range.

    Args:
        scores:       Predicted scores; lower = more favourable.
        labels:       True binary labels.
        min_fraction: Lower FPR integration bound.
        max_fraction: Upper FPR integration bound.

    Returns:
        logAUC value in [0, 1].
    """
    n = len(scores)
    if n == 0:
        return 0.0
    pairs = sorted(zip(scores, labels), key=lambda x: x[0])
    n_actives = sum(1 for lbl in labels if lbl)
    n_decoys = n - n_actives
    if n_actives == 0 or n_decoys == 0:
        return 0.0

    tpr_vals = []
    fpr_vals = []
    tp = fp = 0
    for _, lbl in pairs:
        if lbl:
            tp += 1
        else:
            fp += 1
        tpr_vals.append(tp / n_actives)
        fpr_vals.append(fp / n_decoys)

    # Integrate via trapezoidal rule on log-FPR axis
    log_min = math.log(min_fraction)
    log_max = math.log(max_fraction)
    integral = 0.0
    prev_fpr = None
    prev_tpr = None
    for fpr, tpr in zip(fpr_vals, tpr_vals):
        if fpr < min_fraction:
            prev_fpr, prev_tpr = fpr, tpr
            continue
        if fpr > max_fraction:
            break
        if prev_fpr is not None and prev_fpr >= min_fraction:
            d_log = math.log(fpr) - math.log(prev_fpr)
            integral += d_log * (tpr + prev_tpr) / 2.0
        prev_fpr, prev_tpr = fpr, tpr

    return integral / (log_max - log_min)


# ---------------------------------------------------------------------------
# Scoring power (affinity correlation)
# ---------------------------------------------------------------------------


def scoring_power(
    predicted: Sequence[float],
    experimental: Sequence[float],
) -> dict[str, float]:
    """Scoring power: Pearson r and RMSE vs experimental affinities.

    Args:
        predicted:    Predicted binding free energies (kcal/mol).
        experimental: Experimental binding free energies (kcal/mol).

    Returns:
        Dict with keys ``pearson_r``, ``spearman_r``, ``rmse``, ``mae``.
    """
    pred = np.asarray(predicted, dtype=float)
    exp = np.asarray(experimental, dtype=float)
    if len(pred) < 2:
        return {"pearson_r": float("nan"), "spearman_r": float("nan"),
                "rmse": float("nan"), "mae": float("nan")}

    # Pearson r
    corr_matrix = np.corrcoef(pred, exp)
    pearson_r = float(corr_matrix[0, 1])

    # Spearman r (rank-based)
    pred_rank = np.argsort(np.argsort(pred)).astype(float)
    exp_rank = np.argsort(np.argsort(exp)).astype(float)
    spearman_r = float(np.corrcoef(pred_rank, exp_rank)[0, 1])

    # RMSE and MAE
    residuals = pred - exp
    rmse = float(np.sqrt(np.mean(residuals**2)))
    mae = float(np.mean(np.abs(residuals)))

    return {
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        "rmse": rmse,
        "mae": mae,
    }


# ---------------------------------------------------------------------------
# Docking power
# ---------------------------------------------------------------------------


def docking_power(
    poses: Sequence[PoseScore],
    rmsd_threshold: float = 2.0,
    top_n: int = 1,
    n_targets: Optional[int] = None,
) -> float:
    """Fraction of targets where the top-N poses contain a near-native pose.

    A target is a "success" if *any* of its top-``top_n`` ranked poses (by
    ``total_score``) has ``0.0 <= rmsd < rmsd_threshold``.  Sentinel RMSD
    values (-1 = not computed, 999 = no pose emitted) are never successes.

    The denominator is the number of targets *attempted*, not the number that
    produced a usable pose: a target that timed out or emitted nothing is a
    failed docking, not an excluded sample.  Pass ``n_targets`` to pin it to
    the full dataset size (e.g. 85 for Astex Diverse) so that targets missing
    from ``poses`` entirely still count against the rate.

    Args:
        poses:          All poses across all targets.
        rmsd_threshold: RMSD cut-off for near-native (Å).
        top_n:          Number of top poses considered.
        n_targets:      Fixed denominator; defaults to the number of distinct
                        targets present in ``poses``.

    Returns:
        Success rate in [0, 1].
    """
    from collections import defaultdict

    by_target: dict[str, list[PoseScore]] = defaultdict(list)
    for p in poses:
        by_target[p.target_id].append(p)

    n_success = 0
    for target_poses in by_target.values():
        # Rank over every emitted pose. Filtering the sentinels out first would
        # promote a worse-scoring pose into the top-N and inflate the rate.
        ranked = sorted(target_poses, key=lambda p: p.total_score)[:top_n]
        if any(0.0 <= p.rmsd < rmsd_threshold for p in ranked):
            n_success += 1

    denom = n_targets if n_targets is not None else len(by_target)
    return n_success / denom if denom > 0 else 0.0


# ---------------------------------------------------------------------------
# Target specificity Z-score
# ---------------------------------------------------------------------------


def target_specificity_zscore(
    target_scores: Sequence[float],
    background_scores: Sequence[float],
) -> float:
    """Z-score of binder scores relative to a random background.

    Z = (μ_target − μ_background) / σ_background

    For docking scores where lower is better, a negative Z-score indicates
    better-than-background binding — callers may negate the result to
    interpret it as an enrichment signal.

    Args:
        target_scores:     Scores of known binders.
        background_scores: Scores of random/decoy compounds.

    Returns:
        Z-score (float; NaN if background has zero variance).
    """
    t = np.asarray(target_scores, dtype=float)
    b = np.asarray(background_scores, dtype=float)
    std_b = float(np.std(b))
    if std_b == 0.0:
        return float("nan")
    return float((np.mean(t) - np.mean(b)) / std_b)


# ---------------------------------------------------------------------------
# Hit rate in top-N
# ---------------------------------------------------------------------------


def hit_rate_top_n(
    scores: Sequence[float],
    labels: Sequence[bool],
    n: int,
) -> float:
    """Fraction of true actives in the top-N ranked compounds.

    Args:
        scores: Predicted scores; lower = more favourable.
        labels: True binary labels.
        n:      Number of top compounds to consider.

    Returns:
        Hit rate in [0, 1].
    """
    pairs = sorted(zip(scores, labels), key=lambda x: x[0])
    top = pairs[:n]
    if not top:
        return 0.0
    hits = sum(1 for _, lbl in top if lbl)
    return hits / len(top)


# ---------------------------------------------------------------------------
# Bootstrap confidence interval
# ---------------------------------------------------------------------------


def bootstrap_ci(
    metric_fn: Callable[[list], float],
    data: list,
    n_resamples: int = 10_000,
    confidence: float = 0.95,
    rng_seed: int = 42,
) -> Tuple[float, float]:
    """Non-parametric bootstrap confidence interval for any metric.

    Resamples ``data`` with replacement ``n_resamples`` times, applies
    ``metric_fn`` to each resample, and returns the percentile CI.

    Args:
        metric_fn:   Function that maps a list → scalar float.
        data:        Input data (list of any objects metric_fn can handle).
        n_resamples: Number of bootstrap resamples.
        confidence:  Desired coverage (default 0.95 → 95% CI).
        rng_seed:    Random seed for reproducibility.

    Returns:
        ``(lower, upper)`` bounds at the requested confidence level.
    """
    rng = random.Random(rng_seed)
    n = len(data)
    stats: list[float] = []
    for _ in range(n_resamples):
        resample = [data[rng.randrange(n)] for _ in range(n)]
        try:
            val = metric_fn(resample)
            if math.isfinite(val):
                stats.append(val)
        except Exception:
            pass

    if not stats:
        return (float("nan"), float("nan"))

    stats.sort()
    alpha = 1.0 - confidence
    lo_idx = int(math.floor(alpha / 2 * len(stats)))
    hi_idx = int(math.ceil((1 - alpha / 2) * len(stats))) - 1
    hi_idx = min(hi_idx, len(stats) - 1)
    return (stats[lo_idx], stats[hi_idx])


# ---------------------------------------------------------------------------
# Convenience: compute all standard metrics from a pose list
# ---------------------------------------------------------------------------


def compute_all_metrics(
    poses: Sequence[PoseScore],
    requested: Optional[Sequence[str]] = None,
    bootstrap: bool = False,
    n_resamples: int = 5_000,
    n_targets: Optional[int] = None,
) -> dict[str, float]:
    """Compute a standard suite of benchmark metrics from a pose list.

    Args:
        poses:      All poses across all targets.
        requested:  Metric names to compute; ``None`` = compute all.
        bootstrap:  Whether to compute 95% bootstrap CIs (slow).
        n_resamples: Bootstrap sample count when ``bootstrap=True``.
        n_targets:  Fixed docking-power denominator (full dataset size).

    Returns:
        Dict of ``{metric_name: value}`` (and ``{metric_name + "_ci_lo/hi"}``
        when bootstrapping).

    Note:
        Callers running with synthetic/dry-run poses must not treat
        ``docking_power_*`` values as real docking success rates.  The
        DatasetRunner strips those keys when ``dry_run=True``.
    """
    _all = {
        "entropy_rescue_rate",
        "docking_power_top1",
        "docking_power_top3",
        "mean_rmsd",
        "median_rmsd",
        "ef_1pct",
        "ef_5pct",
        "log_auc",
        "scoring_power_pearson_r",
        "scoring_power_rmse",
        "scoring_power_spearman_r",
        "hit_rate_top10",
    }
    to_compute = set(requested) if requested else _all

    results: dict[str, float] = {}

    if "entropy_rescue_rate" in to_compute:
        results["entropy_rescue_rate"] = entropy_rescue_rate(poses)

    if "docking_power_top1" in to_compute:
        results["docking_power_top1"] = docking_power(
            poses, top_n=1, n_targets=n_targets)

    if "docking_power_top3" in to_compute:
        results["docking_power_top3"] = docking_power(
            poses, top_n=3, n_targets=n_targets)

    # Pose-accuracy aggregates: per-target best-pose RMSD, then mean/median
    # across targets — mirrors benchmark.BenchmarkSummary's median-of-best-pose
    # definition.  Sentinel RMSDs (< 0 = not computed, 999 = no pose) are
    # excluded, so a target that produced no usable pose contributes NOTHING to
    # the aggregate.  This is a deliberate asymmetry with docking_power(), whose
    # denominator counts an empty target as a FAILURE: a success *rate* is
    # defined over attempts, but you cannot average a distance that does not
    # exist.  The risk of that asymmetry is that a mostly-failed multi-target
    # run reports a clean mean over its few survivors — so the coverage count is
    # emitted alongside (``mean_rmsd_n_targets``) to put the denominator in the
    # report rather than leave it inferred.  When no target yields a valid
    # best-pose RMSD the metric is left ABSENT so the #326 completeness gate
    # reports it unmeasured rather than emitting a misleading 0.0.
    if "mean_rmsd" in to_compute or "median_rmsd" in to_compute:
        from collections import defaultdict as _defaultdict
        _rmsd_by_target: dict[str, list[float]] = _defaultdict(list)
        for p in poses:
            # Exclude both sentinels: -1 (not computed) and 999 (no pose
            # emitted).  Averaging 999 would silently poison the aggregate,
            # unlike docking_power where the rmsd < 2.0 threshold neutralises it.
            if 0.0 <= p.rmsd < 999.0:
                _rmsd_by_target[p.target_id].append(p.rmsd)
        best_rmsds = sorted(min(v) for v in _rmsd_by_target.values() if v)
        if best_rmsds:
            n = len(best_rmsds)
            # Coverage: number of targets the aggregate is actually over, so a
            # reader sees "3.16 over 5 targets" rather than a clean-looking mean.
            # Informational only — not a declared baseline, so it is neither
            # gated by completeness nor direction-checked by check_regressions.
            results["mean_rmsd_n_targets"] = float(n)
            if "mean_rmsd" in to_compute:
                results["mean_rmsd"] = sum(best_rmsds) / n
            if "median_rmsd" in to_compute:
                mid = n // 2
                results["median_rmsd"] = (
                    best_rmsds[mid] if n % 2 == 1
                    else (best_rmsds[mid - 1] + best_rmsds[mid]) / 2.0
                )

    # Virtual-screening metrics require active/decoy labels and scores
    all_scores = [p.total_score for p in poses]
    all_labels = [p.is_active for p in poses]

    if any(lbl for lbl in all_labels):
        if "ef_1pct" in to_compute:
            results["ef_1pct"] = enrichment_factor(all_scores, all_labels, 0.01)
        if "ef_5pct" in to_compute:
            results["ef_5pct"] = enrichment_factor(all_scores, all_labels, 0.05)
        if "log_auc" in to_compute:
            results["log_auc"] = log_auc(all_scores, all_labels)
        if "hit_rate_top10" in to_compute:
            results["hit_rate_top10"] = hit_rate_top_n(all_scores, all_labels, 10)

    # Scoring power requires experimental affinities
    exp_pairs = [(p.total_score, p.exp_affinity) for p in poses
                 if p.exp_affinity is not None]
    if len(exp_pairs) >= 2:
        pred_vals = [x[0] for x in exp_pairs]
        exp_vals = [x[1] for x in exp_pairs]
        sp = scoring_power(pred_vals, exp_vals)
        if "scoring_power_pearson_r" in to_compute:
            results["scoring_power_pearson_r"] = sp["pearson_r"]
        if "scoring_power_rmse" in to_compute:
            results["scoring_power_rmse"] = sp["rmse"]
        if "scoring_power_spearman_r" in to_compute:
            results["scoring_power_spearman_r"] = sp["spearman_r"]

    return results
