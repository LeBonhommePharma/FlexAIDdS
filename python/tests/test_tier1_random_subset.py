"""Tier-1 target selection: fixed vs seeded-random.

The gate used to take ``targets[:N]`` unconditionally, so on astex_diverse only
1gpk/1mq6 were ever exercised and the other 83 codes were never docked.
``tier1_selection: random`` rotates the draw while keeping any single run exactly
reproducible from its logged seed.

These tests pin both halves of that contract: the rotation (different seeds give
different draws) and, more importantly, the reproducibility (same seed always
gives the same draw) -- a random gate that cannot be replayed would make every
failure un-investigable.
"""

from __future__ import annotations

import pytest

from flexaidds.dataset_runner.runner import DatasetConfig

CODES = [f"t{i:03d}" for i in range(85)]


def _cfg(**kw) -> DatasetConfig:
    base = dict(slug="s", name="n", description="d", targets=list(CODES),
                tier1_subset_size=2)
    base.update(kw)
    return DatasetConfig(**base)


def test_fixed_selection_is_unchanged_default():
    """Default stays first-N so existing datasets keep their behaviour."""
    cfg = _cfg()
    assert cfg.tier1_selection == "fixed"
    assert cfg.tier1_targets() == CODES[:2]


def test_random_selection_is_reproducible_for_a_seed(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    cfg = _cfg(tier1_selection="random", tier1_seed=12345)
    first = cfg.tier1_targets()
    # Repeated calls inside one run must not redraw -- callers invoke this more
    # than once (runner.py and benchmarks/nextgen/runner.py both do).
    assert cfg.tier1_targets() == first
    assert _cfg(tier1_selection="random", tier1_seed=12345).tier1_targets() == first


def test_random_selection_rotates_across_seeds(monkeypatch):
    """The whole point: the draw must not be constant across seeds."""
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    draws = {
        tuple(_cfg(tier1_selection="random", tier1_seed=s).tier1_targets())
        for s in range(40)
    }
    assert len(draws) > 1, "seeded sampling produced a constant subset"


def test_random_draw_is_valid_subset(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    for s in range(25):
        got = _cfg(tier1_selection="random", tier1_seed=s, tier1_subset_size=5).tier1_targets()
        assert len(got) == 5
        assert len(set(got)) == 5, "sampled with replacement"
        assert set(got) <= set(CODES)
        # dataset order preserved within the draw, so run order is stable
        assert got == sorted(got, key=CODES.index)


def test_env_seed_overrides_yaml_seed(monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", "999")
    from_env = _cfg(tier1_selection="random", tier1_seed=1).tier1_targets()
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED")
    assert from_env == _cfg(tier1_selection="random", tier1_seed=999).tier1_targets()


def test_bad_env_seed_fails_loudly(monkeypatch):
    """A mistyped seed must not silently fall back to a different draw."""
    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", "not-a-number")
    with pytest.raises(ValueError, match="FLEXAIDDS_TIER1_SEED"):
        _cfg(tier1_selection="random").tier1_targets()


def test_unknown_selection_mode_rejected(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    with pytest.raises(ValueError, match="tier1_selection"):
        _cfg(tier1_selection="shuffle").tier1_targets()


def test_subset_larger_than_target_list_is_clamped(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    cfg = _cfg(tier1_selection="random", tier1_seed=7, tier1_subset_size=500)
    assert len(cfg.tier1_targets()) == len(CODES)


def test_astex_diverse_yaml_uses_random_selection():
    """Both on-disk copies of the dataset must agree (they have drifted before)."""
    import pathlib

    import yaml

    root = pathlib.Path(__file__).resolve().parents[2]
    copies = [
        root / "benchmarks/datasets/astex_diverse.yaml",
        root / "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml",
    ]
    seen = []
    for p in copies:
        if not p.exists():
            pytest.skip(f"{p} not present")
        d = yaml.safe_load(p.read_text())
        seen.append((d.get("tier1_selection"), d.get("tier1_subset_size")))
    # 4 = the free-coverage point: one batch of 4 workers, same wall as 2.
    assert seen[0] == ("random", 4)
    assert seen[0] == seen[1], f"astex_diverse copies disagree: {seen}"


# --------------------------------------------------------------------------
# Tier-2 subsetting + CI time-budget sizing
# --------------------------------------------------------------------------

def test_tier2_defaults_to_every_target():
    """size 0 means all -- tier-2's historical behaviour must be the default."""
    cfg = _cfg()
    assert cfg.tier2_subset_size == 0
    assert cfg.tier_targets(2) == CODES


def test_tier2_random_subset_is_seeded_and_reproducible(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER2_SEED", raising=False)
    a = _cfg(tier2_selection="random", tier2_subset_size=12, tier2_seed=5)
    b = _cfg(tier2_selection="random", tier2_subset_size=12, tier2_seed=5)
    assert len(a.tier_targets(2)) == 12
    assert a.tier_targets(2) == b.tier_targets(2)
    assert a.tier_targets(2) != _cfg(
        tier2_selection="random", tier2_subset_size=12, tier2_seed=6
    ).tier_targets(2)


def test_tier_seeds_are_independent(monkeypatch):
    """TIER1 env must not silently steer the tier-2 draw."""
    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", "111")
    monkeypatch.delenv("FLEXAIDDS_TIER2_SEED", raising=False)
    cfg = _cfg(tier2_selection="random", tier2_subset_size=6, tier2_seed=222)
    assert cfg.resolve_tier_seed(1) == 111
    assert cfg.resolve_tier_seed(2) == 222


def test_random_with_size_at_or_above_total_returns_all(monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_TIER2_SEED", raising=False)
    cfg = _cfg(tier2_selection="random", tier2_subset_size=len(CODES) + 10)
    assert cfg.tier_targets(2) == CODES


def test_wall_estimate_matches_batch_model(monkeypatch):
    """ceil(n/workers)*per_target -- the model used to size against the CI cap."""
    monkeypatch.delenv("FLEXAIDDS_TIER1_SEED", raising=False)
    cfg = _cfg(tier1_subset_size=4, tier1_selection="fixed")
    # 4 targets on 4 workers is one batch: same wall as 2 targets would be.
    assert cfg.estimate_tier_wall_minutes(1, 4, per_target_minutes=26.9) == pytest.approx(26.9)
    assert cfg.estimate_tier_wall_minutes(1, 2, per_target_minutes=26.9) == pytest.approx(53.8)


def test_full_tier2_would_exceed_the_ci_cap(monkeypatch):
    """Pins the reason tier-2 is subset: the full set does not fit in 360 min."""
    monkeypatch.delenv("FLEXAIDDS_TIER2_SEED", raising=False)
    full = _cfg()  # tier2_subset_size=0 -> all 85
    assert full.estimate_tier_wall_minutes(2, 4, per_target_minutes=26.9) > 360
    sized = _cfg(tier2_selection="random", tier2_subset_size=12, tier2_seed=1)
    assert sized.estimate_tier_wall_minutes(2, 4, per_target_minutes=26.9) <= 360


def test_astex_yaml_tier_sizes_fit_their_ci_budgets():
    """The shipped sizes must actually fit the caps in the workflows."""
    import pathlib

    import yaml
    from flexaidds.dataset_runner.runner import DatasetConfig

    root = pathlib.Path(__file__).resolve().parents[2]
    for rel in ("benchmarks/datasets/astex_diverse.yaml",
                "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml"):
        p = root / rel
        if not p.exists():
            pytest.skip(f"{rel} not present")
        cfg = DatasetConfig.from_yaml(p)
        assert cfg.tier1_selection == "random" and cfg.tier2_selection == "random"
        # tier-1 docking step cap is 105 min; tier-2 job cap is 360 min.
        assert cfg.estimate_tier_wall_minutes(1, 4) <= 105
        assert cfg.estimate_tier_wall_minutes(2, 4) <= 360
