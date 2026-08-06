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
    assert seen[0] == ("random", 2)
    assert seen[0] == seen[1], f"astex_diverse copies disagree: {seen}"
