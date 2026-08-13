"""Guard: regenerable A/B scratch dirs must not sit at the repo root."""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import check_repo_hygiene as hygiene  # noqa: E402


def test_check_root_artifacts_clean(tmp_path: Path):
    (tmp_path / "LIB").mkdir()
    (tmp_path / "scripts").mkdir()
    assert hygiene.check_root_artifacts(tmp_path) == []


def test_check_root_artifacts_flags_ab_mac(tmp_path: Path):
    (tmp_path / "ab_mac_20260806T133329").mkdir()
    errors = hygiene.check_root_artifacts(tmp_path)
    assert len(errors) == 1
    assert "ab_mac_20260806T133329/" in errors[0]


def test_check_root_artifacts_flags_ab_prefix(tmp_path: Path):
    (tmp_path / "ab_smoke_test").mkdir()
    errors = hygiene.check_root_artifacts(tmp_path)
    assert len(errors) == 1
    assert "ab_smoke_test/" in errors[0]


def test_check_root_artifacts_flags_wt_pre(tmp_path: Path):
    (tmp_path / "wt_pre_cpu").mkdir()
    errors = hygiene.check_root_artifacts(tmp_path)
    assert any("wt_pre_cpu/" in e for e in errors)
