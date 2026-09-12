"""Synthetic regressions for saved-pool coverage; no docking or campaign reads."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def load_script(name):
    spec = importlib.util.spec_from_file_location(
        f"_pose_budget_test_{name}", SCRIPTS / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def elector(monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPTS))
    return load_script("offline_elector")


def pose_file(prefix, suffix, cf=-1.0):
    path = Path(f"{prefix}_{suffix}.pdb")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"REMARK CF= {cf}\nREMARK Frequency: 1\n"
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    )
    return path


def write_ranks(prefix, count):
    for rank in range(count):
        pose_file(prefix, rank, cf=-float(rank + 1))


@pytest.fixture
def scoring_provider(tmp_path):
    """A fixed deterministic stand-in for the lane's external RMSD provider.

    The numeric CF in each synthetic file encodes its oracle value. This
    exercises real CLI/provider/pool wiring without loading RDKit or spyRMSD.
    """
    path = tmp_path / "scoring" / "score_electors.py"
    path.parent.mkdir()
    path.write_text(textwrap.dedent("""\
        from pathlib import Path

        class Reference:
            def __init__(self, path):
                assert Path(path).is_file()
                self.n = 1

        def pool_tiers(poses, ref, score_fn, top_n=10):
            assert ref.n == 1
            if not poses:
                return {'oracle': None, 'oracle_naive': None,
                        'top_n': None, 'top_n_naive': None, 'n_pool': 0}
            ordered = sorted(poses, key=score_fn)
            oracle = min(1000.0 + p['cf'] for p in poses)
            top = min(1000.0 + p['cf'] for p in ordered[:top_n])
            return {'oracle': oracle, 'oracle_naive': oracle,
                    'top_n': top, 'top_n_naive': top, 'n_pool': len(poses)}
        """))
    return path


def test_explicit_budget_reaches_a_decisive_pose_after_rank_49(elector, tmp_path):
    prefix = tmp_path / "TEST"
    write_ranks(prefix, 251)

    old_heads, old_stats = elector.enumerate_heads(str(prefix))
    wide_heads, wide_stats = elector.enumerate_heads(str(prefix), budget=250)
    assert len(old_heads) == 50
    assert old_heads[-1][1] == 49
    assert len(wide_heads) == 250
    assert wide_heads[-1][1] == 249
    assert old_stats["cf_truncated"] is True
    assert wide_stats["cf_truncated"] is True

    old_pool, _ = elector.build_pool([str(prefix)])
    wide_pool, _ = elector.build_pool([str(prefix)], budget=250)
    old_winner, _ = elector.elect_mincf(old_pool)
    wide_winner, _ = elector.elect_mincf(wide_pool)
    assert Path(old_winner["path"]).name == "TEST_49.pdb"
    assert Path(wide_winner["path"]).name == "TEST_249.pdb"


def test_budget_is_per_restart_and_reaches_build_pool(elector, tmp_path):
    prefixes = [str(tmp_path / f"r{restart}" / "TEST") for restart in range(3)]
    for prefix in prefixes:
        write_ranks(prefix, 80)

    default_pool, default_stats = elector.build_pool(prefixes)
    wide_pool, wide_stats = elector.build_pool(prefixes, budget=250)
    assert len(default_pool) == 150
    assert len(wide_pool) == 240
    assert [p["restart"] for p in wide_pool] == [0] * 80 + [1] * 80 + [2] * 80
    assert [s["budget"] for s in default_stats] == [50, 50, 50]
    assert [s["budget"] for s in wide_stats] == [250, 250, 250]


def test_default_50_does_not_inherit_live_runner_environment(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(SCRIPTS))
    monkeypatch.setenv("FLEXAIDDS_MAX_RESULTS", "250")
    module = load_script("offline_elector")
    prefix = tmp_path / "TEST"
    write_ranks(prefix, 60)
    heads, _ = module.enumerate_heads(str(prefix))
    pool, _ = module.build_pool([str(prefix)])
    assert len(heads) == len(pool) == 50


def test_cf_and_fastoptics_share_budget_and_preserve_order(elector, tmp_path):
    prefix = tmp_path / "TEST"
    for suffix in ["2_2", "10_0", "0", "2_0", "2"]:
        pose_file(prefix, suffix)
    for invalid in ["INI", "2_0_extra", "-1", "2_-1", "two_1"]:
        pose_file(prefix, invalid, cf=-999.0)
    pose_file(tmp_path / "OTHER", "0", cf=-999.0)

    heads, stats = elector.enumerate_heads(str(prefix), budget=4)
    assert [Path(p).name for p, _, _ in heads] == [
        "TEST_0.pdb", "TEST_2.pdb", "TEST_2_0.pdb", "TEST_2_2.pdb"
    ]
    assert [(rank, min_pts) for _, rank, min_pts in heads] == [
        (0, -1), (2, -1), (0, 2), (2, 2)
    ]
    assert stats["cf_found"] == 2
    assert stats["fo_found"] == 3
    assert stats["fo_kept"] == 2
    assert stats["fo_truncated"] is True


def test_ini_and_published_elected_copy_are_not_restart_heads(elector, tmp_path):
    cell = tmp_path / "CELL"
    target_dir = cell / "TEST"
    prefixes = elector.cell_prefixes(str(target_dir), "TEST", n_restarts=3)
    for prefix in prefixes:
        pose_file(prefix, "0", cf=-1.0)
        pose_file(prefix, "INI", cf=-999.0)
    published = pose_file(cell / "TEST", "0", cf=-1000.0)

    pool, _ = elector.build_pool(prefixes, budget=250)
    assert len(pool) == 3
    assert all(p["cf"] == -1.0 for p in pool)
    assert str(published) not in [p["path"] for p in pool]


@pytest.mark.parametrize("budget", [0, -1, 5001, True, False, 50.0, "50", None])
@pytest.mark.parametrize("api", ["enumerate_heads", "build_pool"])
def test_invalid_budgets_fail_even_on_empty_input(elector, tmp_path, budget, api):
    arg = str(tmp_path / "TEST") if api == "enumerate_heads" else []
    with pytest.raises((TypeError, ValueError)):
        getattr(elector, api)(arg, budget=budget)


@pytest.mark.parametrize("budget", [1, 5000])
def test_budget_range_endpoints_are_accepted(elector, tmp_path, budget):
    prefix = tmp_path / "TEST"
    pose_file(prefix, "0")
    pool, stats = elector.build_pool([str(prefix)], budget=budget)
    assert len(pool) == 1
    assert stats[0]["budget"] == budget


@pytest.mark.parametrize("bad_headers", [False, True])
def test_cli_scores_the_same_files_at_both_budgets(tmp_path, scoring_provider, bad_headers):
    prefix = tmp_path / "poses" / "TEST"
    write_ranks(prefix, 80)
    # Deliberate parse failures: raw enumeration must not be relabeled as
    # successfully parsed or scored coverage.
    if bad_headers:
        pose_file(prefix, "2", cf="nan")
        pose_file(prefix, "3", cf="not-a-number")
    reference = tmp_path / "reference.sdf"
    reference.write_text("synthetic provider fixture\n")
    environment = dict(os.environ, FLEXAIDDS_MAX_RESULTS="5000")
    reports = {}
    for budget in [50, 250]:
        output = tmp_path / f"budget_{budget}.json"
        command = [
            sys.executable, str(SCRIPTS / "score_emitted_pools.py"),
            "--prefix", str(prefix), "--reference", str(reference),
            "--scoring-code", str(scoring_provider),
            "--pose-budget", str(budget), "--out", str(output),
        ]
        result = subprocess.run(command, capture_output=True, text=True, env=environment)
        assert result.returncode == (2 if bad_headers else 0), result.stdout + result.stderr
        reports[budget] = json.loads(output.read_text())

    assert reports[50]["pose_budget"] == 50
    assert reports[250]["pose_budget"] == 250
    assert reports[50]["n_enumerated"] == 50
    assert reports[250]["n_enumerated"] == 80
    failures = 2 if bad_headers else 0
    assert reports[50]["n_parsed"] == reports[50]["n_pool_scored"] == 50 - failures
    assert reports[250]["n_parsed"] == reports[250]["n_pool_scored"] == 80 - failures
    for report in reports.values():
        assert report["n_parse_failed"] == failures
        assert report["status"] == ("parse_failures" if bad_headers else "scored")
        assert report["n_enumerated"] == report["n_parsed"] + report["n_parse_failed"]
        assert len(report["enum_stats"][0]["parse_errors"]) == failures
    assert reports[50]["oracle_spyrmsd_A"] == 950.0
    assert reports[250]["oracle_spyrmsd_A"] == 920.0
    provider_hash = hashlib.sha256(scoring_provider.read_bytes()).hexdigest()
    assert reports[50]["score_code_sha256"] == provider_hash
    assert reports[250]["score_code_sha256"] == provider_hash
    assert reports[250]["enum_stats"][0]["budget"] == 250
    assert not (scoring_provider.parent / "__pycache__").exists()


def test_cli_refuses_to_overwrite_existing_receipt(tmp_path, scoring_provider):
    prefix = tmp_path / "TEST"
    pose_file(prefix, "0")
    reference = tmp_path / "reference.sdf"
    reference.write_text("synthetic provider fixture\n")
    output = tmp_path / "receipt.json"
    original = b'{"existing": "preserve exactly"}\n'
    output.write_bytes(original)
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / "score_emitted_pools.py"),
         "--prefix", str(prefix), "--reference", str(reference),
         "--scoring-code", str(scoring_provider), "--out", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "already exists" in result.stderr
    assert output.read_bytes() == original


def test_adapter_rejects_duplicate_restart_prefixes(monkeypatch, tmp_path, scoring_provider):
    monkeypatch.syspath_prepend(str(SCRIPTS))
    adapter = load_script("score_emitted_pools")
    prefix = tmp_path / "TEST"
    reference = tmp_path / "reference.sdf"
    reference.write_text("synthetic provider fixture\n")
    with pytest.raises(ValueError, match="unique"):
        adapter.score_saved_pool(
            [str(prefix), str(prefix.parent / "." / prefix.name)],
            reference, scoring_provider, pose_budget=250,
        )


def test_adapter_rejects_provider_with_missing_scoring_api(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(SCRIPTS))
    adapter = load_script("score_emitted_pools")
    provider = tmp_path / "incomplete_provider.py"
    provider.write_text("class Reference: pass\n")
    with pytest.raises(ValueError, match="Reference and pool_tiers"):
        adapter.load_scoring_code(provider)


@pytest.mark.parametrize("case", ["missing_all", "missing_one_restart", "provider_drops_all"])
def test_cli_does_not_report_incomplete_or_empty_scoring_as_success(
    tmp_path, scoring_provider, case
):
    prefixes = [tmp_path / f"r{restart}" / "TEST" for restart in range(3)]
    if case == "missing_one_restart":
        pose_file(prefixes[0], "0", cf=-2.0)
        pose_file(prefixes[2], "0", cf=-3.0)
    elif case == "provider_drops_all":
        for prefix in prefixes:
            pose_file(prefix, "0")
        # Mirror a legitimate scorer's empty result after filtering or a
        # ligand-matching failure; CF-header parsing itself still succeeds.
        with scoring_provider.open("a") as handle:
            handle.write(textwrap.dedent("""\

                def pool_tiers(poses, ref, score_fn, top_n=10):
                    assert len(poses) == 3
                    return {'oracle': None, 'oracle_naive': None,
                            'top_n': None, 'top_n_naive': None, 'n_pool': 0}
                """))
    reference = tmp_path / "reference.sdf"
    reference.write_text("synthetic provider fixture\n")
    output = tmp_path / "incomplete_receipt.json"
    command = [sys.executable, str(SCRIPTS / "score_emitted_pools.py")]
    for prefix in prefixes:
        command.extend(["--prefix", str(prefix)])
    command.extend([
        "--reference", str(reference), "--scoring-code", str(scoring_provider),
        "--pose-budget", "250", "--out", str(output),
    ])
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 2, result.stdout + result.stderr
    report = json.loads(output.read_text())
    assert report["n_parse_failed"] == 0
    assert report["n_enumerated"] == report["n_parsed"]
    if case == "missing_all":
        assert report["status"] == "missing_or_empty_prefixes"
        assert report["empty_prefixes"] == [str(p.resolve()) for p in prefixes]
        assert report["n_enumerated"] == report["n_pool_scored"] == 0
        assert report["oracle_spyrmsd_A"] is None
    elif case == "missing_one_restart":
        assert report["status"] == "missing_or_empty_prefixes"
        assert report["empty_prefixes"] == [str(prefixes[1].resolve())]
        assert report["n_enumerated"] == report["n_pool_scored"] == 2
        assert report["oracle_spyrmsd_A"] == 997.0
        assert [s["enumerated"] for s in report["enum_stats"]] == [1, 0, 1]
    else:
        assert report["status"] == "empty_pool"
        assert report["empty_prefixes"] == []
        assert report["n_enumerated"] == 3
        assert report["n_pool_scored"] == 0
        assert report["n_not_in_scored_pool"] == 3
        assert report["oracle_spyrmsd_A"] is None
