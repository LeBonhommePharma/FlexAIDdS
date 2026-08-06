"""Early-termination tracking: the GA's exit path must reach the artifact.

The runner captures engine stdout but used to scan it only for ``[GRAND]`` and
throw the rest away.  Every early-exit path announces itself there, so a run that
stopped at 14% of its configured generation budget landed in the results JSON
looking exactly like one that ran the full 2000 generations -- same schema, same
fields, no trace of the truncation.  Its mean RMSD then gets compared against
full-budget baselines as though the two came from the same amount of search.

These tests pin three things:

* every exit message the engine can actually print is recognised, with its
  governing numbers extracted (the regexes are checked against the literal
  ``printf`` formats in ``LIB/gaboom.cpp``, not against paraphrases);
* the governing knobs are recorded from the environment the *subprocess* saw,
  including the one that matters -- ``FLEXAIDDS_ADAPTIVE_GENERATIONS`` is not
  covered by ``FLEXAIDDS_NO_SEC`` (gaboom.cpp:927), so the runner's
  "full budget is always consumed" guarantee does not hold when it is set;
* absence of a record is never synthesised into "ran clean" -- a dry run and a
  full-budget run must not produce the same evidence.
"""

from __future__ import annotations

import pathlib
import re

import pytest

from flexaidds.dataset_runner.runner import (
    DatasetConfig,
    DatasetResult,
    DatasetRunner,
    TargetResult,
    early_exit_parameters,
    parse_early_termination,
)

REPO = pathlib.Path(__file__).resolve().parents[2]
GABOOM = REPO / "LIB" / "gaboom.cpp"

# Verbatim renderings of the five printf() sites in LIB/gaboom.cpp, with the
# format specifiers filled in. If the engine changes its wording these strings
# must change with it -- that is the point of testing against them.
ADAPTIVE = (
    "[P5-ADAPTIVE-GEN] GA converged: best-CF plateau for 50 gens "
    "(best_CF=-123.4500) at gen 279 — early stop (max_generations=2000)\n"
)
H_PLATEAU = (
    "Early exit at gen 812: H plateau < 0.0010 "
    "(H_now=0.863000 nats, delta=0.000400 nats)\n"
)
CF_STAGNANT = (
    "GA terminated: CF stagnant for 300 gens (best_CF=-98.7600) "
    "with gene-space collapsed\n"
)
ENTROPY = "GA terminated early by entropy convergence\n"
FITNESS = "GA terminated early by fitness stagnation\n"


# ---------------------------------------------------------------------------
# parse_early_termination
# ---------------------------------------------------------------------------


def test_clean_run_is_not_flagged():
    out = parse_early_termination("Generation:     1\nGeneration:  2000\ndone\n")
    assert out["terminated_early"] is False
    assert out["reason"] is None
    assert out["last_generation"] == 2000


def test_empty_stdout_is_not_flagged():
    assert parse_early_termination("")["terminated_early"] is False
    assert parse_early_termination(None)["terminated_early"] is False  # type: ignore[arg-type]


def test_adaptive_plateau_is_parsed():
    out = parse_early_termination("Generation:   279\n" + ADAPTIVE)
    assert out["terminated_early"] is True
    assert out["reason"] == "adaptive_cf_plateau"
    # The numbers are the whole reason to record the event: "stopped at 279 of
    # 2000" is actionable, "stopped early" is not.
    assert out["generation"] == 279
    assert out["max_generations"] == 2000
    assert out["plateau"] == 50
    assert out["best_cf"] == pytest.approx(-123.45)


def test_h_plateau_is_parsed():
    out = parse_early_termination(H_PLATEAU)
    assert out["reason"] == "h_plateau"
    assert out["generation"] == 812
    assert out["eps"] == pytest.approx(0.001)
    assert out["h_now"] == pytest.approx(0.863)
    assert out["delta"] == pytest.approx(0.0004)


def test_cf_stagnation_is_parsed():
    out = parse_early_termination(CF_STAGNANT)
    assert out["reason"] == "cf_stagnant_gene_collapsed"
    assert out["plateau"] == 300
    assert out["best_cf"] == pytest.approx(-98.76)


@pytest.mark.parametrize(
    "text,reason",
    [(ENTROPY, "entropy_convergence"), (FITNESS, "fitness_stagnation")],
)
def test_bare_termination_messages_are_parsed(text, reason):
    out = parse_early_termination(text)
    assert out["terminated_early"] is True
    assert out["reason"] == reason


def test_last_generation_takes_the_final_occurrence():
    """Progress lines repeat; the truncation point is the last one seen."""
    text = "".join(f"Generation: {i:5d}\n" for i in (1, 2, 3, 154))
    assert parse_early_termination(text)["last_generation"] == 154


def test_termination_is_found_amid_ordinary_engine_chatter():
    noisy = (
        "Generation:   279\n"
        "   7 (   2046.00      10.00 )  cf=-111058.302 fitnes=    2.627\n"
        + ADAPTIVE
        + "[GRAND] log_Z=12.5\n"
    )
    assert parse_early_termination(noisy)["reason"] == "adaptive_cf_plateau"


@pytest.mark.skipif(not GABOOM.is_file(), reason="engine source not present")
@pytest.mark.parametrize(
    "needle",
    [
        "[P5-ADAPTIVE-GEN] GA converged: best-CF plateau for %d ",
        "Early exit at gen %d: H plateau < %.4f ",
        "GA terminated: CF stagnant for %d gens (best_CF=%.4f) ",
        "GA terminated early by entropy convergence",
        "GA terminated early by fitness stagnation",
    ],
)
def test_patterns_track_the_real_engine_messages(needle):
    """Guard against the parser silently going deaf.

    A reworded printf() would make every regex above miss, and a missed exit
    reads as a clean full-budget run -- failing open, in the one direction that
    matters. Asserting the literal format strings still exist makes that a test
    failure instead of a quiet loss of evidence.
    """
    assert needle in GABOOM.read_text(errors="replace")


# ---------------------------------------------------------------------------
# early_exit_parameters
# ---------------------------------------------------------------------------


def test_parameters_come_from_the_subprocess_env_not_the_parent():
    """The engine sees ``sub_env``; recording os.environ would misreport it."""
    got = early_exit_parameters({"FLEXAIDDS_NO_SEC": "1", "FLEXAIDDS_BENCHMARK": "1"})
    assert got["FLEXAIDDS_NO_SEC"] == "1"
    assert got["FLEXAIDDS_BENCHMARK"] == "1"
    assert got["FLEXAIDDS_ADAPTIVE_GENERATIONS"] is None
    assert got["adaptive_bypasses_no_sec"] is False


def test_adaptive_generations_is_flagged_as_bypassing_no_sec():
    """The finding this tracking exists for.

    ``FLEXAIDDS_NO_SEC=1`` is what the runner relies on to guarantee the full
    generation budget, and it does cover the stagnation / entropy / SEC /
    H-plateau exits. It does NOT cover the adaptive best-CF plateau exit, which
    is guarded only by ``ag_patience > 0``. So NO_SEC=1 and ADAPTIVE=50 together
    still produce a truncated run, and the artifact has to say so.
    """
    got = early_exit_parameters(
        {"FLEXAIDDS_NO_SEC": "1", "FLEXAIDDS_ADAPTIVE_GENERATIONS": "50",
         "FLEXAIDDS_ADAPTIVE_EPS": "1.0"}
    )
    assert got["adaptive_bypasses_no_sec"] is True
    assert got["FLEXAIDDS_ADAPTIVE_GENERATIONS"] == "50"
    assert got["FLEXAIDDS_ADAPTIVE_EPS"] == "1.0"


@pytest.mark.parametrize("value,expected", [("0", False), ("", False), ("  ", False),
                                            ("1", True), ("50", True)])
def test_adaptive_disabled_values_do_not_raise_the_flag(value, expected):
    got = early_exit_parameters({"FLEXAIDDS_ADAPTIVE_GENERATIONS": value})
    assert got["adaptive_bypasses_no_sec"] is expected


@pytest.mark.skipif(not GABOOM.is_file(), reason="engine source not present")
def test_recorded_compiled_thresholds_match_the_engine():
    """The compiled-in constants are not configurable, so the artifact carries
    them as literals. Literals drift; pin them to the source."""
    src = GABOOM.read_text(errors="replace")
    params = early_exit_parameters({})

    def _const(pattern: str) -> str:
        m = re.search(pattern, src)
        assert m, f"constant not found in gaboom.cpp: {pattern}"
        return m.group(1)

    assert int(_const(r"STAGNATION_LIMIT\s*=\s*(\d+)")) == params["compiled_stagnation_limit"]
    assert int(_const(r"kHPlateauWindow\s*=\s*(\d+)")) == params["compiled_h_plateau_window"]
    assert float(_const(r"kHPlateauEps\s*=\s*([\d.eE+-]+)")) == pytest.approx(
        params["compiled_h_plateau_eps_nats"]
    )

    ga_const = REPO / "LIB" / "ga_constants.h"
    if ga_const.is_file():
        m = re.search(r"GA_DEFAULT_ENTROPY_CHECK_INTERVAL\s*=?\s*(\d+)", ga_const.read_text())
        if m:
            assert int(m.group(1)) == params["compiled_entropy_check_interval"]


# ---------------------------------------------------------------------------
# Runner plumbing: _run_flexaid records, run_dataset consumes
# ---------------------------------------------------------------------------


def _runner(tmp_path) -> DatasetRunner:
    return DatasetRunner(results_dir=tmp_path / "results", dry_run=True, n_workers=1)


def test_no_record_is_not_synthesised_into_a_clean_run(tmp_path):
    """An entry with no evidence must stay empty.

    Dry runs and exec failures never reach the stdout parse. Filling in a
    default "terminated_early: False" there would assert something the run never
    established -- the same conflation, one layer up.
    """
    r = _runner(tmp_path)
    got = r._early_termination_for("1gpk", "holo")
    assert got == {"early_termination": {}, "early_exit_params": {}}


def test_recorded_evidence_reaches_the_target_result_kwargs(tmp_path):
    r = _runner(tmp_path)
    key = "1gpk/holo"
    r._entry_early_exit_params[key] = early_exit_parameters({"FLEXAIDDS_NO_SEC": "1"})
    r._record_early_termination(key, "lig1", parse_early_termination(ADAPTIVE))

    got = r._early_termination_for("1gpk", "holo")
    assert got["early_termination"]["terminated_early"] is True
    assert got["early_termination"]["reason"] == "adaptive_cf_plateau"
    assert got["early_exit_params"]["FLEXAIDDS_NO_SEC"] == "1"
    # TargetResult accepts them as-is (the call site splats this dict).
    tr = TargetResult(target_id="1gpk", structural_state="holo", poses=[], **got)
    assert tr.early_termination["generation"] == 279


def test_one_truncated_ligand_marks_the_whole_entry(tmp_path):
    """Scores are pooled across an entry's ligands, so one truncated ligand is
    enough to make the entry's numbers incomparable."""
    r = _runner(tmp_path)
    key = "1gpk/holo"
    r._record_early_termination(key, "ligA", parse_early_termination("Generation:  2000\n"))
    r._record_early_termination(key, "ligB", parse_early_termination(ADAPTIVE))

    rec = r._entry_early_termination[key]
    assert rec["terminated_early"] is True
    assert rec["reason"] == "adaptive_cf_plateau"
    # ...and which ligand it was stays recoverable.
    assert rec["per_ligand"]["ligA"]["terminated_early"] is False
    assert rec["per_ligand"]["ligB"]["reason"] == "adaptive_cf_plateau"


def test_entry_record_keeps_the_furthest_generation_seen(tmp_path):
    r = _runner(tmp_path)
    key = "t/holo"
    r._record_early_termination(key, "a", parse_early_termination("Generation:   154\n"))
    r._record_early_termination(key, "b", parse_early_termination("Generation:  1900\n"))
    r._record_early_termination(key, "c", parse_early_termination("Generation:   300\n"))
    assert r._entry_early_termination[key]["last_generation"] == 1900


def test_entries_are_tracked_independently(tmp_path):
    r = _runner(tmp_path)
    r._record_early_termination("a/holo", "l", parse_early_termination(ADAPTIVE))
    r._record_early_termination("b/holo", "l", parse_early_termination("Generation: 2000\n"))
    assert r._early_termination_for("a", "holo")["early_termination"]["terminated_early"]
    assert not r._early_termination_for("b", "holo")["early_termination"]["terminated_early"]


# ---------------------------------------------------------------------------
# Persistence -- the evidence has to survive the round trip
# ---------------------------------------------------------------------------


def test_evidence_survives_save_and_resume(tmp_path):
    """--resume reloads entries from JSON. Dropping the record there would put
    the truncation back out of sight for exactly the runs that were interrupted.
    """
    r = _runner(tmp_path)
    cfg = DatasetConfig(slug="s", name="n", description="d", targets=["1gpk"])
    tr = TargetResult(
        target_id="1gpk",
        structural_state="holo",
        poses=[],
        early_termination=parse_early_termination(ADAPTIVE),
        early_exit_params=early_exit_parameters(
            {"FLEXAIDDS_NO_SEC": "1", "FLEXAIDDS_ADAPTIVE_GENERATIONS": "50"}
        ),
    )
    r._save_target_result(tr, cfg, tier=1)
    back = r._load_target_result(cfg, 1, "1gpk", "holo")

    assert back is not None
    assert back.early_termination["reason"] == "adaptive_cf_plateau"
    assert back.early_termination["generation"] == 279
    assert back.early_exit_params["adaptive_bypasses_no_sec"] is True


def test_dataset_report_exposes_the_rollup():
    dr = DatasetResult(
        config=DatasetConfig(slug="s", name="n", description="d", targets=[]), tier=1
    )
    dr.early_terminations = {"1gpk/holo": parse_early_termination(ADAPTIVE)}
    dr.early_exit_params = early_exit_parameters({"FLEXAIDDS_ADAPTIVE_GENERATIONS": "50"})
    payload = dr.to_dict()
    assert payload["early_terminations"]["1gpk/holo"]["reason"] == "adaptive_cf_plateau"
    assert payload["early_exit_params"]["adaptive_bypasses_no_sec"] is True


def test_dataset_report_rollup_defaults_to_empty():
    dr = DatasetResult(
        config=DatasetConfig(slug="s", name="n", description="d", targets=[]), tier=1
    )
    payload = dr.to_dict()
    assert payload["early_terminations"] == {}
    assert payload["early_exit_params"] == {}
