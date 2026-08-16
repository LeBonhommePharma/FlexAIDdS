"""Lock the coverage.yml gate: denominator, ratchet, and comment status contract.

The job gate lives in the Calculate bash step. The PR comment must never
invent a threshold or report PASS/FAIL when it cannot judge.
"""
from __future__ import annotations

import math
import os
import re
import stat
import subprocess
import textwrap
from pathlib import Path

import pytest

try:
    import yaml
except ImportError:
    yaml = None

REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / ".github" / "workflows" / "coverage.yml"


def _workflow() -> str:
    return WORKFLOW.read_text(encoding="utf-8")


def _block_after(text: str, heading: str, key: str) -> str:
    """Return the unindented `|` scalar under `key:` after a step heading."""
    hpos = text.find(heading)
    assert hpos >= 0, f"missing heading {heading!r}"
    kpos = text.find(f"{key}: |", hpos)
    assert kpos >= 0, f"missing {key}: | after {heading!r}"
    after = text[kpos + len(f"{key}: |") :]
    lines = after.splitlines()
    if lines and lines[0] == "":
        lines = lines[1:]
    body: list[str] = []
    indent = None
    for line in lines:
        if not line.strip():
            body.append("")
            continue
        leading = len(line) - len(line.lstrip(" "))
        if indent is None:
            indent = leading
        if leading < indent:
            break
        body.append(line[indent:])
    return "\n".join(body).rstrip() + "\n"


def _synthetic_info(path: Path, hit: int, found: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        textwrap.dedent(
            f"""\
            TN:synthetic
            SF:LIB/Vcontacts.cpp
            DA:1,1
            LH:{hit}
            LF:{found}
            end_of_record
            """
        ),
        encoding="utf-8",
    )


def _run_calculate(tmp_path: Path, info: Path | None) -> subprocess.CompletedProcess[str]:
    script = _block_after(_workflow(), "Calculate coverage percentage", "run")
    work = tmp_path / "work"
    work.mkdir()
    if info is not None:
        dest = work / "build" / "coverage" / "core.info"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(info.read_bytes())
    env_file = tmp_path / "github.env"
    summary = tmp_path / "summary.md"
    env_file.write_text("", encoding="utf-8")
    summary.write_text("", encoding="utf-8")
    runner = work / "calculate.sh"
    runner.write_text("#!/usr/bin/env bash\n" + script, encoding="utf-8")
    runner.chmod(runner.stat().st_mode | stat.S_IEXEC)
    return subprocess.run(
        ["bash", str(runner)],
        cwd=work,
        env={
            **os.environ,
            "GITHUB_ENV": str(env_file),
            "GITHUB_STEP_SUMMARY": str(summary),
        },
        capture_output=True,
        text=True,
        timeout=15,
        check=False,
    )


def _parse_github_env(tmp_path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (tmp_path / "github.env").read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k] = v
    return out


def comment_status(coverage: str | None, threshold: str | None) -> tuple[str, str]:
    """Python port of the coverage.yml github-script status contract."""
    try:
        threshold_n = (
            float("nan")
            if threshold is None or threshold == ""
            else float(threshold)
        )
    except ValueError:
        threshold_n = float("nan")
    available = coverage is not None and coverage != "" and coverage != "unavailable"
    try:
        pct = float("nan") if coverage is None or coverage == "" else float(coverage)
    except ValueError:
        pct = float("nan")
    # Match coverage.yml github-script: Number.isFinite rejects NaN and ±Infinity.
    # math.isnan lets inf through, which would false-PASS comment_status("inf", "45").
    t_ok = math.isfinite(threshold_n)
    pct_ok = math.isfinite(pct)
    can_judge = available and pct_ok and t_ok
    passed = can_judge and pct >= threshold_n
    if t_ok:
        shown = int(threshold_n) if threshold_n == int(threshold_n) else threshold_n
        target = f">= {shown}%"
    else:
        target = ">= ?%"
    status = "UNAVAILABLE" if not can_judge else ("PASS" if passed else "FAIL")
    return status, target


def test_yaml_parses_and_paths_filter_skips_github():
    if yaml is None:
        pytest.skip("pyyaml not installed")
    data = yaml.safe_load(_workflow())
    assert data["name"] == "Code Coverage"
    # PyYAML 1.1 treats the key `on` as boolean True.
    on = data.get("on", data.get(True))
    pr_paths = on["pull_request"]["paths"]
    assert ".github/**" not in pr_paths
    assert "LIB/**" in pr_paths


def test_denominator_includes_dataset_runner_excludes_only_top():
    text = _workflow()
    assert not re.search(r"--exclude\s+'?\*/DatasetRunner\.cpp'?", text)
    assert not re.search(r"--remove[\s\S]{0,120}'?\*/DatasetRunner\.cpp'?", text)
    assert "'*/top.cpp'" in text
    assert "--exclude '*/benchmarks/*'" not in text
    assert "tests/test_dataset_runner.cpp" in text


def test_dataset_runner_gtest_target_still_registered():
    """Ratio gate is inverted for under-covered files: deleting this target
    would drop DatasetRunner.cpp LF/LH records and raise the percentage.
    Keep the executable and ctest name in CMakeLists.txt.
    """
    cmake = (REPO / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "add_executable(test_dataset_runner" in cmake
    assert "tests/test_dataset_runner.cpp" in cmake
    assert "add_test(NAME DatasetRunnerTests COMMAND test_dataset_runner)" in cmake
    assert (REPO / "tests" / "test_dataset_runner.cpp").is_file()


def test_single_threshold_source_is_45():
    text = _workflow()
    calc = _block_after(text, "Calculate coverage percentage", "run")
    comment = _block_after(text, "Comment PR with coverage", "script")
    assert re.search(r"^THRESHOLD=45\s*$", calc, re.M)
    assert "echo \"THRESHOLD=${THRESHOLD}\" >> \"$GITHUB_ENV\"" in calc
    assert "pct >= 50" not in comment
    assert ">= 50%" not in comment
    assert "Number(process.env.THRESHOLD)" not in comment


def test_calculate_script_bash_n(tmp_path: Path):
    script = _block_after(_workflow(), "Calculate coverage percentage", "run")
    path = tmp_path / "calculate.sh"
    path.write_text("#!/usr/bin/env bash\n" + script, encoding="utf-8")
    proc = subprocess.run(
        ["bash", "-n", str(path)],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_comment_js_treats_empty_threshold_as_unknown():
    comment = _block_after(_workflow(), "Comment PR with coverage", "script")
    assert "rawT === undefined || rawT === ''" in comment
    assert "canJudge" in comment
    assert "!canJudge ? 'UNAVAILABLE'" in comment


def test_calculate_pass_at_live_capture(tmp_path: Path):
    info = tmp_path / "core.info"
    _synthetic_info(info, hit=21804, found=44046)
    proc = _run_calculate(tmp_path, info)
    assert proc.returncode == 0, proc.stderr + proc.stdout
    env = _parse_github_env(tmp_path)
    assert env["THRESHOLD"] == "45"
    assert env["COVERAGE"] == "49.5"
    assert "meets threshold 45%" in proc.stdout


def test_calculate_fail_below_ratchet(tmp_path: Path):
    info = tmp_path / "core.info"
    _synthetic_info(info, hit=4400, found=10000)
    proc = _run_calculate(tmp_path, info)
    assert proc.returncode == 1, proc.stdout
    env = _parse_github_env(tmp_path)
    assert env["THRESHOLD"] == "45"
    assert env["COVERAGE"] == "44.0"
    assert "below threshold 45%" in proc.stdout


def test_calculate_empty_tracefile_unavailable_still_exports_threshold(
    tmp_path: Path,
):
    proc = _run_calculate(tmp_path, info=None)
    assert proc.returncode == 0, proc.stderr
    env = _parse_github_env(tmp_path)
    assert env["THRESHOLD"] == "45"
    assert env["COVERAGE"] == "unavailable"


@pytest.mark.parametrize(
    "coverage,threshold,status,target",
    [
        ("49.5", "45", "PASS", ">= 45%"),
        ("47.4", "45", "PASS", ">= 45%"),
        ("45.0", "45", "PASS", ">= 45%"),
        ("44.0", "45", "FAIL", ">= 45%"),
        ("unavailable", "45", "UNAVAILABLE", ">= 45%"),
        (None, None, "UNAVAILABLE", ">= ?%"),
        ("47.4", None, "UNAVAILABLE", ">= ?%"),
        ("47.4", "", "UNAVAILABLE", ">= ?%"),
        ("inf", "45", "UNAVAILABLE", ">= 45%"),
        ("Infinity", "45", "UNAVAILABLE", ">= 45%"),
        ("-inf", "45", "UNAVAILABLE", ">= 45%"),
        ("49.5", "inf", "UNAVAILABLE", ">= ?%"),
        ("49.5", "Infinity", "UNAVAILABLE", ">= ?%"),
    ],
)
def test_comment_status_contract(coverage, threshold, status, target):
    got_status, got_target = comment_status(coverage, threshold)
    assert got_status == status
    assert got_target == target
