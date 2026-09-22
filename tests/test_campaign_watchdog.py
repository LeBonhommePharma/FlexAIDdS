"""Guards for scripts/campaign_watchdog.py and scripts/campaign_preflight.sh.

Every fixture in this file mirrors a shape measured on the real tree at
/Users/lp.more/flexaidds_results/campaigns/wall_paired_85_seed2 -- the nested
DatasetRunner layout, the rN restart directories, the real healthy stderr
lines and the real fatal line from the 2026-09-20 TMPDIR incident. Fixtures
are built under pytest's tmp_path. Nothing here reads or writes the live
campaign tree.

The detectors are only worth their exit codes if they also stay QUIET, so
roughly half of these tests assert silence on shapes that superficially look
like failure -- a flat file count, a single hard target, a recovered burst.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import campaign_watchdog as W  # noqa: E402

PREFLIGHT = REPO_ROOT / "scripts" / "campaign_preflight.sh"

# --- real text, copied from the live tree -----------------------------------

# From a healthy cell's stderr.log. None of this is fatal, and the [WARN]
# line recurs across most targets -- it is the single most dangerous
# false-positive source for the signature detector.
HEALTHY_STDERR = """\
[DIRECT-IC] GPA topology atoms=90001,90013,90002 local=1,13,2
[GRID-CACHE] Saved 4453 grid points -> /out/1JD0/1JD0/1JD0.rrg
[NATIVE-SEED-RMSD] /cache/1JD0/1JD0_ligand.sdf round-trip RMSD = 2.97 A (5 genes, 13 atoms, gene0=0.000)
[FRAME_CHART] status=warn rmsd=2.965 strict=0 warn_A=1.0 strict_A=0.1
[WARN] NATIVE-SEED-RMSD = 2.97 A exceeds 1.0 A threshold.
"""

# The 2026-09-20 fatal line, verbatim apart from the workspace UUID.
FATAL_TMPDIR = (
    'Fatal error: filesystem error: in temp_directory_path: path '
    '"/Users/x/.claude-science/orgs/AAAA/workspaces/BBBB/.tmp" '
    "is not a directory: Not a directory\n"
)

RESULT_HEADER = "pdb_id,num_poses,best_score,wall_time_s,docking_exit_code,docking_completed,success\n"


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def make_cell(
    root: Path,
    target: str,
    *,
    arm: str = "c1_on",
    poses: int = 6,
    wall_s: float = 1800.0,
    exit_code: int = 0,
    stderr_text: str = HEALTHY_STDERR,
    restarts: int = 2,
    mtime: float | None = None,
    nested: bool = True,
    write_result_csv: bool = True,
    quarantine_suffix: str = "",
    stdout_bytes: int = 4096,
) -> Path:
    """Build one target cell in the real nested DatasetRunner shape.

    <root>/<arm>/<TARGET><suffix>/          outer runner dir (result.csv, run.log)
    <root>/<arm>/<TARGET><suffix>/<TARGET>/ inner dock cell (stdout/stderr/poses)
    <root>/<arm>/<TARGET><suffix>/<TARGET>/rN/   restart dirs
    """
    outer = root / arm / (target + quarantine_suffix)
    outer.mkdir(parents=True, exist_ok=True)
    (outer / "run.log").write_text("runner\n", encoding="utf-8")
    if nested:
        (outer / "result.csv").write_text(
            RESULT_HEADER + f"{target},{poses},-1.0,{wall_s},{exit_code},1,0\n",
            encoding="utf-8",
        )
        cell = outer / target
        cell.mkdir(parents=True, exist_ok=True)
    else:
        cell = outer

    cell.mkdir(parents=True, exist_ok=True)
    (cell / "stdout.log").write_text("x" * stdout_bytes, encoding="utf-8")
    (cell / "stderr.log").write_text(stderr_text, encoding="utf-8")
    (cell / "dock_config.json").write_text("{}", encoding="utf-8")

    if write_result_csv:
        (cell / "result.csv").write_text(
            RESULT_HEADER + f"{target},{poses},-1.0,{wall_s},{exit_code},1,0\n",
            encoding="utf-8",
        )

    # Pose files, plus the input frame that must NOT count as a pose.
    (cell / f"{target}_INI.pdb").write_text("INI\n", encoding="utf-8")
    for i in range(poses):
        (cell / f"{target}_{i}.pdb").write_text("ATOM\n", encoding="utf-8")
    if poses:
        (cell / "elected_pose.pdb").write_text("ATOM\n", encoding="utf-8")

    for r in range(1, restarts + 1):
        rd = cell / f"r{r}"
        rd.mkdir(exist_ok=True)
        (rd / "stdout.log").write_text("x" * 128, encoding="utf-8")
        (rd / "stderr.log").write_text(stderr_text, encoding="utf-8")

    if mtime is not None:
        for path in sorted(outer.rglob("*"), reverse=True):
            os.utime(path, (mtime, mtime))
        os.utime(outer, (mtime, mtime))
    return cell


def build_campaign(root: Path, specs: list[dict], base_time: float | None = None) -> Path:
    """Build a campaign whose cells complete in the order given.

    Completion order is set explicitly via mtime because the watchdog's
    "consecutive" logic is defined on completion order, and filesystem
    creation order is not a reliable proxy for it.
    """
    base_time = time.time() - 3600 if base_time is None else base_time
    for i, spec in enumerate(specs):
        spec = dict(spec)
        spec.setdefault("mtime", base_time + i * 60)
        make_cell(root, spec.pop("target"), **spec)
    return root


def fresh(root: Path) -> None:
    """Stamp the whole tree as just-modified, so the stall detector is quiet."""
    now = time.time()
    for path in sorted(root.rglob("*"), reverse=True):
        os.utime(path, (now, now))
    os.utime(root, (now, now))


def run_watchdog(root: Path, *args: str) -> tuple[int, str]:
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "campaign_watchdog.py"), str(root), *args],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr


def snapshot(root: Path) -> dict:
    """(size, mtime_ns) for every path. Deliberately excludes atime, which a
    read legitimately changes."""
    out = {}
    for path in sorted(root.rglob("*")):
        st = path.stat()
        out[str(path)] = (st.st_size, st.st_mtime_ns)
    return out


# ===========================================================================
# Cell discovery -- the nesting facts that made the first run miscount
# ===========================================================================


def test_nested_runner_layout_counts_each_target_once(tmp_path):
    """Outer runner dir and inner dock dir are ONE cell, not two.

    Before the deepest-wins rule the watchdog reported 267 cells for a
    97-target arm because both directories carry a result.csv.
    """
    build_campaign(tmp_path, [{"target": "1ABC"}, {"target": "2DEF"}])
    assert len(W.discover_cells(tmp_path)) == 2


def test_restart_dirs_are_not_cells(tmp_path):
    """rN directories carry stdout/stderr but belong to the cell above."""
    build_campaign(tmp_path, [{"target": "1ABC", "restarts": 4}])
    cells = W.discover_cells(tmp_path)
    assert len(cells) == 1
    assert not any(W.RESTART_DIR_RE.match(c.name) for c in cells)


def test_unnested_cell_still_discovered(tmp_path):
    """A cell that died before the inner dir existed is still seen."""
    make_cell(tmp_path, "1ABC", nested=False, poses=0, restarts=0)
    assert len(W.discover_cells(tmp_path)) == 1


# ===========================================================================
# Pose counting
# ===========================================================================


def test_ini_pdb_is_not_a_pose(tmp_path):
    """<target>_INI.pdb is the input frame and exists even in a dead cell."""
    cell = make_cell(tmp_path, "1ABC", poses=0, restarts=0, write_result_csv=False)
    assert (cell / "1ABC_INI.pdb").exists()
    count, source = W.count_poses(cell, None)
    assert (count, source) == (0, "pose_files")


def test_pose_count_falls_back_to_files(tmp_path):
    """No result.csv -> count poses on disk rather than skipping the cell."""
    cell = make_cell(tmp_path, "1ABC", poses=5, restarts=0, write_result_csv=False)
    count, source = W.count_poses(cell, None)
    assert source == "pose_files"
    assert count == 6  # 5 numbered poses + elected_pose.pdb


# ===========================================================================
# Detector: consecutive zero-pose cells
# ===========================================================================


def test_healthy_campaign_is_quiet(tmp_path):
    build_campaign(tmp_path, [{"target": t} for t in ("1ABC", "2DEF", "3GHI", "4JKL")])
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out
    assert "no detector fired" in out


def test_zero_pose_run_fires(tmp_path):
    """Three consecutive poseless cells at the tail -> systemic."""
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC"},
            {"target": "2DEF", "poses": 0, "exit_code": 1, "wall_s": 1.09},
            {"target": "3GHI", "poses": 0, "exit_code": 1, "wall_s": 1.09},
            {"target": "4JKL", "poses": 0, "exit_code": 1, "wall_s": 1.09},
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_SYSTEMIC, out
    assert "zero_pose_run" in out
    assert "3 consecutive cells" in out


def test_zero_pose_run_below_threshold_quiet(tmp_path):
    """Two poseless cells with a default threshold of three -> quiet.

    Hard targets exist. A watchdog that fires on the second failure would be
    switched off within a day.
    """
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC"},
            {"target": "2DEF", "poses": 0, "exit_code": 1},
            {"target": "3GHI", "poses": 0, "exit_code": 1},
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out


def test_recovered_burst_does_not_fire(tmp_path):
    """Poseless run followed by a producing cell is history, not an alarm."""
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC", "poses": 0, "exit_code": 1},
            {"target": "2DEF", "poses": 0, "exit_code": 1},
            {"target": "3GHI", "poses": 0, "exit_code": 1},
            {"target": "4JKL", "poses": 400},
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out
    assert "recovered" in out


# ===========================================================================
# Detector: recurring fatal signature
# ===========================================================================


def test_healthy_stderr_contains_no_fatal_line(tmp_path):
    """The real healthy stderr must not register as fatal.

    [WARN] NATIVE-SEED-RMSD recurs at nearly every target. If the classifier
    admitted it, every campaign ever run would report systemic failure.
    """
    for line in HEALTHY_STDERR.splitlines():
        assert not W.is_fatal_line(line), line


# The engine-hardening lane added LIB/fs_safe.h, whose detail::warn emits
#     [FS] <op> failed on '<path>' (<error>); <consequence>
# every time a filesystem call degrades instead of aborting. That line is the
# engine RECOVERING -- the exact outcome the 2026-09-20 fix was written to
# produce. Once it ships, a healthy campaign emits these at many targets, so a
# fatal classifier that matched on "failed" would convert the fix into a
# fleet-wide false alarm: same signature, many distinct targets, systemic.
FS_SAFE_RECOVERY = [
    "[FS] temp_directory_path failed on '/Users/x/.claude-science/orgs/A/workspaces/B/.tmp'"
    " (Not a directory); falling back to /tmp",
    "[FS] create_directories failed on '/out/1ABC' (Permission denied); continuing without scratch",
    "[FS] remove failed on '/tmp/probe.tmp' (No such file or directory); ignoring",
]


def test_fs_safe_recovery_lines_are_not_fatal():
    """A degraded-but-alive filesystem call must not read as a process death."""
    for line in FS_SAFE_RECOVERY:
        assert not W.is_fatal_line(line), line


def test_recovered_campaign_with_fs_warnings_stays_healthy(tmp_path):
    """End to end: every cell recovers and reports it, verdict is healthy."""
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC", "stderr_text": HEALTHY_STDERR + FS_SAFE_RECOVERY[0] + "\n"},
            {"target": "2DEF", "stderr_text": HEALTHY_STDERR + FS_SAFE_RECOVERY[0] + "\n"},
            {"target": "3GHI", "stderr_text": HEALTHY_STDERR + FS_SAFE_RECOVERY[0] + "\n"},
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out


def test_recurring_fatal_across_targets_fires(tmp_path):
    """The exact 2026-09-20 shape: one environment fault at many targets."""
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC"},
            {"target": "2DEF", "poses": 0, "exit_code": 1, "stderr_text": FATAL_TMPDIR},
            {"target": "3GHI", "poses": 0, "exit_code": 1, "stderr_text": FATAL_TMPDIR},
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_SYSTEMIC, out
    assert "recurring_fatal" in out
    assert "2 distinct targets" in out


def test_fatal_at_one_target_is_not_systemic(tmp_path):
    """One target, four restarts, four identical fatal lines -> NOT systemic.

    Counting lines instead of distinct targets would turn every retried hard
    target into a fleet-wide alarm.
    """
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC"},
            {"target": "2DEF"},
            {
                "target": "3GHI",
                "poses": 0,
                "exit_code": 1,
                "stderr_text": FATAL_TMPDIR,
                "restarts": 4,
            },
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out
    assert "recurring_fatal" not in out


def test_signature_groups_across_differing_paths():
    """Same fault, different path operand -> one signature."""
    a = 'Fatal error: filesystem error: in temp_directory_path: path "/a/b/.tmp" is not a directory'
    b = 'Fatal error: filesystem error: in temp_directory_path: path "/x/y/z/.tmp" is not a directory'
    assert W.normalise_signature(a) == W.normalise_signature(b)


def test_signature_keeps_distinct_causes_distinct():
    """Over-normalising would merge unrelated faults into one alarm."""
    a = 'Fatal error: filesystem error: in temp_directory_path: path "/a" is not a directory'
    b = "terminate called after throwing an instance of 'std::bad_alloc'"
    assert W.normalise_signature(a) != W.normalise_signature(b)


def test_since_minutes_bounds_fatal_evidence(tmp_path):
    """Old damage can be excluded so a fixed campaign stops latching."""
    build_campaign(
        tmp_path,
        [
            {"target": "2DEF", "poses": 0, "exit_code": 1, "stderr_text": FATAL_TMPDIR},
            {"target": "3GHI", "poses": 0, "exit_code": 1, "stderr_text": FATAL_TMPDIR},
        ],
    )
    old = time.time() - 6 * 3600
    for path in sorted(tmp_path.rglob("*"), reverse=True):
        os.utime(path, (old, old))

    unbounded = W.observe(tmp_path)
    assert W.detect_recurring_fatal(unbounded.cells, 2)

    bounded = W.observe(tmp_path, since=time.time() - 600)
    assert W.detect_recurring_fatal(bounded.cells, 2) == []


# ===========================================================================
# Detector: stall -- bytes and mtime, never file count
# ===========================================================================


def test_flat_file_count_with_growing_bytes_is_not_a_stall(tmp_path):
    """THE non-false-positive test.

    During docking the engine appends to an EXISTING stdout.log. The file
    count is flat for minutes while the run is perfectly healthy. A
    count-based liveness check calls that wedged -- an error actually made
    against this campaign. This test pins all three facts at once: count
    unchanged, bytes grown, verdict healthy.
    """
    build_campaign(tmp_path, [{"target": "1ABC"}, {"target": "2DEF"}])
    fresh(tmp_path)
    before = W.observe(tmp_path)

    # Append to an existing log. Create nothing; delete nothing.
    log = next(tmp_path.rglob("stdout.log"))
    with log.open("a", encoding="utf-8") as handle:
        handle.write("y" * 50_000)
    os.utime(log, None)

    after = W.observe(tmp_path)

    assert after.file_count == before.file_count, "no file was created"
    assert after.total_bytes > before.total_bytes, "bytes must have grown"
    assert W.detect_stall(after, W.DEFAULT_STALL_FACTOR, W.DEFAULT_STALL_FLOOR_S) == []

    # Also exercise the two-observation path, where a count-based detector is
    # most tempting: it sees an unchanged file count between snapshots and
    # concludes "wedged". Bytes and mtime both say otherwise, and they win.
    previous = {
        "total_bytes": before.total_bytes,
        "max_mtime": before.max_mtime,
        "file_count": before.file_count,
    }
    assert W.detect_stall(after, W.DEFAULT_STALL_FACTOR, W.DEFAULT_STALL_FLOOR_S,
                          previous=previous) == []

    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out


def test_stall_fires_when_nothing_has_changed(tmp_path):
    build_campaign(tmp_path, [{"target": "1ABC"}, {"target": "2DEF"}])
    old = time.time() - 6 * 3600
    for path in sorted(tmp_path.rglob("*"), reverse=True):
        os.utime(path, (old, old))
    os.utime(tmp_path, (old, old))

    obs = W.observe(tmp_path)
    findings = W.detect_stall(obs, W.DEFAULT_STALL_FACTOR, W.DEFAULT_STALL_FLOOR_S)
    assert findings and findings[0]["detector"] == "stall"


def test_stall_budget_ignores_poseless_wall_times(tmp_path):
    """The dead cells of 2026-09-20 each "took" 1.09 s.

    Asserted on an ALL-POSELESS tree, which is the only shape where the
    exclusion changes anything: max() is already robust to short outliers, so
    on a mixed campaign including the corpses is a no-op and a mixed-tree
    assertion would pass against a detector that had no exclusion at all.
    Here every wall time is ~1 s, so a detector that trusted them would set a
    ~3 s stall budget and fire on the first healthy pause.
    """
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC", "poses": 0, "exit_code": 1, "wall_s": 1.09},
            {"target": "2DEF", "poses": 0, "exit_code": 1, "wall_s": 1.08},
            {"target": "3GHI", "poses": 0, "exit_code": 1, "wall_s": 1.09},
        ],
    )
    obs = W.observe(tmp_path)
    estimate = W.observed_cell_seconds(obs.cells)
    assert estimate is None or estimate > 10.0, (
        f"budget estimate {estimate} came from poseless wall times"
    )

    # And the floor must keep the detector quiet on a freshly-written tree.
    fresh(tmp_path)
    assert W.detect_stall(W.observe(tmp_path), W.DEFAULT_STALL_FACTOR,
                          W.DEFAULT_STALL_FLOOR_S) == []


def test_mixed_tree_budget_is_set_by_the_slowest_target(tmp_path):
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC", "wall_s": 1800.0},
            {"target": "2DEF", "poses": 0, "exit_code": 1, "wall_s": 1.09},
        ],
    )
    obs = W.observe(tmp_path)
    assert W.observed_cell_seconds(obs.cells) == pytest.approx(1800.0)


def test_growth_between_observations_suppresses_stall(tmp_path):
    """An old mtime plus grown bytes is a slow run, not a dead one."""
    build_campaign(tmp_path, [{"target": "1ABC"}])
    old = time.time() - 6 * 3600
    for path in sorted(tmp_path.rglob("*"), reverse=True):
        os.utime(path, (old, old))
    obs = W.observe(tmp_path)
    previous = {"total_bytes": obs.total_bytes - 1000, "max_mtime": obs.max_mtime - 10}
    assert W.detect_stall(obs, W.DEFAULT_STALL_FACTOR, W.DEFAULT_STALL_FLOOR_S,
                          previous=previous) == []


# ===========================================================================
# Vacuity -- exit 3, and emphatically not exit 0
# ===========================================================================


def test_empty_campaign_is_vacuous_not_healthy(tmp_path):
    empty = tmp_path / "campaign"
    empty.mkdir()
    code, out = run_watchdog(empty)
    assert code == W.EXIT_VACUOUS, out
    assert code != W.EXIT_HEALTHY
    assert "NOT a clean bill of health" in out


def test_missing_root_is_vacuous(tmp_path):
    code, out = run_watchdog(tmp_path / "nope")
    assert code == W.EXIT_VACUOUS, out


def test_all_quarantined_is_vacuous(tmp_path):
    """Every cell set aside by a human leaves nothing live to judge."""
    make_cell(tmp_path, "1ABC", quarantine_suffix=".TMPDIRFAIL_20260920T210340Z")
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_VACUOUS, out


def test_quarantined_cells_excluded_by_default(tmp_path):
    """A triaged historical failure must not re-alarm on every poll."""
    build_campaign(
        tmp_path,
        [
            {"target": "1ABC"},
            {
                "target": "2DEF",
                "poses": 0,
                "exit_code": 1,
                "stderr_text": FATAL_TMPDIR,
                "quarantine_suffix": ".TMPDIRFAIL_20260920T210340Z",
            },
            {
                "target": "3GHI",
                "poses": 0,
                "exit_code": 1,
                "stderr_text": FATAL_TMPDIR,
                "quarantine_suffix": ".TMPDIRFAIL_20260920T210340Z",
            },
        ],
    )
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path)
    assert code == W.EXIT_HEALTHY, out

    code2, out2 = run_watchdog(tmp_path, "--include-quarantined")
    assert code2 == W.EXIT_SYSTEMIC, out2


# ===========================================================================
# The observer must not perturb what it measures
# ===========================================================================


def test_watchdog_does_not_modify_the_tree(tmp_path):
    """Sizes and mtimes identical before and after a full run.

    Checked on a STATIC fixture: the same check against the live campaign is
    meaningless, because the campaign itself is writing while the scan runs.
    """
    build_campaign(tmp_path, [{"target": "1ABC"}, {"target": "2DEF", "poses": 0}])
    fresh(tmp_path)
    before = snapshot(tmp_path)
    run_watchdog(tmp_path)
    assert snapshot(tmp_path) == before


def test_state_file_inside_campaign_root_is_refused(tmp_path):
    """A state file in the tree would corrupt the byte total it records."""
    build_campaign(tmp_path, [{"target": "1ABC"}])
    fresh(tmp_path)
    code, out = run_watchdog(tmp_path, "--state-file", str(tmp_path / "wd.json"))
    assert code == W.EXIT_USAGE, out
    assert "refusing to write state file" in out
    assert not (tmp_path / "wd.json").exists()


def test_state_file_outside_root_is_written(tmp_path):
    root = tmp_path / "campaign"
    root.mkdir()
    build_campaign(root, [{"target": "1ABC"}])
    fresh(root)
    state = tmp_path / "state" / "wd.json"
    code, out = run_watchdog(root, "--state-file", str(state))
    assert code == W.EXIT_HEALTHY, out
    assert state.is_file()


# ===========================================================================
# campaign_preflight.sh
# ===========================================================================


@pytest.fixture
def durable_tmp():
    """A scratch directory that is NOT inside an agent session workspace.

    pytest's own tmp_path is unusable for the pass-expecting preflight tests:
    on this machine TMPDIR points at the agent session workspace, so tmp_path
    resolves under .claude-science/.../workspaces/... and the ephemerality
    check correctly condemns it. The check firing on pytest's fixtures is the
    detector working, not a bug -- but a test that wants to isolate one
    failure needs ground the check approves of.
    """
    import shutil
    import tempfile

    base = tempfile.mkdtemp(prefix="flexaidds_preflight_", dir="/tmp")
    try:
        yield Path(base)
    finally:
        # Created by this test moments ago, under /tmp, outside any research
        # tree. Nothing here is ever a research artefact.
        shutil.rmtree(base, ignore_errors=True)


DURABLE_TMPDIR = "/tmp"


def make_engine(tmp_path: Path, *, missing: str = "") -> Path:
    """An engine binary with its data files beside it."""
    bindir = tmp_path / "build_lane"
    bindir.mkdir(parents=True, exist_ok=True)
    engine = bindir / "FlexAIDdS"
    engine.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    engine.chmod(0o755)
    for name in ("MC_st0r5.2_6.dat", "AMINO.def", "NUCLEOTIDES.def"):
        if name == missing:
            continue
        (bindir / name).write_text("data\n", encoding="utf-8")
    return engine


def make_cache(tmp_path: Path, n: int = 4, stray: bool = False) -> Path:
    cache = tmp_path / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    for i in range(n):
        (cache / f"{i:04d}").mkdir(exist_ok=True)
    if stray:
        (cache / "README.txt").write_text("notes\n", encoding="utf-8")
    return cache


def run_preflight(*args: str, env: dict | None = None) -> tuple[int, str]:
    full_env = dict(os.environ)
    # Overwrite rather than setdefault: this session's real TMPDIR points
    # into the agent workspace, which the ephemerality check condemns (and
    # should). Tests that want to exercise that path set TMPDIR explicitly.
    full_env["TMPDIR"] = DURABLE_TMPDIR
    for key in ("TMP", "TEMP", "FLEXAIDDS_TMPDIR"):
        full_env.pop(key, None)
    if env:
        full_env.update(env)
    proc = subprocess.run(
        ["/bin/bash", str(PREFLIGHT), *args],
        capture_output=True,
        text=True,
        env=full_env,
    )
    return proc.returncode, proc.stdout + proc.stderr


def test_preflight_passes_on_a_sound_setup(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    cache = make_cache(tmp_path, n=4)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--cache", str(cache),
        "--expect-targets", "4",
        "--min-free-gb", "0",
    )
    assert code == 0, out
    assert "OK to launch" in out


def test_preflight_rejects_ephemeral_tmpdir(durable_tmp):
    tmp_path = durable_tmp
    """The 2026-09-20 root cause, caught before launch."""
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    ephemeral = "/Users/x/.claude-science/orgs/AAA/workspaces/BBB/.tmp"
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--min-free-gb", "0",
        env={"TMPDIR": ephemeral},
    )
    assert code == 1, out
    assert "ephemeral:TMPDIR" in out
    assert "REFUSING TO LAUNCH" in out


def test_preflight_accepts_durable_tmpdir_containing_workspaces(durable_tmp):
    tmp_path = durable_tmp
    """Only BOTH markers condemn a path.

    A durable directory that merely contains the word "workspaces" is not an
    agent session workspace, and failing it would train the operator to pass
    --require-all=off and ignore the whole gate.
    """
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--min-free-gb", "0",
        env={"TMPDIR": "/Volumes/big/workspaces/scratch"},
    )
    assert code == 0, out


def test_preflight_rejects_missing_engine(durable_tmp):
    tmp_path = durable_tmp
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(tmp_path / "nope" / "FlexAIDdS"),
        "--output-root", str(out_root),
        "--min-free-gb", "0",
    )
    assert code == 1, out
    assert "binary not found" in out


def test_preflight_rejects_non_executable_engine(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    engine.chmod(0o644)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--min-free-gb", "0",
    )
    assert code == 1, out
    assert "not executable" in out


def test_preflight_rejects_missing_data_file(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path, missing="AMINO.def")
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--min-free-gb", "0",
    )
    assert code == 1, out
    assert "AMINO.def" in out


def test_preflight_rejects_cache_count_mismatch(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    cache = make_cache(tmp_path, n=3)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--cache", str(cache),
        "--expect-targets", "85",
        "--min-free-gb", "0",
    )
    assert code == 1, out
    assert "cache holds 3 target dirs, expected 85" in out


def test_preflight_rejects_cache_strays(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    cache = make_cache(tmp_path, n=4, stray=True)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--cache", str(cache),
        "--expect-targets", "4",
        "--min-free-gb", "0",
    )
    assert code == 1, out
    assert "cache-strays" in out
    assert "README.txt" in out


def test_preflight_rejects_insufficient_disk(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", str(engine),
        "--output-root", str(out_root),
        "--min-free-gb", "999999999",
    )
    assert code == 1, out
    assert "floor is 999999999 GB" in out


def test_preflight_resolves_nested_launcher_variables(durable_tmp):
    tmp_path = durable_tmp
    """THE bug this script was written against.

    A preflight that sed-substituted a launcher's text expanded $R but not
    $BB, tested the literal string "$BB/FlexAIDdS", and reported the engine
    ABSENT while the engine was present. Recovering variables by eval'ing
    assignments in order resolves the chain R -> BB -> ENGINE correctly.
    """
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    launcher = tmp_path / "launch.sh"
    launcher.write_text(
        "#!/usr/bin/env bash\n"
        f"R={tmp_path}\n"
        "BB=$R/build_lane\n"
        "ENGINE=$BB/FlexAIDdS\n"
        f"OUTROOT={out_root}\n"
        'echo "launching"\n',
        encoding="utf-8",
    )
    code, out = run_preflight("--launcher", str(launcher), "--min-free-gb", "0")
    assert code == 0, out
    assert str(engine) in out
    assert "$BB" not in out, "an unexpanded variable reached a check"


def test_preflight_resolves_exported_launcher_variables(durable_tmp):
    """`export VAR=...` lines must be honoured, not silently dropped.

    They were. In a bash `case`, [A-Za-z_]*=* matches "export FOO=bar" (the
    leading * swallows "xport FOO"), so the plain-assignment branch won, the
    "export " prefix survived into the variable name, the identifier check
    saw a space and dropped the line -- without counting it. Launchers that
    export their paths got a preflight that silently checked nothing. The
    unparsable counter is asserted too, because an uncounted drop is how this
    stayed invisible.
    """
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    launcher = tmp_path / "launch.sh"
    launcher.write_text(
        "#!/usr/bin/env bash\n"
        f"export R={tmp_path}\n"
        "export BB=$R/build_lane\n"
        "export ENGINE=$BB/FlexAIDdS\n"
        f"export OUTROOT={out_root}\n",
        encoding="utf-8",
    )
    code, out = run_preflight("--launcher", str(launcher), "--min-free-gb", "0")
    assert "4 assignments eval'd" in out, out
    assert "0 unparsable" in out, out
    assert code == 0, out
    assert str(engine) in out


def test_preflight_reports_unparsable_launcher_lines(durable_tmp):
    """A line the loader cannot use is COUNTED, never dropped in silence.

    The fixture must be a line that actually reaches the identifier check:
    "my-var=1" matches the assignment pattern (the leading * swallows
    "y-var") but its name holds a hyphen. A line with no "=" at all never
    gets that far and would leave the counter untested.
    """
    tmp_path = durable_tmp
    make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    launcher = tmp_path / "launch.sh"
    launcher.write_text(
        f"export R={tmp_path}\n"
        "export BB=$R/build_lane\n"
        "export ENGINE=$BB/FlexAIDdS\n"
        f"export OUTROOT={out_root}\n"
        "my-var=1\n"
        "not an assignment at all\n",
        encoding="utf-8",
    )
    code, out = run_preflight("--launcher", str(launcher), "--min-free-gb", "0")
    assert "4 assignments eval'd" in out, out
    assert "1 unparsable" in out, out
    assert code == 0, out


def test_preflight_expands_a_variable_passed_in_a_path_argument(durable_tmp):
    """The sed bug in its purest form.

    A path argument that still contains $BB must be EXPANDED before it is
    tested. The broken preflight stat'd the literal string "$BB/FlexAIDdS"
    and reported the engine absent while it sat right there. This exercises
    resolve_path directly, which the launcher test does not: by the time
    load_launcher_vars has run, ENGINE is already a concrete path.
    """
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    code, out = run_preflight(
        "--engine", "$BB/FlexAIDdS",
        "--output-root", str(out_root),
        "--min-free-gb", "0",
        env={"BB": str(tmp_path / "build_lane")},
    )
    assert code == 0, out
    assert str(engine) in out
    assert "$BB" not in out, "an unexpanded variable reached a check"


def test_preflight_does_not_execute_command_substitution(durable_tmp):
    tmp_path = durable_tmp
    """A preflight must not run code it finds in a launcher."""
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    canary = tmp_path / "canary"
    launcher = tmp_path / "launch.sh"
    launcher.write_text(
        f"R={tmp_path}\n"
        "BB=$R/build_lane\n"
        "ENGINE=$BB/FlexAIDdS\n"
        f"OUTROOT={out_root}\n"
        f"STAMP=$(touch {canary})\n",
        encoding="utf-8",
    )
    code, out = run_preflight("--launcher", str(launcher), "--min-free-gb", "0")
    assert not canary.exists(), "preflight executed command substitution from a launcher"
    assert "1 skipped" in out
    assert code == 0, out


def test_preflight_env_file_is_sourced(durable_tmp):
    tmp_path = durable_tmp
    engine = make_engine(tmp_path)
    out_root = tmp_path / "out"
    out_root.mkdir()
    env_file = tmp_path / "env.sh"
    env_file.write_text(
        f"BB={tmp_path}/build_lane\n"
        "ENGINE=$BB/FlexAIDdS\n"
        f"OUTROOT={out_root}\n",
        encoding="utf-8",
    )
    code, out = run_preflight("--env-file", str(env_file), "--min-free-gb", "0")
    assert code == 0, out
    assert str(engine) in out
