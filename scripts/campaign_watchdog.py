#!/usr/bin/env python3
"""Read-only systemic-failure detector for a running docking campaign.

WHY THIS EXISTS
---------------
On 2026-09-20 nine cells of an 85-target Astex campaign died at startup:
$TMPDIR pointed into an agent session workspace that its harness swept after
a few hours idle, and std::filesystem::temp_directory_path() (the throwing
overload) aborted the process. The engine wrote zero poses and exited
non-zero; the driver logged only "ERROR: Incomplete docking", which reads
like a hard target rather than a dead environment.

The engine defect is fixed (LIB/temp_dir.{h,cpp}). This file addresses the
*second* property that made the incident expensive: the campaign kept walking
the roster for six hours, producing empty cell after empty cell, and nobody
knew until a human opened a stderr.log. Nothing watched the run.

WHAT THIS IS, AND IS NOT
------------------------
This is an OBSERVER. It opens files for reading, and does nothing else. It
never kills a process, never writes into the campaign tree, never renames or
deletes anything. It reports and sets an exit code; a human or a driver
decides what to do. That restraint is deliberate: a watchdog that can act on
a false positive is a new way to lose a campaign, and this project has lost
work to confident automation before.

EXIT CODES
----------
    0   healthy      -- cells observed, no detector fired
    1   systemic     -- at least one detector fired
    3   vacuous      -- nothing to observe (missing/empty root, no cells)
    2   usage        -- bad arguments (argparse convention)

Exit 3 exists because "I found no problems" and "I found nothing" are
different answers and must not share a code. An empty campaign directory is
the shape a campaign has when its launcher died before the first cell -- the
single most important case to not report as healthy.

THE MEASUREMENT RULE THAT COST A DAY
------------------------------------
Liveness is judged by BYTES and MTIME, never by FILE COUNT. During docking
the engine appends to an existing stdout.log (2.3 MB in a healthy cell) and
rewrites in place; the file count in a live cell is flat for minutes at a
time while the run is perfectly healthy. Counting files and calling a
stationary count "wedged" is a false positive that was actually made against
this campaign. Every liveness path in this file reads st_mtime and st_size.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Defaults. Every one is overridable from the CLI; none is a path.
# ---------------------------------------------------------------------------

DEFAULT_ZERO_POSE_RUN = 3          # consecutive poseless cells that mean "stop"
DEFAULT_FATAL_TARGETS = 2          # distinct targets sharing one fatal signature
DEFAULT_STALL_FACTOR = 3.0         # multiples of observed per-cell time
DEFAULT_STALL_FLOOR_S = 900.0      # never call a stall before this, whatever the data

# A cell is a directory holding at least one of these.
CELL_MARKERS = ("result.csv", "stdout.log", "stderr.log")

# Restart subdirectories (r1, r2, ...) carry their own stdout/stderr but are
# not cells; they belong to the cell above them.
RESTART_DIR_RE = re.compile(r"^r\d+$")

# Human quarantine markers. This project sets a bad cell aside by RENAMING it
# with a reason and a UTC stamp (never by deleting it), e.g.
#   2GBP.TMPDIRFAIL_20260920T210340Z
#   1GPK.CONTAMINATED_20260920T062650Z
# A quarantined cell is an already-triaged historical failure. Counting it
# against a live run would make a fresh watchdog scream about nine failures a
# human already dealt with, so it is excluded from the verdict by default and
# reported separately. --include-quarantined restores it to the verdict.
QUARANTINE_RE = re.compile(r"\.[A-Z][A-Z0-9_]*_\d{8}T\d{6}Z$")

# Pose files are <stem>_<n>.pdb. The regex deliberately requires digits so
# that <stem>_INI.pdb -- the input frame, present even in a cell that docked
# nothing -- does not count as a pose.
POSE_FILE_RE = re.compile(r"_\d+\.pdb$")
ELECTED_POSE = "elected_pose.pdb"


# ---------------------------------------------------------------------------
# Fatal-line classification
# ---------------------------------------------------------------------------
#
# The signature detector keys on FATAL lines only, and this list is the whole
# reason it does not fire on a healthy campaign. A healthy cell's stderr.log
# is not empty -- it carries routine per-target diagnostics:
#
#   [DIRECT-IC] GPA topology atoms=90001,90013,90002 local=1,13,2
#   [GRID-CACHE] Saved 4453 grid points -> .../1JD0.rrg
#   [NATIVE-SEED-RMSD] ... round-trip RMSD = 2.97 A ...
#   [WARN] NATIVE-SEED-RMSD = 2.97 A exceeds 1.0 A threshold.
#
# That [WARN] recurs across most targets and normalises to a constant string.
# A detector that grouped *any* recurring stderr line would flag every healthy
# campaign ever run. Matching is therefore restricted to process-death
# signatures: the runtime aborting, not the science complaining.
#
# Note what is NOT here: bare "ERROR:". The driver's own "ERROR: Incomplete
# docking" is the misleading message that masked the incident -- it is emitted
# for genuinely hard targets too, so promoting it to fatal would trade a
# missed detection for a noisy one. Incomplete docking is caught instead by
# the zero-pose detector, which keys on the outcome rather than the wording.

FATAL_PATTERNS: Tuple[re.Pattern, ...] = (
    re.compile(r"^\s*Fatal error:"),                    # the 2026-09-20 signature
    re.compile(r"terminate called"),                    # uncaught C++ exception
    re.compile(r"libc\+\+abi"),
    re.compile(r"^\s*what\(\):"),
    re.compile(r"Segmentation fault"),
    re.compile(r"Bus error"),
    re.compile(r"Abort trap"),
    re.compile(r"std::bad_alloc"),
    re.compile(r"Assertion failed"),
    re.compile(r"stack smashing detected"),
    re.compile(r"error while loading shared libraries"),
    re.compile(r"cannot open shared object file"),
    re.compile(r"^\s*Killed\b"),                        # OOM killer
    re.compile(r"dyld\[\d+\]:"),                        # macOS loader failure
    re.compile(r"Sanitizer:"),                          # ASan/UBSan reports
)

# Volatile substrings are erased before grouping, so that one environmental
# cause yields one signature. Order matters: paths before numbers, or a path
# containing digits would be shredded into an unrecognisable stem.
_NORMALISERS: Tuple[Tuple[re.Pattern, str], ...] = (
    (re.compile(r'"[^"]*"'), '"<S>"'),                  # quoted operands (paths)
    (re.compile(r"'[^']*'"), "'<S>'"),
    (re.compile(r"(?<![\w/])/[^\s:,;)]+"), "<PATH>"),   # bare absolute paths
    (re.compile(r"0x[0-9a-fA-F]+"), "<ADDR>"),
    (re.compile(r"\b[0-9a-f]{16,}\b"), "<HEX>"),
    (re.compile(r"\b\d+(?:\.\d+)?\b"), "<N>"),
    (re.compile(r"\s+"), " "),
)

SIGNATURE_MAX_CHARS = 160


def is_fatal_line(line: str) -> bool:
    """True when a stderr line indicates the process died, not that it warned."""
    return any(p.search(line) for p in FATAL_PATTERNS)


def normalise_signature(line: str) -> str:
    """Collapse a fatal line to a cause-stable key.

    Erases paths, addresses and numbers so that the same failure at different
    targets groups together, while leaving the error vocabulary intact so that
    two different failures do not. Over-normalising would merge distinct
    causes into one alarm; that is why only volatile token classes are
    replaced and the remaining words are kept verbatim.
    """
    out = line.strip()
    for pattern, repl in _NORMALISERS:
        out = pattern.sub(repl, out)
    return out.strip()[:SIGNATURE_MAX_CHARS]


# ---------------------------------------------------------------------------
# Tree model
# ---------------------------------------------------------------------------


@dataclass
class Cell:
    """One target's output directory."""

    path: Path
    target: str
    quarantined: bool
    pose_count: int
    pose_source: str                 # "result.csv" | "pose_files"
    wall_time_s: Optional[float]
    docking_exit_code: Optional[int]
    mtime: float                     # completion time, best available
    total_bytes: int
    fatal_signatures: List[str] = field(default_factory=list)

    @property
    def poseless(self) -> bool:
        return self.pose_count == 0


@dataclass
class Observation:
    """Everything one read-only pass over the tree established."""

    root: Path
    cells: List[Cell]
    quarantined: List[Cell]
    max_mtime: float
    total_bytes: int
    file_count: int
    scanned_at: float


def _is_quarantined(path: Path, root: Path) -> bool:
    """True when any component between root and path carries a quarantine stamp."""
    try:
        rel = path.relative_to(root)
    except ValueError:
        rel = path
    return any(QUARANTINE_RE.search(part) for part in rel.parts)


def _read_result_csv(path: Path) -> Optional[Dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle):
                return row
    except (OSError, csv.Error):
        return None
    return None


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    text = value.strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", ""}:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = value.strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", ""}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def count_poses(cell_dir: Path, row: Optional[Dict[str, str]]) -> Tuple[int, str]:
    """Poses this cell produced, and where the number came from.

    result.csv's num_poses is the driver's own accounting and is preferred.
    When it is missing or unparseable the pose files on disk are counted, so
    that a cell whose driver died before writing result.csv is still assessed
    on evidence rather than skipped.
    """
    if row is not None:
        n = _to_int(row.get("num_poses"))
        if n is not None:
            return n, "result.csv"

    count = 0
    try:
        for entry in os.scandir(cell_dir):
            if not entry.is_file(follow_symlinks=False):
                continue
            if POSE_FILE_RE.search(entry.name) or entry.name == ELECTED_POSE:
                count += 1
    except OSError:
        return 0, "pose_files"
    return count, "pose_files"


def _scan_stderr_files(
    cell_dir: Path, max_bytes: int, since: Optional[float] = None
) -> List[str]:
    """Distinct normalised fatal signatures in this cell's stderr logs.

    Reads the cell's own stderr.log plus each restart subdirectory's, because
    an environment failure kills every restart identically and the restart
    logs are corroborating evidence for the same cell.

    `since` is an epoch cutoff: logs last modified before it are ignored.
    Without it the detector answers "has this campaign been damaged", which
    latches forever once damage exists -- running against the live seed2
    campaign, stderr written at 15:58 kept the alarm on hours after the
    operator had moved on. With it the detector answers "is this campaign
    failing NOW", which is what a watchdog polled on a timer needs. Neither
    is wrong; they are different questions, so the bound is explicit and
    unbounded remains the default (never silently forget evidence of damage).
    """
    signatures: List[str] = []
    seen = set()
    candidates = [cell_dir / "stderr.log"]
    try:
        for entry in os.scandir(cell_dir):
            if entry.is_dir(follow_symlinks=False) and RESTART_DIR_RE.match(entry.name):
                candidates.append(Path(entry.path) / "stderr.log")
    except OSError:
        pass

    for log in candidates:
        try:
            if not log.is_file():
                continue
            stat_result = log.stat()
            if since is not None and stat_result.st_mtime < since:
                continue
            size = stat_result.st_size
            with log.open("r", encoding="utf-8", errors="replace") as handle:
                if size > max_bytes:
                    handle.seek(size - max_bytes)
                    handle.readline()          # drop the partial line
                for line in handle:
                    if is_fatal_line(line):
                        sig = normalise_signature(line)
                        if sig and sig not in seen:
                            seen.add(sig)
                            signatures.append(sig)
        except OSError:
            continue
    return signatures


def _dir_stats(path: Path) -> Tuple[float, int, int]:
    """(max mtime, total bytes, file count) for a subtree. Read-only."""
    max_mtime = 0.0
    total = 0
    count = 0
    for dirpath, _dirnames, filenames in os.walk(path, onerror=lambda _e: None):
        try:
            st = os.stat(dirpath)
            max_mtime = max(max_mtime, st.st_mtime)
        except OSError:
            pass
        for name in filenames:
            try:
                st = os.stat(os.path.join(dirpath, name))
            except OSError:
                continue
            count += 1
            total += st.st_size
            if st.st_mtime > max_mtime:
                max_mtime = st.st_mtime
    return max_mtime, total, count


def discover_cells(root: Path) -> List[Path]:
    """Directories that look like a target's output cell.

    Two nesting facts about a real campaign tree drive this, and getting
    either wrong miscounts every cell in the run:

    1. A DatasetRunner target produces TWO nested marker-bearing directories.
       The outer <arm>/<TARGET>/ holds the runner's own result.csv, run.log
       and aggregate CSVs; the inner <arm>/<TARGET>/<TARGET>/ holds the
       engine's stdout.log, stderr.log and poses. Both match a naive marker
       test, which double-counts every target -- measured as 267 "cells" for
       a 97-target arm before this rule existed. The DEEPEST marker-bearing
       directory in a branch is the real cell, because that is where the
       engine actually wrote. When only the outer directory exists (the
       engine died before producing anything) that outer directory is itself
       the deepest, so the cell is still seen rather than dropped.

    2. Restart subdirectories rN/ carry their own stdout.log and stderr.log
       and are DEEPER than the cell. They must be pruned before the
       deepest-wins rule runs, or a four-restart cell would be reported as
       four cells. They are not lost: their stderr is folded into the parent
       cell's evidence by _scan_stderr_files.
    """
    candidates: List[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, onerror=lambda _e: None):
        dirnames[:] = [d for d in dirnames if not RESTART_DIR_RE.match(d)]
        if set(filenames).intersection(CELL_MARKERS):
            candidates.append(Path(dirpath))

    # Deepest wins: drop any candidate that is a strict ancestor of another.
    as_strings = {str(c) + os.sep for c in candidates}
    deepest = [
        c
        for c in candidates
        if not any(other.startswith(str(c) + os.sep) for other in as_strings if other != str(c) + os.sep)
    ]
    return sorted(deepest)


def observe(
    root: Path,
    max_stderr_bytes: int = 262_144,
    since: Optional[float] = None,
) -> Observation:
    """One read-only pass. Opens files for reading; writes nothing, anywhere."""
    scanned_at = time.time()
    cells: List[Cell] = []
    quarantined: List[Cell] = []

    for cell_dir in discover_cells(root):
        row = _read_result_csv(cell_dir / "result.csv")
        poses, pose_source = count_poses(cell_dir, row)
        mtime, total_bytes, _ = _dir_stats(cell_dir)

        result_csv = cell_dir / "result.csv"
        if result_csv.is_file():
            try:
                mtime = result_csv.stat().st_mtime
            except OSError:
                pass

        target = (row or {}).get("pdb_id") or cell_dir.name
        cell = Cell(
            path=cell_dir,
            target=str(target).strip() or cell_dir.name,
            quarantined=_is_quarantined(cell_dir, root),
            pose_count=poses,
            pose_source=pose_source,
            wall_time_s=_to_float((row or {}).get("wall_time_s")),
            docking_exit_code=_to_int((row or {}).get("docking_exit_code")),
            mtime=mtime,
            total_bytes=total_bytes,
            fatal_signatures=_scan_stderr_files(cell_dir, max_stderr_bytes, since=since),
        )
        (quarantined if cell.quarantined else cells).append(cell)

    max_mtime, total_bytes, file_count = _dir_stats(root)
    cells.sort(key=lambda c: (c.mtime, str(c.path)))
    quarantined.sort(key=lambda c: (c.mtime, str(c.path)))
    return Observation(
        root=root,
        cells=cells,
        quarantined=quarantined,
        max_mtime=max_mtime,
        total_bytes=total_bytes,
        file_count=file_count,
        scanned_at=scanned_at,
    )


# ---------------------------------------------------------------------------
# Detectors. Each returns a list of finding dicts; empty means quiet.
# ---------------------------------------------------------------------------


def detect_zero_pose_run(cells: Sequence[Cell], threshold: int) -> List[dict]:
    """Consecutive poseless cells at the END of completion order.

    Cells are ordered by completion time (result.csv mtime), and only the
    TRAILING run counts toward the alarm. That distinction is the difference
    between "the campaign is dead right now" and "three targets failed at
    04:00 and it recovered" -- the second is history and must not stop a run
    that is currently producing poses. A historical burst is still reported,
    as context, without firing.
    """
    findings: List[dict] = []
    if not cells:
        return findings

    trailing: List[Cell] = []
    for cell in reversed(cells):
        if cell.poseless:
            trailing.append(cell)
        else:
            break
    trailing.reverse()

    longest = 0
    run = 0
    for cell in cells:
        run = run + 1 if cell.poseless else 0
        longest = max(longest, run)

    if len(trailing) >= threshold:
        findings.append(
            {
                "detector": "zero_pose_run",
                "severity": "systemic",
                "count": len(trailing),
                "threshold": threshold,
                "targets": [c.target for c in trailing],
                "pose_sources": sorted({c.pose_source for c in trailing}),
                "message": (
                    f"{len(trailing)} consecutive cells finished with zero poses "
                    f"(threshold {threshold}): {', '.join(c.target for c in trailing)}"
                ),
            }
        )
    elif longest >= threshold:
        findings.append(
            {
                "detector": "zero_pose_run",
                "severity": "info",
                "count": longest,
                "threshold": threshold,
                "message": (
                    f"a run of {longest} poseless cells occurred earlier but the "
                    f"campaign recovered; latest cell produced poses"
                ),
            }
        )
    return findings


def detect_recurring_fatal(cells: Sequence[Cell], min_targets: int) -> List[dict]:
    """One fatal signature appearing at several distinct targets.

    A hard target fails in its own way. An environment failure fails
    identically everywhere, which is exactly what the 2026-09-20 TMPDIR sweep
    looked like: the same temp_directory_path abort in nine different cells.
    Counting DISTINCT TARGETS rather than lines is what makes this specific --
    one target retried four times produces four identical lines and must not
    look like a fleet-wide problem.
    """
    by_signature: Dict[str, set] = {}
    for cell in cells:
        for sig in cell.fatal_signatures:
            by_signature.setdefault(sig, set()).add(cell.target)

    findings = []
    for sig, targets in sorted(by_signature.items(), key=lambda kv: -len(kv[1])):
        if len(targets) >= min_targets:
            listed = sorted(targets)
            findings.append(
                {
                    "detector": "recurring_fatal",
                    "severity": "systemic",
                    "signature": sig,
                    "n_targets": len(listed),
                    "min_targets": min_targets,
                    "targets": listed,
                    "message": (
                        f"fatal signature recurs at {len(listed)} distinct targets "
                        f"({', '.join(listed[:8])}"
                        f"{', ...' if len(listed) > 8 else ''}): {sig}"
                    ),
                }
            )
    return findings


def observed_cell_seconds(cells: Sequence[Cell]) -> Optional[float]:
    """Slowest credible per-cell time, from the campaign's own record.

    Prefers wall_time_s as the driver measured it, and takes the MAXIMUM, so
    the budget is set by the slowest target rather than an average that a
    fast target could drag down.

    Poseless cells are excluded. Note precisely what that does and does not
    buy, because the obvious rationale is wrong: max() is already robust to
    short outliers, so on a mixed campaign dropping the 1.09 s corpses
    changes nothing (measured: 1800.0 either way). The exclusion earns its
    place in exactly one shape -- the incident shape, where EVERY cell is
    poseless. Then the surviving wall times are all ~1 s, max() returns ~1 s,
    and a budget built on it would call a healthy pause a stall. Excluding
    them empties the list and falls through to inter-completion gaps, which
    describe the roster's real cadence. DEFAULT_STALL_FLOOR_S is the
    backstop if that is unavailable too.
    """
    walls = [c.wall_time_s for c in cells if c.wall_time_s and not c.poseless]
    if walls:
        return max(walls)

    stamps = sorted(c.mtime for c in cells if c.mtime > 0)
    gaps = [b - a for a, b in zip(stamps, stamps[1:]) if b > a]
    return max(gaps) if gaps else None


def detect_stall(
    obs: Observation,
    factor: float,
    floor_s: float,
    now: Optional[float] = None,
    previous: Optional[dict] = None,
) -> List[dict]:
    """Nothing in the tree has changed for longer than a cell takes.

    Judged on mtime and bytes -- never on file count. A live cell appends to
    an existing stdout.log for minutes without creating a single file, so a
    count-based liveness check reports a healthy run as wedged. That error was
    actually made against this campaign; the byte/mtime rule is the fix.

    With --state-file a previous pass's byte total is compared against this
    one, which upgrades the evidence from "nothing looks recent" to "nothing
    grew between two observations".
    """
    now = time.time() if now is None else now
    budget = max(floor_s, (observed_cell_seconds(obs.cells) or 0.0) * factor)
    idle = now - obs.max_mtime if obs.max_mtime > 0 else None

    if idle is None:
        return []

    grew = None
    if previous is not None:
        prev_bytes = previous.get("total_bytes")
        prev_mtime = previous.get("max_mtime")
        if isinstance(prev_bytes, (int, float)) and isinstance(prev_mtime, (int, float)):
            grew = obs.total_bytes != prev_bytes or obs.max_mtime > prev_mtime

    if idle <= budget:
        return []
    if grew:
        return []

    detail = (
        f"no file in the tree has been modified for {idle / 60.0:.1f} min "
        f"(budget {budget / 60.0:.1f} min)"
    )
    if grew is False:
        detail += "; total bytes unchanged since the previous observation"

    return [
        {
            "detector": "stall",
            "severity": "systemic",
            "idle_seconds": round(idle, 1),
            "budget_seconds": round(budget, 1),
            "total_bytes": obs.total_bytes,
            "bytes_grew_since_previous": grew,
            "message": detail,
        }
    ]


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

EXIT_HEALTHY = 0
EXIT_SYSTEMIC = 1
EXIT_USAGE = 2
EXIT_VACUOUS = 3


def assess(
    obs: Observation,
    zero_pose_run: int = DEFAULT_ZERO_POSE_RUN,
    fatal_targets: int = DEFAULT_FATAL_TARGETS,
    stall_factor: float = DEFAULT_STALL_FACTOR,
    stall_floor_s: float = DEFAULT_STALL_FLOOR_S,
    now: Optional[float] = None,
    previous: Optional[dict] = None,
    check_stall: bool = True,
) -> Tuple[int, List[dict]]:
    """Run every detector. Returns (exit code, findings)."""
    if not obs.cells and not obs.quarantined:
        return EXIT_VACUOUS, [
            {
                "detector": "vacuity",
                "severity": "vacuous",
                "message": (
                    f"no observable cells under {obs.root} "
                    f"({obs.file_count} files, {obs.total_bytes} bytes): "
                    "nothing to judge -- this is NOT a clean bill of health"
                ),
            }
        ]
    if not obs.cells:
        return EXIT_VACUOUS, [
            {
                "detector": "vacuity",
                "severity": "vacuous",
                "message": (
                    f"every cell under {obs.root} is quarantined "
                    f"({len(obs.quarantined)} of them): no live cell to judge"
                ),
            }
        ]

    findings: List[dict] = []
    findings += detect_zero_pose_run(obs.cells, zero_pose_run)
    findings += detect_recurring_fatal(obs.cells, fatal_targets)
    if check_stall:
        findings += detect_stall(obs, stall_factor, stall_floor_s, now=now, previous=previous)

    systemic = [f for f in findings if f.get("severity") == "systemic"]
    return (EXIT_SYSTEMIC if systemic else EXIT_HEALTHY), findings


def render(obs: Observation, code: int, findings: Sequence[dict]) -> str:
    """Loud, greppable, human-first. The incident was missed by reading past
    a quiet line, so the verdict is a banner and never a suffix."""
    banner = {
        EXIT_HEALTHY: "HEALTHY",
        EXIT_SYSTEMIC: "SYSTEMIC FAILURE",
        EXIT_VACUOUS: "VACUOUS -- NOTHING OBSERVED",
    }.get(code, "UNKNOWN")

    poseless = sum(1 for c in obs.cells if c.poseless)
    lines = [
        "=" * 72,
        f"CAMPAIGN WATCHDOG: {banner}",
        "=" * 72,
        f"root            : {obs.root}",
        f"cells observed  : {len(obs.cells)} live, {len(obs.quarantined)} quarantined",
        f"poseless cells  : {poseless}",
        f"tree            : {obs.file_count} files, {obs.total_bytes} bytes",
    ]
    if obs.max_mtime > 0:
        # Measured against now, not against scan start: on a live tree the
        # campaign keeps writing while the walk runs, so a file can be newer
        # than the moment the scan began and a scan-relative age goes
        # negative. A negative "minutes ago" is a reporting bug that makes a
        # healthy, actively-writing run look nonsensical.
        age = max(0.0, time.time() - obs.max_mtime)
        lines.append(f"newest mtime    : {age / 60.0:.1f} min ago")
    lines.append("")

    if not findings:
        lines.append("no detector fired.")
    for finding in findings:
        tag = {"systemic": "!! SYSTEMIC", "info": "-- info", "vacuous": "?? VACUOUS"}.get(
            finding.get("severity", ""), "--"
        )
        lines.append(f"{tag} [{finding['detector']}] {finding['message']}")

    lines.append("")
    lines.append(f"exit {code} ({banner})")
    if code == EXIT_SYSTEMIC:
        lines.append(
            "This tool changed nothing. Stop the campaign yourself if the "
            "finding is real."
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _load_state(path: Optional[Path]) -> Optional[dict]:
    if path is None or not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _save_state(path: Optional[Path], obs: Observation, root: Path) -> None:
    """Persist byte/mtime totals for the next pass.

    Refuses to write anywhere inside the campaign root. The whole value of
    this tool is that it cannot perturb what it measures, and a state file
    dropped into the tree would both violate that and corrupt the very byte
    total it records.
    """
    if path is None:
        return
    resolved = path.resolve()
    root_resolved = root.resolve()
    if resolved == root_resolved or root_resolved in resolved.parents:
        raise ValueError(
            f"refusing to write state file inside the campaign root: {resolved}"
        )
    resolved.parent.mkdir(parents=True, exist_ok=True)
    resolved.write_text(
        json.dumps(
            {
                "root": str(root_resolved),
                "scanned_at": obs.scanned_at,
                "max_mtime": obs.max_mtime,
                "total_bytes": obs.total_bytes,
                "file_count": obs.file_count,
                "n_cells": len(obs.cells),
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read-only systemic-failure detector for a docking campaign.",
        epilog="Exit: 0 healthy, 1 systemic failure, 2 usage, 3 vacuous.",
    )
    parser.add_argument("root", type=Path, help="campaign output root")
    parser.add_argument(
        "--zero-pose-run", type=int, default=DEFAULT_ZERO_POSE_RUN,
        help=f"consecutive poseless cells that mean systemic (default {DEFAULT_ZERO_POSE_RUN})",
    )
    parser.add_argument(
        "--fatal-targets", type=int, default=DEFAULT_FATAL_TARGETS,
        help=f"distinct targets sharing a fatal signature (default {DEFAULT_FATAL_TARGETS})",
    )
    parser.add_argument(
        "--stall-factor", type=float, default=DEFAULT_STALL_FACTOR,
        help=f"multiples of observed per-cell time before stall (default {DEFAULT_STALL_FACTOR})",
    )
    parser.add_argument(
        "--stall-floor-seconds", type=float, default=DEFAULT_STALL_FLOOR_S,
        help=f"never report a stall sooner than this (default {DEFAULT_STALL_FLOOR_S})",
    )
    parser.add_argument(
        "--no-stall-check", action="store_true",
        help="skip the stall detector (for offline analysis of a finished tree)",
    )
    parser.add_argument(
        "--include-quarantined", action="store_true",
        help="count human-quarantined cells toward the verdict",
    )
    parser.add_argument(
        "--since-minutes", type=float, default=None,
        help=(
            "ignore stderr logs older than this many minutes; bounds the fatal-"
            "signature detector to recent damage (default: unbounded)"
        ),
    )
    parser.add_argument(
        "--state-file", type=Path, default=None,
        help="record byte/mtime totals here for growth comparison (never inside root)",
    )
    parser.add_argument("--json", action="store_true", help="emit findings as JSON")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    root: Path = args.root

    if not root.exists():
        message = f"campaign root does not exist: {root}"
        if args.json:
            print(json.dumps({"exit": EXIT_VACUOUS, "findings": [
                {"detector": "vacuity", "severity": "vacuous", "message": message}]}, indent=2))
        else:
            print(f"CAMPAIGN WATCHDOG: VACUOUS -- NOTHING OBSERVED\n{message}\nexit {EXIT_VACUOUS}")
        return EXIT_VACUOUS
    if not root.is_dir():
        print(f"campaign root is not a directory: {root}", file=sys.stderr)
        return EXIT_USAGE

    since = None
    if args.since_minutes is not None:
        since = time.time() - args.since_minutes * 60.0

    obs = observe(root, since=since)
    if args.include_quarantined and obs.quarantined:
        obs.cells = sorted(obs.cells + obs.quarantined, key=lambda c: (c.mtime, str(c.path)))
        obs.quarantined = []

    previous = _load_state(args.state_file)
    code, findings = assess(
        obs,
        zero_pose_run=args.zero_pose_run,
        fatal_targets=args.fatal_targets,
        stall_factor=args.stall_factor,
        stall_floor_s=args.stall_floor_seconds,
        previous=previous,
        check_stall=not args.no_stall_check,
    )

    try:
        _save_state(args.state_file, obs, root)
    except ValueError as exc:
        print(f"{exc}", file=sys.stderr)
        return EXIT_USAGE

    if args.json:
        print(json.dumps(
            {
                "exit": code,
                "root": str(root),
                "n_cells": len(obs.cells),
                "n_quarantined": len(obs.quarantined),
                "total_bytes": obs.total_bytes,
                "file_count": obs.file_count,
                "max_mtime": obs.max_mtime,
                "findings": list(findings),
            },
            indent=2,
        ))
    else:
        print(render(obs, code, findings))
    return code


if __name__ == "__main__":
    sys.exit(main())
