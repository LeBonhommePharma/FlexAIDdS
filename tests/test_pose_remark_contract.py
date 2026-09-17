"""Golden-record contract for the pose REMARK block.

WHY THIS EXISTS
---------------
The pose REMARK block is written by FOUR hand-maintained code paths that have
drifted apart. Measured 2026-09-16 over LIB/ (76 distinct fields, 197 format
literals):

    common core emitted by all four backends      14 fields
    emitted by SOME but not all                   56 fields
    fields duplicated at >1 site inside one file  34 fields

The drift is not cosmetic. Two concrete consequences were measured, not inferred:

  1. DSVIB.* (the project's dS_vib receipt -- status, S_complex_minus_apo,
     S_ligand_free, S_apo, dS_vib, minus_T_dS_vib) is emitted by cluster.cpp
     ONLY. Selecting the DensityPeak or FastOPTICS backend therefore drops the
     entropy provenance from every pose, silently, with no error and no missing
     file. A banked campaign would simply have no entropy receipt.

  2. `rmsd_raw = %.5f` / `rmsd_sym = %.5f` exist in BindingMode.cpp and
     cluster.cpp and in neither DensityPeak_Cluster.cpp nor FOPTICS.cpp.

This suite pins the record format so a field disappearing from one backend is a
TEST FAILURE rather than a silent gap discovered after a campaign is banked.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
It does NOT unify clustering SEMANTICS. The pose-limit cap means a structurally
different thing in each backend BY DESIGN -- CF bounds cluster CONSTRUCTION,
DensityPeak bounds EMISSION only with its creation-side cap intentionally
disabled -- and that is conserved deliberately and guarded separately by
tests/test_dp_partition_invariant.py. Unify the WRITER, never the semantics.

NON-VACUITY
-----------
Every assertion here must be observed failing before it is trusted. The
documented way to do that: delete one DSVIB line from cluster.cpp, or one
`CF.com` snprintf from DensityPeak_Cluster.cpp, run this file, see it go red,
then restore. A guard never seen failing is not a guard.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

LIB = Path(__file__).resolve().parents[1] / "LIB"

# The four code paths that compose a pose REMARK block.
POSE_WRITERS = (
    "cluster.cpp",                 # CF clustering
    "DensityPeak_Cluster.cpp",     # DP clustering
    "FOPTICS.cpp",                 # FastOPTICS ordering
    "BindingMode.cpp",             # binding-mode / OPTICS centroid emission
)

# Measured 2026-09-16: emitted by ALL FOUR writers. This is the floor, not a wish.
COMMON_CORE = (
    "CF",
    "CF.app",
    "CF.com",
    "CF.con",
    "CF.gist",
    "CF.hbond",
    "CF.sas",
    "CF.wal",
    "RMSD to ref. structure (no symmetry correction)",
    "RMSD to ref. structure     (symmetry corrected)",
    "Residue has an overall SAS of",
    "optimizable residue",
    "optimized structure",
)

# The dS_vib receipt. Currently cluster.cpp ONLY -- see test_dsvib_reaches_every_writer.
DSVIB_FIELDS = (
    "DSVIB.status",
    "DSVIB.S_complex_minus_apo",
    "DSVIB.S_ligand_free",
    "DSVIB.S_apo",
    "DSVIB.dS_vib",
    "DSVIB.minus_T_dS_vib",
)


def _src(name: str) -> str:
    p = LIB / name
    assert p.is_file(), f"missing source file {p}"
    return p.read_text(errors="replace")


def _remark_literals(name: str) -> list[str]:
    """Every "REMARK ..." format literal in one file."""
    return re.findall(r'"REMARK[^"]*"', _src(name))


@pytest.fixture(scope="module")
def literals() -> dict[str, list[str]]:
    return {f: _remark_literals(f) for f in POSE_WRITERS}


def test_sources_are_present_and_nonempty(literals):
    """NON-VACUITY GATE. Run before believing any assertion below.

    If a filename is wrong or a file is empty, every `in` check would pass or
    fail for the wrong reason. Assert we are actually reading REMARK-bearing
    source before asserting anything about its contents.
    """
    for f in POSE_WRITERS:
        n = len(literals[f])
        # FastOPTICS_cluster.cpp legitimately emits 0 and is NOT in POSE_WRITERS
        # precisely for that reason: it delegates to BindingMode::output_Population.
        assert n > 0, f"{f} yielded 0 REMARK literals -- wrong path or empty read"
    assert sum(len(v) for v in literals.values()) > 100, (
        "implausibly few REMARK literals across all writers; the reader is wrong"
    )


@pytest.mark.parametrize("field", COMMON_CORE)
def test_common_core_reaches_every_pose_writer(field, literals):
    """A field in the measured common core must not vanish from one backend.

    This is the guard the DSVIB gap shows we needed: a pose record silently
    missing a field on one code path is indistinguishable, downstream, from a
    field that was never meant to be there.
    """
    missing = [f for f in POSE_WRITERS
               if not any(field in lit for lit in literals[f])]
    assert not missing, (
        f"REMARK field {field!r} is emitted by some pose writers but is MISSING "
        f"from {missing}. Either add it there or move the whole block behind the "
        f"shared writer -- do not leave the record format backend-dependent."
    )


def test_rmsd_pair_is_emitted_together(literals):
    """The two RMSD lines are a PAIR. One without the other is a parse hazard.

    scripts/parse_flexaid_arm_results.py prefers the (symmetry corrected) line
    and falls back to (no symmetry correction). A writer emitting only the
    fallback silently changes which quantity every downstream table reports.
    """
    sym = "RMSD to ref. structure     (symmetry corrected)"
    nosym = "RMSD to ref. structure (no symmetry correction)"
    for f in POSE_WRITERS:
        has_sym = any(sym in lit for lit in literals[f])
        has_nosym = any(nosym in lit for lit in literals[f])
        assert has_sym == has_nosym, (
            f"{f} emits only one of the RMSD pair (sym={has_sym}, nosym={has_nosym}). "
            f"The parser's preference order makes a half-pair a silent metric switch."
        )


def test_symmetry_corrected_line_is_the_hungarian_call(literals):
    """Pin the MEANING of the (symmetry corrected) line, not just its presence.

    Measured at every site -- DensityPeak_Cluster.cpp:682/686,
    BindingMode.cpp:984/990 and :1154/1160, cluster.cpp:892/899,
    FOPTICS.cpp:441/445 -- the pattern is:

        bool Hungarian = false;
        ... calc_rmsd(..., Hungarian)   -> "(no symmetry correction)"
        Hungarian = true;
        ... calc_rmsd(..., Hungarian)   -> "(symmetry corrected)"

    and calc_rmsd.cpp dispatches Hungarian==true to calc_Hungarian_RMSD. So the
    "(symmetry corrected)" line is the engine's HUNGARIAN-ASSIGNMENT RMSD -- NOT
    a graph-automorphism symmetry correction, and NOT serial.

    This matters beyond naming: Hungarian assignment is OVER-PERMISSIVE relative
    to a graph-automorphism metric, so any success flag derived from this field
    is an UPPER BOUND. scripts/aggregate_claim_metrics.py labels the same column
    `serial_legacy_top1`, which is wrong in the dangerous direction.

    If someone changes which quantity feeds that line, this test must fail.
    """
    for f in POSE_WRITERS:
        s = _src(f)
        if "(symmetry corrected)" not in s:
            continue
        for m in re.finditer(r'"REMARK[^"]*\(symmetry corrected\)[^"]*"', s):
            window = s[max(0, m.start() - 900):m.start()]
            assert "Hungarian = true" in window or "Hungarian=true" in window, (
                f"{f}: a '(symmetry corrected)' REMARK is emitted without a preceding "
                f"`Hungarian = true`. Either the quantity changed or the flip moved -- "
                f"both invalidate every downstream label built on this field."
            )


def test_calc_rmsd_double_call_order_is_preserved(literals):
    """calc_rmsd() IS NOT PURE -- gaboom.cpp documents that it writes state.

    Every RMSD-REMARK site calls it TWICE in immediate succession with the flag
    flipped between calls. A refactor that memoizes, reorders, or hoists those
    calls can change the second value. This test exists so that a future shared
    writer cannot quietly 'optimize' the double call away.
    """
    gab = _src("gaboom.cpp")
    assert "calc_rmsd() is not a pure function" in gab, (
        "The non-purity note in gaboom.cpp is gone. Either calc_rmsd was made pure "
        "(then delete this test and say so in the commit) or the warning was lost -- "
        "the second case is how a refactor silently changes RMSD values."
    )
    for f in POSE_WRITERS:
        s = _src(f)
        if "(symmetry corrected)" not in s:
            continue
        # the raw call must precede the sym call in every emitting region
        raw_positions = [m.start() for m in
                         re.finditer(r"\(no symmetry correction\)", s)]
        sym_positions = [m.start() for m in
                         re.finditer(r"\(symmetry corrected\)", s)]
        assert raw_positions and sym_positions, f"{f}: incomplete RMSD pair"
        for sp in sym_positions:
            assert any(rp < sp for rp in raw_positions), (
                f"{f}: a '(symmetry corrected)' emission precedes every raw emission. "
                f"calc_rmsd is stateful; the raw-then-sym order is load-bearing."
            )


@pytest.mark.parametrize("field", DSVIB_FIELDS)
def test_dsvib_reaches_every_writer(field, literals):
    missing = [f for f in POSE_WRITERS
               if not any(field in lit for lit in literals[f])]
    assert not missing, f"{field} missing from {missing}"


def test_no_field_is_emitted_twice_within_one_file(literals):
    offenders = {}
    for f in POSE_WRITERS:
        seen = {}
        for lit in literals[f]:
            key = re.sub(r"%[-0-9.+ #]*[a-zA-Z]+", "%", lit)
            seen[key] = seen.get(key, 0) + 1
        dupes = {k: n for k, n in seen.items() if n > 1}
        if dupes:
            offenders[f] = dupes
    assert not offenders, f"duplicated REMARK emissions within a file: {offenders}"
