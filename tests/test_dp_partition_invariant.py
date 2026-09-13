"""Guard the clustering cap SEMANTICS of each backend.

WHY THIS EXISTS
---------------
The three clustering backends apply FA->max_results at DIFFERENT STAGES, and the
difference is structural rather than parametric -- the same number means a
different thing depending on which one runs:

    CF  cluster.cpp             cap bounds CONSTRUCTION.  The cluster-building
                                loop BREAKS at the cap, so clusters past it never
                                exist and raising the cap CHANGES THE PARTITION.

    DP  DensityPeak_Cluster.cpp cap bounds EMISSION only.  `nResults = nClusters`
                                builds the complete partition; max_results is
                                applied in the output loops.  Raising the cap
                                reveals more of a FIXED set.

    FO  FastOPTICS_cluster.cpp  max_results reaches only output_Population().

DP's emission-only behaviour is DELIBERATE AND TO BE CONSERVED.  The creation-side
cap is present but commented out, immediately above the live line, so the decision
stays legible.

Re-enabling it would import CF's confound into DP -- pose counts stop being
monotone in search effort, and a scoring-function change moves both the partition
and the point at which it is truncated -- while changing nothing visible in any
output.  That is a regression no benchmark result would flag, which is exactly why
it needs a test rather than a comment.

These are SOURCE-LEVEL assertions.  They cannot prove runtime behaviour; they
prove the code still says what it was designed to say.  A behavioural test would
need a full DP run and is not a substitute for reading these lines.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DP = ROOT / "LIB" / "DensityPeak_Cluster.cpp"
CF = ROOT / "LIB" / "cluster.cpp"
FO = ROOT / "LIB" / "FastOPTICS_cluster.cpp"


def _src(p):
    assert p.exists(), f"{p} missing -- the invariant cannot be checked"
    return p.read_text(encoding="utf-8", errors="replace")


# --- DP: the invariant being conserved --------------------------------------

def test_dp_builds_the_complete_partition():
    """nResults = nClusters must be LIVE: DP clusters the whole population."""
    s = _src(DP)
    live = re.search(r"^\s*nResults\s*=\s*nClusters\s*;", s, re.M)
    assert live, "DP no longer assigns nResults = nClusters -- the complete-partition invariant is broken"


def test_dp_creation_side_cap_stays_disabled():
    """The creation-side cap must remain commented out, and must not reappear live."""
    s = _src(DP)
    for line in s.splitlines():
        stripped = line.strip()
        if stripped.startswith("//"):
            continue
        assert not re.match(r"if\s*\(\s*nClusters\s*<\s*FA->max_results\s*\)", stripped), (
            "The creation-side cap has been re-enabled in DensityPeak_Cluster.cpp. "
            "This is the DELIBERATE design decision documented above the line: DP builds "
            "the complete partition and truncates only at emission. Re-enabling this "
            "imports CF's construction-bound confound into DP invisibly."
        )


def test_dp_applies_the_cap_at_emission():
    """max_results must still bound the OUTPUT loops -- emission-side, not creation-side."""
    s = _src(DP)
    n = len(re.findall(r"i\s*<\s*nResults\s*&&\s*i\s*<\s*FA->max_results", s))
    assert n >= 2, f"expected >=2 emission-side guards in DP, found {n}"


def test_dp_design_intent_is_documented_in_source():
    """The rationale must travel with the code, not only in a test file."""
    s = _src(DP)
    assert "DESIGN INVARIANT" in s, "the design-invariant comment block is gone from DensityPeak_Cluster.cpp"


# --- CF: the contrast that makes DP's choice meaningful ----------------------

def test_cf_still_bounds_construction():
    """Documents the asymmetry deliberately. If CF is ever unified onto emission
    bounding, this test should FAIL and be updated in the same commit -- so the
    change is recorded rather than absorbed silently."""
    s = _src(CF)
    brk = re.search(r"if\s*\(\s*num_of_clusters\s*==\s*num_of_results\s*\)\s*\{\s*break\s*;\s*\}", s)
    assert brk, (
        "cluster.cpp no longer breaks the cluster-building loop at max_results. "
        "If CF was intentionally moved to emission bounding, update this test and "
        "the DP invariant comment together -- the two backends' semantics are a "
        "documented pair, not independent facts."
    )


# --- FO: recorded as UNVERIFIED, deliberately -------------------------------

@pytest.mark.xfail(reason="FO cap semantics inferred from the call site only; "
                          "output_Population's body has not been read",
                   strict=False)
def test_fo_cap_semantics_are_verified():
    """FO passes max_results to output_Population() and nowhere else, which LOOKS
    emission-side. That is an inference from one call site, not a verified fact.
    This test is xfail on purpose: it records an open question instead of letting
    an unread assumption pass as established. Read output_Population, then make
    this a real assertion."""
    s = _src(FO)
    n = len(re.findall(r"FA->max_results", s))
    assert n == 1, f"expected exactly 1 max_results reference in FO, found {n}"
    pytest.fail("output_Population body unread -- semantics not established")
