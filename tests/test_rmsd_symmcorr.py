"""Contract tests for the symmetry-corrected claim metric (#365).

What is pinned here:
  * DIRECTION. symmcorr <= serial, always. Symmetry correction minimises RMSD
    over the ligand's graph automorphisms and the identity mapping is one of
    them, so the corrected value can never exceed the ordered one. This is why
    the correction is safe: it moves targets fail -> PASS only, never the
    reverse, and the serial gate is a conservative lower bound rather than a
    wrong answer.
  * ORDERING. hungarian <= symmcorr <= serial. Element-only Hungarian
    minimises over ALL same-element bijections, a superset of the chemically
    valid automorphisms, so it is over-permissive — which is why it is banned
    from the claim path rather than merely deprecated.
  * FAIL-CLOSED. No silent substitution of a weaker metric, ever.
  * SAME-POSE. A symmcorr value may only be joined to a claim row that
    describes the same pose (pose_sha256).

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


symm = _load("rmsd_symmcorr", "scripts/rmsd_symmcorr.py")
agg = _load("aggregate_claim_metrics", "scripts/aggregate_claim_metrics.py")

import numpy as np

# ── dependency contract: this module must FAIL, never skip ───────────────────
# `scripts/rmsd_symmcorr.py` fails closed: when spyrmsd is absent it raises
# SpyrmsdUnavailable rather than substitute a weaker metric. The tests for that
# contract previously opened with `pytest.importorskip("spyrmsd")`, which fails
# OPEN — on a machine without spyrmsd this module collected ZERO tests and the
# run was reported as success. A gate that reports success without executing is
# precisely the failure mode the code under test exists to prevent, so the test
# suite was contradicting its own subject.
#
# Why a fixture and not a bare `import spyrmsd` at module level: an unguarded
# import raises at collection time, and pytest aborts the WHOLE session on a
# collection error (`!!! Interrupted: 1 error during collection !!!`). One
# missing optional dependency anywhere would then take the entire suite down.
# That trades a silent pass for a noisy suite, and someone would rightly revert
# it. Instead this module ALWAYS imports and ALWAYS collects; if the dependency
# is missing, every test in it fails loudly and locally, and every other test
# file — including ones that skip deliberately — is completely unaffected.
_MISSING_DEPS: list[str] = []
try:
    import spyrmsd as _spyrmsd_pkg  # noqa: F401
except ImportError as _exc:  # pragma: no cover - environment-dependent
    _MISSING_DEPS.append(f"spyrmsd ({_exc})")

_DEP_HINT = (
    "Install the pinned version with:  pip install 'spyrmsd==0.9.0'\n"
    "This repo has no requirements-test.txt; test dependencies are declared in\n"
    "the pip-install steps of .github/workflows/ci.yml, where spyrmsd==0.9.0 is\n"
    "pinned for the job that runs this file."
)


@pytest.fixture(autouse=True)
def _require_symmcorr_deps():
    """Turn a missing load-bearing dependency into a FAILURE, not a skip."""
    if _MISSING_DEPS:
        pytest.fail(
            "REQUIRED test dependency missing: "
            + "; ".join(_MISSING_DEPS)
            + ".\ntests/test_rmsd_symmcorr.py gates a fail-closed metric and "
            "must never be skipped — a skipped gate reports success without "
            "executing.\n" + _DEP_HINT,
            pytrace=False,
        )


def test_symmcorr_dependencies_are_installed():
    """Tripwire: this file must never again collect zero tests.

    If the guard above is ever reverted to `importorskip`, the module stops
    collecting and its silence looks identical to success. This test exists so
    that there is always at least one collected item whose failure is visible.
    """
    assert not _MISSING_DEPS, "missing: " + "; ".join(_MISSING_DEPS)


# ── fixtures: a molecule with a real automorphism ────────────────────────────
# O=C=O. The two oxygens are graph-equivalent, so swapping them is an
# automorphism: a pose with the oxygens exchanged is the SAME structure, and a
# correct metric must score it 0.

SDF = """CO2
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.1600    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1600    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  2  0
M  END
$$$$
"""


def _pdb_line(serial: int, name: str, element: str, x: float, y: float, z: float) -> str:
    ln = list(" " * 80)
    ln[0:6] = list("HETATM")
    ln[6:11] = list(f"{serial:5d}")
    ln[12:16] = list(f"{name:<4}")
    ln[17:20] = list("LIG")
    ln[21] = "A"
    ln[22:26] = list(f"{900:4d}")
    ln[30:38] = list(f"{x:8.3f}")
    ln[38:46] = list(f"{y:8.3f}")
    ln[46:54] = list(f"{z:8.3f}")
    ln[76:78] = list(f"{element:>2}")
    return "".join(ln).rstrip()


def _pose(coords) -> bytes:
    names = ("O1", "C1", "O2")
    els = ("O", "C", "O")
    body = "\n".join(
        _pdb_line(90001 + i, names[i % 3], els[i % 3], *c) for i, c in enumerate(coords)
    )
    return (body + "\nEND\n").encode()


@pytest.fixture()
def crystal(tmp_path: Path) -> str:
    p = tmp_path / "LIG_ligand.sdf"
    p.write_text(SDF)
    return str(p)


SWAPPED = [(1.16, 0.0, 0.0), (0.0, 0.0, 0.0), (-1.16, 0.0, 0.0)]
IDENTICAL = [(-1.16, 0.0, 0.0), (0.0, 0.0, 0.0), (1.16, 0.0, 0.0)]


# ── direction ────────────────────────────────────────────────────────────────


def test_automorphic_swap_scores_zero(crystal):
    """The defect in one test: serial calls an identical structure a failure."""
    pose = _pose(SWAPPED)
    res = symm.symmcorr_rmsd(crystal, pose)
    assert res["status"] == "ok"
    assert res["rmsd"] == pytest.approx(0.0, abs=1e-6)

    serial = symm.serial_rmsd(crystal, pose)
    assert serial > 1.0  # ordered mapping sees a 2.32 A displacement on each O
    assert res["rmsd"] < serial


def test_identity_pose_agrees_with_serial(crystal):
    pose = _pose(IDENTICAL)
    assert symm.symmcorr_rmsd(crystal, pose)["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert symm.serial_rmsd(crystal, pose) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("jitter", [0.0, 0.1, 0.4, 1.3])
def test_symmcorr_never_exceeds_serial(crystal, jitter):
    """The invariant the whole correction rests on, over perturbed poses."""
    rng = np.random.default_rng(20260829)
    for _ in range(25):
        base = np.array(SWAPPED if rng.random() < 0.5 else IDENTICAL, dtype=float)
        pose = _pose(base + rng.normal(0.0, jitter, base.shape))
        s = symm.symmcorr_rmsd(crystal, pose)
        if s["status"] != "ok":
            continue
        serial = symm.serial_rmsd(crystal, pose)
        assert s["rmsd"] <= serial + 1e-9, (s["rmsd"], serial)


def test_hungarian_never_exceeds_symmcorr(crystal):
    """Why rmsd_hungarian stays banned: it is a LOWER bound, not a variant.

    Element-only assignment minimises over all same-element bijections, a
    superset of the graph automorphisms, so it can only ever score at or below
    the chemically valid answer. The ordering must never invert.
    """
    scipy_opt = pytest.importorskip("scipy.optimize")
    pose_raw = _pose(SWAPPED)
    cx_all, an_all, _ = symm.parse_sdf(crystal)
    heavy = an_all > 1
    cx = cx_all[heavy]
    an = an_all[heavy]
    px = symm.pose_ligand_coords(pose_raw)

    asg = np.full(len(cx), -1)
    for el in set(an.tolist()):
        idx = np.where(an == el)[0]
        d = np.linalg.norm(cx[idx][:, None] - px[idx][None], axis=2)
        r, c = scipy_opt.linear_sum_assignment(d)
        for a, b in zip(r, c):
            asg[idx[a]] = idx[b]
    diff = cx - px[asg]
    hung = float(np.sqrt((diff * diff).sum(1).mean()))

    spy = symm.symmcorr_rmsd(crystal, pose_raw)["rmsd"]
    serial = symm.serial_rmsd(crystal, pose_raw)
    assert hung <= spy + 1e-9
    assert spy <= serial + 1e-9


# ── fail-closed ──────────────────────────────────────────────────────────────


def test_atom_count_mismatch_fails_closed(crystal):
    pose = _pose(SWAPPED[:2])
    res = symm.symmcorr_rmsd(crystal, pose)
    assert res["rmsd"] is None
    assert res["status"].startswith("atom_count_mismatch")


def test_no_ligand_serial_fails_closed(crystal):
    """Atoms below the serial>=90000 invariant are receptor, not ligand."""
    body = "\n".join(
        _pdb_line(i + 1, "C1", "C", *c) for i, c in enumerate(SWAPPED)
    ).encode()
    res = symm.symmcorr_rmsd(crystal, body)
    assert res["rmsd"] is None
    assert res["status"] == "no_ligand_serial_ge_90000"


def test_unavailable_spyrmsd_raises_rather_than_degrading(monkeypatch):
    """Absent spyrmsd must abort, never fall back to a weaker metric."""

    def boom():
        raise symm.SymmCorrUnavailable("simulated")

    monkeypatch.setattr(symm, "_spyrmsd", boom)
    with pytest.raises(symm.SymmCorrUnavailable):
        symm.symmcorr_rmsd("ignored", b"")


# ── aggregator wiring ────────────────────────────────────────────────────────


def test_symmcorr_supersedes_engine_serial_success_flag():
    """A stale engine success_rmsd=0 must not veto a corrected pass."""
    row = {
        "pdb_id": "1TZ8",
        "rmsd_to_crystal": "6.8360",
        "rmsd_symmcorr": "1.0871",
        "success_rmsd": "0",
        "seed_echo": "0",
    }
    val, metric = agg.elected_rmsd_labelled(row)
    assert metric == "symmcorr"
    assert val == pytest.approx(1.0871)
    assert agg.is_s1(row) is True


def test_serial_fallback_is_conservative_and_labelled():
    row = {
        "pdb_id": "1TZ8",
        "rmsd_to_crystal": "6.8360",
        "success_rmsd": "0",
        "seed_echo": "0",
    }
    val, metric = agg.elected_rmsd_labelled(row)
    assert metric == "serial"
    assert val == pytest.approx(6.8360)
    assert agg.is_s1(row) is False


def test_hungarian_is_never_consulted():
    row = {
        "pdb_id": "1TZ8",
        "rmsd_to_crystal": "6.8360",
        "rmsd_hungarian": "0.9755",
        "seed_echo": "0",
    }
    assert agg.is_s1(row) is False


def test_join_refuses_pose_sha_mismatch():
    rows = [{"pdb_id": "1TZ8", "pose_sha256": "aaaa"}]
    sidecar = {
        "1TZ8": {"pdb_id": "1TZ8", "rmsd_symmcorr": "1.0871", "pose_sha256": "bbbb"}
    }
    prov = agg.join_symmcorr(rows, sidecar)
    assert prov["joined"] == 0
    assert prov["refused_sha_mismatch"] == ["1TZ8"]
    assert agg.SYMMCORR_COL not in rows[0]


def test_join_accepts_matching_pose_sha():
    sha = hashlib.sha256(b"pose").hexdigest()
    rows = [{"pdb_id": "1TZ8", "pose_sha256": sha}]
    sidecar = {
        "1TZ8": {"pdb_id": "1TZ8", "rmsd_symmcorr": "1.0871", "pose_sha256": sha}
    }
    prov = agg.join_symmcorr(rows, sidecar)
    assert prov["joined"] == 1
    assert rows[0][agg.SYMMCORR_COL] == "1.0871"


def test_sidecar_loader_drops_non_ok_rows(tmp_path: Path):
    p = tmp_path / "s.csv"
    p.write_text(
        "pdb_id,rmsd_symmcorr,status,pose_sha256\n"
        "1AAA,1.0,ok,x\n"
        "1BBB,,no_pose_artifact,y\n"
        "1CCC,9.9,spyrmsd_error:ValueError,z\n"
    )
    loaded = agg.load_symmcorr_sidecar(p)
    assert set(loaded) == {"1AAA"}
