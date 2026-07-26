"""Contract: DUMP_POP on modern AUTONOMOUS path sets refstructure via RMSDST.

The prior instrument failure: FLEXAIDDS_DUMP_POP=1 alone never wrote .pop.tsv
because FA->refstructure stayed 0 on the JSON/direct path (no classic RMSDST).
These tests pin the wiring that fixes it — they exercise the shipped source
text so a revert breaks CI, not only a manual dock.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_native_score_exports_dump_pop_loader():
    h = (ROOT / "LIB" / "native_score.h").read_text(encoding="utf-8")
    cpp = (ROOT / "LIB" / "native_score.cpp").read_text(encoding="utf-8")
    assert "bool load_dump_pop_refstructure(" in h
    assert "load_dump_pop_refstructure" in cpp
    assert "FLEXAIDDS_DUMP_POP" in cpp
    assert "FLEXAIDDS_RMSDST" in cpp
    # Must set refstructure and fill coor_ref; must not overwrite coor[]
    assert "FA->refstructure = 1" in cpp or "FA->refstructure=1" in cpp
    assert "coor_ref" in cpp
    # Audit-only: do not assign atoms[i].coor from crystal in this function
    # (score_native_pose may; dump loader must only touch coor_ref).
    fn = cpp.split("bool load_dump_pop_refstructure")[1].split("\nbool ")[0]
    # crude: no "atoms[i].coor[" assignments in dump loader body
    for line in fn.splitlines():
        s = line.strip()
        if s.startswith("//"):
            continue
        assert "atoms[i].coor[" not in s, line


def test_top_calls_loader_before_native_score():
    top = (ROOT / "LIB" / "top.cpp").read_text(encoding="utf-8")
    assert "load_dump_pop_refstructure(FA, atoms, residue)" in top
    i_dump = top.index("load_dump_pop_refstructure")
    i_native = top.index("Native-pose CF diagnostic")
    assert i_dump < i_native


def test_dataset_runner_injects_rmsdst_for_dump_pop_autonomous():
    dr = (ROOT / "LIB" / "DatasetRunner.cpp").read_text(encoding="utf-8")
    assert "dump_pop_on" in dr
    assert "FLEXAIDDS_DUMP_POP" in dr
    # AUTONOMOUS-safe branch must inject RMSDST without SCORE_NATIVE
    assert "injecting FLEXAIDDS_RMSDST for .pop.tsv audit" in dr
    # Ensure SCORE_NATIVE is not forced in the dump-only branch
    dump_branch = dr.split("else if (dump_pop_on")[1].split("cmd << \"OMP_NUM_THREADS")[0]
    assert "FLEXAIDDS_RMSDST=" in dump_branch
    assert "FLEXAIDDS_SCORE_NATIVE" not in dump_branch


def test_cluster_dump_still_gated_on_refstructure_and_env():
    cl = (ROOT / "LIB" / "cluster.cpp").read_text(encoding="utf-8")
    assert "FLEXAIDDS_DUMP_POP" in cl
    assert ".pop.tsv" in cl
    # Gate remains dual: refstructure AND env
    assert "FA->refstructure == 1" in cl
