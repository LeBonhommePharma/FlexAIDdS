"""#509 item 1: cluster.cpp must reset Hungarian per emitted pose.

The cluster() result loop used a function-scoped `bool Hungarian = false`
and set `Hungarian = true` after the serial REMARK. It never reset the flag
inside `for (j = 0; j < num_of_results; ++j)`, so ranks ≥ 1 labeled BOTH
RMSD REMARK lines as Hungarian (the "no symmetry correction" line carried
the symmetry-corrected value).

BindingMode.cpp already scopes false→raw then true→sym per pose. This is a
source-contract test so the sticky-flag regression fails in CI without a
full engine rebuild. Do not treat it as a docking-success claim.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLUSTER = ROOT / "LIB" / "cluster.cpp"
BINDING_MODE = ROOT / "LIB" / "BindingMode.cpp"


def _strip_line_comments(text: str) -> str:
    return re.sub(r"//.*?$", "", text, flags=re.M)


def _brace_block(text: str, open_brace: int) -> str:
    assert text[open_brace] == "{", "open_brace must point at '{'"
    depth = 0
    for i, ch in enumerate(text[open_brace:], open_brace):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return text[open_brace : i + 1]
    raise AssertionError("unbalanced braces")


def _emission_loop(text: str) -> str:
    needle = "for(j=0;j<num_of_results;++j)"
    start = text.find(needle)
    assert start != -1, "cluster.cpp result loop not found"
    open_brace = text.find("{", start)
    assert open_brace != -1
    return _brace_block(text, open_brace)


def _refstructure_rmsd_block(loop: str) -> str:
    marker = "no symmetry correction"
    idx = loop.find(marker)
    assert idx != -1, "serial RMSD REMARK missing from emission loop"
    guard = loop.rfind("if(FA->refstructure == 1)", 0, idx)
    assert guard != -1, "refstructure guard missing before serial RMSD REMARK"
    open_brace = loop.find("{", guard)
    assert open_brace != -1
    return _brace_block(loop, open_brace)


def test_cluster_emission_loop_resets_hungarian_per_pose():
    """The result loop must assign Hungarian=false before serial calc_rmsd."""
    text = CLUSTER.read_text(encoding="utf-8")
    loop = _strip_line_comments(_emission_loop(text))
    resets = list(re.finditer(r"Hungarian\s*=\s*false\s*;", loop))
    assert resets, (
        "cluster.cpp result loop must reset Hungarian=false per pose (#509); "
        "otherwise ranks ≥1 inherit Hungarian=true and both REMARK lines are "
        "symmetry-corrected"
    )
    true_assign = re.search(r"Hungarian\s*=\s*true\s*;", loop)
    assert true_assign, "symmetry-corrected path must still set Hungarian=true"
    first_reset = resets[0].start()
    assert first_reset < true_assign.start()


def test_cluster_serial_rmsd_uses_hungarian_false_then_true():
    """Mirror BindingMode: false→raw REMARK, then true→sym REMARK, per pose."""
    text = CLUSTER.read_text(encoding="utf-8")
    block = _strip_line_comments(_refstructure_rmsd_block(_emission_loop(text)))
    raw = block.find("no symmetry correction")
    sym = block.find("(symmetry corrected)")
    assert 0 <= raw < sym
    i_false = re.search(r"Hungarian\s*=\s*false\s*;", block)
    i_true = re.search(r"Hungarian\s*=\s*true\s*;", block)
    assert i_false, "refstructure RMSD block must reset Hungarian=false before serial"
    assert i_true, "refstructure RMSD block must set Hungarian=true before symmetry"
    calcs = [m.start() for m in re.finditer(r"calc_rmsd\s*\(", block)]
    assert len(calcs) >= 2, "expected serial + Hungarian calc_rmsd pair"
    assert i_false.start() < calcs[0] < i_true.start() < calcs[1]
    assert i_false.start() < raw < i_true.start() < sym


def test_bindingmode_still_scopes_hungarian_per_pose():
    """#509: cluster.cpp should keep matching this BindingMode pattern."""
    text = BINDING_MODE.read_text(encoding="utf-8")
    locals_false = re.findall(r"bool\s+Hungarian\s*=\s*false\s*;", text)
    assert len(locals_false) >= 2, (
        "BindingMode.cpp must keep per-pose `bool Hungarian = false` "
        "(false→raw then true→sym)"
    )
