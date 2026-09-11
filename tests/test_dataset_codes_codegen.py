#!/usr/bin/env python3
"""The engine's PDB code lists must equal what their canonical sources produce.

Why this test exists
--------------------
``LIB/DatasetRunner.cpp`` held the Astex Diverse 85 as string literals, and that
literal list -- not any YAML -- was what ``--dataset astex_diverse`` executed. It
differed from the published set on 13 of 85 targets for months. Six files carried
the same 85 codes and nothing compared them, so a gate that diffed the YAMLs
against each other would have passed throughout.

The fix is generation, not comparison: ``scripts/gen_dataset_codes.py`` emits
``LIB/generated/dataset_codes.inc`` from ONE declared canonical source per dataset.
This test is the enforcement half -- it regenerates and diffs, so a hand-edit of the
generated header, or an edit to a canonical source without regenerating, fails.

Three separate obligations, three tests:

1. ``test_generated_file_matches_canonical_sources`` -- the committed generated file
   is what the generator produces NOW. Catches hand-edits and stale regeneration.
2. ``test_generated_lists_are_byte_identical_to_the_hardcoded_baseline`` -- the
   generated lists still equal the literal lists that were in DatasetRunner.cpp
   before the codegen landed. This is the acceptance criterion: a codegen that
   silently changes the executed set is strictly worse than the hardcoding. When a
   roster is switched ON PURPOSE, this test is the thing that must be updated in
   the same commit, by the person switching it.
3. ``test_cpp_no_longer_hardcodes_the_lists`` -- the literals are actually gone from
   the C++, so there is no second copy to drift.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "gen_dataset_codes.py"
GENERATED = ROOT / "LIB" / "generated" / "dataset_codes.inc"
BASELINE = ROOT / "tests" / "fixtures" / "dataset_codes" / "hardcoded_baseline.json"
CPP = ROOT / "LIB" / "DatasetRunner.cpp"

_SPEC = importlib.util.spec_from_file_location("gen_dataset_codes", SCRIPT)
gen = importlib.util.module_from_spec(_SPEC)
sys.modules["gen_dataset_codes"] = gen
_SPEC.loader.exec_module(gen)


def _codes_sha(codes):
    return hashlib.sha256(("\n".join(codes) + "\n").encode()).hexdigest()


def test_registry_is_not_empty():
    """A generator with an empty registry would pass every other test vacuously."""
    assert len(gen.REGISTRY) >= 3, "registry holds %d lists" % len(gen.REGISTRY)
    assert {e.slug for e in gen.REGISTRY} >= {"astex_diverse", "hap2", "casf2016"}


def test_generated_file_exists_and_carries_the_banner():
    assert GENERATED.is_file(), "%s is missing; run scripts/gen_dataset_codes.py" % GENERATED
    text = GENERATED.read_text(encoding="utf-8")
    assert text.startswith(gen.BANNER_MARK), "generated file must open with the DO-NOT-EDIT banner"
    for e in gen.REGISTRY:
        assert e.source in text, "banner does not name the canonical source for %s" % e.slug
        assert gen.sha256_file(ROOT / e.source) in text, (
            "banner does not carry the CURRENT sha256 of %s -- regenerate" % e.source)


def test_generated_file_matches_canonical_sources():
    """Regenerate and diff. Fails with both sha256 values, per the brief."""
    want = gen.render()
    have = GENERATED.read_text(encoding="utf-8")
    if have != want:
        pytest.fail(
            "LIB/generated/dataset_codes.inc does not match what its canonical sources "
            "produce.\n  committed sha256: %s\n  generated sha256: %s\n"
            "Run scripts/gen_dataset_codes.py and review the diff: a change here changes "
            "the set the engine executes."
            % (hashlib.sha256(have.encode()).hexdigest(),
               hashlib.sha256(want.encode()).hexdigest()))


def test_check_mode_returns_zero():
    """The CLI face CI calls must agree with the in-process comparison."""
    r = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, "rc=%d\nstdout:\n%s\nstderr:\n%s" % (
        r.returncode, r.stdout, r.stderr)


def test_check_mode_fails_on_a_tampered_generated_file(tmp_path):
    """Discrimination control. A --check that cannot fail proves nothing."""
    have = GENERATED.read_text(encoding="utf-8")
    backup = have
    try:
        GENERATED.write_text(have.replace('"1G9V"', '"9ZZZ"', 1), encoding="utf-8")
        r = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                           cwd=str(ROOT), capture_output=True, text=True, timeout=180)
        assert r.returncode == 1, "tampered generated file did NOT fail --check (rc=%d)" % r.returncode
        assert "committed sha256" in r.stdout and "generated sha256" in r.stdout, (
            "--check failure must print both sha256 values; got:\n%s" % r.stdout)
    finally:
        GENERATED.write_text(backup, encoding="utf-8")


def test_generated_lists_are_byte_identical_to_the_hardcoded_baseline():
    """The acceptance criterion for the codegen: nothing the engine runs changed."""
    assert BASELINE.is_file(), "baseline fixture missing: %s" % BASELINE
    base = json.loads(BASELINE.read_text(encoding="utf-8"))["lists"]
    assert base, "baseline fixture holds no lists -- this test would be vacuous"
    checked = 0
    for e in gen.REGISTRY:
        assert e.slug in base, (
            "%s has no hardcoded baseline. If this is a NEW list, add its baseline in the "
            "same commit so a later change to it is detectable." % e.slug)
        want = base[e.slug]["codes"]
        got = gen.read_codes(e)
        assert got == want, (
            "%s: generated list differs from the pre-codegen hardcoded list.\n"
            "  generated n=%d sha256=%s\n  baseline  n=%d sha256=%s\n"
            "  only generated: %s\n  only baseline:  %s\n"
            "This CHANGES THE SET THE ENGINE EXECUTES. If deliberate, update "
            "tests/fixtures/dataset_codes/hardcoded_baseline.json in the same commit."
            % (e.slug, len(got), _codes_sha(got), len(want), _codes_sha(want),
               sorted(set(got) - set(want))[:8], sorted(set(want) - set(got))[:8]))
        assert _codes_sha(got) == base[e.slug]["codes_sha256"]
        checked += len(got)
    assert checked >= 400, "only %d codes compared -- denominator too small to mean anything" % checked


def test_cpp_no_longer_hardcodes_the_lists():
    """No second copy of the lists in the C++, or the codegen bought nothing."""
    text = CPP.read_text(encoding="utf-8")
    for e in gen.REGISTRY:
        fn = e.cpp_function.split("::")[-1]
        m = re.search(r"std::vector<std::string> DatasetRunner::%s\(\)\s*\{" % re.escape(fn), text)
        assert m, "%s not found in LIB/DatasetRunner.cpp" % fn
        seg = text[m.end():]
        depth = 1
        for i, ch in enumerate(seg):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    seg = seg[:i]
                    break
        body = "\n".join(l for l in seg.splitlines() if not l.strip().startswith("//"))
        literals = re.findall(r'"([0-9][A-Za-z][A-Za-z0-9]{2})"', body)
        assert not literals, (
            "%s still holds %d hardcoded PDB literals (%s...) -- it must return the "
            "generated list" % (fn, len(literals), literals[:4]))
        assert 'dataset_codes("%s")' % e.slug in body, (
            "%s does not return the generated list for slug %r" % (fn, e.slug))


def test_every_list_declares_its_provenance_honestly():
    """A circular source must say so. The v1 Astex manifest was circular and silent."""
    allowed = {"ROSTER_CSV", "DEFINITION_YAML", "EXTRACTED_FROM_CODE"}
    for e in gen.REGISTRY:
        assert e.provenance in allowed, "%s: unknown provenance %r" % (e.slug, e.provenance)
        if e.provenance == "EXTRACTED_FROM_CODE":
            assert e.primary_source == "UNVERIFIED", (
                "%s was extracted from the code, so its primary source is not verified by "
                "anything; it must be labelled UNVERIFIED" % e.slug)
            assert "CIRCULAR" in e.note.upper(), (
                "%s: a code-extracted source must say it is circular in its note" % e.slug)


def test_generator_refuses_an_unknown_slug():
    """The engine must not run a list nobody declared."""
    with pytest.raises(SystemExit):
        gen.read_codes(gen.CodeList(
            slug="ghost", cpp_function="DatasetRunner::ghost_codes",
            source="benchmarks/datasets/does_not_exist.csv", reader="csv:pdb_id",
            expected_n=1, provenance="ROSTER_CSV", primary_source="UNVERIFIED", note=""))


def test_generator_refuses_a_source_whose_count_changed(tmp_path):
    """A roster that grew or shrank must not be emitted silently."""
    e = next(x for x in gen.REGISTRY if x.slug == "hap2")
    with pytest.raises(SystemExit) as ei:
        gen.read_codes(e._replace(expected_n=e.expected_n + 1))
    assert "expected_n" in str(ei.value) or "registry declares" in str(ei.value)
