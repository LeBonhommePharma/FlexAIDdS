#!/usr/bin/env python3
"""A run must record the roster it executed, and that record must be checkable.

Pytest face of scripts/check_run_receipt_codes.py.

Why this exists
---------------
RUN_RECEIPT.json recorded matrix hash, binary hash, GA knobs and a full
ProtocolConfig snapshot -- everything about HOW a run was configured and nothing
about WHAT set of targets it enumerated. A roster that was 72/85 of the published
Astex Diverse Set ran for months, and no receipt it wrote could have contradicted
the claim that it was the published set.

Two halves have to hold together, and a test that only exercises the Python half
would pass while the C++ emitted nothing:

  * the C++ emitter (LIB/RunReceipt.cpp) writes the keys the validator requires,
    asserted here against the source so the contract cannot drift silently;
  * the validator passes a truthful receipt and fails every way of lying, proved
    by its own --selftest arms.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_run_receipt_codes.py"
RECEIPT_CPP = ROOT / "LIB" / "RunReceipt.cpp"
RECEIPT_H = ROOT / "LIB" / "RunReceipt.h"

_SPEC = importlib.util.spec_from_file_location("check_run_receipt_codes", SCRIPT)
CRC = importlib.util.module_from_spec(_SPEC)
sys.modules["check_run_receipt_codes"] = CRC
_SPEC.loader.exec_module(CRC)

# Keys the validator reads out of a receipt. Every one must be emitted by the C++.
_CONTRACT = (
    "executed_codes", "executed_codes_count", "executed_codes_sha256",
    "executed_codes_file", "declared_codes_count", "dataset_source_path",
    "dataset_source_sha256", "dataset_source_provenance", "dataset_source_primary",
)


def test_script_exists():
    assert SCRIPT.is_file(), "scripts/check_run_receipt_codes.py is missing"


def test_cpp_emits_every_key_the_validator_reads():
    """The two halves must agree, or the validator checks fields nobody writes."""
    src = RECEIPT_CPP.read_text(encoding="utf-8")
    emitted = set(re.findall(r'\\"([a-z_0-9]+)\\":', src))
    assert emitted, "no JSON keys found in LIB/RunReceipt.cpp -- the regex moved"
    missing = [k for k in _CONTRACT if k not in emitted]
    assert not missing, (
        "LIB/RunReceipt.cpp does not emit %s, but scripts/check_run_receipt_codes.py "
        "reads them. A validator that checks fields nobody writes is vacuous." % missing)


def test_schema_version_was_bumped_for_the_roster_fields():
    """The bump is what distinguishes a pre-roster run from a stripped receipt."""
    h = RECEIPT_H.read_text(encoding="utf-8")
    m = re.search(r"kRunReceiptSchemaVersion\s*=\s*(\d+)", h)
    assert m, "kRunReceiptSchemaVersion not found"
    assert int(m.group(1)) >= CRC.MIN_SCHEMA, (
        "receipt schema is %s but the validator requires >= %d"
        % (m.group(1), CRC.MIN_SCHEMA))


def test_selftest_is_non_vacuous():
    r = subprocess.run([sys.executable, str(SCRIPT), "--selftest"],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, "selftest rc=%d\n%s\n%s" % (r.returncode, r.stdout, r.stderr)
    assert "NON-VACUOUS" in r.stdout
    m = re.search(r"performed (\d+) field comparisons", r.stdout)
    assert m and int(m.group(1)) > 0, (
        "selftest reported zero field comparisons; the arms prove nothing:\n%s" % r.stdout)


def _write_run(tmp_path: Path, body: dict) -> Path:
    d = tmp_path / "run"
    d.mkdir(exist_ok=True)
    (d / "RUN_RECEIPT.json").write_text(json.dumps(body, indent=2), encoding="utf-8")
    if body.get("executed_codes"):
        (d / "executed_codes.txt").write_text(
            "\n".join(body["executed_codes"]) + "\n", encoding="utf-8")
    return d


def _truthful(src_rel="benchmarks/astex_diverse/astex_diverse_set.csv") -> dict:
    codes = CRC.read_source_codes(src_rel, ROOT)
    assert codes, "cannot read %s" % src_rel
    return CRC._synth_receipt(codes, src_rel, ROOT)


def test_truthful_receipt_passes_with_a_nonzero_denominator(tmp_path):
    d = _write_run(tmp_path, _truthful())
    findings, info = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    assert not findings, [f._asdict() for f in findings]
    assert int(info["checks"]) >= 5, "only %s comparisons performed" % info["checks"]


def test_a_code_the_source_does_not_declare_fails(tmp_path):
    body = _truthful()
    body["executed_codes"] = body["executed_codes"] + ["9ZZZ"]
    body["executed_codes_count"] = len(body["executed_codes"])
    body["executed_codes_sha256"] = CRC.sha256_text("\n".join(body["executed_codes"]) + "\n")
    d = _write_run(tmp_path, body)
    findings, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    kinds = {f.kind for f in findings}
    assert "ROSTER_MISMATCH" in kinds, kinds
    assert any(CRC.is_failure(k) for k in kinds)


def test_a_tampered_digest_fails(tmp_path):
    body = _truthful()
    body["executed_codes_sha256"] = "0" * 64
    d = _write_run(tmp_path, body)
    findings, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    assert "SIDECAR_SHA_MISMATCH" in {f.kind for f in findings}


def test_a_sidecar_edited_after_the_run_fails(tmp_path):
    """The digest is over a file on disk, so editing the file must be detectable."""
    body = _truthful()
    d = _write_run(tmp_path, body)
    (d / "executed_codes.txt").write_text("1G9V\n", encoding="utf-8")
    findings, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    assert "SIDECAR_SHA_MISMATCH" in {f.kind for f in findings}


def test_a_pre_schema2_receipt_is_uncheckable_not_passing(tmp_path):
    """Absence of evidence must not read as evidence of correctness."""
    d = _write_run(tmp_path, {"schema_version": 1, "dataset": "astex_diverse"})
    findings, info = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    assert {f.kind for f in findings} == {"UNCHECKABLE_SCHEMA"}
    assert int(info["checks"]) == 0
    assert not CRC.is_failure("UNCHECKABLE_SCHEMA")
    strict, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT, strict=True)
    assert CRC.is_failure(strict[0].kind), "--strict must turn it into a failure"


def test_a_schema2_receipt_with_the_roster_stripped_fails(tmp_path):
    d = _write_run(tmp_path, {"schema_version": 2, "dataset": "astex_diverse"})
    findings, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    assert "MISSING_FIELD" in {f.kind for f in findings}


def test_a_short_roster_is_reported_with_both_numbers(tmp_path):
    """A silently short roster is how a 72/85 set passes for 85."""
    body = _truthful()
    body["executed_codes"] = body["executed_codes"][:10]
    body["executed_codes_count"] = 10
    body["executed_codes_sha256"] = CRC.sha256_text("\n".join(body["executed_codes"]) + "\n")
    d = _write_run(tmp_path, body)
    findings, _ = CRC.check_receipt(d / "RUN_RECEIPT.json", ROOT)
    inc = [f for f in findings if f.kind == "INCOMPLETE_RUN"]
    assert inc, {f.kind for f in findings}
    assert "10 of 85" in inc[0].detail, inc[0].detail


def test_every_generated_dataset_source_is_readable_by_the_validator():
    """The validator must be able to read every source the codegen declares, or a
    real receipt for that dataset would report SOURCE_SHA_MISMATCH spuriously."""
    gen_spec = importlib.util.spec_from_file_location(
        "gen_dataset_codes", ROOT / "scripts" / "gen_dataset_codes.py")
    gen = importlib.util.module_from_spec(gen_spec)
    gen_spec.loader.exec_module(gen)
    n = 0
    for e in gen.REGISTRY:
        codes = CRC.read_source_codes(e.source, ROOT)
        assert codes, "validator cannot read the declared source for %s (%s)" % (e.slug, e.source)
        assert len(codes) == e.expected_n, (
            "%s: validator read %d codes from %s, registry declares %d"
            % (e.slug, len(codes), e.source, e.expected_n))
        n += len(codes)
    assert n >= 400, "only %d codes reachable through the validator's readers" % n
