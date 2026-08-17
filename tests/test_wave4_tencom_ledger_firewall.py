#!/usr/bin/env python3
"""Wave 4: tENCoM λ must stay ledger-only (inert on election)."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _fn_body(text: str, name: str) -> str:
    token = f"{name}("
    start = text.find(token)
    assert start >= 0, f"missing {name}"
    brace = text.find("{", start)
    assert brace >= 0
    depth = 0
    for i, ch in enumerate(text[brace:], brace):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return text[brace : i + 1]
    raise AssertionError(f"unclosed {name}")


def test_ledger_flag_defaults_off():
    header = (ROOT / "LIB" / "tencom_ledger.h").read_text(encoding="utf-8")
    assert 'env_bool("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", false)' in header
    flags = (ROOT / "LIB" / "flexaidds_flags.cpp").read_text(encoding="utf-8")
    assert "FLEXAIDDS_LEDGER_TENCOM_LAMBDA" in flags


def test_compute_energy_does_not_consume_tencom_lambda():
    text = (ROOT / "LIB" / "BindingMode.cpp").read_text(encoding="utf-8")
    body = _fn_body(text, "BindingMode::compute_energy")
    assert "collect_tencom_lambda" not in body
    assert "tencom_lambda_ledger" not in body
    assert "lambda_min" not in body
    assert "format_tencom_lambda_remark" not in body


def test_compute_vibrational_correction_stays_fail_closed_zero():
    text = (ROOT / "LIB" / "BindingMode.cpp").read_text(encoding="utf-8")
    body = _fn_body(text, "BindingMode::compute_vibrational_correction")
    assert "collect_tencom_lambda" not in body
    assert "vib_correction_cache_ = 0.0" in body


def test_cluster_acf_does_not_add_tencom_lambda():
    text = (ROOT / "LIB" / "cluster.cpp").read_text(encoding="utf-8")
    assert "Clus_ACF" in text
    for line in text.splitlines():
        if "Clus_ACF[" in line and "=" in line and "tencom" in line.lower():
            raise AssertionError(f"tENCoM λ fed into cluster ACF ranking: {line}")


def test_remark_writers_tag_inert_on_election():
    bm = (ROOT / "LIB" / "BindingMode.cpp").read_text(encoding="utf-8")
    cl = (ROOT / "LIB" / "cluster.cpp").read_text(encoding="utf-8")
    assert "append_tencom_lambda_ledger_remark" in bm
    assert "format_tencom_lambda_remark" in cl
    header = (ROOT / "LIB" / "tencom_ledger.h").read_text(encoding="utf-8")
    assert "inert_on_election=%d" in header
    assert "kTencomLambdaInertOnElection" in header
