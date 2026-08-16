from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.validate_thermo_claims import validate_repository


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _make_claim_surface(tmp_path: Path, audits: list[dict[str, object]]) -> tuple[Path, Path, Path]:
    registry = tmp_path / "docs" / "entropy-help" / "audits" / "audits.json"
    site = tmp_path / "site" / "entropy-help"
    site.mkdir(parents=True)
    (site / "index.html").write_text(
        "<html><body><p>Planned audit framework; no result claimed.</p></body></html>\n",
        encoding="utf-8",
    )
    _write_json(registry, {"version": 1, "audits": audits})
    return tmp_path, registry, site


def _deposit_provenance_backed_record(registry: Path, audit_id: str) -> None:
    """Write one genuinely provenance-backed record plus its three artifacts.

    Used by the loophole tests: the registry must contain a *valid* record so
    that ``claim_ids`` is non-empty, which is precisely the condition the old
    validator let unrelated claim lines ride on.
    """

    audit_dir = registry.parent
    report_json = audit_dir / f"{audit_id}.json"
    report_markdown = audit_dir / f"{audit_id}.md"
    provenance = audit_dir / f"{audit_id}.provenance.json"

    _write_json(
        report_json,
        {"audit_id": audit_id, "record_status": "PUBLISHED", "signature": None},
    )
    report_markdown.write_text(
        f"# {audit_id}\n\nResult summary backed by the adjacent provenance record.\n",
        encoding="utf-8",
    )
    _write_json(
        provenance,
        {
            "git_sha": hashlib.sha1(b"source revision").hexdigest(),
            "timestamp": "2026-08-08T12:00:00Z",
            "runner": "fixture-runner --input deposited-ensemble.json",
            "binary_sha256": hashlib.sha256(b"binary bytes").hexdigest(),
            "raw_ensemble_digest": hashlib.sha256(b"ensemble bytes").hexdigest(),
        },
    )
    _write_json(
        registry,
        {
            "version": 1,
            "audits": [
                {
                    "id": audit_id,
                    "status": "PUBLISHED_REPRODUCIBLE",
                    "artifacts": {
                        "report_json": report_json.name,
                        "report_markdown": report_markdown.name,
                        "provenance": provenance.name,
                    },
                }
            ],
        },
    )


def test_planned_unverified_record_needs_no_result_artifacts(tmp_path: Path) -> None:
    root, _, _ = _make_claim_surface(
        tmp_path,
        [
            {
                "id": "AUD-PLAN-001",
                "target": "candidate complex",
                "status": "PLANNED_UNVERIFIED",
                "claim": "No result claimed.",
            }
        ],
    )

    assert validate_repository(root) == []


def test_published_claim_requires_all_three_provenance_artifacts(tmp_path: Path) -> None:
    root, _, _ = _make_claim_surface(
        tmp_path,
        [
            {
                "id": "AUD-CLAIM-001",
                "target": "claimed complex",
                "status": "PUBLISHED_REPRODUCIBLE",
            }
        ],
    )

    errors = validate_repository(root)

    for field in ("report_json", "report_markdown", "provenance"):
        assert any(f"requires artifact '{field}'" in error for error in errors)


def test_missing_linked_audit_artifact_fails_even_for_planned_record(tmp_path: Path) -> None:
    root, _, site = _make_claim_surface(
        tmp_path,
        [
            {
                "id": "AUD-PLAN-002",
                "status": "PLANNED_UNVERIFIED",
                "json_url": "AUD-PLAN-002.json",
            }
        ],
    )
    (site / "index.html").write_text(
        '<html><body><a href="AUD-PLAN-002.json">planned artifact</a></body></html>\n',
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("linked artifact 'json_url' is missing" in error for error in errors)
    assert any("linked audit artifact is missing" in error for error in errors)


def test_placeholder_digest_and_signature_fail(tmp_path: Path) -> None:
    root, registry, _ = _make_claim_surface(tmp_path, [])
    _write_json(
        registry.parent.parent / "fake-example.json",
        {
            "record_status": "EXAMPLE_UNVERIFIED",
            "raw_ensemble_digest": "sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
            "signature": {
                "method": "manual-review + future-gpg",
                "auditor": "example",
            },
        },
    )

    errors = validate_repository(root)

    assert any("placeholder/fake digest" in error for error in errors)
    assert any("signature object has no signature payload" in error for error in errors)


def test_public_completion_language_requires_claim_record(tmp_path: Path) -> None:
    root, _, site = _make_claim_surface(tmp_path, [])
    (site / "index.html").write_text(
        "<html><body><p>The AUD-TEST audit is published and reproducible.</p></body></html>\n",
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("completion/publication claim has no provenance-backed" in error for error in errors)


def test_provenance_backed_published_fixture_passes(tmp_path: Path) -> None:
    audit_id = "AUD-VALID-001"
    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, audit_id)
    (site / "index.html").write_text(
        f"<html><body><p>The {audit_id} audit is published and reproducible.</p></body></html>\n",
        encoding="utf-8",
    )

    assert validate_repository(root) == []


# ---------------------------------------------------------------------------
# Hostile tests: each of these passed against the pre-firewall validator.
# ---------------------------------------------------------------------------


def test_uncited_completion_claim_cannot_borrow_an_unrelated_record(tmp_path: Path) -> None:
    """A claim line that cites no audit ID must not ride on someone else's receipt.

    The old guard tested ``not claim_ids or unsupported_ids``.  With a valid
    record present, ``claim_ids`` was non-empty and an uncited line had an empty
    ``unsupported_ids``, so the guard passed for a claim with zero local evidence.
    """

    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, "AUD-VALID-001")
    (site / "index.html").write_text(
        "<html><body><p>The thermodynamic audit is published and reproducible.</p></body></html>\n",
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("completion/publication claim has no provenance-backed" in e for e in errors)


def test_uncited_quantitative_claim_cannot_borrow_an_unrelated_record(tmp_path: Path) -> None:
    """Same loophole on the quantitative branch: evidence must be local to the claim."""

    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, "AUD-VALID-001")
    (site / "index.html").write_text(
        "<html><body><p>Configurational entropy recovers 91.8% of binding modes.</p></body></html>\n",
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("quantitative thermodynamic result has no provenance-backed" in e for e in errors)


def test_claim_citing_only_an_unknown_audit_id_still_fails(tmp_path: Path) -> None:
    """Citing an ID is not enough — the cited ID has to resolve to a real record."""

    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, "AUD-VALID-001")
    (site / "index.html").write_text(
        "<html><body><p>The AUD-GHOST-999 audit is published and reproducible.</p></body></html>\n",
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("completion/publication claim has no provenance-backed" in e for e in errors)


def test_signature_like_keys_are_not_skipped(tmp_path: Path) -> None:
    """``detached_signature`` / ``pgp_signature`` are signatures, not stray strings.

    The old check was ``if key == "signature"``, so any decorated signature key
    was skipped by the signature branch and — not being a digest either — was
    never validated at all.
    """

    root, registry, _ = _make_claim_surface(tmp_path, [])
    _write_json(
        registry.parent.parent / "detached-signature-example.json",
        {
            "record_status": "EXAMPLE_UNVERIFIED",
            "detached_signature": "REPLACE_ME",
            "pgp_signature": "not-a-real-signature",
        },
    )

    errors = validate_repository(root)

    assert any(
        "detached_signature: placeholder/fake signature is forbidden" in e for e in errors
    )
    assert any(
        "pgp_signature: signature payload is too short to be credible" in e for e in errors
    )


def test_incidental_qualifier_cannot_exempt_a_quantitative_claim(tmp_path: Path) -> None:
    """A stray "if" must not disable claim validation for the whole line.

    ``SAFE_QUALIFIER_RE`` used to contain bare ``\\bif\\b`` / ``\\bnot\\b`` /
    ``\\bwill\\b``, so an incidental function word exempted an asserted result.
    """

    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, "AUD-VALID-001")
    (site / "index.html").write_text(
        "<html><body><p>If the ensemble converges, the binding mode entropy shifts "
        "the result by 3.5 kcal/mol.</p></body></html>\n",
        encoding="utf-8",
    )

    errors = validate_repository(root)

    assert any("quantitative thermodynamic result has no provenance-backed" in e for e in errors)


def test_explicit_nonclaim_marker_still_exempts_a_line(tmp_path: Path) -> None:
    """The tightened qualifier must keep exempting genuine non-claims."""

    root, registry, site = _make_claim_surface(tmp_path, [])
    _deposit_provenance_backed_record(registry, "AUD-VALID-001")
    (site / "index.html").write_text(
        "<html><body><p>Planned example only: no result is claimed for the "
        "binding mode entropy shift of 3.5 kcal/mol.</p></body></html>\n",
        encoding="utf-8",
    )

    assert validate_repository(root) == []


def test_fake_digest_inside_a_list_is_not_invisible(tmp_path: Path) -> None:
    """Scalar members of arrays are inspected, and inherit their field's key."""

    root, registry, _ = _make_claim_surface(tmp_path, [])
    _write_json(
        registry.parent.parent / "array-digests.json",
        {
            "record_status": "EXAMPLE_UNVERIFIED",
            "raw_ensemble_digests": [
                "sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a"
            ],
        },
    )

    errors = validate_repository(root)

    assert any("placeholder/fake digest is forbidden" in e for e in errors)


def test_primary_docs_do_not_claim_unwired_physical_election() -> None:
    root = Path(__file__).resolve().parents[1]
    surfaces = {
        "README.md": (root / "README.md").read_text(encoding="utf-8"),
        "docs/SCORING.md": (root / "docs" / "SCORING.md").read_text(
            encoding="utf-8"
        ),
        "docs/USERGUIDE.md": (root / "docs" / "USERGUIDE.md").read_text(
            encoding="utf-8"
        ),
    }
    forbidden = (
        "full binding free energy reported as `G_bind`",
        "ensuring they can never be elected rank-0",
        "ranks binding modes by Helmholtz free energy",
        "Total ΔG (F + vib corr)",
    )

    for name, text in surfaces.items():
        for phrase in forbidden:
            assert phrase not in text, f"{name} retains unsupported claim: {phrase}"

    ranking = (root / "docs" / "classic_entropy_ranking.md").read_text(
        encoding="utf-8"
    )
    header = (root / "LIB" / "BindingMode.h").read_text(encoding="utf-8")
    assert "Vibrational entropy stays." not in ranking
    assert "fail-closed" in ranking
    assert "Ranking F also adds ENCoM vib correction" not in header
    assert "fail-closes to 0.0 and does not elect" in header


def test_runtime_proxy_outputs_declare_their_domain() -> None:
    root = Path(__file__).resolve().parents[1]
    binding_mode = (root / "LIB" / "BindingMode.cpp").read_text(encoding="utf-8")
    classic_cluster = (root / "LIB" / "cluster.cpp").read_text(encoding="utf-8")
    gaboom = (root / "LIB" / "gaboom.cpp").read_text(encoding="utf-8")
    top = (root / "LIB" / "top.cpp").read_text(encoding="utf-8")

    assert binding_mode.count("thermo_claim_validity = proxy_only") >= 2
    assert "thermo_claim_validity = proxy_only" in classic_cluster
    assert "Physical StatMech ledger" not in classic_cluster
    assert "claim_validity=proxy_only energy_domain=cf_arbitrary_units" in gaboom
    assert "enforced_in_final_election=0" in gaboom
    assert "Post-GA CF-proxy ensemble diagnostics" in top
    assert "Post-GA Ensemble Thermodynamics" not in top
    assert "fail_closed: no eigenvalue channel" in binding_mode
    assert "format_vibrational_diagnostic_remark" in binding_mode
    assert "CF.elec_gist_con_status = gated_inert_on_claim_path" in binding_mode
    assert "p_bind_like" in top
    assert top.count("claim_validity=proxy_only") >= 4
