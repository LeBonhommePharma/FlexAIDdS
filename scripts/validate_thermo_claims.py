#!/usr/bin/env python3
"""Fail-closed validation for the public entropy.help claims surface.

The public registry may contain planned, example, or explicitly unverified
records without result artifacts.  A record described as complete, published,
validated, verified, or reproducible is claim-bearing and must instead point to
three repository-local artifacts:

* a machine-readable JSON report;
* a human-readable Markdown report; and
* a separate JSON provenance record.

This validator does not decide whether a thermodynamic method is scientifically
valid.  It enforces the narrower precondition that public claims have inspectable
artifacts and that cryptographic-looking fields are not placeholders.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import unquote, urlparse


CLAIM_STATUS_WORDS = {
    "complete",
    "completed",
    "published",
    "reproducible",
    "validated",
    "verified",
}
NONCLAIM_STATUS_WORDS = {
    "candidate",
    "draft",
    "example",
    "pending",
    "planned",
    "proposed",
    "unverified",
}
REQUIRED_ARTIFACTS = ("report_json", "report_markdown", "provenance")
REQUIRED_PROVENANCE_FIELDS = (
    "git_sha",
    "timestamp",
    "runner",
    "binary_sha256",
    "raw_ensemble_digest",
)

SHA256_RE = re.compile(r"(?:sha256:)?([0-9a-fA-F]{64})\Z")
GIT_SHA_RE = re.compile(r"[0-9a-fA-F]{40}(?:[0-9a-fA-F]{24})?\Z")
AUDIT_LINK_RE = re.compile(r"(?:AUD-[A-Za-z0-9_-]+|/audits/|audits/).*\.(?:json|md)\Z", re.I)
AUDIT_ID_RE = re.compile(r"AUD-[A-Za-z0-9_-]+", re.I)
HTML_LINK_RE = re.compile(r"\bhref\s*=\s*[\"']([^\"']+)[\"']", re.I)
MARKDOWN_LINK_RE = re.compile(r"\]\(([^)]+)\)")
PUBLIC_CLAIM_RE = re.compile(
    r"\b(?:complete|completed|published|reproducible|validated|verified)\b",
    re.I,
)
CLAIM_CONTEXT_RE = re.compile(r"\b(?:audit|report|result|ledger)\b", re.I)
# The trailing guard has to be (?!\w), not \b: a \b after "%" requires a word
# character immediately after the percent sign, which never happens in prose, so
# the "%" alternative could never match.  Percentages are exactly the unit that
# benchmark success-rate claims are written in, so that alternative being dead
# was the largest hole in the quantitative check.
QUANTITATIVE_CLAIM_RE = re.compile(
    r"(?:\b\d+(?:\.\d+)?\s*(?:%|kcal(?:\s*mol)?(?:\^?-?1)?|Å|angstroms?)(?!\w)"
    r"|\bpearson\s+r\s*[≈=]\s*-?\d)",
    re.I,
)
THERMO_CONTEXT_RE = re.compile(
    r"(?:entropy|F_config|S_config|ΔG|delta\s*G|RMSD|affinit|binding mode)",
    re.I,
)
# A "safe qualifier" exempts a line from claim validation, so it must be an
# explicit non-claim marker or an explicit negation of the claim itself.  Bare
# function words (no / not / if / when / will / must / requires / until) used to
# live here, which meant any incidental "if" or "not" anywhere on a line silently
# disabled the quantitative-claim check for that line.  Only constructions that
# cannot occur by accident in an asserted result are accepted.
SAFE_QUALIFIER_RE = re.compile(
    r"(?:"
    r"\bnever\b|\bplanned\b|\bexample\b|\bexamples\b|\bunverified\b|\bdraft\b|"
    r"\bcandidate\b|\bproposed\b|\bpending\b|\btemplate\b|\bplaceholder\b|"
    r"\billustrative\b|\bhypothetical\b|\bnotional\b|\bTBD\b|"
    r"\bnot (?:yet|a|an|the|claimed|verified|validated|validation|published|"
    r"reproducible|provenance|physical|calibrated|consumed|enforced|wired|"
    r"interpretable|evidence|receipted)\b|"
    r"\bno (?:result|results|claim|claims|receipt|receipts|provenance|"
    r"calibration|physical|completed|complete|validated|verified|published|"
    r"reproducible)\b|"
    r"\bis not\b|\bare not\b|\bwas not\b|\bwere not\b|\bdoes not\b|\bdo not\b|"
    r"\bcannot\b|\bmust not\b|\bwill not\b|\bmay not\b|"
    r"\bwithout (?:a |an |any )?(?:receipt|receipts|provenance|calibration|"
    r"artifact|artifacts|evidence|verification)\b|"
    r"\bonly after\b|"
    r"\brequires? (?:a |an |the )?(?:receipt|receipts|provenance|calibration|"
    r"artifact|artifacts|evidence|matched)\b"
    r")",
    re.I,
)
PLACEHOLDER_RE = re.compile(
    r"(?:\.\.\.|\b(?:placeholder|fake|dummy|example|future|todo|tbd|replace[-_ ]?me)\b|"
    r"^[<\[].*[>\]]$)",
    re.I,
)

KNOWN_FAKE_DIGESTS = {
    # SHA-256 of empty content: useful as a sentinel, not evidence of a report.
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    # Historical entropy.help example filler, never derived from an ensemble.
    "3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
}
KNOWN_FAKE_GIT_SHAS = {
    "a1b2c3d4e5f6789012345678901234567890abcd",
}


def _words(value: Any) -> set[str]:
    return set(re.findall(r"[a-z]+", str(value).lower()))


def _record_status(record: dict[str, Any]) -> str:
    values = [
        record.get("status", ""),
        record.get("record_status", ""),
        record.get("verdict", ""),
    ]
    summary = record.get("audit_summary")
    if isinstance(summary, dict):
        values.append(summary.get("verdict", ""))
    return " ".join(str(value) for value in values if value is not None)


def _is_claim_bearing(record: dict[str, Any]) -> bool:
    for flag in ("complete", "completed", "published", "reproducible", "validated", "verified"):
        if record.get(flag) is True:
            return True

    words = _words(_record_status(record))
    if words & CLAIM_STATUS_WORDS:
        return True
    if words & NONCLAIM_STATUS_WORDS:
        return False

    # Legacy records with results or verdicts but no explicit status are claims.
    result_keys = {
        "f_config_kcal_mol",
        "s_config_kcal_mol_K",
        "total_sampled",
        "audit_summary",
    }
    return bool(result_keys & record.keys())


def _looks_placeholder(value: str) -> bool:
    stripped = value.strip()
    if not stripped or PLACEHOLDER_RE.search(stripped):
        return True
    compact = stripped.lower().removeprefix("sha256:")
    if compact in KNOWN_FAKE_DIGESTS or compact in KNOWN_FAKE_GIT_SHAS:
        return True
    if len(compact) >= 32 and len(set(compact)) <= 2:
        return True
    return False


def _repo_path_from_url(value: str, root: Path, base: Path) -> tuple[Path | None, str | None]:
    """Resolve a linked audit artifact to a repository-local path."""

    parsed = urlparse(value)
    if parsed.scheme in {"http", "https"}:
        path_parts = [unquote(part) for part in parsed.path.split("/") if part]
        if parsed.netloc == "github.com" and "blob" in path_parts:
            blob_index = path_parts.index("blob")
            if len(path_parts) > blob_index + 2:
                relative = Path(*path_parts[blob_index + 2 :])
                return (root / relative).resolve(), None
        if parsed.netloc == "raw.githubusercontent.com" and len(path_parts) >= 4:
            relative = Path(*path_parts[3:])
            return (root / relative).resolve(), None
        return None, "external audit links are not on-disk provenance"

    if parsed.scheme or parsed.netloc:
        return None, f"unsupported link scheme: {parsed.scheme or 'network path'}"

    raw_path = Path(unquote(parsed.path))
    if raw_path.is_absolute():
        root_relative = Path(*raw_path.parts[1:])
        if root_relative.parts and root_relative.parts[0] in {"docs", "site"}:
            candidate = (root / root_relative).resolve()
        else:
            return None, "absolute artifact paths are not portable"

    elif raw_path.parts and raw_path.parts[0] in {"docs", "site"}:
        candidate = (root / raw_path).resolve()
    else:
        candidate = (base / raw_path).resolve()

    try:
        candidate.relative_to(root)
    except ValueError:
        return None, "artifact path escapes repository root"
    return candidate, None


def _iter_json_values(value: Any, prefix: tuple[str, ...] = ()) -> Iterable[tuple[tuple[str, ...], Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            path = (*prefix, str(key))
            yield path, child
            yield from _iter_json_values(child, path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            path = (*prefix, str(index))
            # Scalar list members must be yielded too, otherwise a fake digest or
            # signature hidden inside an array is never inspected.
            yield path, child
            yield from _iter_json_values(child, path)


def _classifying_key(key_path: tuple[str, ...]) -> str:
    """Return the nearest named key, skipping list indices.

    ``{"digests": ["sha256:..."]}`` yields the member under path
    ``digests.0``; the member must still be classified as a digest, so numeric
    path components are transparent for classification purposes.
    """

    for part in reversed(key_path):
        if not part.isdigit():
            return part.lower()
    return key_path[-1].lower() if key_path else ""


def _validate_signature(value: Any, location: str, errors: list[str]) -> None:
    if value is None:
        return
    if isinstance(value, list):
        # Members are visited individually and inherit this classifying key.
        return
    if isinstance(value, dict):
        candidates = [
            value.get("value"),
            value.get("armored"),
            value.get("detached"),
            value.get("signature"),
        ]
        candidates = [candidate for candidate in candidates if isinstance(candidate, str)]
        if not candidates:
            errors.append(f"{location}: signature object has no signature payload")
            return
        for candidate in candidates:
            _validate_signature(candidate, location, errors)
        return
    if not isinstance(value, str):
        errors.append(f"{location}: signature must be a string, object, or null")
        return
    if _looks_placeholder(value):
        errors.append(f"{location}: placeholder/fake signature is forbidden")
        return
    if len(value.strip()) < 64:
        errors.append(f"{location}: signature payload is too short to be credible")
        return
    if "BEGIN PGP SIGNATURE" in value and "END PGP SIGNATURE" not in value:
        errors.append(f"{location}: incomplete armored PGP signature")


def _validate_crypto_fields(document: Any, source: Path, errors: list[str]) -> None:
    for key_path, value in _iter_json_values(document):
        key = _classifying_key(key_path)
        location = f"{source}:{'.'.join(key_path)}"

        # Any signature-like key, not just the exact leaf name "signature":
        # detached_signature / pgp_signature / signature_armored are signatures
        # and are not digests either, so exact equality left them unvalidated.
        if "signature" in key:
            _validate_signature(value, location, errors)
            continue

        if key == "git_sha" and value is not None:
            if isinstance(value, (dict, list)):
                continue
            if not isinstance(value, str) or _looks_placeholder(value):
                errors.append(f"{location}: placeholder/fake Git object ID is forbidden")
            elif GIT_SHA_RE.fullmatch(value.strip()) is None:
                errors.append(f"{location}: malformed Git object ID")
            continue

        is_digest = "digest" in key or key.endswith("sha256") or key == "content_hash"
        if not is_digest or value is None:
            continue
        if isinstance(value, (dict, list)):
            # Container members are inspected individually via _iter_json_values.
            continue
        if not isinstance(value, str):
            errors.append(f"{location}: digest must be a string or null")
            continue
        if _looks_placeholder(value):
            errors.append(f"{location}: placeholder/fake digest is forbidden")
            continue
        if SHA256_RE.fullmatch(value.strip()) is None:
            errors.append(f"{location}: expected 64 hexadecimal SHA-256 characters")


def _load_json(path: Path, errors: list[str]) -> Any | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        errors.append(f"{path}: invalid or unreadable JSON: {exc}")
        return None


def _validate_provenance(path: Path, audit_id: str, errors: list[str]) -> None:
    provenance = _load_json(path, errors)
    if not isinstance(provenance, dict):
        if provenance is not None:
            errors.append(f"{path}: provenance must be a JSON object")
        return

    for field in REQUIRED_PROVENANCE_FIELDS:
        value = provenance.get(field)
        if value is None or (isinstance(value, str) and _looks_placeholder(value)):
            errors.append(f"{audit_id}: provenance file {path} lacks real field '{field}'")

    git_sha = provenance.get("git_sha")
    if isinstance(git_sha, str) and GIT_SHA_RE.fullmatch(git_sha.strip()) is None:
        errors.append(f"{audit_id}: provenance git_sha is not a full Git object ID")

    for field in ("binary_sha256", "raw_ensemble_digest"):
        digest = provenance.get(field)
        if isinstance(digest, str) and SHA256_RE.fullmatch(digest.strip()) is None:
            errors.append(f"{audit_id}: provenance {field} is not a SHA-256 digest")


def _validate_registry(root: Path, registry_path: Path, errors: list[str]) -> set[str]:
    registry = _load_json(registry_path, errors)
    if registry is None:
        return set()
    if not isinstance(registry, dict) or not isinstance(registry.get("audits"), list):
        errors.append(f"{registry_path}: expected a top-level 'audits' list")
        return set()

    claim_ids: set[str] = set()
    for index, record in enumerate(registry["audits"]):
        if not isinstance(record, dict):
            errors.append(f"{registry_path}: audits[{index}] must be an object")
            continue
        audit_id = str(record.get("id") or record.get("audit_id") or f"audits[{index}]")
        claim_bearing = _is_claim_bearing(record)
        if claim_bearing:
            claim_ids.add(audit_id.upper())

        artifacts = record.get("artifacts")
        if not isinstance(artifacts, dict):
            artifacts = {}

        legacy_artifacts = {
            "report_json": record.get("json_url"),
            "report_markdown": record.get("md_url"),
            "provenance": record.get("provenance_url"),
        }

        resolved_artifacts: dict[str, Path] = {}
        for field in REQUIRED_ARTIFACTS:
            raw_value = artifacts.get(field, legacy_artifacts[field])
            if claim_bearing and not isinstance(raw_value, str):
                errors.append(f"{audit_id}: claim requires artifact '{field}'")
                continue
            if not isinstance(raw_value, str):
                continue
            resolved, reason = _repo_path_from_url(raw_value, root, registry_path.parent)
            if reason is not None or resolved is None:
                errors.append(f"{audit_id}: artifact '{field}' {reason}")
            elif not resolved.is_file():
                errors.append(f"{audit_id}: linked artifact '{field}' is missing: {resolved}")
            else:
                resolved_artifacts[field] = resolved

        # Optional legacy links are still promises for planned records. Claimed
        # records were already resolved above through the compatibility mapping.
        for field in ("json_url", "md_url", "provenance_url"):
            if claim_bearing:
                continue
            raw_value = record.get(field)
            if not isinstance(raw_value, str):
                continue
            resolved, reason = _repo_path_from_url(raw_value, root, registry_path.parent)
            if reason is not None or resolved is None:
                errors.append(f"{audit_id}: linked '{field}' {reason}")
            elif not resolved.is_file():
                errors.append(f"{audit_id}: linked artifact '{field}' is missing: {resolved}")

        report_json_path = resolved_artifacts.get("report_json")
        if claim_bearing and report_json_path is not None:
            report = _load_json(report_json_path, errors)
            if not isinstance(report, dict):
                if report is not None:
                    errors.append(f"{audit_id}: report_json must contain a JSON object")
            else:
                report_id = report.get("audit_id") or report.get("id")
                if report_id != audit_id:
                    errors.append(
                        f"{audit_id}: report_json identifies {report_id!r}, expected {audit_id!r}"
                    )

        provenance_path = resolved_artifacts.get("provenance")
        if claim_bearing and provenance_path is not None:
            _validate_provenance(provenance_path, audit_id, errors)

    return claim_ids


def _extract_links(path: Path) -> Iterable[str]:
    text = path.read_text(encoding="utf-8")
    yield from HTML_LINK_RE.findall(text)
    yield from MARKDOWN_LINK_RE.findall(text)


def _validate_public_links(root: Path, files: Iterable[Path], errors: list[str]) -> None:
    for source in files:
        for link in _extract_links(source):
            clean_link = link.split("#", 1)[0].split("?", 1)[0]
            if not AUDIT_LINK_RE.search(clean_link):
                continue
            resolved, reason = _repo_path_from_url(clean_link, root, source.parent)
            if reason is not None or resolved is None:
                errors.append(f"{source}: linked audit artifact {link!r} {reason}")
            elif not resolved.is_file():
                errors.append(f"{source}: linked audit artifact is missing: {resolved}")


def _visible_line(raw_line: str) -> str:
    without_tags = re.sub(r"<[^>]+>", " ", raw_line)
    return re.sub(r"\s+", " ", without_tags).strip()


def _validate_public_language(
    files: Iterable[Path], claim_ids: set[str], errors: list[str]
) -> None:
    for source in files:
        for line_number, raw_line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
            line = _visible_line(raw_line)
            if not line or SAFE_QUALIFIER_RE.search(line):
                continue
            referenced_ids = {value.upper() for value in AUDIT_ID_RE.findall(line)}
            unsupported_ids = referenced_ids - claim_ids
            # Evidence must be LOCAL to the claim.  Testing only
            # ``not claim_ids or unsupported_ids`` let a claim-bearing line that
            # cites no audit ID at all pass, because the empty citation set has
            # no unsupported members — the line rode on some unrelated valid
            # record elsewhere in the registry.  A claim-bearing line must cite
            # at least one ID itself, and every ID it cites must resolve to a
            # provenance-backed registry record.
            locally_unsupported = not referenced_ids or bool(unsupported_ids)
            if PUBLIC_CLAIM_RE.search(line) and CLAIM_CONTEXT_RE.search(line):
                if locally_unsupported:
                    errors.append(
                        f"{source}:{line_number}: completion/publication claim has no provenance-backed registry record"
                    )
            if QUANTITATIVE_CLAIM_RE.search(line) and THERMO_CONTEXT_RE.search(line):
                if locally_unsupported:
                    errors.append(
                        f"{source}:{line_number}: quantitative thermodynamic result has no provenance-backed registry record"
                    )


def _validate_text_crypto(files: Iterable[Path], errors: list[str]) -> None:
    digest_re = re.compile(r"sha256:\s*([^\s<]+)", re.I)
    pgp_re = re.compile(
        r"-----BEGIN PGP SIGNATURE-----(.*?)-----END PGP SIGNATURE-----",
        re.I | re.S,
    )
    for source in files:
        text = source.read_text(encoding="utf-8")
        for value in digest_re.findall(text):
            if _looks_placeholder(value) or SHA256_RE.fullmatch(f"sha256:{value}") is None:
                errors.append(f"{source}: placeholder/fake SHA-256 text is forbidden")
        for match in pgp_re.finditer(text):
            body = match.group(1).strip()
            if _looks_placeholder(body) or len(body) < 64:
                errors.append(f"{source}: placeholder/fake armored signature is forbidden")


def validate_repository(root: Path) -> list[str]:
    """Return all entropy.help claim-validation errors under *root*."""

    root = root.resolve()
    docs_dir = root / "docs" / "entropy-help"
    site_dir = root / "site" / "entropy-help"
    registry_path = docs_dir / "audits" / "audits.json"
    errors: list[str] = []

    for required in (docs_dir, site_dir, registry_path):
        if not required.exists():
            errors.append(f"required claims surface is missing: {required}")
    if errors:
        return errors

    claim_ids = _validate_registry(root, registry_path, errors)

    public_files = sorted(site_dir.rglob("*.html"))
    public_files.extend(sorted(docs_dir.rglob("*.md")))
    public_files.extend(sorted(docs_dir.rglob("*.txt")))
    # The repository landing page and the primary scoring/user documentation are
    # public claim surfaces too: they are the first thing a reader sees and they
    # historically carried unreceipted benchmark and kcal/mol language.  Scanning
    # only docs/entropy-help left the loudest surface unguarded.
    for extra in (
        root / "README.md",
        root / "docs" / "SCORING.md",
        root / "docs" / "USERGUIDE.md",
    ):
        if extra.is_file():
            public_files.append(extra)

    _validate_public_links(root, public_files, errors)
    _validate_public_language(public_files, claim_ids, errors)
    _validate_text_crypto(public_files, errors)

    for json_path in sorted(docs_dir.rglob("*.json")):
        document = _load_json(json_path, errors)
        if document is not None:
            _validate_crypto_fields(document, json_path, errors)

    return errors


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate entropy.help public claims against on-disk audit provenance."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="repository root (default: parent of scripts/)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    errors = validate_repository(args.root)
    if errors:
        print("THERMODYNAMIC CLAIM VALIDATION FAILED", file=sys.stderr)
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"{len(errors)} error(s)", file=sys.stderr)
        return 1

    print("THERMODYNAMIC CLAIM VALIDATION PASSED")
    print("Public audit records are either explicitly non-claiming or provenance-backed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
