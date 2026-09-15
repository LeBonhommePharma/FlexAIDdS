"""Tests for scripts/safe_archive_to_icloud.py keep-local pack+SHA-256 verify.

Drives the shipped pack_and_verify entry (not a reimplementation). Local APFS
only — no CloudDocs walks.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import icloud_safe_io as iso  # noqa: E402
import safe_archive_to_icloud as sai  # noqa: E402


def _write_tree(root: Path) -> None:
    (root / "nested").mkdir(parents=True)
    (root / "a.txt").write_text("alpha\n", encoding="utf-8")
    (root / "nested" / "b.bin").write_bytes(b"\x00\x01\x02" * 50)
    (root / "nested" / "empty").write_bytes(b"")
    (root / "pose.pdb").write_text("ATOM      1  CA  ALA A   1\n", encoding="utf-8")


def test_pack_and_verify_roundtrip_sha256(tmp_path: Path):
    src = tmp_path / "campaign"
    _write_tree(src)
    archive = tmp_path / "campaign.tar.zst"
    extract = tmp_path / "extract"
    report = sai.pack_and_verify(src, archive, extract)
    assert report.ok is True
    assert report.n_source == 4
    assert report.n_matched == 4
    assert report.missing == []
    assert report.mismatched == []
    assert archive.is_file()
    # Every relative path's digest matches the source file (shipped compare).
    src_map = {rel: digest for rel, digest in report.matched}
    for rel, digest in src_map.items():
        expected = hashlib.sha256((src / rel).read_bytes()).hexdigest()
        assert digest == expected
        extracted = extract / rel
        assert extracted.is_file()
        assert hashlib.sha256(extracted.read_bytes()).hexdigest() == expected
    # Source payloads untouched
    assert (src / "a.txt").read_text(encoding="utf-8") == "alpha\n"


def test_pack_and_verify_fails_on_digest_mismatch(tmp_path: Path):
    src = tmp_path / "campaign"
    _write_tree(src)
    archive = tmp_path / "campaign.tar.zst"
    extract = tmp_path / "extract"
    sai.pack_directory(src, archive)
    (src / "a.txt").write_text("tampered\n", encoding="utf-8")
    report = sai.verify_archive(src, archive, extract)
    assert report.ok is False
    assert "a.txt" in report.mismatched
    assert (src / "a.txt").read_text(encoding="utf-8") == "tampered\n"


def test_pack_and_verify_fails_if_source_file_missing_from_archive(tmp_path: Path):
    src = tmp_path / "campaign"
    _write_tree(src)
    archive = tmp_path / "campaign.tar.zst"
    extract = tmp_path / "extract"
    sai.pack_directory(src, archive)
    extra = src / "late.txt"
    extra.write_text("not in archive\n", encoding="utf-8")
    report = sai.verify_archive(src, archive, extract)
    assert report.ok is False
    assert "late.txt" in report.missing


def test_refuses_clouddocs_source(tmp_path: Path):
    cloud = (
        Path.home()
        / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/nope"
    )
    with pytest.raises(ValueError, match="CloudDocs"):
        sai.iter_source_files(cloud)
    with pytest.raises(ValueError, match="CloudDocs"):
        sai.pack_directory(cloud, tmp_path / "x.tar.zst")


def test_keep_local_does_not_delete_source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    src = tmp_path / "campaign"
    _write_tree(src)
    dest = tmp_path / "archived"
    pack_dir = tmp_path / "packs"
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "local_root"))
    report = sai.archive_keep_local(
        src,
        dest,
        local_pack_dir=pack_dir,
        extract_parent=tmp_path / "ex",
        copy_icloud=True,
        keep_local_pack=True,
    )
    assert report.ok is True
    assert (src / "a.txt").is_file()
    assert (src / "nested" / "b.bin").is_file()
    packs = list(pack_dir.glob("*.tar.zst"))
    assert len(packs) == 1
    assert (dest / packs[0].name).is_file()


def test_cli_keep_local_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    src = tmp_path / "campaign"
    _write_tree(src)
    dest = tmp_path / "archived"
    pack_dir = tmp_path / "packs"
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "lr"))
    rc = sai._cli(
        [
            "--source",
            str(src),
            "--dest",
            str(dest),
            "--keep-local",
            "--local-pack-dir",
            str(pack_dir),
            "--extract-parent",
            str(tmp_path / "ex"),
        ]
    )
    assert rc == 0
    assert (src / "pose.pdb").is_file()
    assert list(pack_dir.glob("*.tar.zst"))
    assert list(dest.glob("*.tar.zst"))
    assert list(dest.glob("*.sha256"))
    assert list(dest.glob("*.manifest.json"))
