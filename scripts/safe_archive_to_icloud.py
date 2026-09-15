#!/usr/bin/env python3
"""Keep-local pack + SHA-256 verify, then optional timeout-bounded iCloud copy.

Contract
--------
* Pack and verify entirely on local APFS (never ``rglob`` / ``find`` CloudDocs).
* Round-trip: decompress/extract and compare SHA-256 of every file by relative
  path against the source tree. Fail closed on missing or mismatched digest.
* ``--keep-local`` never deletes, moves, or rewrites campaign payloads.
  Eviction is not implemented (safer-than-safe).
* iCloud receives *packed archives + sidecar manifests only*, copied one file
  at a time via ``icloud_safe_io.safe_copy_file``.

Callers: ``scripts/reproduce_astex85.sh``, ``compare_astex_2015_fair_vs_full_ds.py``.

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import icloud_safe_io as iso  # noqa: E402

ICLOUD_STANDARD = "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"

# Directories we never pack when walking a whole results root (Mach-O / scratch).
SKIP_DIR_NAMES = frozenset(
    {
        "_safe_archives",
        "pins",
        "three_engine_entropy_q1",
        "__pycache__",
        ".git",
    }
)
SKIP_DIR_PREFIXES = ("build_", "archive_scratch_")


@dataclass
class VerifyReport:
    ok: bool
    n_source: int
    n_matched: int
    missing: List[str] = field(default_factory=list)
    mismatched: List[str] = field(default_factory=list)
    extra: List[str] = field(default_factory=list)
    archive: str = ""
    source: str = ""
    seconds: float = 0.0
    matched: List[Tuple[str, str]] = field(default_factory=list)  # rel, sha256

    def sidecar_sha256_text(self) -> str:
        lines = [f"{digest}  {rel}" for rel, digest in sorted(self.matched)]
        return "\n".join(lines) + ("\n" if lines else "")


def get_icloud_base() -> Path:
    """Return iCloud base dir for FlexAIDdS benchmarks.

    Honors ``FLEXAIDDS_ICLOUD``. If that path's name is not
    ``FlexAIDdS_benchmarks``, append it. Falls back to the standard CloudDocs
    container. Does not walk the tree.
    """
    env = os.environ.get("FLEXAIDDS_ICLOUD", "").strip()
    if env:
        p = Path(env).expanduser()
        if p.name != "FlexAIDdS_benchmarks":
            p = p / "FlexAIDdS_benchmarks"
        return p
    return Path.home() / ICLOUD_STANDARD


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_source_files(
    root: Path,
    *,
    skip_names: Iterable[str] = (),
    direct_files_only: bool = False,
) -> List[str]:
    """Sorted relative POSIX paths of regular files under *root* (local APFS).

    Raises ``ValueError`` if *root* is under CloudDocs. Does not follow
    symlinks. Skips nested CloudDocs path components if they appear.
    """
    root = Path(root)
    if iso.is_clouddocs(root):
        raise ValueError(f"refusing to walk CloudDocs source: {root}")
    skip = set(skip_names)
    rels: List[str] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        if "Mobile Documents" in dirpath:
            dirnames[:] = []
            continue
        if direct_files_only:
            dirnames[:] = []
        else:
            dirnames[:] = [
                d
                for d in dirnames
                if d not in skip
                and d not in SKIP_DIR_NAMES
                and not d.startswith(SKIP_DIR_PREFIXES)
                and "Mobile Documents" not in d
            ]
        for fn in filenames:
            full = Path(dirpath) / fn
            try:
                full.lstat()
            except OSError:
                continue
            if os.path.islink(full) or not os.path.isfile(full):
                continue
            rel = full.relative_to(root).as_posix()
            rels.append(rel)
    rels.sort()
    return rels


def _open_tar(path: Path, mode: str) -> tarfile.TarFile:
    # Python 3.14 supports w:zst / r:zst; r:* auto-detects gz/bz2/xz/zst.
    return tarfile.open(path, mode)


def pack_directory(
    src: Path,
    archive: Path,
    *,
    direct_files_only: bool = False,
) -> Path:
    """Write a compressed tar of regular files under *src*. Returns *archive*.

    Format is ``.tar.zst`` when the suffix contains ``zst``, else gzip.
    Source must be local APFS. Never deletes *src*.
    """
    src = Path(src)
    archive = Path(archive)
    if iso.is_clouddocs(src):
        raise ValueError(f"refusing to pack CloudDocs source: {src}")
    if iso.is_clouddocs(archive):
        raise ValueError(f"refusing to write archive onto CloudDocs: {archive}")
    rels = iter_source_files(src, direct_files_only=direct_files_only)
    archive.parent.mkdir(parents=True, exist_ok=True)
    tmp = archive.with_name(archive.name + ".part")
    if tmp.exists():
        tmp.unlink()
    suffix = archive.name
    if suffix.endswith(".tar.zst") or suffix.endswith(".tzst"):
        mode = "w:zst"
    elif suffix.endswith(".tar.gz") or suffix.endswith(".tgz"):
        mode = "w:gz"
    else:
        mode = "w:zst"
    with _open_tar(tmp, mode) as tf:
        for rel in rels:
            tf.add(src / rel, arcname=rel, recursive=False)
    tmp.replace(archive)
    return archive


def verify_archive(
    src: Path,
    archive: Path,
    extract_dir: Path,
    *,
    direct_files_only: bool = False,
) -> VerifyReport:
    """Extract *archive* into *extract_dir* and compare SHA-256 by relative path.

    Fails (``ok=False``) if any source file is missing from the extract or the
    digest differs. Extra files in the extract are also recorded as failure.
    """
    src = Path(src)
    archive = Path(archive)
    extract_dir = Path(extract_dir)
    t0 = time.time()
    if iso.is_clouddocs(src) or iso.is_clouddocs(archive) or iso.is_clouddocs(extract_dir):
        raise ValueError("verify_archive requires local APFS source, archive, and extract dir")
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)
    with _open_tar(archive, "r:*") as tf:
        try:
            tf.extractall(extract_dir, filter="data")
        except TypeError:
            tf.extractall(extract_dir)

    source_rels = iter_source_files(src, direct_files_only=direct_files_only)
    extract_rels = iter_source_files(extract_dir)
    src_set = set(source_rels)
    ext_set = set(extract_rels)
    missing = sorted(src_set - ext_set)
    extra = sorted(ext_set - src_set)
    mismatched: List[str] = []
    matched: List[Tuple[str, str]] = []
    for rel in source_rels:
        if rel in missing:
            continue
        src_h = sha256_file(src / rel)
        ext_h = sha256_file(extract_dir / rel)
        if src_h != ext_h:
            mismatched.append(rel)
        else:
            matched.append((rel, src_h))
    ok = not missing and not mismatched and not extra and len(matched) == len(source_rels)
    return VerifyReport(
        ok=ok,
        n_source=len(source_rels),
        n_matched=len(matched),
        missing=missing,
        mismatched=mismatched,
        extra=extra,
        archive=str(archive),
        source=str(src),
        seconds=time.time() - t0,
        matched=matched,
    )


def pack_and_verify(
    src: Path,
    archive: Path,
    extract_dir: Path,
    *,
    direct_files_only: bool = False,
) -> VerifyReport:
    """Shipped entry: pack *src* to *archive*, extract, SHA-256 compare.

    On verify failure the partial archive is removed. Source is never touched.
    """
    src = Path(src)
    archive = Path(archive)
    pack_directory(src, archive, direct_files_only=direct_files_only)
    report = verify_archive(
        src, archive, extract_dir, direct_files_only=direct_files_only
    )
    if not report.ok:
        try:
            archive.unlink(missing_ok=True)
        except TypeError:
            if archive.exists():
                archive.unlink()
        part = archive.with_name(archive.name + ".part")
        if part.exists():
            part.unlink()
    return report


def write_sidecars(archive: Path, report: VerifyReport) -> Tuple[Path, Path]:
    sha_path = archive.with_suffix(archive.suffix + ".sha256")
    if archive.name.endswith(".tar.zst"):
        sha_path = Path(str(archive)[: -len(".tar.zst")] + ".sha256")
    elif archive.name.endswith(".tar.gz"):
        sha_path = Path(str(archive)[: -len(".tar.gz")] + ".sha256")
    man_path = Path(str(sha_path)[: -len(".sha256")] + ".manifest.json")
    sha_path.write_text(report.sidecar_sha256_text(), encoding="utf-8")
    payload = {
        "ok": report.ok,
        "source": report.source,
        "archive": report.archive,
        "n_source": report.n_source,
        "n_matched": report.n_matched,
        "missing": report.missing,
        "mismatched": report.mismatched,
        "extra": report.extra,
        "seconds": report.seconds,
        "keep_local": True,
    }
    man_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return sha_path, man_path


def dir_nbytes(root: Path) -> int:
    total = 0
    for rel in iter_source_files(root):
        try:
            total += (root / rel).lstat().st_size
        except OSError:
            pass
    return total


def _skip_dir(name: str) -> bool:
    if name in SKIP_DIR_NAMES:
        return True
    return name.startswith(SKIP_DIR_PREFIXES)


@dataclass(frozen=True)
class PackUnit:
    path: Path
    direct_files_only: bool = False


def plan_pack_units(root: Path, *, max_bytes: int) -> List[PackUnit]:
    """Split *root* into local subtrees that fit in *max_bytes* for extract-verify."""
    root = Path(root)
    if iso.is_clouddocs(root):
        raise ValueError(f"refusing to plan CloudDocs: {root}")

    def _walk(p: Path, depth: int) -> List[PackUnit]:
        if _skip_dir(p.name) and p != root:
            return []
        nbytes = dir_nbytes(p)
        if nbytes <= max_bytes or depth >= 6:
            return [PackUnit(p)] if nbytes > 0 else []
        try:
            entries = sorted(p.iterdir(), key=lambda x: x.name)
        except OSError:
            return [PackUnit(p)] if nbytes > 0 else []
        subdirs = [
            e
            for e in entries
            if e.is_dir() and not e.is_symlink() and not _skip_dir(e.name)
        ]
        files = [e for e in entries if e.is_file() and not e.is_symlink()]
        if not subdirs:
            return [PackUnit(p)]
        kids: List[PackUnit] = []
        for d in subdirs:
            kids.extend(_walk(d, depth + 1))
        if files:
            kids.append(PackUnit(p, direct_files_only=True))
        return kids

    units = _walk(root, 0)
    seen = set()
    out: List[PackUnit] = []
    for u in units:
        key = (str(u.path), u.direct_files_only)
        if key in seen:
            continue
        seen.add(key)
        out.append(u)
    return out


def _archive_name_for(src: Path, root: Path, *, direct_files_only: bool = False) -> str:
    try:
        rel = src.relative_to(root).as_posix()
    except ValueError:
        rel = src.name
    if rel in (".", ""):
        rel = src.name
    safe = rel.replace("/", "__").replace(" ", "_")
    if direct_files_only:
        safe = f"{safe}__rootfiles"
    return f"{safe}.tar.zst"


def archive_keep_local(
    src: Path,
    dest_dir: Path,
    *,
    local_pack_dir: Optional[Path] = None,
    extract_parent: Optional[Path] = None,
    copy_icloud: bool = True,
    keep_local_pack: bool = True,
    direct_files_only: bool = False,
    name_root: Optional[Path] = None,
) -> VerifyReport:
    """Pack+verify on local APFS, write sidecars, optionally copy pack to dest.

    *dest_dir* may be CloudDocs; the pack itself is always written locally first.
    """
    src = Path(src)
    dest_dir = Path(dest_dir)
    if iso.is_clouddocs(src):
        raise ValueError(f"refusing CloudDocs source: {src}")
    local_pack_dir = Path(local_pack_dir or (iso.local_root() / "_safe_archives"))
    local_pack_dir.mkdir(parents=True, exist_ok=True)
    archive = local_pack_dir / _archive_name_for(
        src, name_root or src.parent, direct_files_only=direct_files_only
    )
    extract_parent = Path(extract_parent or local_pack_dir / "_extract_tmp")
    extract_dir = extract_parent / (archive.stem + "_rt")
    sha_guess = Path(str(archive)[: -len(".tar.zst")] + ".sha256") if archive.name.endswith(".tar.zst") else archive.with_suffix(".sha256")
    man_guess = Path(str(sha_guess)[: -len(".sha256")] + ".manifest.json")
    if archive.is_file() and sha_guess.is_file() and man_guess.is_file():
        try:
            prior = json.loads(man_guess.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            prior = {}
        if prior.get("ok") is True:
            report = VerifyReport(
                ok=True,
                n_source=int(prior.get("n_source") or 0),
                n_matched=int(prior.get("n_matched") or 0),
                archive=str(archive),
                source=str(src),
            )
            sha_path, man_path = sha_guess, man_guess
            # Fall through to iCloud copy of existing verified pack.
            if copy_icloud:
                if iso.is_clouddocs(dest_dir):
                    if not iso.safe_mkdir(dest_dir, timeout_s=30.0):
                        report.ok = False
                        report.extra.append("ICLOUD_MKDIR_TIMEOUT")
                        return report
                else:
                    dest_dir.mkdir(parents=True, exist_ok=True)
                for local_f in (archive, sha_path, man_path):
                    dest_f = dest_dir / local_f.name
                    timeout = max(120.0, local_f.stat().st_size / (2 * 1024 * 1024) + 60.0)
                    n = iso.safe_copy_file(local_f, dest_f, timeout_s=timeout)
                    if n is None:
                        report.ok = False
                        report.extra.append(f"ICLOUD_COPY_TIMEOUT:{local_f.name}")
                        return report
            return report
    report = pack_and_verify(
        src, archive, extract_dir, direct_files_only=direct_files_only
    )
    shutil.rmtree(extract_dir, ignore_errors=True)
    if not report.ok:
        return report
    sha_path, man_path = write_sidecars(archive, report)
    if copy_icloud:
        if iso.is_clouddocs(dest_dir):
            if not iso.safe_mkdir(dest_dir, timeout_s=30.0):
                report.ok = False
                report.extra.append("ICLOUD_MKDIR_TIMEOUT")
                return report
        else:
            dest_dir.mkdir(parents=True, exist_ok=True)
        for local_f in (archive, sha_path, man_path):
            dest_f = dest_dir / local_f.name
            timeout = max(120.0, local_f.stat().st_size / (2 * 1024 * 1024) + 60.0)
            n = iso.safe_copy_file(local_f, dest_f, timeout_s=timeout)
            if n is None:
                report.ok = False
                report.extra.append(f"ICLOUD_COPY_TIMEOUT:{local_f.name}")
                return report
        if not keep_local_pack:
            # Only the packed bytes — never the source tree.
            for p in (archive,):
                try:
                    p.unlink()
                except OSError:
                    pass
    return report


def _cli(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Keep-local tar.zst pack with SHA-256 round-trip verify; "
        "optional timeout-bounded copy of packs to iCloud."
    )
    p.add_argument("--source", required=True, help="Local APFS directory to pack")
    p.add_argument(
        "--dest",
        required=True,
        help="Directory for verified archives (may be iCloud archived/)",
    )
    p.add_argument(
        "--keep-local",
        action="store_true",
        help="Never delete source payloads (the only supported mode; implied).",
    )
    p.add_argument(
        "--local-pack-dir",
        default="",
        help="Local APFS directory for pack+sidecars (default: $FLEXAIDDS_LOCAL_ROOT/_safe_archives)",
    )
    p.add_argument(
        "--children",
        action="store_true",
        help="Pack each immediate child of --source as its own archive",
    )
    p.add_argument(
        "--max-unit-bytes",
        type=int,
        default=3 * 1024 * 1024 * 1024,
        help="Split subtrees larger than this so extract-verify fits on disk",
    )
    p.add_argument(
        "--no-icloud-copy",
        action="store_true",
        help="Pack+verify locally only (still writes dest if dest is local)",
    )
    p.add_argument(
        "--extract-parent",
        default="",
        help="Local scratch parent for extract-verify (deleted after each unit)",
    )
    args = p.parse_args(argv)

    src = Path(args.source).expanduser()
    dest = Path(args.dest).expanduser()
    if not args.keep_local:
        print(
            "[INFO] --keep-local not passed; eviction is not implemented. Source will not be deleted.",
            file=sys.stderr,
        )
    if not src.is_dir():
        print(f"FATAL: source is not a local directory: {src}", file=sys.stderr)
        return 2
    if iso.is_clouddocs(src):
        print(f"FATAL: source is CloudDocs: {src}", file=sys.stderr)
        return 2

    local_pack = Path(args.local_pack_dir).expanduser() if args.local_pack_dir else (
        iso.local_root() / "_safe_archives"
    )
    extract_parent = (
        Path(args.extract_parent).expanduser()
        if args.extract_parent
        else (local_pack / "_extract_tmp")
    )

    if args.children:
        units: List[PackUnit] = []
        for child in sorted(src.iterdir(), key=lambda x: x.name):
            if child.is_symlink():
                continue
            if child.is_dir() and _skip_dir(child.name):
                print(f"[SKIP] {child.name} (binaries/scratch)")
                continue
            if child.is_dir():
                units.extend(plan_pack_units(child, max_bytes=args.max_unit_bytes))
            elif child.is_file():
                continue
        root_files = [
            e for e in src.iterdir() if e.is_file() and not e.is_symlink()
        ]
        if root_files:
            units.append(PackUnit(src, direct_files_only=True))
    else:
        units = plan_pack_units(src, max_bytes=args.max_unit_bytes)

    copy_icloud = not args.no_icloud_copy
    n_ok = 0
    n_bad = 0
    n_skip = 0
    for unit in units:
        if not unit.path.is_dir():
            n_skip += 1
            continue
        nfiles = len(
            iter_source_files(unit.path, direct_files_only=unit.direct_files_only)
        )
        if nfiles == 0:
            n_skip += 1
            continue
        print(
            f"[PACK] {unit.path}  files={nfiles}  direct_only={unit.direct_files_only}"
        )
        report = archive_keep_local(
            unit.path,
            dest,
            local_pack_dir=local_pack,
            extract_parent=extract_parent,
            copy_icloud=copy_icloud,
            keep_local_pack=True,
            direct_files_only=unit.direct_files_only,
            name_root=src,
        )
        if report.ok:
            n_ok += 1
            print(
                f"  VERIFIED matched={report.n_matched}/{report.n_source} "
                f"in {report.seconds:.1f}s  archive={report.archive}"
            )
        else:
            n_bad += 1
            print(
                f"  VERIFY-FAIL missing={report.missing[:5]} mismatched={report.mismatched[:5]} "
                f"extra={report.extra[:5]} — SOURCE UNTOUCHED"
            )
    print(f"done verified={n_ok} failed={n_bad} skipped={n_skip}")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(_cli())
