"""Bonhomme Fleet: resumable DatasetRunner orchestration with iCloud archival.

The controller uses only the Python standard library. Docking remains in the
version-pinned ``benchmark_datasets`` and ``FlexAIDdS`` binaries.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = 1
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
SAFE_ENV_KEYS = {
    "EVAL_SCALE_DIHEDRAL",
    "FLEXAIDDS_HVIB",
    "FLEXAIDDS_NATIVE_SEED_FRAC",
    "FLEXAIDDS_PARALLEL_RESTARTS",
    "FLEXAIDDS_POSEBUST_BACKEND",
    "FLEXAIDDS_POSEBUSTERS_BIN",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
    "FLEXAIDDS_RESTARTS",
    "FLEXAIDDS_SEED_ELITISM",
    "FLEXAIDDS_TENCOM_BIN",
    "FLEXAIDDS_THERMO",
    "FLEXAIDDS_VCT_R0",
    "OMP_NUM_THREADS",
    "SHARING_ALPHA",
}
SECRET_KEY_PATTERN = re.compile(r"(?:TOKEN|SECRET|PASSWORD|CREDENTIAL|API_KEY)", re.I)


class FleetError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _parse_time(value: str) -> float:
    return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()


def _canonical_json(data: Any) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    try:
        fd = os.open(str(path), os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


def atomic_write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{uuid.uuid4().hex}")
    payload = _canonical_json(data)
    with temporary.open("xb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    _fsync_directory(path.parent)


def exclusive_write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(_canonical_json(data))
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        try:
            path.unlink()
        except OSError:
            pass
        raise
    _fsync_directory(path.parent)


def read_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise FleetError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise FleetError(f"expected JSON object in {path}")
    return data


def confined_path(root: Path, *parts: str) -> Path:
    base = root.resolve()
    candidate = base.joinpath(*parts).resolve(strict=False)
    try:
        common = Path(os.path.commonpath((str(base), str(candidate))))
    except ValueError as exc:
        raise FleetError("path is outside campaign root") from exc
    if common != base:
        raise FleetError("path is outside campaign root")
    return candidate


def validate_id(value: str, label: str) -> str:
    if not ID_PATTERN.fullmatch(value):
        raise FleetError(f"invalid {label}: {value!r}")
    return value


def parse_codes(path: Path) -> List[str]:
    codes: List[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].replace(",", " ")
        codes.extend(token.upper() for token in line.split())
    if not codes:
        raise FleetError(f"no target codes found in {path}")
    if len(set(codes)) != len(codes):
        raise FleetError("target code list contains duplicates")
    for code in codes:
        validate_id(code, "target code")
    return codes


def parse_environment(items: Sequence[str]) -> Dict[str, str]:
    environment: Dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise FleetError(f"environment override must be KEY=VALUE: {item!r}")
        key, value = item.split("=", 1)
        if SECRET_KEY_PATTERN.search(key) or key not in SAFE_ENV_KEYS:
            raise FleetError(f"environment key is not allowed in Fleet manifests: {key}")
        environment[key] = value
    return dict(sorted(environment.items()))


def _manifest_paths(campaign: Path) -> Tuple[Path, Path]:
    return campaign / "manifest.json", campaign / "manifest.sha256"


def load_manifest(campaign: Path) -> Tuple[Dict[str, Any], str]:
    manifest_path, hash_path = _manifest_paths(campaign)
    manifest = read_json(manifest_path)
    try:
        expected = hash_path.read_text(encoding="ascii").strip()
    except OSError as exc:
        raise FleetError(f"cannot read manifest hash: {exc}") from exc
    actual = sha256_file(manifest_path)
    if actual != expected:
        raise FleetError(f"campaign manifest hash mismatch: expected {expected}, got {actual}")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise FleetError("unsupported Fleet manifest schema")
    return manifest, actual


def plan_campaign(args: argparse.Namespace) -> Dict[str, Any]:
    campaign = Path(args.campaign).expanduser().resolve()
    if campaign.exists() and any(campaign.iterdir()):
        raise FleetError(f"campaign directory is not empty: {campaign}")

    campaign_id = validate_id(args.campaign_id or campaign.name, "campaign ID")
    runner = Path(args.runner).expanduser().resolve()
    engine = Path(args.engine).expanduser().resolve() if args.engine else runner.with_name("FlexAIDdS")
    for label, binary in (("runner", runner), ("engine", engine)):
        if not binary.is_file():
            raise FleetError(f"{label} binary does not exist: {binary}")
    posebusters = Path(args.posebusters_bin).expanduser().resolve() if args.posebusters_bin else None
    if posebusters is None:
        env_bust = os.environ.get("FLEXAIDDS_POSEBUSTERS_BIN")
        discovered = env_bust or shutil.which("bust")
        posebusters = Path(discovered).expanduser().resolve() if discovered else None
    if posebusters is None or not posebusters.is_file() or not os.access(posebusters, os.X_OK):
        raise FleetError("official PoseBusters 'bust' executable is required; use --posebusters-bin")

    codes = parse_codes(Path(args.codes_file).expanduser().resolve())
    if args.chunks < 1 or args.chunks > len(codes):
        raise FleetError("--chunks must be between 1 and the number of target codes")
    shards = [codes[index::args.chunks] for index in range(args.chunks)]
    chunks = [
        {"id": f"chunk-{index:04d}", "index": index, "codes": shard}
        for index, shard in enumerate(shards)
    ]

    compute_root = Path(args.compute_root).expanduser().resolve()
    environment = parse_environment(args.env)
    required_environment = {
        "FLEXAIDDS_HVIB": "1",
        "FLEXAIDDS_NATIVE_SEED_FRAC": "0",
        "FLEXAIDDS_POSEBUST_BACKEND": "bust_cli",
        "FLEXAIDDS_POSEBUSTERS_BIN": str(posebusters),
        "FLEXAIDDS_SEED_ELITISM": "0",
    }
    for key, required_value in required_environment.items():
        if key in environment and environment[key] != required_value:
            raise FleetError(f"production Fleet requires {key}={required_value}")
        environment[key] = required_value

    if args.threads < 1 or args.omp_threads < 1:
        raise FleetError("thread counts must be positive")
    if args.ga_population < 1 or args.ga_generations < 1:
        raise FleetError("GA population and generations must be positive")
    if args.lease_seconds < 30 or args.max_attempts < 1 or args.min_free_gb < 0:
        raise FleetError("lease must be >=30 s, attempts >=1, and free-space guard non-negative")

    manifest: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "type": "flexaidds-fleet-campaign",
        "campaign_id": campaign_id,
        "created_at": utc_now(),
        "dataset": args.dataset or args.benchmark,
        "benchmark": args.benchmark,
        "runner": {"path": str(runner), "sha256": sha256_file(runner)},
        "engine": {"path": str(engine), "sha256": sha256_file(engine)},
        "validators": {
            "posebusters": {
                "backend": "bust_cli",
                "path": str(posebusters),
                "sha256": sha256_file(posebusters),
            },
            "tencom_eigen": {"backend": "dataset_runner_builtin", "required": True},
        },
        "protocol": {
            "mode": args.mode,
            "threads": args.threads,
            "omp_threads": args.omp_threads,
            "ga_population": args.ga_population,
            "ga_generations": args.ga_generations,
            "temperature_k": args.temperature,
            "grid_spacing_a": args.grid_spacing,
            "job_timeout_s": args.job_timeout_seconds,
            "clustering": args.clustering,
            "gpu_backend": args.gpu,
            "cache": str(Path(args.cache).expanduser().resolve()) if args.cache else None,
        },
        "environment": dict(sorted(environment.items())),
        "compute_root": str(compute_root),
        "lease_seconds": args.lease_seconds,
        "max_attempts": args.max_attempts,
        "min_free_bytes": int(args.min_free_gb * 1024**3),
        "chunks": chunks,
    }
    if manifest["protocol"]["mode"] == "oracle-ceiling":
        raise FleetError("oracle-ceiling is diagnostic-only and is not accepted by Fleet production planning")

    campaign.mkdir(parents=True, exist_ok=True)
    for directory in ("chunks", "attempts", "heartbeats", "artifacts", "quarantine", "aggregate"):
        (campaign / directory).mkdir()
    manifest_path, hash_path = _manifest_paths(campaign)
    exclusive_write_json(manifest_path, manifest)
    manifest_hash = sha256_file(manifest_path)
    with hash_path.open("x", encoding="ascii") as stream:
        stream.write(manifest_hash + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    _fsync_directory(campaign)
    return {"campaign": str(campaign), "manifest_sha256": manifest_hash, "chunks": len(chunks)}


def _chunk_dir(campaign: Path, chunk_id: str) -> Path:
    validate_id(chunk_id, "chunk ID")
    path = confined_path(campaign, "chunks", chunk_id)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _attempt_files(campaign: Path, chunk_id: str) -> List[Path]:
    directory = confined_path(campaign, "attempts", chunk_id)
    return sorted(directory.glob("*.json")) if directory.is_dir() else []


def _heartbeat_path(campaign: Path, claim: Mapping[str, Any]) -> Path:
    return confined_path(
        campaign, "heartbeats", str(claim["chunk_id"]), f"{claim['lease_token']}.json"
    )


def lease_owned(campaign: Path, claim: Mapping[str, Any]) -> bool:
    claim_path = _chunk_dir(campaign, str(claim["chunk_id"])) / "claim.json"
    try:
        current = read_json(claim_path)
    except FleetError:
        return False
    return current.get("lease_token") == claim.get("lease_token") and current.get("lease_epoch") == claim.get("lease_epoch")


def write_heartbeat(campaign: Path, claim: Mapping[str, Any]) -> bool:
    if not lease_owned(campaign, claim):
        return False
    atomic_write_json(_heartbeat_path(campaign, claim), {
        "campaign_id": claim["campaign_id"],
        "chunk_id": claim["chunk_id"],
        "attempt_id": claim["attempt_id"],
        "lease_token": claim["lease_token"],
        "lease_epoch": claim["lease_epoch"],
        "worker_id": claim["worker_id"],
        "timestamp": utc_now(),
    })
    return lease_owned(campaign, claim)


def claim_is_stale(campaign: Path, claim: Mapping[str, Any], lease_seconds: int) -> bool:
    heartbeat = _heartbeat_path(campaign, claim)
    timestamp = str(claim.get("claimed_at", "1970-01-01T00:00:00Z"))
    if heartbeat.is_file():
        try:
            timestamp = str(read_json(heartbeat).get("timestamp", timestamp))
        except FleetError:
            pass
    try:
        age = time.time() - _parse_time(timestamp)
    except (TypeError, ValueError):
        return True
    return age > lease_seconds


def _new_claim(manifest: Mapping[str, Any], chunk_id: str, worker_id: str, epoch: int) -> Dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": manifest["campaign_id"],
        "chunk_id": chunk_id,
        "attempt_id": f"attempt-{uuid.uuid4().hex}",
        "lease_token": uuid.uuid4().hex,
        "lease_epoch": epoch,
        "worker_id": worker_id,
        "pid": os.getpid(),
        "claimed_at": utc_now(),
    }


def _acquire_chunk_lock(
    campaign: Path, manifest: Mapping[str, Any], chunk_id: str, worker_id: str
) -> Optional[Dict[str, Any]]:
    lock_path = _chunk_dir(campaign, chunk_id) / "lease.lock"
    lock = {"token": uuid.uuid4().hex, "worker_id": worker_id, "timestamp": utc_now()}
    try:
        exclusive_write_json(lock_path, lock)
        return lock
    except FileExistsError:
        try:
            existing = read_json(lock_path)
            age = time.time() - _parse_time(str(existing.get("timestamp", "")))
            if age <= int(manifest["lease_seconds"]):
                return None
            lock_path.unlink()
            exclusive_write_json(lock_path, lock)
            return lock
        except (FleetError, FileExistsError, OSError, TypeError, ValueError):
            return None


def _release_chunk_lock(campaign: Path, chunk_id: str, lock: Mapping[str, Any]) -> None:
    lock_path = _chunk_dir(campaign, chunk_id) / "lease.lock"
    try:
        if read_json(lock_path).get("token") == lock.get("token"):
            lock_path.unlink()
    except (FleetError, OSError):
        pass


def _reclaim_stale(campaign: Path, manifest: Mapping[str, Any], chunk_id: str, worker_id: str) -> Optional[Dict[str, Any]]:
    chunk_dir = _chunk_dir(campaign, chunk_id)
    claim_path = chunk_dir / "claim.json"
    lock = _acquire_chunk_lock(campaign, manifest, chunk_id, worker_id)
    if lock is None:
        return None

    try:
        try:
            old_claim = read_json(claim_path)
        except FleetError:
            return None
        if not claim_is_stale(campaign, old_claim, int(manifest["lease_seconds"])):
            return None
        claim = _new_claim(manifest, chunk_id, worker_id, int(old_claim.get("lease_epoch", 0)) + 1)
        atomic_write_json(claim_path, claim)
        write_heartbeat(campaign, claim)
        return claim
    finally:
        _release_chunk_lock(campaign, chunk_id, lock)


def claim_next(campaign: Path, manifest: Mapping[str, Any], worker_id: str) -> Optional[Tuple[Dict[str, Any], Dict[str, Any]]]:
    validate_id(worker_id, "worker ID")
    for chunk in manifest["chunks"]:
        chunk_id = str(chunk["id"])
        chunk_dir = _chunk_dir(campaign, chunk_id)
        if (chunk_dir / "result.json").is_file():
            continue
        if len(_attempt_files(campaign, chunk_id)) >= int(manifest["max_attempts"]):
            continue
        claim = _new_claim(manifest, chunk_id, worker_id, 1)
        try:
            exclusive_write_json(chunk_dir / "claim.json", claim)
            write_heartbeat(campaign, claim)
            return dict(chunk), claim
        except FileExistsError:
            try:
                old_claim = read_json(chunk_dir / "claim.json")
            except FleetError:
                continue
            if claim_is_stale(campaign, old_claim, int(manifest["lease_seconds"])):
                reclaimed = _reclaim_stale(campaign, manifest, chunk_id, worker_id)
                if reclaimed is not None:
                    return dict(chunk), reclaimed
    return None


def release_claim(
    campaign: Path, manifest: Mapping[str, Any], claim: Mapping[str, Any]
) -> None:
    lock = _acquire_chunk_lock(
        campaign, manifest, str(claim["chunk_id"]), str(claim["worker_id"])
    )
    if lock is None:
        return
    claim_path = _chunk_dir(campaign, str(claim["chunk_id"])) / "claim.json"
    try:
        if lease_owned(campaign, claim):
            claim_path.unlink()
    finally:
        _release_chunk_lock(campaign, str(claim["chunk_id"]), lock)


def verify_binary_pins(manifest: Mapping[str, Any]) -> None:
    for label in ("runner", "engine"):
        path = Path(str(manifest[label]["path"]))
        if not path.is_file():
            raise FleetError(f"pinned {label} binary is missing: {path}")
        actual = sha256_file(path)
        if actual != manifest[label]["sha256"]:
            raise FleetError(f"pinned {label} binary changed: expected {manifest[label]['sha256']}, got {actual}")
    posebusters = manifest["validators"]["posebusters"]
    pb_path = Path(str(posebusters["path"]))
    if not pb_path.is_file() or not os.access(pb_path, os.X_OK):
        raise FleetError(f"pinned PoseBusters executable is missing: {pb_path}")
    actual_pb = sha256_file(pb_path)
    if actual_pb != posebusters["sha256"]:
        raise FleetError(
            f"pinned PoseBusters executable changed: expected {posebusters['sha256']}, got {actual_pb}"
        )


def summary_matches_chunk(
    summary: Mapping[str, Any], manifest: Mapping[str, Any], manifest_hash: str,
    chunk: Mapping[str, Any], claim: Mapping[str, Any],
) -> bool:
    targets = summary.get("targets")
    if not isinstance(targets, list) or not all(isinstance(target, dict) for target in targets):
        return False
    expected_codes = [str(code).upper() for code in chunk["codes"]]
    actual_codes = [str(target.get("pdb_id", "")).upper() for target in targets]
    target_summary = summary.get("summary")
    target_count_matches = (
        isinstance(target_summary, dict)
        and target_summary.get("target_count") == len(expected_codes)
    )
    return all((
        summary.get("schema_version") == SCHEMA_VERSION,
        summary.get("type") == "flexaidds-fleet-chunk-result",
        summary.get("campaign_id") == manifest["campaign_id"],
        summary.get("chunk_id") == chunk["id"],
        summary.get("attempt_id") == claim["attempt_id"],
        summary.get("worker_id") == claim["worker_id"],
        summary.get("dataset") == manifest["benchmark"],
        target_count_matches,
        len(actual_codes) == len(expected_codes),
        len(set(actual_codes)) == len(actual_codes),
        set(actual_codes) == set(expected_codes),
        summary.get("provenance", {}).get("manifest_sha256") == manifest_hash,
        summary.get("provenance", {}).get("runner_sha256") == manifest["runner"]["sha256"],
        summary.get("provenance", {}).get("engine_sha256") == manifest["engine"]["sha256"],
    ))


def build_runner_command(
    manifest: Mapping[str, Any], manifest_hash: str, chunk: Mapping[str, Any],
    claim: Mapping[str, Any], local_root: Path,
) -> List[str]:
    protocol = manifest["protocol"]
    command = [
        str(manifest["runner"]["path"]),
        "--fleet",
        "--campaign-id", str(manifest["campaign_id"]),
        "--chunk-id", str(chunk["id"]),
        "--attempt-id", str(claim["attempt_id"]),
        "--worker-id", str(claim["worker_id"]),
        "--manifest-sha256", manifest_hash,
        "--runner-sha256", str(manifest["runner"]["sha256"]),
        "--engine-sha256", str(manifest["engine"]["sha256"]),
        "--output-json", str(local_root / "fleet_chunk_result.json"),
        "--benchmark", str(manifest["benchmark"]),
        "--output", str(local_root / "outputs"),
        "--only-codes", ",".join(chunk["codes"]),
        "--mode", str(protocol["mode"]),
        "--threads", str(protocol["threads"]),
        "--omp-threads", str(protocol["omp_threads"]),
        "--ga-population", str(protocol["ga_population"]),
        "--ga-generations", str(protocol["ga_generations"]),
        "--temperature", str(protocol["temperature_k"]),
        "--grid-spacing", str(protocol["grid_spacing_a"]),
        "--job-timeout-seconds", str(protocol["job_timeout_s"]),
        "--clustering", str(protocol["clustering"]),
    ]
    if protocol.get("cache"):
        command.extend(("--cache", str(protocol["cache"])))
    if protocol.get("gpu_backend"):
        command.extend(("--gpu", str(protocol["gpu_backend"])))
    return command


def copy_tree_verified(source: Path, destination: Path) -> str:
    if destination.exists():
        raise FleetError(f"immutable archive already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{uuid.uuid4().hex}")
    source_hash = tree_sha256(source)
    shutil.copytree(source, temporary, copy_function=shutil.copy2)
    copied_hash = tree_sha256(temporary)
    if copied_hash != source_hash:
        shutil.rmtree(temporary, ignore_errors=True)
        raise FleetError(f"archive verification failed: source {source_hash}, copy {copied_hash}")
    os.replace(temporary, destination)
    _fsync_directory(destination.parent)
    final_hash = tree_sha256(destination)
    if final_hash != source_hash:
        raise FleetError(f"archive changed after publication: expected {source_hash}, got {final_hash}")
    return source_hash


def _attempt_record_path(campaign: Path, claim: Mapping[str, Any]) -> Path:
    return confined_path(campaign, "attempts", str(claim["chunk_id"]), f"{claim['attempt_id']}.json")


def _quarantine(campaign: Path, claim: Mapping[str, Any], record: Mapping[str, Any]) -> None:
    path = confined_path(campaign, "quarantine", str(claim["chunk_id"]), f"{claim['attempt_id']}.json")
    exclusive_write_json(path, dict(record))


def run_claim(
    campaign: Path, manifest: Mapping[str, Any], manifest_hash: str,
    chunk: Mapping[str, Any], claim: Mapping[str, Any], keep_local: bool,
) -> bool:
    verify_binary_pins(manifest)
    compute_root = Path(str(manifest["compute_root"]))
    compute_root.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(compute_root).free
    if free_bytes < int(manifest["min_free_bytes"]):
        raise FleetError(f"local free-space guard failed: {free_bytes} bytes available")

    local_root = compute_root / str(manifest["campaign_id"]) / str(chunk["id"]) / str(claim["attempt_id"])
    local_root.mkdir(parents=True, exist_ok=False)
    (local_root / "outputs").mkdir()
    command = build_runner_command(manifest, manifest_hash, chunk, claim, local_root)
    environment = os.environ.copy()
    environment.update({str(key): str(value) for key, value in manifest["environment"].items()})
    environment["FLEXAIDDS_BINARY"] = str(manifest["engine"]["path"])

    started = utc_now()
    return_code: Optional[int] = None
    lease_lost = False
    launch_error: Optional[str] = None
    process = None
    stdout_path = local_root / "runner.stdout.log"
    stderr_path = local_root / "runner.stderr.log"
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            try:
                process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=environment)
            except OSError as exc:
                launch_error = str(exc)
            else:
                heartbeat_interval = max(1.0, min(10.0, int(manifest["lease_seconds"]) / 3.0))
                while True:
                    if not write_heartbeat(campaign, claim):
                        lease_lost = True
                        process.terminate()
                        try:
                            process.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            process.kill()
                        break
                    try:
                        return_code = process.wait(timeout=heartbeat_interval)
                        break
                    except subprocess.TimeoutExpired:
                        continue
                if return_code is None:
                    return_code = process.wait()
    except KeyboardInterrupt:
        if process is not None and process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
        raise

    summary_path = local_root / "fleet_chunk_result.json"
    summary: Optional[Dict[str, Any]] = None
    if summary_path.is_file():
        try:
            summary = read_json(summary_path)
        except FleetError:
            summary = None

    record: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": manifest["campaign_id"],
        "chunk_id": chunk["id"],
        "attempt_id": claim["attempt_id"],
        "worker_id": claim["worker_id"],
        "lease_epoch": claim["lease_epoch"],
        "manifest_sha256": manifest_hash,
        "started_at": started,
        "finished_at": utc_now(),
        "return_code": return_code,
        "command": command,
        "local_root": str(local_root),
        "launch_error": launch_error,
    }

    if lease_lost or not lease_owned(campaign, claim):
        record["status"] = "quarantined_lease_lost"
        _quarantine(campaign, claim, record)
        return False

    valid_summary = summary is not None and summary_matches_chunk(
        summary, manifest, manifest_hash, chunk, claim
    )
    if return_code != 0 or not valid_summary:
        record["status"] = "failed"
        record["error"] = launch_error or "runner failed or emitted an invalid Fleet summary"
        record["stderr_path"] = str(stderr_path)
        exclusive_write_json(_attempt_record_path(campaign, claim), record)
        release_claim(campaign, manifest, claim)
        return False

    write_heartbeat(campaign, claim)
    archive = confined_path(campaign, "artifacts", str(chunk["id"]), str(claim["attempt_id"]))
    archive_hash = copy_tree_verified(local_root, archive)
    write_heartbeat(campaign, claim)
    record.update({
        "status": "completed",
        "archive": str(archive.relative_to(campaign)),
        "archive_tree_sha256": archive_hash,
        "summary": summary,
    })
    exclusive_write_json(_attempt_record_path(campaign, claim), record)

    lock = _acquire_chunk_lock(
        campaign, manifest, str(chunk["id"]), str(claim["worker_id"])
    )
    if lock is None:
        record["status"] = "quarantined_publish_lock"
        _quarantine(campaign, claim, record)
        return False
    try:
        if not lease_owned(campaign, claim):
            record["status"] = "quarantined_after_archive"
            _quarantine(campaign, claim, record)
            return False
        accepted = {
            "schema_version": SCHEMA_VERSION,
            "status": "completed",
            "campaign_id": manifest["campaign_id"],
            "chunk_id": chunk["id"],
            "attempt_id": claim["attempt_id"],
            "worker_id": claim["worker_id"],
            "lease_epoch": claim["lease_epoch"],
            "manifest_sha256": manifest_hash,
            "accepted_at": utc_now(),
            "archive": str(archive.relative_to(campaign)),
            "archive_tree_sha256": archive_hash,
            "summary": summary,
        }
        exclusive_write_json(_chunk_dir(campaign, str(chunk["id"])) / "result.json", accepted)
        claim_path = _chunk_dir(campaign, str(chunk["id"])) / "claim.json"
        if lease_owned(campaign, claim):
            claim_path.unlink()
    finally:
        _release_chunk_lock(campaign, str(chunk["id"]), lock)
    if not keep_local:
        shutil.rmtree(local_root, ignore_errors=True)
    return True


def campaign_snapshot(campaign: Path, manifest: Mapping[str, Any], write: bool = True) -> Dict[str, Any]:
    states = {"completed": 0, "running": 0, "stale": 0, "failed": 0, "pending": 0}
    active_chunks: List[Dict[str, Any]] = []
    devices: Dict[str, Dict[str, Any]] = {}
    for chunk in manifest["chunks"]:
        chunk_id = str(chunk["id"])
        chunk_dir = _chunk_dir(campaign, chunk_id)
        if (chunk_dir / "result.json").is_file():
            states["completed"] += 1
            continue
        attempts = len(_attempt_files(campaign, chunk_id))
        claim_path = chunk_dir / "claim.json"
        if claim_path.is_file():
            try:
                claim = read_json(claim_path)
                stale = claim_is_stale(campaign, claim, int(manifest["lease_seconds"]))
                state = "stale" if stale else "running"
                states[state] += 1
                active_chunks.append({
                    "id": chunk_id,
                    "status": state,
                    "claimedBy": claim.get("worker_id"),
                    "attemptID": claim.get("attempt_id"),
                    "leaseEpoch": claim.get("lease_epoch"),
                })
                if not stale:
                    worker = str(claim.get("worker_id", "unknown"))
                    devices[worker] = {"deviceID": worker, "model": worker, "thermalState": "unknown"}
                continue
            except FleetError:
                states["stale"] += 1
                continue
        if attempts >= int(manifest["max_attempts"]):
            states["failed"] += 1
        else:
            states["pending"] += 1

    old_revision = 0
    status_path = campaign / "status.json"
    if status_path.is_file():
        try:
            old_revision = int(read_json(status_path).get("revision", 0))
        except FleetError:
            pass
    total = len(manifest["chunks"])
    snapshot: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": manifest["campaign_id"],
        "dataset": manifest["dataset"],
        "timestamp": utc_now(),
        "revision": old_revision + 1,
        "states": states,
        "devices": list(devices.values()),
        "activeChunks": active_chunks,
        "completedChunks": states["completed"],
        "totalChunks": total,
        "failedChunks": states["failed"],
        "orphanedChunks": states["stale"],
        "metrics": {
            "jobID": manifest["campaign_id"],
            "totalChunks": total,
            "completedChunks": states["completed"],
            "failedChunks": states["failed"],
            "orphanedChunks": states["stale"],
            "activeDevices": len(devices),
            "totalTFLOPS": 0.0,
            "meanChunkTimeSeconds": None,
            "estimatedRemainingSeconds": None,
            "timestamp": utc_now(),
        },
    }
    if write:
        atomic_write_json(status_path, snapshot)
    return snapshot


def run_worker(args: argparse.Namespace) -> int:
    campaign = Path(args.campaign).expanduser().resolve()
    manifest, manifest_hash = load_manifest(campaign)
    worker_id = validate_id(args.worker_id or f"{socket.gethostname()}-{os.getpid()}", "worker ID")
    attempted = 0
    while True:
        claimed = claim_next(campaign, manifest, worker_id)
        if claimed is None:
            campaign_snapshot(campaign, manifest)
            return 0
        chunk, claim = claimed
        attempted += 1
        try:
            run_claim(campaign, manifest, manifest_hash, chunk, claim, args.keep_local)
        except BaseException:
            if lease_owned(campaign, claim):
                release_claim(campaign, manifest, claim)
            raise
        finally:
            campaign_snapshot(campaign, manifest)
        if args.once or attempted >= args.max_chunks:
            return 0


AGGREGATE_FIELDS = [
    "campaign_id", "chunk_id", "attempt_id", "pdb_id", "execution_completed",
    "rmsd_hungarian_a", "success_rmsd", "posebusters_ran", "posebusters_pass",
    "success_pb", "tencom_status", "eigen_status", "validators_complete",
    "protocol_claim_eligible", "claim_ready", "seed_echo", "pose_source",
    "pose_sha256", "rmsd_pose_sha256", "posebusters_pose_sha256",
    "tencom_pose_sha256", "elected_pose_path", "elected_pose_source",
    "elected_restart", "elected_cluster", "elected_cf", "best_cluster_rmsd_a",
]


def aggregate_campaign(campaign: Path, manifest: Mapping[str, Any], manifest_hash: str) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    completed_chunks = 0
    for chunk in manifest["chunks"]:
        result_path = _chunk_dir(campaign, str(chunk["id"])) / "result.json"
        if not result_path.is_file():
            continue
        result = read_json(result_path)
        if result.get("manifest_sha256") != manifest_hash:
            raise FleetError(f"chunk result has wrong manifest hash: {result_path}")
        completed_chunks += 1
        summary = result.get("summary", {})
        for target in summary.get("targets", []):
            row = {field: target.get(field) for field in AGGREGATE_FIELDS}
            row.update({
                "campaign_id": manifest["campaign_id"],
                "chunk_id": chunk["id"],
                "attempt_id": result["attempt_id"],
            })
            rows.append(row)

    aggregate_dir = campaign / "aggregate"
    aggregate_dir.mkdir(exist_ok=True)
    csv_path = aggregate_dir / "targets.csv"
    temporary = csv_path.with_name(f".{csv_path.name}.tmp-{uuid.uuid4().hex}")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=AGGREGATE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, csv_path)

    def count(field: str) -> int:
        return sum(value.get(field) is True for value in rows)

    total = len(rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": manifest["campaign_id"],
        "manifest_sha256": manifest_hash,
        "created_at": utc_now(),
        "complete": completed_chunks == len(manifest["chunks"]),
        "chunks_completed": completed_chunks,
        "chunks_total": len(manifest["chunks"]),
        "target_count": total,
        "execution_completed": count("execution_completed"),
        "success_rmsd": count("success_rmsd"),
        "success_pb": count("success_pb"),
        "validators_complete": count("validators_complete"),
        "claim_ready": count("claim_ready"),
        "success_pb_rate": count("success_pb") / total if total else 0.0,
        "claim_ready_rate": count("claim_ready") / total if total else 0.0,
        "targets_csv": str(csv_path),
    }
    atomic_write_json(aggregate_dir / "summary.json", summary)
    return summary


def _print_json(data: Any) -> None:
    print(json.dumps(data, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="flexaidds-fleet", description="Bonhomme Fleet benchmark orchestrator")
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan = subparsers.add_parser("plan", help="create an immutable campaign manifest")
    plan.add_argument("campaign")
    plan.add_argument("--campaign-id")
    plan.add_argument("--runner", required=True)
    plan.add_argument("--engine")
    plan.add_argument("--posebusters-bin")
    plan.add_argument("--benchmark", required=True)
    plan.add_argument("--dataset")
    plan.add_argument("--mode", required=True, choices=("defined-cleft-redock", "autonomous"))
    plan.add_argument("--codes-file", required=True)
    plan.add_argument("--chunks", type=int, default=1)
    plan.add_argument("--compute-root", default=os.environ.get("FLEXAIDDS_FLEET_COMPUTE_ROOT", "/private/tmp/flexaidds_fleet"))
    plan.add_argument("--cache")
    plan.add_argument("--threads", type=int, default=1)
    plan.add_argument("--omp-threads", type=int, default=4)
    plan.add_argument("--ga-population", type=int, default=1000)
    plan.add_argument("--ga-generations", type=int, default=6000)
    plan.add_argument("--temperature", type=float, default=298.0)
    plan.add_argument("--grid-spacing", type=float, default=0.375)
    plan.add_argument("--job-timeout-seconds", type=int, default=10800)
    plan.add_argument("--clustering", default="CF")
    plan.add_argument("--gpu", choices=("metal", "cuda"))
    plan.add_argument("--lease-seconds", type=int, default=300)
    plan.add_argument("--max-attempts", type=int, default=3)
    plan.add_argument("--min-free-gb", type=float, default=5.0)
    plan.add_argument("--env", action="append", default=[])

    for name in ("worker", "run", "resume"):
        worker = subparsers.add_parser(name, help="claim and execute campaign chunks")
        worker.add_argument("campaign")
        worker.add_argument("--worker-id")
        worker.add_argument("--once", action="store_true")
        worker.add_argument("--max-chunks", type=int, default=sys.maxsize)
        worker.add_argument("--keep-local", action="store_true")

    status = subparsers.add_parser("status", help="show current campaign state")
    status.add_argument("campaign")

    aggregate = subparsers.add_parser("aggregate", help="build strict target-level CSV and summary")
    aggregate.add_argument("campaign")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "plan":
            _print_json(plan_campaign(args))
        elif args.command in ("worker", "run", "resume"):
            return run_worker(args)
        elif args.command == "status":
            campaign = Path(args.campaign).expanduser().resolve()
            manifest, _ = load_manifest(campaign)
            _print_json(campaign_snapshot(campaign, manifest))
        elif args.command == "aggregate":
            campaign = Path(args.campaign).expanduser().resolve()
            manifest, manifest_hash = load_manifest(campaign)
            _print_json(aggregate_campaign(campaign, manifest, manifest_hash))
    except (FleetError, FileExistsError, OSError, ValueError) as exc:
        print(f"fleet error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
