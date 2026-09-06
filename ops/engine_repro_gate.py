#!/usr/bin/env python3
"""Isolated engine receipts and exact scientific-output comparisons.

run-one creates a NEW output directory, pins inputs/environment/build metadata,
runs one explicitly selected binary, and admits only complete outputs. compare
rechecks artifact hashes and shapes, then compares ranks and grid ordering. Only
FLEXAID.commit and FLEXAID.dirty values on their dedicated provenance REMARK may
be normalized; every other emitted byte, including the seed, is retained.

All gate comparisons require generation-zero receipts. Explicit observer-off
transparency runs are non-gate diagnostics and cannot pass a gate comparison. METHODOLOGY.md sections 1-2 define gate interpretation.
Canonical FLEXAID_ and FLEXAIDDS_ scientific controls are pinned internally;
they cannot be supplied or overridden through --env.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import signal
import struct
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

RUN_SCHEMA = "flexaidds.engine_repro_run.v1"
GEN0_SCHEMA = "flexaidds.ga_population_receipt.v1"
SEED = "12345"
EXPECTED_POPULATION = 1000
RANK_COUNT = 10
PROVENANCE = re.compile(rb"^(REMARK FLEXAID\.commit=)([^ \r\n]+)( FLEXAID\.dirty=)([012])( FLEXAID\.seed=([0-9]+)[^\r\n]*)(\r?\n|$)", re.MULTILINE)
CF_LINE = re.compile(rb"^REMARK CF=\s*([^\s\r\n]+)\s*$", re.MULTILINE)
BITS64 = re.compile(r"0x[0-9a-fA-F]{16}\Z")
BITS32 = re.compile(r"0x[0-9a-fA-F]{8}\Z")
CF_FIELDS = {"com_bits", "con_bits", "wal_bits", "sas_bits", "elec_bits", "gist_bits",
             "hbond_bits", "gist_desolv_bits", "metal_coord_bits", "h_rep_bits", "entropy_bits",
             "pb_clash_bits", "totsas_bits", "nor_bits", "rclash"}
RECORD_FIELDS = {"index", "status", "genes", "cf", "evalue_bits", "app_evalue_bits", "fitnes_bits",
                 "boltzmann_weight_bits", "free_energy_bits", "ring_phases_bits", "ring_six", "ring_five"}
PROTECTED_ENV = {"FLEXAID_SEED", "OMP_NUM_THREADS", "FLEXAIDDS_NO_SEC", "FLEXAIDDS_RESTARTS",
                 "FLEXAIDDS_DATA_DIR", "FLEXAIDDS_PARALLEL_REPRODUCE", "FLEXAIDDS_GEN0_RECEIPT", "FLEXAIDDS_CACHE_READONLY",
                 "FLEXAIDDS_VORONOI_KEYED_JITTER", "FLEXAIDDS_VCF_DIAG"}
SCIENCE_ENV_PREFIXES = ("FLEXAID_", "FLEXAIDDS_")
CANONICAL_SCIENCE_ENV = {"FLEXAID_SEED", "FLEXAIDDS_NO_SEC", "FLEXAIDDS_RESTARTS",
                         "FLEXAIDDS_CACHE_READONLY", "FLEXAIDDS_VCF_DIAG", "FLEXAIDDS_DATA_DIR",
                         "FLEXAIDDS_PARALLEL_REPRODUCE", "FLEXAIDDS_GEN0_RECEIPT"}
SCORING_FAILURE_MARKERS = (
    b"[VCF-DIAG]",
    b"Atom OptRes pointer is not owned by the master OptRes array",
    b"Atom/OptRes workspace sizes changed after binding capture",
)


class GateError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pin_file(path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    if not path.is_file() or path.stat().st_size == 0:
        raise GateError(f"missing or empty file: {path}")
    if "Mobile Documents" in path.parts or "com~apple~CloudDocs" in path.parts:
        raise GateError("stage all inputs locally; CloudDocs is not a compute filesystem")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def json_read(path: Path) -> dict:
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise GateError(f"duplicate JSON key {key!r} in {path}")
            result[key] = value
        return result
    def finite_float(text):
        value = float(text)
        if not math.isfinite(value):
            raise GateError(f"nonfinite JSON number in {path}")
        return value
    value = json.loads(path.read_text(), object_pairs_hook=unique_pairs, parse_float=finite_float,
                       parse_constant=lambda value: (_ for _ in ()).throw(GateError(f"nonfinite JSON token: {value}")))
    if not isinstance(value, dict):
        raise GateError(f"JSON object required: {path}")
    return value


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def checked_git(source: Path, *args: str) -> bytes:
    result = subprocess.run(["git", "-C", str(source), *args], capture_output=True, timeout=30)
    if result.returncode:
        raise GateError(f"cannot capture source metadata: git {' '.join(args)}: {result.stderr.decode(errors='replace')}")
    return result.stdout


def source_build_metadata(source: Path, build: Path) -> dict:
    source, build = source.resolve(strict=True), build.resolve(strict=True)
    cache = build / "CMakeCache.txt"
    cache_pin = pin_file(cache)
    entries = {}
    for line in cache.read_text().splitlines():
        if line and not line.startswith(("#", "//")) and "=" in line and ":" in line.split("=", 1)[0]:
            key, value = line.split("=", 1)
            entries[key.split(":", 1)[0]] = value
    if not entries.get("CMAKE_HOME_DIRECTORY") or Path(entries["CMAKE_HOME_DIRECTORY"]).resolve() != source:
        raise GateError("CMake cache source directory does not match --source-dir")
    generated = [build / name for name in ("build.ninja", "Makefile") if (build / name).is_file()]
    if not generated:
        raise GateError("build metadata requires a generated build.ninja or Makefile")
    diff = checked_git(source, "diff", "--binary", "HEAD")
    untracked = []
    for raw in checked_git(source, "ls-files", "--others", "--exclude-standard", "-z").split(b"\0"):
        if raw:
            path = source / os.fsdecode(raw)
            if path.is_file():
                untracked.append({"relative_path": os.fsdecode(raw), "sha256": sha256(path), "bytes": path.stat().st_size})
    files = [cache_pin, *(pin_file(path) for path in generated)]
    if entries.get("FLEXAIDS_USE_OPENMP", "").upper() not in ("ON", "1", "TRUE", "YES"):
        raise GateError("merge validation requires FLEXAIDS_USE_OPENMP=ON")
    compile_flags = []
    commands_path = build / "compile_commands.json"
    if commands_path.is_file():
        commands = json.loads(commands_path.read_text())
        for entry in commands:
            command = entry.get("command", "") or " ".join(entry.get("arguments", []))
            if Path(entry.get("file", "")).name == "gaboom.cpp" and "CMakeFiles/FlexAIDdS.dir/" in command:
                compile_flags.append(command)
    if not compile_flags:
        flags_path = build / "CMakeFiles/FlexAIDdS.dir/flags.make"
        if flags_path.is_file():
            compile_flags.append(flags_path.read_text())
            files.append(pin_file(flags_path))
    if not compile_flags or any("-fopenmp" not in text or re.search(r"[-/]D\s*FLEXAID_DETERMINISTIC(?:\b|=)", text) for text in compile_flags):
        raise GateError("actual FlexAIDdS/gaboom compile flags must enable OpenMP and omit FLEXAID_DETERMINISTIC")
    for name in ("compile_commands.json", "CMakeFiles/TargetDirectories.txt"):
        if (build / name).is_file():
            files.append(pin_file(build / name))
    return {"source_dir": str(source), "source_commit": checked_git(source, "rev-parse", "HEAD").decode().strip(),
            "source_status": checked_git(source, "status", "--porcelain=v1", "--untracked-files=all").decode(),
            "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(), "untracked_files": untracked,
            "build_dir": str(build), "build_files": files, "openmp_compile_flags": compile_flags,
            "cmake": {key: value for key, value in entries.items() if key.startswith(("CMAKE_", "BUILD_", "USE_", "ENABLE_", "FLEXAID", "OpenMP_"))},
            "evidence_limit": "Source/build metadata captured; this receipt does not independently prove compiler execution."}


def validate_config(path: Path) -> dict:
    config = json_read(path)
    ga = config.get("ga", {})
    reference = config.get("reference_ligand", {})
    if ga.get("num_generations") != 2000 or ga.get("num_chromosomes") != 1000:
        raise GateError("merge gate requires explicit ga.num_generations=2000 and ga.num_chromosomes=1000")
    if reference.get("pose_seed_enabled") is not False or reference.get("seed_fraction") != 0:
        raise GateError("merge gate requires explicit disabled native-pose seeding and zero seed_fraction")
    if config.get("grid_file"):
        raise GateError("grid_file reuse is forbidden: this gate must generate and compare new grids")
    return config


def normalize_pose(path: Path) -> tuple[bytes, dict]:
    data = path.read_bytes()
    matches = list(PROVENANCE.finditer(data))
    if len(matches) != 1 or matches[0].group(6).decode() != SEED:
        raise GateError(f"{path.name}: exactly one recognized commit/dirty/seed={SEED} REMARK required")
    cf = list(CF_LINE.finditer(data))
    if len(cf) != 1:
        raise GateError(f"{path.name}: exactly one elected CF REMARK required")
    try:
        if not math.isfinite(float(cf[0].group(1))):
            raise ValueError()
    except ValueError:
        raise GateError(f"{path.name}: nonfinite or invalid elected CF")
    if not data.splitlines() or data.splitlines()[-1] != b"END":
        raise GateError(f"{path.name}: missing final END record")
    atoms = [line for line in data.splitlines() if line.startswith((b"ATOM  ", b"HETATM"))]
    if not atoms:
        raise GateError(f"{path.name}: no coordinate records")
    for line in atoms:
        try:
            if not all(math.isfinite(float(line[start:start + 8])) for start in (30, 38, 46)):
                raise ValueError()
        except ValueError:
            raise GateError(f"{path.name}: malformed/nonfinite coordinate")
    normalized = PROVENANCE.sub(lambda m: m.group(1) + b"<commit>" + m.group(3) + b"<dirty>" + m.group(5) + m.group(7), data)
    return normalized, {"cf_text": cf[0].group(1).decode(), "coordinate_record_count": len(atoms),
                        "commit_remark": matches[0].group(2).decode(), "dirty_remark": matches[0].group(4).decode(),
                        "seed_remark": matches[0].group(6).decode()}


def validate_grid(path: Path) -> dict:
    data = path.read_bytes()
    if len(data) < 12:
        raise GateError("missing/truncated RRG header")
    for endian in ("<", ">"):
        magic, version, count = struct.unpack(endian + "IIi", data[:12])
        if magic == 0x56435400:
            break
    else:
        raise GateError("unrecognized RRG magic")
    if version != 1 or count <= 0 or len(data) != 12 + count * 32:
        raise GateError("invalid RRG version/count/record width or trailing bytes")
    for offset in range(12, len(data), 32):
        values = struct.unpack(endian + "ii6f", data[offset:offset + 32])
        if not all(math.isfinite(value) for value in values[2:]):
            raise GateError("RRG contains nonfinite coordinate/internal-coordinate values")
    return {"version": version, "point_count": count, "record_bytes": 32,
            "byte_order": "little" if endian == "<" else "big",
            "ordered_records_sha256": hashlib.sha256(data[12:]).hexdigest()}


def _integer(value: Any, low: int, high: int, field: str) -> None:
    if type(value) is not int or not low <= value <= high:
        raise GateError(f"invalid integer {field}")


def _bits(value: Any, width: int, field: str) -> None:
    pattern = BITS64 if width == 64 else BITS32
    if not isinstance(value, str) or not pattern.fullmatch(value):
        raise GateError(f"invalid IEEE bit string for {field}")
    # Nonfinite values cannot establish a usable initial scoring population.
    number = struct.unpack(">d" if width == 64 else ">f", bytes.fromhex(value[2:]))[0]
    if not math.isfinite(number):
        raise GateError(f"nonfinite generation-zero {field}")


def validate_execution(value: Any, expected_threads: int) -> dict:
    """Require witnessed, complete RANDOM initialization on the requested team."""
    if (not isinstance(value, dict)
            or set(value) != {"openmp_compiled", "deterministic_compile", "evaluation_batches"}
            or value["openmp_compiled"] is not True or value["deterministic_compile"] is not False):
        raise GateError("generation-zero execution requires OpenMP without deterministic compilation")
    if type(expected_threads) is not int or expected_threads not in (1, 4):
        raise GateError("invalid expected OpenMP team size")
    batches = value["evaluation_batches"]
    if not isinstance(batches, list) or len(batches) != 1:
        raise GateError("canonical RANDOM initialization requires one witnessed evaluation batch")
    batch = batches[0]
    if (not isinstance(batch, dict)
            or set(batch) != {"region", "population_count", "popoffset", "workspace_slots", "workers"}
            or batch["region"] != "populate_chromosomes"):
        raise GateError("invalid initial-population evaluation batch")
    _integer(batch["population_count"], EXPECTED_POPULATION, EXPECTED_POPULATION, "batch population_count")
    _integer(batch["popoffset"], 0, 0, "batch popoffset")
    _integer(batch["workspace_slots"], expected_threads, 100000, "batch workspace_slots")
    workers = batch["workers"]
    if not isinstance(workers, list) or len(workers) != expected_threads:
        raise GateError("every requested OpenMP worker must witness completed evaluations")
    ids, evaluated = set(), 0
    for worker in workers:
        if not isinstance(worker, dict) or set(worker) != {"worker_id", "team_size", "evaluated_chromosomes"}:
            raise GateError("invalid observed worker record")
        _integer(worker["worker_id"], 0, expected_threads - 1, "worker_id")
        _integer(worker["team_size"], expected_threads, expected_threads, "observed team_size")
        _integer(worker["evaluated_chromosomes"], 1, EXPECTED_POPULATION, "evaluated_chromosomes")
        if worker["worker_id"] in ids:
            raise GateError("duplicate observed worker id")
        ids.add(worker["worker_id"])
        evaluated += worker["evaluated_chromosomes"]
    if evaluated != EXPECTED_POPULATION:
        raise GateError("witnessed initial evaluations must total configured 1000 chromosomes")
    return {**value, "observed_threads": expected_threads, "evaluated_chromosomes": evaluated}


def validate_gen0(path: Path, expected_threads: int = 1) -> dict:
    value = json_read(path)
    if set(value) != {"schema", "boundary", "generation", "population_count", "n_genes", "seed", "execution", "records", "complete"}:
        raise GateError("unexpected or missing generation-zero top-level fields")
    if (value.get("schema") != GEN0_SCHEMA or value.get("boundary") != "initial_population_complete_before_reproduction"
            or type(value.get("generation")) is not int or value.get("generation") != 0
            or value.get("complete") is not True or value.get("seed") != SEED):
        raise GateError("invalid generation-zero schema, completion, boundary, generation or seed")
    count, n_genes = value.get("population_count"), value.get("n_genes")
    if count != EXPECTED_POPULATION or type(count) is not int:
        raise GateError("generation-zero population must match configured 1000 chromosomes")
    execution = validate_execution(value["execution"], expected_threads)
    _integer(n_genes, 1, 100000, "n_genes")
    records = value.get("records")
    if not isinstance(records, list) or len(records) != count:
        raise GateError("generation-zero population count mismatch")
    indices, canonical = set(), []
    for record in records:
        if not isinstance(record, dict) or set(record) != RECORD_FIELDS:
            raise GateError("incomplete generation-zero record")
        index = record["index"]
        _integer(index, 0, count - 1, "record index")
        if index in indices:
            raise GateError("duplicate generation-zero record index")
        indices.add(index)
        if type(record["status"]) is not int or record["status"] != ord("n"):
            raise GateError("generation-zero record is not evaluated (status must be n)")
        if not isinstance(record["genes"], list) or len(record["genes"]) != n_genes:
            raise GateError("generation-zero gene count mismatch")
        for gene in record["genes"]:
            if not isinstance(gene, dict) or set(gene) != {"to_int32", "to_ic_bits"}:
                raise GateError("incomplete generation-zero gene")
            _integer(gene["to_int32"], -(2**31), 2**31 - 1, "gene.to_int32")
            _bits(gene["to_ic_bits"], 64, "gene.to_ic_bits")
        cf = record["cf"]
        if not isinstance(cf, dict) or set(cf) != CF_FIELDS:
            raise GateError("incomplete generation-zero CF components")
        for key in CF_FIELDS - {"rclash"}:
            _bits(cf[key], 64, key)
        _integer(cf["rclash"], -(2**31), 2**31 - 1, "rclash")
        for key in ("evalue_bits", "app_evalue_bits", "fitnes_bits", "boltzmann_weight_bits", "free_energy_bits"):
            _bits(record[key], 64, key)
        for key in ("ring_phases_bits", "ring_six", "ring_five"):
            if not isinstance(record[key], list) or len(record[key]) != 16:
                raise GateError(f"invalid {key}: schema v1 requires MAX_RING_FLEX=16 entries")
        for bits in record["ring_phases_bits"]:
            _bits(bits, 32, "ring phase")
        for key in ("ring_six", "ring_five"):
            for integer in record[key]:
                _integer(integer, 0, 255, key)
        canonical.append(json.dumps({k: v for k, v in record.items() if k != "index"}, sort_keys=True, separators=(",", ":"), allow_nan=False))
    # Sort whole gene/component pairs. Never sort each channel independently or
    # discard duplicates: either operation could disguise a broken pairing.
    canonical.sort()
    payload = "\n".join(canonical).encode()
    return {"population_count": count, "n_genes": n_genes,
            "order_independent_records_sha256": hashlib.sha256(payload).hexdigest(),
            "distinct_record_count": len(Counter(canonical)), "boundary": value["boundary"],
            "execution": execution}


def inspect_outputs(out: Path, require_gen0: bool, expected_threads: int = 1) -> dict:
    found = sorted(path.name for path in out.glob("d_*.pdb") if re.fullmatch(r"d_[0-9]+\.pdb", path.name))
    expected = sorted(f"d_{rank}.pdb" for rank in range(RANK_COUNT))
    if found != expected:
        raise GateError(f"exactly ranks 0..9 required; found {found}")
    poses = []
    for rank in range(RANK_COUNT):
        path = out / f"d_{rank}.pdb"
        if path.is_symlink():
            raise GateError("output artifacts must not be symlinks")
        normalized, detail = normalize_pose(path)
        poses.append({"rank": rank, "file": path.name, **pin_file(path), **detail,
                      "scientific_sha256": hashlib.sha256(normalized).hexdigest()})
    grid = out / "d.rrg"
    if grid.is_symlink():
        raise GateError("output artifacts must not be symlinks")
    grid_info = {"file": grid.name, **pin_file(grid), **validate_grid(grid)}
    receipt = out / "gen0.json"
    if require_gen0:
        if receipt.is_symlink():
            raise GateError("output artifacts must not be symlinks")
        gen0 = {"file": receipt.name, **pin_file(receipt), **validate_gen0(receipt, expected_threads)}
    else:
        gen0 = {"status": "observer_off_transparency_non_gate"}
        if receipt.exists():
            raise GateError("unexpected generation-zero file when observer was not requested")
    return {"poses": poses, "grid": grid_info, "gen0": gen0}


def verify_pins(pins: dict) -> None:
    for group in ("inputs", "runtime_data"):
        for record in pins[group].values():
            if pin_file(Path(record["path"])) != record:
                raise GateError(f"input changed after pinning: {record['path']}")


def pin_run_log(path: Path) -> dict:
    # Empty stdout/stderr is allowed, but a missing or redirected log is not.
    if path.is_symlink() or not path.is_file():
        raise GateError("run.log must be an existing regular file, not a symlink")
    return {"path": str(path.resolve(strict=True)), "bytes": path.stat().st_size,
            "sha256": sha256(path)}


def validate_run_log(path: Path) -> None:
    """Reject skipped scoring and ownership failures even after child exit zero."""
    carry = b""
    overlap = max(map(len, SCORING_FAILURE_MARKERS)) - 1
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            data = carry + block
            for marker in SCORING_FAILURE_MARKERS:
                if marker in data:
                    raise GateError("run.log reports skipped scoring or an ownership invariant failure: "
                                    + marker.decode("ascii"))
            carry = data[-overlap:]


def run_one(args: argparse.Namespace) -> int:
    args.require_gen0 = not args.observer_off_transparency
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    manifest_path = out / "run_manifest.json"
    manifest = {"schema": RUN_SCHEMA, "label": args.label, "run_dir": str(out), "status": "preflight_failed",
                "started_unix_ns": time.time_ns(), "require_gen0": args.require_gen0, "gate_eligible": args.require_gen0,
                "normalization": "only commit and dirty values in the dedicated FLEXAID provenance REMARK",
                "errors": []}
    try:
        if not math.isfinite(args.timeout) or args.timeout <= 0:
            raise GateError("timeout must be finite and positive")
        if not args.label.strip() or any(c in args.label for c in "\r\n"):
            raise GateError("run label must be nonempty and single-line")
        inputs = {name: pin_file(getattr(args, name)) for name in ("engine", "receptor", "ligand", "config")}
        if not os.access(inputs["engine"]["path"], os.X_OK):
            raise GateError("engine is not executable")
        validate_config(Path(inputs["config"]["path"]))
        data_dir = args.data_dir.resolve(strict=True)
        matrix = data_dir / "MC_st0r5.2_6.dat"
        runtime_data = {matrix.name: pin_file(matrix)}
        for path in sorted(data_dir.iterdir()):
            if path.is_file() and path.suffix.lower() in {".dat", ".def", ".tbl"}:
                runtime_data[path.name] = pin_file(path)
        metadata = source_build_metadata(args.source_dir, args.build_dir)
        (out / "tmp").mkdir()
        environment = {"PATH": os.environ.get("PATH", os.defpath), "HOME": os.environ.get("HOME", str(out)),
                       "TMPDIR": str(out / "tmp"), "LANG": "C", "LC_ALL": "C", "TZ": "UTC",
                       "FLEXAID_SEED": SEED, "FLEXAIDDS_NO_SEC": "1", "FLEXAIDDS_RESTARTS": "1",
                       "FLEXAIDDS_CACHE_READONLY": "1", "FLEXAIDDS_VCF_DIAG": "1",
                       "FLEXAIDDS_DATA_DIR": str(data_dir), "OMP_NUM_THREADS": str(args.omp_threads),
                       "OMP_DYNAMIC": "FALSE", "FLEXAIDDS_PARALLEL_REPRODUCE": "1" if args.parallel_reproduce == "on" else "0"}
        for entry in args.env:
            key, separator, value = entry.partition("=")
            if (not separator or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key) or key in PROTECTED_ENV
                    or key.startswith(SCIENCE_ENV_PREFIXES + ("OMP_", "KMP_", "GOMP_")) or "DETERMINISTIC" in key):
                raise GateError(f"invalid/protected explicit environment override: {entry!r}")
            if key in environment:
                raise GateError(f"duplicate explicit environment key: {key}")
            environment[key] = value
        if args.require_gen0:
            environment["FLEXAIDDS_GEN0_RECEIPT"] = str(out / "gen0.json")
        argv = [inputs["engine"]["path"], inputs["receptor"]["path"], inputs["ligand"]["path"],
                "-c", inputs["config"]["path"], "-o", str(out / "d"), "--data-dir", str(data_dir)]
        manifest.update(inputs=inputs, runtime_data=runtime_data, source_build=metadata, argv=argv,
                        environment=environment, status="running")
        write_json(manifest_path, manifest)
        start = time.monotonic()
        with (out / "run.log").open("wb") as log:
            process = subprocess.Popen(argv, cwd=out, env=environment, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                manifest["exit_code"] = process.wait(timeout=args.timeout)
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
                # Terminate only the process group launched for this run, including
                # child tools. Never leave timed-out engines running untracked.
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                manifest["exit_code"] = process.returncode
                manifest["timed_out"] = isinstance(exc, subprocess.TimeoutExpired)
                raise GateError("engine timeout" if manifest["timed_out"] else "run interrupted")
        manifest["wall_seconds"] = time.monotonic() - start
        manifest["log"] = pin_run_log(out / "run.log")
        validate_run_log(out / "run.log")
        if manifest["exit_code"] != 0:
            raise GateError(f"engine exited {manifest['exit_code']}")
        verify_pins(manifest)
        manifest["outputs"] = inspect_outputs(out, args.require_gen0, args.omp_threads)
        manifest["status"] = "complete"
    except (GateError, OSError, ValueError, subprocess.SubprocessError) as exc:
        manifest["status"] = "failed"
        manifest["errors"].append(str(exc))
    manifest["observed_artifacts"] = []
    for path in sorted(out.iterdir()):
        if path.is_file() and not path.is_symlink() and path.name != "run_manifest.json":
            manifest["observed_artifacts"].append({"file": path.name, "bytes": path.stat().st_size,
                                                   "sha256": sha256(path)})
    manifest["finished_unix_ns"] = time.time_ns()
    write_json(manifest_path, manifest)
    print(json.dumps({"manifest": str(manifest_path), "status": manifest["status"], "errors": manifest["errors"]}))
    return 0 if manifest["status"] == "complete" else 1


def load_run(path: Path) -> dict:
    path = path / "run_manifest.json" if path.is_dir() else path
    manifest = json_read(path)
    out = path.parent.resolve()
    if manifest.get("schema") != RUN_SCHEMA or manifest.get("run_dir") != str(out):
        raise GateError("invalid run schema or manifest/output-directory identity")
    if manifest.get("status") != "complete" or manifest.get("exit_code") != 0 or manifest.get("timed_out"):
        raise GateError(f"run did not complete successfully: {out}")
    if manifest.get("gate_eligible") is not True or manifest.get("require_gen0") is not True:
        raise GateError("observer-off transparency runs cannot pass a gate: generation-zero evidence required")
    unsupported = sorted(key for key in manifest["environment"]
                         if key.startswith(SCIENCE_ENV_PREFIXES) and key not in CANONICAL_SCIENCE_ENV)
    if unsupported:
        raise GateError("unsupported scientific environment keys in saved run: " + ", ".join(unsupported))
    if manifest["environment"].get("FLEXAIDDS_VCF_DIAG") != "1":
        raise GateError("VCF diagnostic guard must be explicitly enabled for a gate")
    log_pin = pin_run_log(out / "run.log")
    validate_run_log(out / "run.log")
    if log_pin != manifest.get("log"):
        raise GateError(f"run.log changed after receipt: {out}")
    verify_pins(manifest)
    current = inspect_outputs(out, manifest.get("require_gen0") is True, int(manifest["environment"]["OMP_NUM_THREADS"]))
    if current != manifest.get("outputs"):
        raise GateError(f"output artifacts changed after receipt: {out}")
    return manifest


def comparison(args: argparse.Namespace) -> int:
    report = {"schema": "flexaidds.engine_repro_comparison.v1", "kind": args.kind, "pass": False, "errors": []}
    try:
        runs = [load_run(path) for path in args.runs]
        if len({run["run_dir"] for run in runs}) != len(runs):
            raise GateError("cannot compare a run with itself")
        report["runs"] = [{"label": run["label"], "run_dir": run["run_dir"], "engine_sha256": run["inputs"]["engine"]["sha256"],
                           "source_commit": run["source_build"]["source_commit"], "outputs": run["outputs"]} for run in runs]
        if args.kind == "parity":
            if len(runs) != 2 or any(run["environment"]["OMP_NUM_THREADS"] != "1" or run["environment"]["FLEXAIDDS_PARALLEL_REPRODUCE"] != "0" for run in runs):
                raise GateError("parity requires exactly two OMP1 runs with parallel reproduction OFF")
            if runs[0]["inputs"]["engine"]["path"] == runs[1]["inputs"]["engine"]["path"]:
                raise GateError("parity requires distinct baseline and candidate engine paths")
            if len({run["require_gen0"] for run in runs}) != 1:
                raise GateError("parity observer requirements must match")
        else:
            if len(runs) != 4 or sorted(run["environment"]["OMP_NUM_THREADS"] for run in runs) != ["1", "1", "4", "4"]:
                raise GateError("determinism requires two OMP1 and two OMP4 runs")
            if any(run["environment"]["FLEXAIDDS_PARALLEL_REPRODUCE"] != "1" or not run["require_gen0"] for run in runs):
                raise GateError("determinism requires reproduction ON and generation-zero receipts in all four runs")
            if len({run["inputs"]["engine"]["sha256"] for run in runs}) != 1:
                raise GateError("determinism requires the same candidate binary")
        reference = runs[0]
        for run in runs[1:]:
            for name in ("receptor", "ligand", "config"):
                if run["inputs"][name]["sha256"] != reference["inputs"][name]["sha256"]:
                    raise GateError(f"{name} differs across runs")
            data = lambda r: {k: v["sha256"] for k, v in r["runtime_data"].items()}
            if data(run) != data(reference):
                raise GateError("runtime data or matrix differs across runs")
            excluded = {"OMP_NUM_THREADS", "FLEXAIDDS_DATA_DIR", "FLEXAIDDS_GEN0_RECEIPT", "TMPDIR"}
            env = lambda r: {k: v for k, v in r["environment"].items() if k not in excluded}
            if env(run) != env(reference):
                raise GateError("scientific environment differs across runs")
            for left, right in zip(reference["outputs"]["poses"], run["outputs"]["poses"]):
                if left["rank"] != right["rank"] or left["scientific_sha256"] != right["scientific_sha256"] or left["cf_text"] != right["cf_text"]:
                    raise GateError(f"rank {left['rank']} scientific payload/CF differs")
            if run["outputs"]["grid"]["sha256"] != reference["outputs"]["grid"]["sha256"]:
                raise GateError("RRG shape or ordered grid payload differs")
            if reference["require_gen0"]:
                for key in ("population_count", "n_genes", "order_independent_records_sha256"):
                    if run["outputs"]["gen0"][key] != reference["outputs"]["gen0"][key]:
                        raise GateError(f"generation-zero exact gene/component population differs: {key}")
        report["raw_pose_bytes_identical"] = all([p["sha256"] for p in run["outputs"]["poses"]] == [p["sha256"] for p in reference["outputs"]["poses"]] for run in runs[1:])
        report["generation_zero"] = "pass" if reference["require_gen0"] else "not_requested_for_parity"
        report["normalization"] = "only FLEXAID.commit and FLEXAID.dirty values; all seed and science bytes retained"
        report["pass"] = True
    except (GateError, OSError, ValueError, KeyError) as exc:
        report["errors"].append(str(exc))
    if args.json:
        write_json(args.json, report)
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run-one")
    for name in ("engine", "source-dir", "build-dir", "receptor", "ligand", "config", "data-dir", "out"):
        run.add_argument("--" + name, type=Path, required=True)
    run.add_argument("--label", required=True)
    run.add_argument("--omp-threads", type=int, choices=(1, 4), required=True)
    run.add_argument("--parallel-reproduce", choices=("off", "on"), required=True)
    run.add_argument("--require-gen0", action="store_true", help="Generation-zero receipt is mandatory for gates (default)")
    run.add_argument("--observer-off-transparency", action="store_true", help="Non-gate diagnostic only; cannot pass compare")
    run.add_argument("--timeout", type=float, default=21600)
    run.add_argument("--env", action="append", default=[],
                     help="Additional nonscientific environment only; FLEXAID_/FLEXAIDDS_ overrides are forbidden")
    compare = commands.add_parser("compare")
    compare.add_argument("--kind", choices=("parity", "determinism"), required=True)
    compare.add_argument("--runs", type=Path, nargs="+", required=True)
    compare.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    try:
        return run_one(args) if args.command == "run-one" else comparison(args)
    except (GateError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
