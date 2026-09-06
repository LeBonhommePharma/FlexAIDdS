"""Synthetic executable/artifact tests; these never invoke FlexAIDdS."""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import struct
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ops/engine_repro_gate.py"
WRAPPER = ROOT / "ops/gate_parity.sh"
spec = importlib.util.spec_from_file_location("engine_repro_gate", SCRIPT)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def bits(value):
    return "0x" + struct.pack(">d", value).hex()


def execution(threads=1):
    return dict(openmp_compiled=True, deterministic_compile=False,
                evaluation_batches=[dict(region="populate_chromosomes", population_count=1000,
                                         popoffset=0, workspace_slots=threads,
                                         workers=[dict(worker_id=i, team_size=threads,
                                                       evaluated_chromosomes=1000 // threads)
                                                  for i in range(threads)])])


def gen0(count=1000, threads=1):
    records = []
    for index in range(count):
        record = dict(index=index, status=110, genes=[dict(to_int32=index, to_ic_bits=bits(float(index)))],
                      cf={key: bits(0.0) for key in gate.CF_FIELDS - {"rclash"}},
                      ring_phases_bits=["0x00000000"] * 16, ring_six=[0] * 16, ring_five=[0] * 16)
        record["cf"].update(rclash=0, com_bits=bits(-float(index + 1)))
        for key in ("evalue_bits", "app_evalue_bits", "fitnes_bits", "boltzmann_weight_bits", "free_energy_bits"):
            record[key] = bits(float(index))
        records.append(record)
    return dict(schema=gate.GEN0_SCHEMA, boundary="initial_population_complete_before_reproduction", generation=0,
                population_count=count, n_genes=1, seed="12345", execution=execution(threads), records=records, complete=True)


def fake_engine(path, mode="normal", commit="a" * 40, dirty="0"):
    path.parent.mkdir(parents=True, exist_ok=True)
    code = f'''#!{sys.executable}
import copy,json,os,struct,sys,time
from pathlib import Path
MODE={mode!r}
COMMIT={commit!r}
DIRTY={dirty!r}
receipt={gen0(1)!r}
template=receipt['records'][0]
receipt['records']=[]
receipt['population_count']=1000
for index in range(1000):
    record=copy.deepcopy(template)
    record['index']=index
    record['genes'][0]={{'to_int32':index,'to_ic_bits':'0x'+struct.pack('>d',float(index)).hex()}}
    record['cf']['com_bits']='0x'+struct.pack('>d',-float(index+1)).hex()
    for key in ['evalue_bits','app_evalue_bits','fitnes_bits','boltzmann_weight_bits','free_energy_bits']:
        record[key]='0x'+struct.pack('>d',float(index)).hex()
    receipt['records'].append(record)
args=sys.argv[1:]
prefix=Path(args[args.index('-o')+1])
if MODE=='timeout':
    time.sleep(30)
if MODE=='no_artifacts':
    print('no output');sys.exit(0)
for rank in range(11 if MODE=='extra_rank' else (9 if MODE=='missing_rank' else 10)):
    seed='12346' if MODE=='wrong_seed' else '12345'
    score='nan' if MODE=='nan_cf' else '-1.00000'
    x=2.0 if MODE=='different_coordinate' else 1.0
    atom=f"HETATM    1  C1  LIG A 900    {{x:8.3f}}{{2.0:8.3f}}{{3.0:8.3f}}  1.00  0.00           C\\n"
    body=f"REMARK FLEXAID.commit={{COMMIT}} FLEXAID.dirty={{DIRTY}} FLEXAID.seed={{seed}}\\nREMARK CF={{score}}\\nREMARK CF.com=-1.00000\\n"
    if MODE=='science_remark': body+='REMARK temperature=299\\n'
    body+=atom+'END\\n'
    if MODE=='no_end': body=body.replace('END\\n','')
    if MODE=='no_cf': body=body.replace(f'REMARK CF={{score}}\\n','')
    if MODE=='no_atoms': body='\\n'.join(line for line in body.splitlines() if not line.startswith('HETATM'))+'\\n'
    if MODE=='unrecognized_provenance': body=body.replace('FLEXAID.dirty=', 'OTHER.dirty=')
    Path(str(prefix)+f'_{{rank}}.pdb').write_text(body)
if MODE!='missing_grid':
    grid=struct.pack('<IIi',0x56435400,1,2)+struct.pack('<ii6f',0,10,1,2,3,4,5,6)+struct.pack('<ii6f',1,11,7,8,9,10,11,12)
    if MODE=='reordered_grid': grid=grid[:12]+grid[44:]+grid[12:44]
    if MODE=='truncated_grid': grid=grid[:-1]
    Path(str(prefix)+'.rrg').write_bytes(grid)
if os.environ.get('FLEXAIDDS_GEN0_RECEIPT') and MODE!='missing_gen0':
    threads=int(os.environ['OMP_NUM_THREADS'])
    receipt['execution']['evaluation_batches'][0]['workspace_slots']=threads
    receipt['execution']['evaluation_batches'][0]['workers']=[{{'worker_id':i,'team_size':threads,'evaluated_chromosomes':1000//threads}} for i in range(threads)]
    if MODE=='serialized_team':
        receipt['execution']['evaluation_batches'][0]['workers']=[{{'worker_id':0,'team_size':1,'evaluated_chromosomes':1000}}]
    if os.environ.get('OMP_NUM_THREADS')=='4':
        receipt['records'].reverse()
        for i,record in enumerate(receipt['records']): record['index']=i
    if MODE=='mismatched_gene_pair' and os.environ.get('OMP_NUM_THREADS')=='4':
        genes=[record['genes'] for record in receipt['records']]
        receipt['records'][0]['genes'],receipt['records'][1]['genes']=genes[1],genes[0]
    if MODE=='incomplete_gen0': receipt['complete']=False
    Path(os.environ['FLEXAIDDS_GEN0_RECEIPT']).write_text(json.dumps(receipt))
print('synthetic fake-executable run completed')
if MODE=='exit_nonzero': sys.exit(7)
'''
    path.write_text(code)
    path.chmod(0o755)
    return path


@pytest.fixture
def setup(tmp_path):
    config = tmp_path / "parity.json"
    config.write_text(json.dumps(dict(ga=dict(num_chromosomes=1000, num_generations=2000),
                                     reference_ligand=dict(pose_seed_enabled=False, seed_fraction=0))))
    receptor = tmp_path / "receptor.pdb"; receptor.write_text("synthetic receptor fixture\n")
    ligand = tmp_path / "ligand.sdf"; ligand.write_text("synthetic ligand fixture\n")
    data = tmp_path / "data"; data.mkdir(); (data / "MC_st0r5.2_6.dat").write_text("synthetic matrix fixture\n")
    build = tmp_path / "build"; build.mkdir()
    (build / "CMakeCache.txt").write_text(f"CMAKE_HOME_DIRECTORY:INTERNAL={ROOT}\nCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\nFLEXAIDS_USE_OPENMP:BOOL=ON\n")
    (build / "Makefile").write_text("# synthetic generated build metadata\n")
    (build / "CMakeFiles/FlexAIDdS.dir").mkdir(parents=True)
    (build / "CMakeFiles/FlexAIDdS.dir/flags.make").write_text("CXX_FLAGS = -O2 -fopenmp\n")
    return dict(root=tmp_path, config=config, receptor=receptor, ligand=ligand, data=data, build=build)


def run_one(setup, label, engine=None, mode="normal", omp=1, reproduce="off", observer=True, timeout=10, extra=()):
    engine = engine or fake_engine(setup["root"] / label / "FlexAIDdS", mode)
    out = setup["root"] / "runs" / label
    args = [sys.executable, "-B", str(SCRIPT), "run-one", "--label", label, "--engine", str(engine),
            "--source-dir", str(ROOT), "--build-dir", str(setup["build"]),
            "--receptor", str(setup["receptor"]), "--ligand", str(setup["ligand"]),
            "--config", str(setup["config"]), "--data-dir", str(setup["data"]), "--out", str(out),
            "--omp-threads", str(omp), "--parallel-reproduce", reproduce, "--timeout", str(timeout), *extra]
    args.append("--require-gen0" if observer else "--observer-off-transparency")
    result = subprocess.run(args, capture_output=True, text=True, timeout=20)
    manifest = json.loads((out / "run_manifest.json").read_text()) if (out / "run_manifest.json").is_file() else None
    return result, manifest, out


def compare(setup, runs, kind="parity"):
    report = setup["root"] / "comparison.json"
    result = subprocess.run([sys.executable, "-B", str(SCRIPT), "compare", "--kind", kind,
                             "--runs", *map(str, runs), "--json", str(report)], capture_output=True, text=True, timeout=20)
    return result, json.loads(report.read_text())


def test_same_binary_basename_has_distinct_runs_and_only_provenance_normalizes(setup):
    a = fake_engine(setup["root"] / "base/FlexAIDdS", commit="a" * 40)
    b = fake_engine(setup["root"] / "candidate/FlexAIDdS", commit="b" * 40, dirty="1")
    ra, ma, pa = run_one(setup, "P0", a)
    rb, mb, pb = run_one(setup, "P1", b)
    assert ra.returncode == rb.returncode == 0, (ra.stdout, rb.stdout)
    result, report = compare(setup, [pa, pb])
    assert result.returncode == 0, report
    assert report["pass"]
    assert not report["raw_pose_bytes_identical"]
    assert report["generation_zero"] == "pass"
    assert ma["inputs"]["engine"]["sha256"] != mb["inputs"]["engine"]["sha256"]
    assert ma["source_build"]["source_commit"]
    assert len(ma["outputs"]["poses"]) == 10
    assert ma["environment"]["FLEXAID_SEED"] == "12345"


@pytest.mark.parametrize("mode", ["missing_rank", "extra_rank", "no_artifacts", "nan_cf", "no_cf", "no_atoms",
                                  "wrong_seed", "unrecognized_provenance", "missing_grid", "truncated_grid", "exit_nonzero", "no_end"])
def test_incomplete_or_failed_fake_engine_never_passes(setup, mode):
    result, manifest, out = run_one(setup, "bad", mode=mode)
    assert result.returncode == 1
    assert manifest["status"] == "failed"
    assert manifest["errors"]
    if mode == "exit_nonzero":
        assert manifest["exit_code"] == 7
        assert len(manifest["observed_artifacts"]) >= 11


def test_reused_output_directory_is_refused_without_overwrite(setup):
    result, manifest, out = run_one(setup, "one")
    assert result.returncode == 0
    before = (out / "run_manifest.json").read_bytes()
    result, _, _ = run_one(setup, "one")
    assert result.returncode == 2
    assert (out / "run_manifest.json").read_bytes() == before


def test_timeout_is_failure_with_exit_receipt(setup):
    result, manifest, _ = run_one(setup, "timeout", mode="timeout", timeout=0.1)
    assert result.returncode == 1
    assert manifest["timed_out"]
    assert manifest["exit_code"] < 0


@pytest.mark.parametrize("mode", ["different_coordinate", "science_remark", "reordered_grid"])
def test_scientific_payload_and_grid_order_differences_fail(setup, mode):
    _, _, left = run_one(setup, "left")
    _, _, right = run_one(setup, "right", mode=mode)
    result, report = compare(setup, [left, right])
    assert result.returncode == 1
    assert not report["pass"]


def test_changed_input_after_run_is_refused(setup):
    _, _, left = run_one(setup, "left")
    _, _, right = run_one(setup, "right")
    setup["ligand"].write_text("changed ligand\n")
    result, report = compare(setup, [left, right])
    assert result.returncode == 1
    assert "input changed" in report["errors"][0]


def test_artifact_changed_after_manifest_is_refused(setup):
    _, _, left = run_one(setup, "left")
    _, _, right = run_one(setup, "right")
    artifact = right / "d_0.pdb"
    artifact.write_bytes(artifact.read_bytes().replace(b"END\n", b"REMARK changed\nEND\n"))
    result, report = compare(setup, [left, right])
    assert result.returncode == 1
    assert "artifacts changed" in report["errors"][0]


def test_cannot_compare_run_with_itself(setup):
    _, _, out = run_one(setup, "same")
    result, report = compare(setup, [out, out])
    assert result.returncode == 1
    assert "itself" in report["errors"][0]


@pytest.mark.parametrize("mode", ["missing_gen0", "incomplete_gen0"])
def test_missing_or_incomplete_requested_gen0_never_passes(setup, mode):
    result, manifest, _ = run_one(setup, "bad", mode=mode, observer=True, reproduce="on")
    assert result.returncode == 1
    assert manifest["status"] == "failed"


def test_full_four_run_matrix_compares_paired_record_multisets(setup):
    engine = fake_engine(setup["root"] / "candidate/FlexAIDdS")
    paths = []
    for label, omp in (("T1a", 1), ("T1b", 1), ("T4a", 4), ("T4b", 4)):
        result, manifest, path = run_one(setup, label, engine, omp=omp, reproduce="on", observer=True)
        assert result.returncode == 0, manifest["errors"]
        paths.append(path)
    result, report = compare(setup, paths, "determinism")
    assert result.returncode == 0, report
    assert report["generation_zero"] == "pass"
    assert gate.sha256(paths[0] / "gen0.json") != gate.sha256(paths[2] / "gen0.json")
    assert report["runs"][2]["outputs"]["gen0"]["execution"]["observed_threads"] == 4


def test_requested_four_threads_cannot_pass_with_a_serialized_actual_team(setup):
    result, manifest, _ = run_one(setup, "serialized", mode="serialized_team", omp=4, reproduce="on")
    assert result.returncode == 1
    assert "every requested OpenMP worker" in manifest["errors"][0]


@pytest.mark.parametrize("mutator", [
    lambda e: e.update(openmp_compiled=False),
    lambda e: e.update(deterministic_compile=True),
    lambda e: e.update(evaluation_batches=[]),
    lambda e: e["evaluation_batches"][0].update(region="calculate_fitness"),
    lambda e: e["evaluation_batches"][0].update(popoffset=1),
    lambda e: e["evaluation_batches"][0].update(population_count=999),
    lambda e: e["evaluation_batches"][0].update(workspace_slots=3),
    lambda e: e["evaluation_batches"][0]["workers"][0].update(evaluated_chromosomes=249),
    lambda e: e["evaluation_batches"][0]["workers"][0].update(evaluated_chromosomes=0),
    lambda e: e["evaluation_batches"][0]["workers"][0].update(worker_id=1),
    lambda e: e["evaluation_batches"][0]["workers"][0].update(team_size=1),
    lambda e: e["evaluation_batches"][0]["workers"][0].update(extra="unknown"),
])
def test_incomplete_or_contradictory_worker_evidence_fails(mutator):
    value = execution(4)
    mutator(value)
    with pytest.raises(gate.GateError):
        gate.validate_execution(value, 4)


def test_worker_scheduling_distribution_does_not_change_scientific_hash(tmp_path):
    first, second = gen0(threads=4), gen0(threads=4)
    workers = second["execution"]["evaluation_batches"][0]["workers"]
    workers[0]["evaluated_chromosomes"] = 249
    workers[1]["evaluated_chromosomes"] = 251
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps(first)); b.write_text(json.dumps(second))
    left, right = gate.validate_gen0(a, 4), gate.validate_gen0(b, 4)
    assert left["execution"] != right["execution"]
    assert left["order_independent_records_sha256"] == right["order_independent_records_sha256"]


def test_gene_component_repairing_cannot_hide_behind_equal_component_sets(setup):
    engine = fake_engine(setup["root"] / "candidate/FlexAIDdS", mode="mismatched_gene_pair")
    paths = []
    for label, omp in (("T1a", 1), ("T1b", 1), ("T4a", 4), ("T4b", 4)):
        result, _, path = run_one(setup, label, engine, omp=omp, reproduce="on", observer=True)
        assert result.returncode == 0
        paths.append(path)
    result, report = compare(setup, paths, "determinism")
    assert result.returncode == 1
    assert "gene/component" in report["errors"][0]


def test_determinism_without_gen0_is_not_a_pass(setup):
    engine = fake_engine(setup["root"] / "candidate/FlexAIDdS")
    paths = []
    for label, omp in (("T1a", 1), ("T1b", 1), ("T4a", 4), ("T4b", 4)):
        _, _, path = run_one(setup, label, engine, omp=omp, reproduce="on", observer=False)
        paths.append(path)
    result, report = compare(setup, paths, "determinism")
    assert result.returncode == 1
    assert "generation-zero evidence required" in report["errors"][0]


@pytest.mark.parametrize("mutator", [
    lambda r: r.update(complete=False),
    lambda r: r.update(population_count=3),
    lambda r: (r.update(population_count=999), r["records"].pop()),
    lambda r: r["records"][0].update(status=111),
    lambda r: r["records"][0].update(unexpected="field"),
    lambda r: r["records"][1].update(index=0),
    lambda r: r["records"][0]["genes"][0].update(to_ic_bits="0x7ff0000000000000"),
    lambda r: r["records"][0]["cf"].pop("wal_bits"),
    lambda r: r["records"][0].update(ring_phases_bits=[]),
])
def test_gen0_schema_and_nonfinite_failures(tmp_path, mutator):
    value = gen0(); mutator(value)
    path = tmp_path / "gen0.json"; path.write_text(json.dumps(value))
    with pytest.raises(gate.GateError): gate.validate_gen0(path)


def test_gen0_duplicate_multiplicity_is_preserved(tmp_path):
    first = gen0()
    second = copy.deepcopy(first)
    second["records"][1] = copy.deepcopy(second["records"][0]); second["records"][1]["index"] = 1
    a = tmp_path / "a.json"; a.write_text(json.dumps(first))
    b = tmp_path / "b.json"; b.write_text(json.dumps(second))
    assert gate.validate_gen0(a)["order_independent_records_sha256"] != gate.validate_gen0(b)["order_independent_records_sha256"]


def test_inherited_science_environment_is_cleared(setup, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_ORACLE_SITE", "unintended")
    monkeypatch.setenv("FLEXAID_SEED", "98765")
    result, manifest, _ = run_one(setup, "clean_env")
    assert result.returncode == 0
    assert "FLEXAIDDS_ORACLE_SITE" not in manifest["environment"]
    assert manifest["environment"]["FLEXAID_SEED"] == "12345"


def test_protected_environment_cannot_override_seed(setup):
    result, manifest, _ = run_one(setup, "override", extra=("--env", "FLEXAID_SEED=2"))
    assert result.returncode == 1
    assert "protected" in manifest["errors"][0]


def test_config_cannot_reuse_grid_or_seed_pose(setup):
    data = json.loads(setup["config"].read_text()); data["grid_file"] = "old.rrg"
    setup["config"].write_text(json.dumps(data))
    result, manifest, _ = run_one(setup, "cached")
    assert result.returncode == 1
    assert "reuse" in manifest["errors"][0]


def test_build_metadata_must_match_declared_source(setup):
    (setup["build"] / "CMakeCache.txt").write_text("CMAKE_HOME_DIRECTORY:INTERNAL=/different/source\n")
    result, manifest, _ = run_one(setup, "wrong_build")
    assert result.returncode == 1
    assert "source directory" in manifest["errors"][0]


def test_wrapper_returns_nonzero_for_failed_engine_and_uses_distinct_names(setup):
    a = fake_engine(setup["root"] / "base/FlexAIDdS")
    b = fake_engine(setup["root"] / "candidate/FlexAIDdS", mode="exit_nonzero")
    out = setup["root"] / "wrapper"
    args = ["bash", str(WRAPPER), str(a), str(b), "--baseline-source", str(ROOT), "--candidate-source", str(ROOT),
            "--baseline-build", str(setup["build"]), "--candidate-build", str(setup["build"]),
            "--receptor", str(setup["receptor"]), "--ligand", str(setup["ligand"]),
            "--config", str(setup["config"]), "--data-dir", str(setup["data"]), "--out", str(out)]
    result = subprocess.run(args, env=dict(os.environ, PYTHON=sys.executable), capture_output=True, text=True, timeout=20)
    assert result.returncode != 0
    assert (out / "P0/run_manifest.json").is_file()
    assert (out / "P1/run_manifest.json").is_file()
    assert not (out / "parity.json").exists()



@pytest.mark.parametrize("timeout", ["nan", "inf", "0", "-1"])
def test_invalid_timeout_rejected_before_start(setup, timeout):
    result, manifest, out = run_one(setup, "bad_timeout", timeout=timeout)
    assert result.returncode == 1
    assert "finite and positive" in manifest["errors"][0]
    assert not (out / "run.log").exists()


def test_missing_cmake_source_is_not_current_directory_by_default(setup):
    (setup["build"] / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\n")
    result, manifest, _ = run_one(setup, "bad_cache")
    assert result.returncode == 1
    assert "source directory" in manifest["errors"][0]



@pytest.mark.parametrize("key", ["OMP_THREAD_LIMIT=1", "FLEXAID_DETERMINISTIC=1", "KMP_HW_SUBSET=1c",
                               "FLEXAIDDS_VORONOI_KEYED_JITTER=1"])
def test_parallel_disabling_overrides_are_rejected(setup, key):
    result, manifest, _ = run_one(setup, "disabled", extra=("--env", key))
    assert result.returncode == 1
    assert "protected" in manifest["errors"][0]


@pytest.mark.parametrize("flags", ["CXX_FLAGS = -O2\n", "CXX_FLAGS = -fopenmp -DFLEXAID_DETERMINISTIC=1\n"])
def test_serial_or_compile_time_deterministic_builds_are_rejected(setup, flags):
    (setup["build"] / "CMakeFiles/FlexAIDdS.dir/flags.make").write_text(flags)
    result, manifest, _ = run_one(setup, "serial_build")
    assert result.returncode == 1
    assert "compile flags" in manifest["errors"][0]


def test_overflowed_json_number_is_not_a_finite_receipt_field(tmp_path):
    path = tmp_path / "overflow.json"
    path.write_text('{"number":1e999}')
    with pytest.raises(gate.GateError, match="nonfinite JSON"):
        gate.json_read(path)
