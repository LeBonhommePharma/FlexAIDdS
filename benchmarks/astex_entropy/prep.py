from __future__ import annotations

import csv
import json
import shutil
import urllib.request
from pathlib import Path
from typing import Any

import yaml

from .io_utils import run_command, write_targets
from .models import TargetRecord
from .tools import prepare_tool_inputs


METAL_RESNAMES = {
    "MG", "MN", "ZN", "CA", "FE", "FE2", "CU", "NA", "K", "CO", "NI", "CL", "BR", "I",
}
WATER_RESNAMES = {"HOH", "WAT", "DOD"}


def _load_yaml_targets(path: Path) -> list[str]:
    data = _load_dataset_yaml(path)
    return [str(x).strip().upper() for x in data.get("targets", []) if str(x).strip()]


def _load_dataset_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Required dataset YAML not found: {path}")
    with path.open() as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Dataset YAML must be a mapping: {path}")
    return data


def _load_native_target_ids(cfg: dict[str, Any]) -> list[str]:
    native_csv = Path(cfg["datasets"]["native_csv"])
    if native_csv.exists():
        with native_csv.open(newline="") as fh:
            return [row["pdb_id"].strip().upper() for row in csv.DictReader(fh) if row.get("pdb_id")]
    return _load_yaml_targets(Path(cfg["datasets"]["native_yaml"]))


def _copy(src: Path, dst: Path, force: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and not force:
        return
    shutil.copy2(src, dst)


def _download_pdb(code: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{code.upper()}.pdb"
    with urllib.request.urlopen(url, timeout=60) as response:
        dst.write_bytes(response.read())


def _local_astex_dir(cfg: dict[str, Any], pdb_id: str) -> Path | None:
    root = Path(cfg["datasets"]["native_data_root"])
    for name in (pdb_id.upper(), pdb_id.lower()):
        candidate = root / name
        if candidate.exists():
            return candidate
    return None


def _filter_receptor_pdb(src: Path, dst: Path, ligand_resname: str, force: bool) -> None:
    if dst.exists() and not force:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    ligand_resname = ligand_resname.upper()
    out: list[str] = []
    for line in src.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if rec == "ATOM":
            out.append(line)
        elif rec == "HETATM" and resname in METAL_RESNAMES and resname not in WATER_RESNAMES:
            out.append(line)
        elif rec in {"TER", "END"}:
            out.append(line)
    if not out or not any(line.startswith("ATOM") for line in out):
        raise RuntimeError(f"No protein atoms found while preparing receptor: {src}")
    if not out[-1].startswith("END"):
        out.append("END")
    dst.write_text("\n".join(out) + "\n")


def _extract_ligand_sdf(
    src_pdb: Path,
    ligand_resname: str,
    ligand_pdb: Path,
    ligand_sdf: Path,
    obabel: str,
    force: bool,
) -> None:
    if ligand_sdf.exists() and not force:
        return
    ligand_resname = ligand_resname.upper()
    lines: list[str] = []
    conect: list[str] = []
    for line in src_pdb.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if rec == "HETATM" and resname == ligand_resname:
            lines.append(line)
        elif rec == "CONECT":
            conect.append(line)
    if not lines:
        raise RuntimeError(f"Ligand {ligand_resname} not found in {src_pdb}")
    ligand_pdb.parent.mkdir(parents=True, exist_ok=True)
    ligand_pdb.write_text("\n".join(lines + conect + ["END"]) + "\n")
    run_command(
        [obabel, "-ipdb", str(ligand_pdb), "-osdf", "-O", str(ligand_sdf), "-h"],
        log_path=ligand_sdf.with_suffix(".obabel.log"),
    )


def prepare_native(
    cfg: dict[str, Any],
    *,
    max_targets: int | None = None,
    force: bool = False,
) -> list[TargetRecord]:
    work = Path(cfg["work_dir"])
    prepared = work / "prepared" / "native"
    native_root = Path(cfg["datasets"]["native_data_root"])
    dataset_yaml = _load_dataset_yaml(Path(cfg["datasets"]["native_yaml"]))
    target_ids = _load_native_target_ids(cfg)
    if max_targets:
        target_ids = target_ids[:max_targets]

    records: list[TargetRecord] = []
    for target_id in target_ids:
        src_dir = native_root / target_id
        receptor = src_dir / f"{target_id}_apo.pdb"
        ligand = src_dir / f"{target_id}_ligand.sdf"
        if not receptor.exists() or not ligand.exists():
            continue

        target_dir = prepared / target_id
        receptor_dst = target_dir / "receptor.pdb"
        ligand_dst = target_dir / "ligand.sdf"
        reference_dst = target_dir / "reference.sdf"
        _copy(receptor, receptor_dst, force)
        _copy(ligand, ligand_dst, force)
        _copy(ligand, reference_dst, force)
        source_complex = src_dir / f"{target_id}.pdb"
        pocket = src_dir / f"{target_id}_binding_site.pdb"
        record = TargetRecord(
            target_id=target_id,
            mode="native",
            receptor_pdb=str(receptor_dst),
            ligand_sdf=str(ligand_dst),
            reference_sdf=str(reference_dst),
            source_complex=str(source_complex) if source_complex.exists() else "",
            pocket_pdb=str(pocket) if pocket.exists() else "",
            cavity_source_pdb=str(source_complex) if source_complex.exists() else "",
            cavity_ligand_sdf=str(reference_dst),
        )
        prepare_tool_inputs(record, cfg, target_dir, force=force)
        records.append(record)

    manifest = work / "manifests" / "native_manifest.csv"
    write_targets(records, manifest)
    _write_manifest_provenance(manifest, cfg["datasets"]["native_yaml"], dataset_yaml, len(records))
    return records


def prepare_non_native(
    cfg: dict[str, Any],
    *,
    max_targets: int | None = None,
    download_missing: bool = False,
    force: bool = False,
) -> list[TargetRecord]:
    work = Path(cfg["work_dir"])
    raw_root = Path(cfg["datasets"]["non_native_data_root"])
    prepared = work / "prepared" / "non_native"
    csv_path = Path(cfg["datasets"]["non_native_csv"])
    obabel = cfg["tools"]["obabel"]
    dataset_yaml = _load_dataset_yaml(Path(cfg["datasets"]["non_native_yaml"]))

    rows: list[dict[str, str]]
    with csv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    records: list[TargetRecord] = []
    skipped: list[dict[str, str]] = []
    for row in rows:
        if max_targets and len(records) >= max_targets:
            break
        target_pdb = row["target_pdb"].strip().lower()
        ligand_pdb = row["ligand_pdb"].strip().lower()
        ligand_id = row["ligand_id"].strip().upper()
        pair_id = f"{target_pdb.upper()}_x_{ligand_pdb.upper()}_{ligand_id}"
        target_raw = raw_root / f"{target_pdb}.pdb"
        ligand_raw = raw_root / f"{ligand_pdb}.pdb"
        local_target_dir = _local_astex_dir(cfg, target_pdb)
        local_ligand_dir = _local_astex_dir(cfg, ligand_pdb)
        local_target_apo = (local_target_dir / f"{target_pdb.upper()}_apo.pdb") if local_target_dir else None
        local_target_complex = (local_target_dir / f"{target_pdb.upper()}.pdb") if local_target_dir else None
        local_target_site = (local_target_dir / f"{target_pdb.upper()}_binding_site.pdb") if local_target_dir else None
        local_target_ligand_sdf = (local_target_dir / f"{target_pdb.upper()}_ligand.sdf") if local_target_dir else None
        local_ligand_sdf = (local_ligand_dir / f"{ligand_pdb.upper()}_ligand.sdf") if local_ligand_dir else None
        local_ligand_complex = (local_ligand_dir / f"{ligand_pdb.upper()}.pdb") if local_ligand_dir else None

        if download_missing:
            try:
                if (not local_target_apo or not local_target_apo.exists()) and not target_raw.exists():
                    _download_pdb(target_pdb, target_raw)
                if (not local_ligand_sdf or not local_ligand_sdf.exists()) and not ligand_raw.exists():
                    _download_pdb(ligand_pdb, ligand_raw)
            except Exception as exc:
                skipped.append({"pair_id": pair_id, "reason": f"download_failed: {exc}"})
                continue
        has_local_target = bool(local_target_apo and local_target_apo.exists())
        has_local_ligand = bool(local_ligand_sdf and local_ligand_sdf.exists())
        if not has_local_target and not target_raw.exists():
            skipped.append({"pair_id": pair_id, "reason": "missing_pdb_file"})
            continue
        if not has_local_ligand and not ligand_raw.exists():
            skipped.append({"pair_id": pair_id, "reason": "missing_ligand_source"})
            continue

        target_dir = prepared / pair_id
        receptor_dst = target_dir / "receptor.pdb"
        ligand_pdb_tmp = target_dir / "ligand.pdb"
        ligand_sdf = target_dir / "ligand.sdf"
        reference_sdf = target_dir / "reference.sdf"
        try:
            if has_local_target:
                _copy(local_target_apo, receptor_dst, force)
            else:
                _filter_receptor_pdb(target_raw, receptor_dst, ligand_id, force)
            if has_local_ligand:
                _copy(local_ligand_sdf, ligand_sdf, force)
            else:
                _extract_ligand_sdf(ligand_raw, ligand_id, ligand_pdb_tmp, ligand_sdf, obabel, force)
            _copy(ligand_sdf, reference_sdf, force)
        except RuntimeError as exc:
            skipped.append({"pair_id": pair_id, "reason": str(exc)})
            continue

        record = TargetRecord(
            target_id=pair_id,
            mode="non_native",
            receptor_pdb=str(receptor_dst),
            ligand_sdf=str(ligand_sdf),
            reference_sdf=str(reference_sdf),
            source_complex=str(local_ligand_complex) if local_ligand_complex and local_ligand_complex.exists() else str(ligand_raw),
            pocket_pdb=str(local_target_site) if local_target_site and local_target_site.exists() else "",
            cavity_source_pdb=str(local_target_complex) if local_target_complex and local_target_complex.exists() else str(target_raw),
            cavity_ligand_sdf=str(local_target_ligand_sdf) if local_target_ligand_sdf and local_target_ligand_sdf.exists() else str(reference_sdf),
            ligand_source_id=ligand_pdb.upper(),
            notes=row.get("target_name", ""),
        )
        try:
            prepare_tool_inputs(record, cfg, target_dir, force=force)
        except RuntimeError as exc:
            skipped.append({"pair_id": pair_id, "reason": str(exc)})
            continue
        records.append(record)

    manifest = work / "manifests" / "non_native_manifest.csv"
    write_targets(records, manifest)
    _write_manifest_provenance(manifest, cfg["datasets"]["non_native_yaml"], dataset_yaml, len(records))
    manifest.with_suffix(".skipped.json").write_text(json.dumps(skipped, indent=2, sort_keys=True) + "\n")
    return records


def prepare_data(
    cfg: dict[str, Any],
    mode: str,
    *,
    max_targets: int | None = None,
    download_missing: bool = False,
    force: bool = False,
) -> dict[str, int]:
    counts: dict[str, int] = {}
    if mode in {"native", "all"}:
        counts["native"] = len(prepare_native(cfg, max_targets=max_targets, force=force))
    if mode in {"non_native", "all"}:
        counts["non_native"] = len(
            prepare_non_native(
                cfg,
                max_targets=max_targets,
                download_missing=download_missing,
                force=force,
            )
        )
    return counts


def _write_manifest_provenance(
    manifest_csv: Path,
    dataset_yaml_path: str | Path,
    dataset_yaml: dict[str, Any],
    n_records: int,
) -> None:
    payload = {
        "manifest_csv": str(manifest_csv),
        "dataset_yaml": str(dataset_yaml_path),
        "dataset_slug": dataset_yaml.get("slug") or dataset_yaml.get("name") or "",
        "dataset_kind": dataset_yaml.get("kind", ""),
        "n_records": n_records,
        "published_baselines": dataset_yaml.get("published_baselines", {}),
    }
    manifest_csv.with_suffix(".provenance.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
