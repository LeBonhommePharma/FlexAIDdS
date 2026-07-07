from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class TargetRecord:
    target_id: str
    mode: str
    receptor_pdb: str
    ligand_sdf: str
    reference_sdf: str
    source_complex: str = ""
    pocket_pdb: str = ""
    cavity_source_pdb: str = ""
    cleft_sphere_pdb: str = ""
    ligand_source_id: str = ""
    notes: str = ""

    @property
    def target_dir_name(self) -> str:
        return self.target_id.replace("/", "_").replace(" ", "_")

    @classmethod
    def from_row(cls, row: dict[str, Any]) -> "TargetRecord":
        fields = {name: str(row.get(name, "") or "") for name in cls.__dataclass_fields__}
        return cls(**fields)

    def to_dict(self) -> dict[str, str]:
        return {k: str(v) for k, v in asdict(self).items()}


@dataclass(frozen=True)
class PoseRecord:
    target_id: str
    mode: str
    tool: str
    pose_id: str
    pose_sdf: str
    receptor_pdb: str
    reference_sdf: str
    raw_score: str = ""
    score_direction: str = "lower"
    source_file: str = ""
    status: str = "ok"

    @classmethod
    def from_row(cls, row: dict[str, Any]) -> "PoseRecord":
        fields = {name: str(row.get(name, "") or "") for name in cls.__dataclass_fields__}
        return cls(**fields)

    def to_dict(self) -> dict[str, str]:
        return {k: str(v) for k, v in asdict(self).items()}


def existing_file(path: str | Path) -> bool:
    return bool(path) and Path(path).is_file()
