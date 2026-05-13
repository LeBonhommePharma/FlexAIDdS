"""Parser for tENCoM vibrational entropy tool output.

Reads PDB REMARK sections and JSON output from the tENCoM tool,
returning FlexModeResult and FlexPopulationResult dataclasses.

Usage:
    from flexaidds.tencom_results import parse_tencom_pdb, parse_tencom_json

    # Parse a single PDB output
    mode = parse_tencom_pdb("tencom_mode_0.pdb")
    print(mode.S_vib, mode.delta_F_vib)

    # Parse JSON results
    population = parse_tencom_json("tencom_results.json")
    for mode in population.modes:
        print(mode.mode_id, mode.delta_S_vib)
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class EigenvalueDiff:
    """Per-mode eigenvalue differential."""
    mode: int
    delta_eigenvalue: float
    overlap: Optional[float] = None


@dataclass
class FlexModeResult:
    """Result for a single flexibility mode (reference or target)."""
    mode_id: int = 0
    mode_type: str = ""             # "reference" or "target"
    source: str = ""                # PDB file path
    S_vib: float = 0.0              # kcal/mol/K
    S_vib_status: str = "model_scale_heuristic"
    eigenvalue_to_omega: float = 1.0
    calibration_label: str = "model-scale"
    calibration_provenance: str = (
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied"
    )
    delta_S_vib: float = 0.0       # kcal/mol/K
    delta_F_vib: float = 0.0       # kcal/mol
    n_modes: int = 0
    n_residues: int = 0
    temperature: float = 300.0      # K
    full_flexibility: bool = True
    bfactors: List[float] = field(default_factory=list)
    delta_bfactors: List[float] = field(default_factory=list)
    per_residue_svib: List[float] = field(default_factory=list)
    per_residue_delta_svib: List[float] = field(default_factory=list)
    eigenvalue_diffs: List[EigenvalueDiff] = field(default_factory=list)
    composition: dict = field(default_factory=dict)
    version: str = ""


@dataclass
class FlexPopulationResult:
    """Collection of FlexModeResults from a tENCoM run."""
    tool: str = ""
    version: str = ""
    temperature: float = 300.0
    full_flexibility: bool = True
    modes: List[FlexModeResult] = field(default_factory=list)

    @property
    def reference(self) -> Optional[FlexModeResult]:
        """Get the reference mode (mode_id=0)."""
        for m in self.modes:
            if m.mode_id == 0:
                return m
        return None

    @property
    def targets(self) -> List[FlexModeResult]:
        """Get all target modes (mode_id > 0)."""
        return [m for m in self.modes if m.mode_id > 0]

    def sorted_by_free_energy(self) -> List[FlexModeResult]:
        """Return targets sorted by delta_F_vib (ascending)."""
        return sorted(self.targets, key=lambda m: m.delta_F_vib)


# ── REMARK key aliases (match standardized KEY=VALUE format) ──────────────

_KEY_ALIASES = {
    "TENCOM_VERSION": "version",
    "TOOL": "tool",
    "MODE_ID": "mode_id",
    "MODE_TYPE": "mode_type",
    "SOURCE": "source",
    "S_VIB": "S_vib",
    "S_VIB_STATUS": "S_vib_status",
    "EIGENVALUE_TO_OMEGA": "eigenvalue_to_omega",
    "CALIBRATION_LABEL": "calibration_label",
    "CALIBRATION_PROVENANCE": "calibration_provenance",
    "DELTA_S_VIB": "delta_S_vib",
    "DELTA_F_VIB": "delta_F_vib",
    "N_MODES": "n_modes",
    "N_RESIDUES": "n_residues",
    "TEMPERATURE": "temperature",
    "FULL_FLEXIBILITY": "full_flexibility",
}


def parse_tencom_pdb(pdb_path: str) -> FlexModeResult:
    """Parse a tENCoM output PDB file and extract REMARK metadata.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        FlexModeResult with all parsed metadata.
    """
    result = FlexModeResult()
    path = Path(pdb_path)

    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    with open(path) as f:
        for line in f:
            if not line.startswith("REMARK"):
                continue
            content = line[6:].strip()

            # Parse KEY=VALUE pairs
            if "=" in content and not content.startswith("EIGENVALUE_DIFF") \
               and not content.startswith("BFACTORS") \
               and not content.startswith("DELTA_BFACTORS") \
               and not content.startswith("PER_RESIDUE") \
               and not content.startswith("COMPOSITION"):
                key, _, val = content.partition("=")
                key = key.strip()
                val = val.strip()
                alias = _KEY_ALIASES.get(key)
                if alias:
                    if alias in ("mode_id", "n_modes", "n_residues"):
                        setattr(result, alias, int(val))
                    elif alias in (
                        "S_vib", "delta_S_vib", "delta_F_vib", "temperature",
                        "eigenvalue_to_omega",
                    ):
                        setattr(result, alias, float(val))
                    elif alias == "full_flexibility":
                        setattr(result, alias, val.upper() == "ON")
                    else:
                        setattr(result, alias, val)

            # Parse EIGENVALUE_DIFF lines
            elif content.startswith("EIGENVALUE_DIFF "):
                parts = content.split()
                diff = EigenvalueDiff(mode=0, delta_eigenvalue=0.0)
                for part in parts[1:]:
                    if part.startswith("MODE="):
                        diff.mode = int(part[5:])
                    elif part.startswith("DELTA_EIG="):
                        diff.delta_eigenvalue = float(part[10:])
                    elif part.startswith("OVERLAP="):
                        diff.overlap = float(part[8:])
                result.eigenvalue_diffs.append(diff)

            # Parse BFACTORS
            elif content.startswith("BFACTORS "):
                result.bfactors = [float(x) for x in content[9:].split()]

            # Parse DELTA_BFACTORS
            elif content.startswith("DELTA_BFACTORS "):
                result.delta_bfactors = [float(x) for x in content[15:].split()]

            # Parse PER_RESIDUE_SVIB
            elif content.startswith("PER_RESIDUE_SVIB "):
                result.per_residue_svib = [float(x) for x in content[17:].split()]

            # Parse PER_RESIDUE_DELTA_SVIB
            elif content.startswith("PER_RESIDUE_DELTA_SVIB "):
                result.per_residue_delta_svib = [float(x) for x in content[23:].split()]

            # Parse COMPOSITION
            elif content.startswith("COMPOSITION"):
                parts = content.split()
                for part in parts[1:]:
                    if "=" in part:
                        k, v = part.split("=", 1)
                        result.composition[k.lower()] = int(v)

    return result


def parse_tencom_json(json_path: str) -> FlexPopulationResult:
    """Parse a tENCoM JSON results file.

    Args:
        json_path: Path to the JSON results file.

    Returns:
        FlexPopulationResult with all modes.
    """
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    with open(path) as f:
        data = json.load(f)

    pop = FlexPopulationResult(
        tool=data.get("tool", ""),
        version=data.get("version", ""),
        temperature=data.get("temperature", 300.0),
        full_flexibility=data.get("full_flexibility", True),
    )

    for m in data.get("modes", []):
        mode = FlexModeResult(
            mode_id=m.get("mode_id", 0),
            mode_type=m.get("type", ""),
            source=m.get("source", ""),
            S_vib=m.get("S_vib", 0.0),
            S_vib_status=m.get("S_vib_status", "model_scale_heuristic"),
            eigenvalue_to_omega=m.get("eigenvalue_to_omega", 1.0),
            calibration_label=m.get("calibration_label", "model-scale"),
            calibration_provenance=m.get(
                "calibration_provenance",
                "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied",
            ),
            delta_S_vib=m.get("delta_S_vib", 0.0),
            delta_F_vib=m.get("delta_F_vib", 0.0),
            n_modes=m.get("n_modes", 0),
            n_residues=m.get("n_residues", 0),
            temperature=pop.temperature,
            bfactors=m.get("bfactors", []),
            delta_bfactors=m.get("delta_bfactors", []),
            per_residue_svib=m.get("per_residue_svib", []),
            per_residue_delta_svib=m.get("per_residue_delta_svib", []),
        )
        comp = m.get("composition", {})
        if comp:
            mode.composition = comp

        for ed in m.get("eigenvalue_diffs", []):
            mode.eigenvalue_diffs.append(EigenvalueDiff(
                mode=ed.get("mode", 0),
                delta_eigenvalue=ed.get("delta", 0.0),
                overlap=ed.get("overlap"),
            ))

        pop.modes.append(mode)

    return pop
