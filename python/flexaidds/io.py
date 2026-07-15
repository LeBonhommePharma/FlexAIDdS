"""PDB I/O and REMARK-header parsing for FlexAID∆S docking output files.

FlexAID∆S writes CF/contact-function scores and (when available) ensemble
thermodynamic ledger fields (Helmholtz F estimate, entropy, heat capacity,
etc.) as ``REMARK`` records in each output PDB file.  CF is a scoring proxy,
not free energy.  This module provides low-level parsers that extract those
values and infer binding-mode and pose-rank identifiers from both REMARK
content and file names.

Also provides general-purpose PDB, MOL2, and FlexAID config file readers
and writers.

Public API:
    - :func:`parse_remark_map` – parse all REMARK lines into a key→value dict.
    - :func:`infer_mode_id` – determine which binding mode a PDB file belongs to.
    - :func:`infer_pose_rank` – determine the pose rank within its mode.
    - :func:`parse_pose_result` – full parse of one PDB file into a
      :class:`~flexaidds.models.PoseResult`.
    - :func:`read_pdb` – parse a PDB file into a PDBStructure.
    - :func:`write_pdb` – write a PDBStructure to a PDB file.
    - :func:`read_flexaid_config` – parse a FlexAID .inp config file.
    - :func:`write_flexaid_config` – write a FlexAID config dictionary.
    - :func:`read_sphere_pdb` – parse a FlexAID sphere PDB file.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
from dataclasses import dataclass, field

from .models import PoseResult

_NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _filename_token_pattern(token: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?<![A-Za-z0-9]){token}[_-]?(\d+)(?![A-Za-z0-9])",
        re.IGNORECASE,
    )


_FILE_MODE_PATTERNS = [
    re.compile(r"binding[_-]?mode[_-]?(\d+)", re.IGNORECASE),
    re.compile(r"\bmode[_-]?(\d+)(?!\d)", re.IGNORECASE),
    re.compile(r"cluster[_-]?(\d+)", re.IGNORECASE),
    re.compile(r"\bbm[_-]?(\d+)(?!\d)", re.IGNORECASE),
]
_FILE_POSE_PATTERNS = [
    _filename_token_pattern("pose"),
    _filename_token_pattern("conformer"),
    _filename_token_pattern("model"),
]


def _normalize_key(raw: str) -> str:
    """Normalise a raw REMARK key to a canonical snake_case identifier.

    Steps:
    1. Strip whitespace and convert to lowercase.
    2. Replace runs of non-alphanumeric characters with ``_``.
    3. Apply a fixed alias table to unify variant spellings
       (e.g. ``dg`` → ``free_energy``, ``cv`` → ``heat_capacity``).

    Args:
        raw: The raw key string extracted from a REMARK line.

    Returns:
        Normalised key string (may differ from the input).
    """
    key = raw.strip().lower()
    key = re.sub(r"[^a-z0-9]+", "_", key)
    key = re.sub(r"_+", "_", key).strip("_")
    aliases = {
        "bindingmode": "binding_mode",
        "binding_mode_id": "binding_mode",
        "modeid": "binding_mode",
        "clusterid": "cluster_id",
        "cluster": "cluster_id",
        "cf_app": "cf_app",
        "cfapp": "cf_app",
        "app_cf": "cf_app",
        "average_cf": "average_cf",
        "delta_g": "free_energy",
        "dg": "free_energy",
        "freeenergy": "free_energy",
        "enthalpy_like": "enthalpy",
        "energy_std": "std_energy",
        "sigma_energy": "std_energy",
        "cv": "heat_capacity",
        "temp": "temperature",
        "rmsd": "rmsd_raw",
        "rmsd_unsym": "rmsd_raw",
        "rmsd_symmetric": "rmsd_sym",
        "sym_rmsd": "rmsd_sym",
    }
    return aliases.get(key, key)


def _coerce_value(raw: str) -> Any:
    """Convert a raw REMARK value string to an appropriate Python type.

    Conversion priority:
    1. Numeric string → ``float`` (or ``int`` if the value is a whole number).
    2. ``"true"``/``"false"`` (case-insensitive) → ``bool``.
    3. Everything else → stripped ``str``.

    Args:
        raw: The raw value string from a REMARK line.

    Returns:
        Coerced Python value.
    """
    value = raw.strip().strip(";,")
    if not value:
        return value
    if _NUMERIC_RE.match(value):
        number = float(value)
        return int(number) if number.is_integer() else number
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    return value


# Classic FlexAID multi-field REMARK (single line with several Key:value pairs)
_FLEXAID_COMPOUND_RE = re.compile(
    r"Binding Mode\s*:\s*(\d+).*?"
    r"Best CF in Binding Mode\s*:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    r"(?:.*?Binding Mode Frequency\s*:\s*(\d+))?",
    re.IGNORECASE,
)

# Production CF clustering line (cluster.cpp):
#   Cluster 4: Rank (top):4 Average CF:-153.81 Frequency:12
_CLUSTER_RANK_RE = re.compile(
    r"Cluster\s+(\d+)\s*:\s*"
    r"(?:Rank\s*\(top\)\s*:\s*(\d+))?"
    r"(?:.*?Average\s+CF\s*:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?))?"
    r"(?:.*?Frequency\s*:\s*(\d+))?",
    re.IGNORECASE,
)

# DensityPeak style: Cluster:N Best CF in Cluster:X ...
_CLUSTER_COLON_RE = re.compile(
    r"Cluster\s*:\s*(\d+).*?"
    r"(?:Best CF in Cluster\s*:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?))?"
    r"(?:.*?Cluster Frequency\s*:\s*(\d+))?",
    re.IGNORECASE,
)

# Secondary Key:value tokens on compound REMARK lines
_COMPOUND_TOKEN_RE = re.compile(
    r"([A-Za-z][A-Za-z0-9 ._\-/()]*?)\s*:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|[A-Za-z0-9_./+-]+)"
)

# Classic numeric RMSD prose lines written by the engine
_RMSD_PROSE_RE = re.compile(
    r"^([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s+RMSD to ref\.\s*structure",
    re.IGNORECASE,
)

# Strip PDB remark class numbers: "REMARK  99  Free Energy: -8.5"
_REMARK_NUMBER_RE = re.compile(r"^\d{1,3}\s+")


def _strip_remark_payload(line: str) -> str:
    """Return payload after REMARK, stripping optional remark class number."""
    if not line.startswith("REMARK"):
        return ""
    payload = line[6:].strip()
    payload = _REMARK_NUMBER_RE.sub("", payload, count=1).strip()
    return payload


def _harvest_compound_tokens(payload: str, remarks: Dict[str, Any]) -> None:
    """Harvest all Key:value tokens from a compound REMARK payload."""
    for tok in _COMPOUND_TOKEN_RE.finditer(payload):
        key = _normalize_key(tok.group(1))
        # Skip tokens whose "key" is pure noise (e.g. bare numbers mis-matched)
        if not key or key[0].isdigit():
            continue
        if key not in remarks:
            remarks[key] = _coerce_value(tok.group(2))


def _merge_flexaid_compound_remark(payload: str, remarks: Dict[str, Any]) -> bool:
    """Extract fields from a classic FlexAID multi-field REMARK line.

    Example payload::

        Binding Mode:1 Best CF in Binding Mode:-5.5 Binding Mode Frequency:3

    Returns True if the compound pattern matched.
    """
    m = _FLEXAID_COMPOUND_RE.search(payload)
    if not m:
        return False
    if "binding_mode" not in remarks:
        remarks["binding_mode"] = int(m.group(1))
    if "cf" not in remarks:
        remarks["cf"] = float(m.group(2))
    if m.group(3) is not None and "frequency" not in remarks:
        remarks["frequency"] = int(m.group(3))
    _harvest_compound_tokens(payload, remarks)
    # Alias best_cf for mode aggregation convenience
    if "best_cf" not in remarks and "cf" in remarks:
        remarks["best_cf"] = remarks["cf"]
    return True


def _merge_cluster_remark(payload: str, remarks: Dict[str, Any]) -> bool:
    """Extract fields from production Cluster REMARK lines.

    Supports::

        Cluster 4: Rank (top):4 Average CF:-153.81 Frequency:12
        Cluster:4 Best CF in Cluster:-5.5 Cluster Frequency:3

    Cluster index maps to ``binding_mode`` / ``cluster_id``.  ``Rank (top)``
    is the GA chromosome index of the representative (stored as
    ``rank_top``), *not* the within-mode pose rank.
    """
    m = _CLUSTER_RANK_RE.search(payload)
    if m:
        cluster_id = int(m.group(1))
        if "binding_mode" not in remarks:
            remarks["binding_mode"] = cluster_id
        if "cluster_id" not in remarks:
            remarks["cluster_id"] = cluster_id
        if m.group(2) is not None and "rank_top" not in remarks:
            remarks["rank_top"] = int(m.group(2))
        if m.group(3) is not None and "average_cf" not in remarks:
            remarks["average_cf"] = float(m.group(3))
        if m.group(4) is not None and "frequency" not in remarks:
            remarks["frequency"] = int(m.group(4))
        # Representative pose for this cluster file
        if "pose_rank" not in remarks:
            remarks["pose_rank"] = 1
        _harvest_compound_tokens(payload, remarks)
        return True

    m2 = _CLUSTER_COLON_RE.search(payload)
    if m2:
        cluster_id = int(m2.group(1))
        if "binding_mode" not in remarks:
            remarks["binding_mode"] = cluster_id
        if "cluster_id" not in remarks:
            remarks["cluster_id"] = cluster_id
        if m2.group(2) is not None and "cf" not in remarks:
            remarks["cf"] = float(m2.group(2))
        if m2.group(3) is not None and "frequency" not in remarks:
            remarks["frequency"] = int(m2.group(3))
        if "pose_rank" not in remarks:
            remarks["pose_rank"] = 1
        _harvest_compound_tokens(payload, remarks)
        if "best_cf" not in remarks and "cf" in remarks:
            remarks["best_cf"] = remarks["cf"]
        return True
    return False


def parse_remark_map(lines: Iterable[str]) -> Dict[str, Any]:
    """Parse ``REMARK`` records from PDB lines into a key→value dictionary.

    Supported REMARK formats:
    - ``REMARK Key = value`` (``=`` or ``:`` delimiter)
    - ``REMARK Key value`` (space delimiter, first-seen wins)
    - ``REMARK  99  Key: value`` (remark class number stripped)
    - Classic FlexAID multi-field lines:
      ``REMARK Binding Mode:N Best CF in Binding Mode:X Binding Mode Frequency:Y``
    - Production CF cluster lines:
      ``REMARK Cluster N: Rank (top):R Average CF:X Frequency:F``
    - Compound thermo lines:
      ``REMARK Free Energy:-8.5 Enthalpy:-10.0 Entropy:0.01``
    - Classic RMSD prose lines (``rmsd_raw`` / ``rmsd_sym``)

    Keys are normalised via :func:`_normalize_key`; values are coerced via
    :func:`_coerce_value`.  Non-REMARK lines are silently ignored.

    Args:
        lines: Iterable of PDB file lines (strings).

    Returns:
        Dictionary mapping normalised keys to coerced values.  If the same
        key appears multiple times the first occurrence wins (except equals/
        colon form which overwrites for simple single-key lines).
    """
    remarks: Dict[str, Any] = {}
    for line in lines:
        if not line.startswith("REMARK"):
            continue
        payload = _strip_remark_payload(line)
        if not payload:
            continue

        # Prefer compound FlexAID / Cluster patterns
        if "Best CF in Binding Mode" in payload or "Binding Mode:" in payload:
            if _merge_flexaid_compound_remark(payload, remarks):
                continue

        if payload.lower().startswith("cluster"):
            if _merge_cluster_remark(payload, remarks):
                continue

        # Classic RMSD prose (must run before space-delimited Key value)
        rmsd_m = _RMSD_PROSE_RE.match(payload)
        if rmsd_m:
            val = float(rmsd_m.group(1))
            pl = payload.lower()
            if "no symmetry" in pl:
                if "rmsd_raw" not in remarks:
                    remarks["rmsd_raw"] = val
            elif "symmetry" in pl:
                if "rmsd_sym" not in remarks:
                    remarks["rmsd_sym"] = val
            else:
                if "rmsd_raw" not in remarks:
                    remarks["rmsd_raw"] = val
            continue

        # Multi-token compound thermo lines: harvest every Key:value
        # e.g. "Free Energy:-8.5 Enthalpy:-10.0 Entropy:0.01"
        colon_tokens = list(_COMPOUND_TOKEN_RE.finditer(payload))
        if len(colon_tokens) >= 2:
            _harvest_compound_tokens(payload, remarks)
            continue

        match = re.match(r"([A-Za-z][A-Za-z0-9 ._\-/]*)\s*(?:=|:)\s*(.+)", payload)
        if match:
            key = _normalize_key(match.group(1))
            # Equals/colon form: last occurrence wins (historical behaviour)
            remarks[key] = _coerce_value(match.group(2))
            continue

        match = re.match(r"([A-Za-z][A-Za-z0-9_\-/]*)\s+(.+)", payload)
        if match:
            key = _normalize_key(match.group(1))
            # Space form: first occurrence wins
            if key not in remarks:
                remarks[key] = _coerce_value(match.group(2))
    return remarks


def infer_mode_id(path: Path, remarks: Dict[str, Any]) -> int:
    """Determine the binding-mode ID for a docking output PDB file.

    Lookup order:
    1. REMARK keys ``binding_mode``, ``mode``, ``cluster_id``, ``cluster``.
    2. Regex patterns applied to the file stem
       (e.g. ``mode3``, ``binding_mode_2``, ``cluster1``).
    3. Fallback: ``1``.

    Args:
        path: Path to the PDB file (used to extract filename patterns).
        remarks: Parsed REMARK map from :func:`parse_remark_map`.

    Returns:
        Integer binding-mode identifier (>= 1).
    """
    for key in ("binding_mode", "mode", "cluster_id", "cluster"):
        value = remarks.get(key)
        if isinstance(value, int):
            return value
        if isinstance(value, float) and value.is_integer():
            return int(value)
    name = path.stem
    for pattern in _FILE_MODE_PATTERNS:
        match = pattern.search(name)
        if match:
            return int(match.group(1))
    return 1


def infer_pose_rank(path: Path, remarks: Dict[str, Any]) -> int:
    """Determine the pose rank within a binding mode for a docking output PDB.

    Lookup order:
    1. REMARK keys ``pose_rank``, ``rank``, ``pose``, ``model``.
    2. Regex patterns applied to the file stem
       (e.g. ``pose5``, ``conformer_2``, ``model3``).
    3. Fallback: ``1``.

    Args:
        path: Path to the PDB file (used to extract filename patterns).
        remarks: Parsed REMARK map from :func:`parse_remark_map`.

    Returns:
        Integer pose rank (>= 1).
    """
    for key in ("pose_rank", "rank", "pose", "model"):
        value = remarks.get(key)
        if isinstance(value, int):
            return value
        if isinstance(value, float) and value.is_integer():
            return int(value)
    name = path.stem
    for pattern in _FILE_POSE_PATTERNS:
        match = pattern.search(name)
        if match:
            return int(match.group(1))
    return 1


def _first_float(remarks: Dict[str, Any], *keys: str) -> Optional[float]:
    """Return the first numeric value found for any of the given REMARK *keys*.

    Args:
        remarks: Parsed REMARK dictionary.
        *keys: One or more candidate key names, tried in order.

    Returns:
        The value cast to ``float``, or ``None`` if no key is present.
    """
    for key in keys:
        value = remarks.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    return None


def parse_pose_result(path: Path) -> PoseResult:
    """Parse one FlexAID∆S output PDB file into a :class:`~flexaidds.models.PoseResult`.

    Reads the file, extracts all REMARK fields, infers the binding-mode ID
    and pose rank, then returns a fully populated (immutable) ``PoseResult``.

    Args:
        path: Absolute or relative path to the PDB output file.

    Returns:
        :class:`~flexaidds.models.PoseResult` with all available fields
        populated from the PDB REMARK section.
    """
    text = path.read_text(encoding="utf-8", errors="ignore")
    remarks = parse_remark_map(text.splitlines())
    mode_id = infer_mode_id(path, remarks)
    pose_rank = infer_pose_rank(path, remarks)

    cf_app = _first_float(remarks, "cf_app", "cfapp", "apparent_cf")
    cf = _first_float(
        remarks,
        "cf",
        "best_cf",
        "best_cf_in_binding_mode",
        "complementarity_function",
        "score",
    )
    if cf is None:
        cf = cf_app

    return PoseResult(
        path=path,
        mode_id=mode_id,
        pose_rank=pose_rank,
        cf=cf,
        cf_app=cf_app,
        rmsd_raw=_first_float(remarks, "rmsd_raw", "rmsd", "rmsd_unsym"),
        rmsd_sym=_first_float(remarks, "rmsd_sym", "rmsd_symmetric", "sym_rmsd"),
        free_energy=_first_float(remarks, "free_energy", "f"),
        enthalpy=_first_float(remarks, "enthalpy", "h", "mean_energy"),
        entropy=_first_float(remarks, "entropy", "s"),
        heat_capacity=_first_float(remarks, "heat_capacity", "cv"),
        std_energy=_first_float(remarks, "std_energy", "sigma_energy"),
        temperature=_first_float(remarks, "temperature", "temp"),
        remarks=remarks,
    )


# ─── General-purpose PDB I/O ─────────────────────────────────────────────────


@dataclass
class Atom:
    """Parsed ATOM/HETATM record from a PDB file."""
    serial: int
    name: str
    altloc: str
    resname: str
    chainid: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    record: str  # 'ATOM' or 'HETATM'

    def __repr__(self) -> str:
        return (
            f"<Atom {self.serial} {self.name} {self.resname} "
            f"{self.chainid}{self.resseq}>"
        )

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


# ── Ion classification ────────────────────────────────────────────────────────
_ION_RESNAMES: frozenset = frozenset({
    "MG", "ZN", "CA", "NA", "K", "FE", "FE2", "FE3",
    "CU", "CU1", "CU2", "MN", "CO", "NI", "CL", "BR",
    "IOD", "LI", "CD", "HG", "PB",
})


def is_ion(atom: "Atom") -> bool:
    """Return True if *atom* is a receptor-bound metal ion or halide.

    Requires ``atom.record == 'HETATM'`` and a recognised single-atom residue
    name (MG, ZN, CA, NA, K, FE, FE2, FE3, CU, CU1, CU2, MN, CO, NI, CL,
    BR, IOD, LI, CD, HG, PB).  Protein Cα atoms (ATOM record, residue ALA
    etc.) are *not* matched even though the atom name may be "CA".
    """
    return atom.record == "HETATM" and atom.resname.strip() in _ION_RESNAMES


@dataclass
class PDBStructure:
    """Parsed PDB structure."""
    atoms: List[Atom] = field(default_factory=list)
    remarks: List[str] = field(default_factory=list)
    title: str = ""

    @property
    def coords(self) -> np.ndarray:
        """All atomic coordinates as (N, 3) array."""
        return np.array([[a.x, a.y, a.z] for a in self.atoms])

    def select_chain(self, chain_id: str) -> "PDBStructure":
        """Return new structure with only atoms from specified chain."""
        filtered = [a for a in self.atoms if a.chainid == chain_id]
        return PDBStructure(atoms=filtered, title=self.title)

    def select_residue(self, resseq: int, chain_id: str = "") -> List[Atom]:
        """Return atoms from a specific residue."""
        return [
            a for a in self.atoms
            if a.resseq == resseq and (not chain_id or a.chainid == chain_id)
        ]

    def get_chain_ids(self) -> List[str]:
        """Return list of unique chain IDs."""
        seen = []
        for a in self.atoms:
            if a.chainid not in seen:
                seen.append(a.chainid)
        return seen

    def __repr__(self) -> str:
        chains = self.get_chain_ids()
        chain_str = ",".join(chains) if chains else "none"
        return (
            f"<PDBStructure atoms={len(self.atoms)} "
            f"chains=[{chain_str}]>"
        )


def read_pdb(path: str) -> PDBStructure:
    """Parse a PDB file into a PDBStructure.

    Args:
        path: Path to PDB file.

    Returns:
        Parsed PDBStructure.
    """
    p = Path(path)
    structure = PDBStructure()

    with open(p) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "TITLE":
                structure.title = line[10:].strip()
            elif rec == "REMARK":
                structure.remarks.append(line[6:].strip())
            elif rec in ("ATOM", "HETATM"):
                try:
                    atom = Atom(
                        record=rec,
                        serial=int(line[6:11]),
                        name=line[12:16].strip(),
                        altloc=line[16].strip(),
                        resname=line[17:20].strip(),
                        chainid=line[21].strip(),
                        resseq=int(line[22:26]),
                        icode=line[26].strip(),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        occupancy=float(line[54:60]) if len(line) > 54 else 1.0,
                        bfactor=float(line[60:66]) if len(line) > 60 else 0.0,
                        element=line[76:78].strip() if len(line) > 76 else "",
                    )
                    structure.atoms.append(atom)
                except (ValueError, IndexError):
                    continue

    return structure


def write_pdb(structure: PDBStructure, path: str) -> None:
    """Write a PDBStructure to a PDB file.

    Args:
        structure: Structure to write.
        path: Output file path.
    """
    with open(path, "w") as fh:
        if structure.title:
            fh.write(f"TITLE     {structure.title}\n")
        for remark in structure.remarks:
            fh.write(f"REMARK{remark}\n")
        for atom in structure.atoms:
            fh.write(
                f"{atom.record:<6}{atom.serial:5d} {atom.name:<4}{atom.altloc:1}"
                f"{atom.resname:<3} {atom.chainid:1}{atom.resseq:4d}{atom.icode:1}   "
                f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}"
                f"          {atom.element:>2}\n"
            )
        fh.write("END\n")


def read_flexaid_config(path: str) -> Dict[str, object]:
    """Parse a FlexAID .inp configuration file into a dictionary.

    The FlexAID config format uses 6-character keywords followed by values.
    Keywords that may appear multiple times (OPTIMZ, FLEXSC) are collected
    into lists.

    Args:
        path: Path to .inp config file.

    Returns:
        Dictionary mapping keyword -> value.
    """
    _list_keys = {"OPTIMZ", "FLEXSC"}
    config: Dict[str, object] = {k: [] for k in _list_keys}

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if len(line) < 6:
                continue
            keyword = line[:6].strip()
            value = line[7:].strip() if len(line) > 7 else ""

            if keyword in _list_keys:
                config[keyword].append(value)  # type: ignore[attr-defined]
            else:
                # Attempt numeric conversion
                try:
                    config[keyword] = int(value)
                except ValueError:
                    try:
                        config[keyword] = float(value)
                    except ValueError:
                        config[keyword] = value if value else True

    return config


def write_flexaid_config(config: Dict[str, object], path: str) -> None:
    """Write a FlexAID configuration dictionary to a .inp file.

    Args:
        config: Dictionary of keyword -> value mappings.
        path: Output file path.
    """
    _list_keys = {"OPTIMZ", "FLEXSC"}
    with open(path, "w") as fh:
        for keyword, value in config.items():
            if keyword in _list_keys and isinstance(value, list):
                for item in value:
                    fh.write(f"{keyword:<6} {item}\n")
            elif value is True:
                fh.write(f"{keyword}\n")
            elif value is not False and value is not None:
                fh.write(f"{keyword:<6} {value}\n")


@dataclass
class SphereRecord:
    """Sphere record from a FlexAID sphere/cleft PDB file."""
    x: float
    y: float
    z: float
    radius: float
    cleft_id: int = 1

    def __repr__(self) -> str:
        return (
            f"<SphereRecord ({self.x:.2f}, {self.y:.2f}, {self.z:.2f}) "
            f"r={self.radius:.2f} cleft={self.cleft_id}>"
        )

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


def read_sphere_pdb(path: str) -> List[SphereRecord]:
    """Parse a FlexAID sphere PDB (B-factor = radius, resSeq = cleft_id).

    Args:
        path: Path to sphere PDB file.

    Returns:
        List of SphereRecord objects.
    """
    spheres: List[SphereRecord] = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                radius = float(line[60:66]) if len(line) > 60 else 1.0
                cleft_id = int(line[22:26]) if len(line) > 22 else 1
                spheres.append(SphereRecord(x=x, y=y, z=z,
                                            radius=radius, cleft_id=cleft_id))
            except (ValueError, IndexError):
                continue
    return spheres
