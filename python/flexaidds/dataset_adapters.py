"""Multi-dataset adapter system for continuous energy matrix training.

Normalizes heterogeneous experimental binding affinity datasets into the
common Complex/ContactPair representation used by train_256x256.py.

Supported datasets:
    - PDBbind (general/refined/core) — Tier 2-3
    - ITC-187 (calorimetric ΔH/TΔS/ΔG) — Tier 1 (gold standard)
    - Binding MOAD — Tier 3
    - BindingDB — Tier 4
    - ChEMBL — Tier 4
    - DUD-E / DEKOIS 2.0 — Validation only (no training)

Dependencies: numpy, scipy (BSD-licensed).  No GPL dependencies.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import re
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import numpy as np

from .energy_matrix import ContactTable, EnergyMatrix, encode_256_type
from .train_256x256 import (
    CONTACT_CUTOFF,
    TEMPERATURE,
    Complex,
    Atom,
    ContactPair,
    enumerate_contacts,
    parse_pdb_atoms,
    parse_mol2_atoms,
    parse_pdbbind_index,
    kB_kcal,
    _MOL2_TYPE_MAP,
    _PDB_ELEMENT_MAP,
)

logger = logging.getLogger(__name__)

# ── constants ────────────────────────────────────────────────────────────────

R_KCAL = kB_kcal  # Gas constant in kcal/(mol·K)
LN10 = math.log(10)

# Reliability tiers and default weights
TIER_WEIGHTS = {
    1: 1.00,   # ITC, PDBbind core
    2: 0.80,   # PDBbind refined
    3: 0.50,   # PDBbind general, MOAD
    4: 0.25,   # BindingDB, ChEMBL
}


# ── affinity normalization ───────────────────────────────────────────────────

def normalize_affinity(
    value: float,
    unit: str,
    temperature: float = TEMPERATURE,
) -> float:
    """Convert heterogeneous affinity measures to ΔG (kcal/mol).

    Args:
        value: Measured value (in the native unit).
        unit: One of 'Kd', 'Ki', 'IC50' (molar), 'pKd', 'pKi', 'pIC50',
              'deltaG' (kcal/mol).
        temperature: Temperature in Kelvin.

    Returns:
        ΔG in kcal/mol (negative = favorable binding).
    """
    unit_lower = unit.lower().replace("_", "").replace("-", "")
    RT = R_KCAL * temperature

    if unit_lower == "deltag":
        return value

    if unit_lower in ("kd", "ki"):
        if value <= 0:
            raise ValueError(f"Non-positive Kd/Ki: {value}")
        return RT * math.log(value)

    if unit_lower == "ic50":
        # Cheng-Prusoff approximation: Ki ≈ IC50 / 2
        if value <= 0:
            raise ValueError(f"Non-positive IC50: {value}")
        ki_approx = value / 2.0
        return RT * math.log(ki_approx)

    if unit_lower in ("pkd", "pki"):
        # pKd = -log10(Kd) → Kd = 10^(-pKd) → ΔG = RT·ln(10^(-pKd))
        return -RT * LN10 * value

    if unit_lower == "pic50":
        pkd_approx = value - math.log10(2)  # Cheng-Prusoff
        return -RT * LN10 * pkd_approx

    raise ValueError(f"Unknown affinity unit: {unit!r}")


# ── dataset metadata ─────────────────────────────────────────────────────────

@dataclass
class DatasetMetadata:
    """Provenance metadata for a dataset."""
    name: str
    version: str = ""
    source_url: str = ""
    n_complexes: int = 0
    reliability_tier: int = 3
    weight: float = 0.5
    date_loaded: str = ""
    checksum: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ── base adapter ─────────────────────────────────────────────────────────────

class DatasetAdapter(ABC):
    """Base class for dataset-specific loaders."""

    @abstractmethod
    def name(self) -> str:
        """Short identifier (e.g. 'pdbbind_refined')."""
        ...

    @abstractmethod
    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        """Load complexes from *data_dir* with contacts enumerated at *cutoff*."""
        ...

    @abstractmethod
    def metadata(self) -> DatasetMetadata:
        """Return provenance metadata for this dataset."""
        ...

    @property
    def is_training_dataset(self) -> bool:
        """True if this dataset provides binding affinities for training."""
        return True


# ── PDBbind adapter ──────────────────────────────────────────────────────────

class PDBbindAdapter(DatasetAdapter):
    """Adapter for PDBbind general/refined/core sets.

    Expected directory layout::

        data_dir/
            INDEX/INDEX_general_PL_data.*
            1a0q/1a0q_protein.pdb
            1a0q/1a0q_ligand.mol2
    """

    def __init__(
        self,
        subset: str = "refined",
        version: str = "v2020",
        tier: int = 2,
        weight: Optional[float] = None,
    ):
        self._subset = subset
        self._version = version
        self._tier = tier
        self._weight = weight if weight is not None else TIER_WEIGHTS.get(tier, 0.5)
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return f"pdbbind_{self._subset}"

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        pdbbind = Path(data_dir)
        if not pdbbind.is_dir():
            raise FileNotFoundError(f"PDBbind directory not found: {data_dir}")

        # Find index file
        affinities: Dict[str, float] = {}
        for pattern in ["INDEX/*PL_data*", "INDEX/*general*", "index/*"]:
            candidates = list(pdbbind.glob(pattern))
            if candidates:
                affinities = parse_pdbbind_index(str(candidates[0]))
                break

        complexes = []
        for subdir in sorted(pdbbind.iterdir()):
            if not subdir.is_dir():
                continue
            code = subdir.name.lower()
            pdb_files = list(subdir.glob("*_protein.pdb")) + list(
                subdir.glob("*_pocket.pdb")
            )
            mol2_files = list(subdir.glob("*_ligand.mol2"))
            if not pdb_files or not mol2_files:
                continue
            try:
                prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
                lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
            except Exception as exc:
                logger.debug("Skipping %s: %s", code, exc)
                continue
            if not contacts:
                continue

            pkd = affinities.get(code, 0.0)
            dg = -R_KCAL * TEMPERATURE * LN10 * pkd if pkd != 0.0 else 0.0

            complexes.append(Complex(
                pdb_code=code,
                protein_atoms=prot_atoms,
                ligand_atoms=lig_atoms,
                contacts=contacts,
                pKd=pkd,
                deltaG=dg,
            ))

        logger.info("PDBbind %s: loaded %d complexes", self._subset, len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(),
            version=self._version,
            n_complexes=len(complexes),
            reliability_tier=self._tier,
            weight=self._weight,
        )
        return complexes

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=self._tier, weight=self._weight,
            )
        return self._meta


# ── ITC-187 adapter ──────────────────────────────────────────────────────────

class ITC187Adapter(DatasetAdapter):
    """Adapter for the ITC-187 calorimetric dataset.

    Expected layout::

        data_dir/
            itc_index.csv          # pdb_code,deltaH,TdeltaS,deltaG,Kd,...
            structures/
                1a0q/
                    1a0q_protein.pdb
                    1a0q_ligand.mol2
    """

    def __init__(self, version: str = "v1.0"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "itc_187"

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        base = Path(data_dir)
        index_path = base / "itc_index.csv"
        struct_dir = base / "structures"

        if not index_path.is_file():
            raise FileNotFoundError(f"ITC index not found: {index_path}")
        if not struct_dir.is_dir():
            raise FileNotFoundError(f"ITC structures not found: {struct_dir}")

        entries = self._parse_index(str(index_path))
        complexes = []

        for entry in entries:
            code = entry["pdb_code"].lower()
            subdir = struct_dir / code
            pdb_files = list(subdir.glob("*_protein.pdb")) + list(
                subdir.glob("*.pdb")
            )
            mol2_files = list(subdir.glob("*_ligand.mol2")) + list(
                subdir.glob("*.mol2")
            )
            if not pdb_files or not mol2_files:
                logger.debug("ITC: missing structure files for %s", code)
                continue

            try:
                prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
                lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
            except Exception as exc:
                logger.debug("ITC: skipping %s: %s", code, exc)
                continue
            if not contacts:
                continue

            dg = entry["deltaG"]
            pkd = -dg / (R_KCAL * TEMPERATURE * LN10) if dg != 0 else 0.0

            complexes.append(Complex(
                pdb_code=code,
                protein_atoms=prot_atoms,
                ligand_atoms=lig_atoms,
                contacts=contacts,
                pKd=pkd,
                deltaG=dg,
            ))

        logger.info("ITC-187: loaded %d complexes", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes),
            reliability_tier=1, weight=1.0,
        )
        return complexes

    @staticmethod
    def _parse_index(index_path: str) -> List[Dict[str, Any]]:
        """Parse ITC index CSV: pdb_code,deltaH,TdeltaS,deltaG,..."""
        entries = []
        with open(index_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    entry = {
                        "pdb_code": row["pdb_code"].strip(),
                        "deltaH": float(row.get("deltaH", 0)),
                        "TdeltaS": float(row.get("TdeltaS", 0)),
                        "deltaG": float(row["deltaG"]),
                    }
                    entries.append(entry)
                except (KeyError, ValueError) as exc:
                    logger.debug("ITC index parse error: %s", exc)
                    continue
        return entries

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=1, weight=1.0,
            )
        return self._meta


# ── Binding MOAD adapter ─────────────────────────────────────────────────────

class BindingMOADAdapter(DatasetAdapter):
    """Adapter for Binding MOAD (Mother of All Databases).

    Expected layout::

        data_dir/
            every.csv              # PDB_ID,Ligand_ID,Binding_Data,...
            structures/
                1a0q/
                    1a0q_protein.pdb
                    1a0q_ligand.mol2
    """

    _AFFINITY_RE = re.compile(
        r"(Kd|Ki|IC50)\s*[=<>~]*\s*([\d.eE+-]+)\s*(nM|uM|mM|M|pM)",
        re.IGNORECASE,
    )

    _UNIT_SCALE = {
        "pm": 1e-12, "nm": 1e-9, "um": 1e-6, "mm": 1e-3, "m": 1.0,
    }

    def __init__(self, version: str = "2024"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "binding_moad"

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        base = Path(data_dir)
        index_path = base / "every.csv"
        struct_dir = base / "structures"

        if not index_path.is_file():
            raise FileNotFoundError(f"MOAD index not found: {index_path}")

        entries = self._parse_index(str(index_path))
        complexes = []

        for entry in entries:
            code = entry["pdb_code"].lower()
            subdir = struct_dir / code
            if not subdir.is_dir():
                continue

            pdb_files = list(subdir.glob("*_protein.pdb")) + list(
                subdir.glob("*.pdb")
            )
            mol2_files = list(subdir.glob("*_ligand.mol2")) + list(
                subdir.glob("*.mol2")
            )
            if not pdb_files or not mol2_files:
                continue

            try:
                prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
                lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
            except Exception as exc:
                logger.debug("MOAD: skipping %s: %s", code, exc)
                continue
            if not contacts:
                continue

            dg = entry["deltaG"]
            pkd = -dg / (R_KCAL * TEMPERATURE * LN10) if dg != 0 else 0.0

            complexes.append(Complex(
                pdb_code=code,
                protein_atoms=prot_atoms,
                ligand_atoms=lig_atoms,
                contacts=contacts,
                pKd=pkd,
                deltaG=dg,
            ))

        logger.info("Binding MOAD: loaded %d complexes", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes),
            reliability_tier=3, weight=TIER_WEIGHTS[3],
        )
        return complexes

    def _parse_index(self, index_path: str) -> List[Dict[str, Any]]:
        """Parse Binding MOAD every.csv."""
        entries = []
        with open(index_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                pdb_code = row.get("PDB_ID", row.get("pdb_id", "")).strip()
                binding_data = row.get("Binding_Data", row.get("binding_data", ""))
                if not pdb_code or not binding_data:
                    continue

                parsed = self._parse_binding_string(binding_data)
                if parsed is None:
                    continue

                measure, value_molar = parsed
                try:
                    dg = normalize_affinity(value_molar, measure)
                except ValueError:
                    continue

                entries.append({"pdb_code": pdb_code, "deltaG": dg})
        return entries

    def _parse_binding_string(
        self, text: str,
    ) -> Optional[Tuple[str, float]]:
        """Extract (measure, value_in_molar) from MOAD binding data string."""
        m = self._AFFINITY_RE.search(text)
        if not m:
            return None
        measure = m.group(1)  # Kd, Ki, or IC50
        try:
            value = float(m.group(2))
        except ValueError:
            return None
        unit = m.group(3).lower()
        scale = self._UNIT_SCALE.get(unit, 1e-9)
        return measure, value * scale

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=3, weight=TIER_WEIGHTS[3],
            )
        return self._meta


# ── BindingDB adapter ────────────────────────────────────────────────────────

class BindingDBAdapter(DatasetAdapter):
    """Adapter for BindingDB TSV export.

    Filters for entries with PDB IDs and converts Ki/Kd/IC50 (nM) to ΔG.

    Expected layout::

        data_dir/
            BindingDB_All.tsv      # Tab-separated export
            structures/            # Optional: PDB + MOL2 per PDB code
                1a0q/1a0q_protein.pdb
                1a0q/1a0q_ligand.mol2
    """

    # Column names in BindingDB TSV export
    _PDB_COL = "PDB ID(s) of Target Chain"
    _KD_COL = "Kd (nM)"
    _KI_COL = "Ki (nM)"
    _IC50_COL = "IC50 (nM)"

    def __init__(self, version: str = "2024"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "bindingdb"

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        base = Path(data_dir)
        tsv_path = base / "BindingDB_All.tsv"
        struct_dir = base / "structures"

        if not tsv_path.is_file():
            # Also try a general .tsv file
            tsv_files = list(base.glob("*.tsv"))
            if tsv_files:
                tsv_path = tsv_files[0]
            else:
                raise FileNotFoundError(f"BindingDB TSV not found in {data_dir}")

        entries = self._parse_tsv(str(tsv_path))
        complexes = []

        for entry in entries:
            code = entry["pdb_code"].lower()
            subdir = struct_dir / code
            if not subdir.is_dir():
                continue

            pdb_files = list(subdir.glob("*_protein.pdb")) + list(
                subdir.glob("*.pdb")
            )
            mol2_files = list(subdir.glob("*_ligand.mol2")) + list(
                subdir.glob("*.mol2")
            )
            if not pdb_files or not mol2_files:
                continue

            try:
                prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
                lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
            except Exception as exc:
                logger.debug("BindingDB: skipping %s: %s", code, exc)
                continue
            if not contacts:
                continue

            complexes.append(Complex(
                pdb_code=code,
                protein_atoms=prot_atoms,
                ligand_atoms=lig_atoms,
                contacts=contacts,
                pKd=entry.get("pKd", 0.0),
                deltaG=entry["deltaG"],
            ))

        logger.info("BindingDB: loaded %d complexes", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes),
            reliability_tier=4, weight=TIER_WEIGHTS[4],
        )
        return complexes

    def _parse_tsv(self, tsv_path: str) -> List[Dict[str, Any]]:
        """Parse BindingDB TSV, preferring Kd > Ki > IC50."""
        entries = []
        seen_codes: set = set()

        with open(tsv_path, newline="", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                pdb_raw = row.get(self._PDB_COL, "")
                if not pdb_raw:
                    continue
                # May contain comma-separated PDB codes
                for pdb_code in pdb_raw.split(","):
                    pdb_code = pdb_code.strip()
                    if len(pdb_code) != 4:
                        continue
                    if pdb_code.lower() in seen_codes:
                        continue

                    # Prefer Kd > Ki > IC50
                    dg = None
                    for col, unit in [
                        (self._KD_COL, "Kd"),
                        (self._KI_COL, "Ki"),
                        (self._IC50_COL, "IC50"),
                    ]:
                        raw = row.get(col, "").strip()
                        if not raw:
                            continue
                        try:
                            val_nm = float(raw)
                            if val_nm <= 0:
                                continue
                            val_molar = val_nm * 1e-9
                            dg = normalize_affinity(val_molar, unit)
                            break
                        except (ValueError, TypeError):
                            continue

                    if dg is None:
                        continue

                    pkd = -dg / (R_KCAL * TEMPERATURE * LN10) if dg != 0 else 0.0
                    seen_codes.add(pdb_code.lower())
                    entries.append({
                        "pdb_code": pdb_code,
                        "deltaG": dg,
                        "pKd": pkd,
                    })

        return entries

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=4, weight=TIER_WEIGHTS[4],
            )
        return self._meta


# ── ChEMBL adapter ───────────────────────────────────────────────────────────

class ChEMBLAdapter(DatasetAdapter):
    """Adapter for ChEMBL bioactivity CSV export.

    Filters for binding assays (assay_type='B') with PDB cross-references
    and pchembl_value (quality-filtered).

    Expected layout::

        data_dir/
            chembl_activities.csv  # ChEMBL export with pchembl_value, PDB columns
            structures/
                1a0q/1a0q_protein.pdb
                1a0q/1a0q_ligand.mol2
    """

    def __init__(self, version: str = "34"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "chembl"

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        base = Path(data_dir)
        csv_path = base / "chembl_activities.csv"
        struct_dir = base / "structures"

        if not csv_path.is_file():
            csv_files = list(base.glob("*.csv"))
            if csv_files:
                csv_path = csv_files[0]
            else:
                raise FileNotFoundError(f"ChEMBL CSV not found in {data_dir}")

        entries = self._parse_csv(str(csv_path))
        complexes = []

        for entry in entries:
            code = entry["pdb_code"].lower()
            subdir = struct_dir / code
            if not subdir.is_dir():
                continue

            pdb_files = list(subdir.glob("*_protein.pdb")) + list(
                subdir.glob("*.pdb")
            )
            mol2_files = list(subdir.glob("*_ligand.mol2")) + list(
                subdir.glob("*.mol2")
            )
            if not pdb_files or not mol2_files:
                continue

            try:
                prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
                lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
            except Exception as exc:
                logger.debug("ChEMBL: skipping %s: %s", code, exc)
                continue
            if not contacts:
                continue

            complexes.append(Complex(
                pdb_code=code,
                protein_atoms=prot_atoms,
                ligand_atoms=lig_atoms,
                contacts=contacts,
                pKd=entry.get("pKd", 0.0),
                deltaG=entry["deltaG"],
            ))

        logger.info("ChEMBL: loaded %d complexes", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes),
            reliability_tier=4, weight=TIER_WEIGHTS[4],
        )
        return complexes

    def _parse_csv(self, csv_path: str) -> List[Dict[str, Any]]:
        """Parse ChEMBL activities CSV, filtering for quality."""
        entries = []
        seen: set = set()

        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                # Quality filters
                assay_type = row.get("assay_type", "").strip()
                if assay_type and assay_type != "B":
                    continue

                pchembl = row.get("pchembl_value", "").strip()
                if not pchembl:
                    continue

                pdb_code = row.get("pdb_id", row.get("PDB_ID", "")).strip()
                if not pdb_code or len(pdb_code) != 4:
                    continue

                if pdb_code.lower() in seen:
                    continue

                try:
                    pchembl_val = float(pchembl)
                except ValueError:
                    continue

                # pchembl_value ≈ -log10(molar) → treat as pKd
                dg = normalize_affinity(pchembl_val, "pKd")
                pkd = pchembl_val

                seen.add(pdb_code.lower())
                entries.append({
                    "pdb_code": pdb_code,
                    "deltaG": dg,
                    "pKd": pkd,
                })

        return entries

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=4, weight=TIER_WEIGHTS[4],
            )
        return self._meta


# ── DUD-E adapter (validation only) ─────────────────────────────────────────

class DUDEAdapter(DatasetAdapter):
    """Adapter for DUD-E virtual screening benchmark.

    Validation only — provides active/decoy labels, no binding affinities.

    Expected layout::

        data_dir/
            targets/
                aa2ar/
                    receptor.pdb
                    actives_final.mol2
                    decoys_final.mol2
    """

    def __init__(self, version: str = "1.0"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "dude"

    @property
    def is_training_dataset(self) -> bool:
        return False

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        """Load DUD-E targets. Complexes have deltaG=0 (no affinity data)."""
        base = Path(data_dir) / "targets"
        if not base.is_dir():
            base = Path(data_dir)

        complexes = []
        for target_dir in sorted(base.iterdir()):
            if not target_dir.is_dir():
                continue
            receptor_files = (
                list(target_dir.glob("receptor.pdb"))
                + list(target_dir.glob("*_receptor.pdb"))
            )
            active_files = list(target_dir.glob("actives*.mol2"))
            decoy_files = list(target_dir.glob("decoys*.mol2"))

            if not receptor_files:
                continue

            # Load a single representative complex per target
            for mol2_path in active_files[:1]:
                try:
                    prot = parse_pdb_atoms(str(receptor_files[0]))
                    lig = parse_mol2_atoms(str(mol2_path))
                    contacts = enumerate_contacts(prot, lig, cutoff)
                except Exception:
                    continue
                if contacts:
                    complexes.append(Complex(
                        pdb_code=target_dir.name,
                        protein_atoms=prot,
                        ligand_atoms=lig,
                        contacts=contacts,
                    ))

        logger.info("DUD-E: loaded %d targets", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes), reliability_tier=0, weight=0.0,
        )
        return complexes

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=0, weight=0.0,
            )
        return self._meta


# ── DEKOIS 2.0 adapter (validation only) ────────────────────────────────────

class DEKOIS2Adapter(DatasetAdapter):
    """Adapter for DEKOIS 2.0 virtual screening benchmark.

    Validation only — similar structure to DUD-E.
    """

    def __init__(self, version: str = "2.0"):
        self._version = version
        self._meta: Optional[DatasetMetadata] = None

    def name(self) -> str:
        return "dekois2"

    @property
    def is_training_dataset(self) -> bool:
        return False

    def load(self, data_dir: str, cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
        base = Path(data_dir) / "targets"
        if not base.is_dir():
            base = Path(data_dir)

        complexes = []
        for target_dir in sorted(base.iterdir()):
            if not target_dir.is_dir():
                continue
            pdb_files = list(target_dir.glob("*.pdb"))
            mol2_files = list(target_dir.glob("actives*.mol2")) + list(
                target_dir.glob("*_active*.mol2")
            )
            if not pdb_files or not mol2_files:
                continue
            try:
                prot = parse_pdb_atoms(str(pdb_files[0]))
                lig = parse_mol2_atoms(str(mol2_files[0]))
                contacts = enumerate_contacts(prot, lig, cutoff)
            except Exception:
                continue
            if contacts:
                complexes.append(Complex(
                    pdb_code=target_dir.name,
                    protein_atoms=prot,
                    ligand_atoms=lig,
                    contacts=contacts,
                ))

        logger.info("DEKOIS 2.0: loaded %d targets", len(complexes))
        self._meta = DatasetMetadata(
            name=self.name(), version=self._version,
            n_complexes=len(complexes), reliability_tier=0, weight=0.0,
        )
        return complexes

    def metadata(self) -> DatasetMetadata:
        if self._meta is None:
            return DatasetMetadata(
                name=self.name(), version=self._version,
                reliability_tier=0, weight=0.0,
            )
        return self._meta


# ── contact table caching ────────────────────────────────────────────────────

def checksum_contact_table(table: ContactTable) -> str:
    """SHA-256 of the serialized contact table (first 16 hex chars)."""
    data = json.dumps({
        "ntypes": table.ntypes,
        "n_structures": table.n_structures,
        "distance_cutoff": table.distance_cutoff,
        "counts_sum": float(table.counts.sum()) if table.counts is not None else 0,
    }, sort_keys=True).encode()
    return f"sha256:{hashlib.sha256(data).hexdigest()[:16]}"


def complexes_to_contact_table(
    complexes: List[Complex],
    ntypes: int = 256,
    cutoff: float = CONTACT_CUTOFF,
) -> ContactTable:
    """Convert a list of Complex objects to a ContactTable."""
    counts = np.zeros((ntypes, ntypes), dtype=np.float64)
    type_totals = np.zeros(ntypes, dtype=np.float64)

    for cpx in complexes:
        for c in cpx.contacts:
            counts[c.type_a, c.type_b] += 1.0
            if c.type_a != c.type_b:
                counts[c.type_b, c.type_a] += 1.0
            type_totals[c.type_a] += 1.0
            type_totals[c.type_b] += 1.0

    return ContactTable(
        ntypes=ntypes,
        counts=counts,
        type_totals=type_totals,
        n_structures=len(complexes),
        distance_cutoff=cutoff,
    )


def get_or_build_contact_table(
    adapter: DatasetAdapter,
    data_dir: str,
    cache_dir: str,
    cutoff: float = CONTACT_CUTOFF,
    force_rebuild: bool = False,
) -> Tuple[ContactTable, List[Complex]]:
    """Load cached ContactTable or build from scratch.

    Returns (ContactTable, complexes) — complexes needed for ridge/LBFGS.
    """
    cache_path = Path(cache_dir) / f"{adapter.name()}.json"
    os.makedirs(cache_dir, exist_ok=True)

    complexes: List[Complex] = []

    if cache_path.exists() and not force_rebuild:
        logger.info("Loading cached contact table: %s", cache_path)
        table = ContactTable.load(str(cache_path))
        # Still need complexes for regression-based training
        complexes = adapter.load(data_dir, cutoff)
        return table, complexes

    complexes = adapter.load(data_dir, cutoff)
    table = complexes_to_contact_table(complexes, cutoff=cutoff)
    table.save(str(cache_path))
    logger.info("Saved contact table cache: %s", cache_path)
    return table, complexes


# ── adapter registry ─────────────────────────────────────────────────────────

ADAPTER_REGISTRY: Dict[str, type] = {
    "pdbbind_core": PDBbindAdapter,
    "pdbbind_refined": PDBbindAdapter,
    "pdbbind_general": PDBbindAdapter,
    "itc_187": ITC187Adapter,
    "binding_moad": BindingMOADAdapter,
    "bindingdb": BindingDBAdapter,
    "chembl": ChEMBLAdapter,
    "dude": DUDEAdapter,
    "dekois2": DEKOIS2Adapter,
}


def create_adapter(name: str, **kwargs: Any) -> DatasetAdapter:
    """Create a dataset adapter by name."""
    if name == "pdbbind_core":
        return PDBbindAdapter(subset="core", tier=1, **kwargs)
    if name == "pdbbind_refined":
        return PDBbindAdapter(subset="refined", tier=2, **kwargs)
    if name == "pdbbind_general":
        return PDBbindAdapter(subset="general", tier=3, **kwargs)
    if name in ADAPTER_REGISTRY:
        cls = ADAPTER_REGISTRY[name]
        return cls(**kwargs)
    raise ValueError(f"Unknown dataset adapter: {name!r}")
