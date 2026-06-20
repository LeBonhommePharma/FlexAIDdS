"""Energy matrix I/O and 256-type projection for FlexAID∆S.

Bridges between:
  - Legacy FlexAID ``.dat`` format (upper-triangle, scalar or density-function entries)
  - 256×256 binary blob format (SHNN magic header, float32)
  - NumPy 2D array representation

The 256-type encoding uses 8 bits: bits 0–5 = base type (64 classes extending
SYBYL — every SYBYL type gets a distinct base, no Solvent fallback),
bit 6 = charge polarity (0 = negative, 1 = positive), bit 7 = H-bond
donor/acceptor flag.  The ``project_to_40()`` method collapses a 256×256
matrix back to the 40-type SYBYL system used by the C++ Voronoi contact
function.

Example:
    >>> mat = EnergyMatrix.from_dat_file("WRK/nrg_mat_BEST_012012.dat")
    >>> print(mat.ntypes, mat.matrix.shape)
    10 (10, 10)
    >>> mat.to_dat_file("/tmp/roundtrip.dat")
"""

import json
import math
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]


# ── constants ────────────────────────────────────────────────────────────────

SHNN_MAGIC = b"SHNN"  # binary blob magic header
SHNN_VERSION = 1
MATRIX_256_SIZE = 256

# SYBYL atom type names (1-indexed in C++, 0-indexed here)
SYBYL_TYPE_NAMES: List[str] = [
    "C.1", "C.2", "C.3", "C.AR", "C.CAT",          # 1-5
    "N.1", "N.2", "N.3", "N.4", "N.AR", "N.AM", "N.PL3",  # 6-12
    "O.2", "O.3", "O.CO2", "O.AR",                  # 13-16
    "S.2", "S.3", "S.O", "S.O2", "S.AR",            # 17-21
    "P.3", "F", "CL", "BR", "I", "SE",              # 22-27
    "MG", "SR", "CU", "MN", "HG", "CD", "NI", "ZN", "CA", "FE", "CO.OH",  # 28-38
    "DUMMY", "SOLVENT",                               # 39-40
]

SYBYL_RADII: Dict[int, float] = {
    1: 1.88, 2: 1.72, 3: 1.88, 4: 1.76, 5: 1.88,
    6: 1.64, 7: 1.64, 8: 1.64, 9: 1.64, 10: 1.64, 11: 1.64, 12: 1.64,
    13: 1.42, 14: 1.46, 15: 1.46, 16: 1.46,
    17: 1.782, 18: 1.782, 19: 1.782, 20: 1.782, 21: 1.782,
    22: 1.871, 23: 1.560, 24: 1.735, 25: 1.978, 26: 2.094, 27: 1.9,
    28: 0.72, 29: 1.18, 30: 1.18, 31: 0.73, 32: 1.02, 33: 0.95,
    34: 0.69, 35: 0.74, 36: 1.00, 37: 0.61, 38: 0.65,
    39: 2.0, 40: 0.0,
}

# Roman numeral labels used in legacy 10-type .dat files
_ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
          "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]


# ── 256-type encoding helpers ────────────────────────────────────────────────

def encode_256_type(base_type: int, donor: bool = False,
                    acceptor: bool = False) -> int:
    # NOTE: Layout is [D:1][A:1][B:6] = donor × acceptor × 64 base classes.
    """Encode a 256-type atom index from components.

    Args:
        base_type:  Base type index (0–31).  Types 0–31 extend SYBYL's 40
                    classes into 32 canonical base types.
        donor:     True if atom can donate an H-bond.
        acceptor:  True if atom can accept an H-bond.

    Returns:
        8-bit integer (0–255).
    """
    base_type = max(0, min(63, int(base_type)))
    return (int(bool(donor)) << 7) | (int(bool(acceptor)) << 6) | base_type


def decode_256_type(code: int) -> Tuple[int, bool, bool]:
    """Decode a 256-type code into (base_type, donor, acceptor)."""
    code = int(code) & 0xFF
    base_type = code & 0x3F
    acceptor = bool((code >> 6) & 0x01)
    donor = bool((code >> 7) & 0x01)
    return base_type, donor, acceptor


# ── 256 → 40 SYBYL projection mapping ───────────────────────────────────────

# Maps each of 32 base types to its SYBYL parent (1-indexed, matching C++).
# Base types 0–31 are a superset of SYBYL types with context-aware splits
# (e.g., C_ar_hetadj, C_pi_bridging).  Types beyond the first 22 collapse
# onto nearby SYBYL parents.
_BASE_TO_SYBYL: List[int] = [
    # 0–4: carbon variants → SYBYL carbon types
    1,   # 0: C_sp  → C.1
    2,   # 1: C_sp2 → C.2
    3,   # 2: C_sp3 → C.3
    4,   # 3: C_ar  → C.AR
    5,   # 4: C_cat → C.CAT
    # 5–11: nitrogen variants
    6,   # 5: N_sp  → N.1
    7,   # 6: N_sp2 → N.2
    8,   # 7: N_sp3 → N.3
    9,   # 8: N_quat → N.4
    10,  # 9: N_ar  → N.AR
    11,  # 10: N_am → N.AM
    12,  # 11: N_pl3 → N.PL3
    # 12–15: oxygen variants
    13,  # 12: O_sp2 → O.2
    14,  # 13: O_sp3 → O.3
    15,  # 14: O_co2 → O.CO2
    16,  # 15: O_ar  → O.AR
    # 16–20: sulfur variants
    17,  # 16: S_sp2 → S.2
    18,  # 17: S_sp3 → S.3
    19,  # 18: S_oxide → S.O
    20,  # 19: S_dioxide → S.O2
    21,  # 20: S_ar → S.AR
    # 21: phosphorus
    22,  # 21: P_sp3 → P.3
    # 22–25: halogens
    23,  # 22: F → F
    24,  # 23: Cl → CL
    25,  # 24: Br → BR
    26,  # 25: I → I
    # 26–27: extended carbon (NATURaL-critical π-system refinements)
    4,   # 26: C_ar_hetadj → C.AR (aromatic C adjacent to heteroatom)
    2,   # 27: C_pi_bridging → C.2 (π-bridging carbon, tryptamine/indole)
    # 28–30: metals
    35,  # 28: Zn → ZN
    36,  # 29: Ca → CA
    37,  # 30: Fe → FE
    # 31: solvent
    40,  # 31: Solvent → SOLVENT
    # 32–41: extended types (no more Solvent fallback)
    27,  # 32: HAL_Se → SE
    28,  # 33: Metal_Mg → MG
    29,  # 34: Metal_Sr → SR
    30,  # 35: Metal_Cu → CU
    31,  # 36: Metal_Mn → MN
    32,  # 37: Metal_Hg → HG
    33,  # 38: Metal_Cd → CD
    34,  # 39: Metal_Ni → NI
    38,  # 40: Metal_Co → CO.OH
    39,  # 41: Dummy → DUMMY
    # 42–63: reserved (→ DUMMY)
    39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
    39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
]


def base_to_sybyl(base_type: int) -> int:
    """Map a base type (0–63) to its SYBYL parent (1–40)."""
    if 0 <= base_type < len(_BASE_TO_SYBYL):
        return _BASE_TO_SYBYL[base_type]
    return 39  # DUMMY fallback


def sybyl_to_base(sybyl_type: int) -> int:
    """Map a SYBYL type (1–40) to the canonical base type (0–63).

    Every SYBYL type maps to a distinct base type — no Solvent fallback.
    Context-aware refinements (C_ar_hetadj, C_pi_bridging) are NOT applied
    here — they require structural context.  Use ``refine_base_type()`` for that.
    """
    # Explicit mapping — every SYBYL type gets its own base type
    _sybyl_to_base_table = {
         1: 0,   2: 1,   3: 2,   4: 3,   5: 4,      # C types
         6: 5,   7: 6,   8: 7,   9: 8,  10: 9,       # N types
        11: 10, 12: 11,
        13: 12, 14: 13, 15: 14, 16: 15,               # O types
        17: 16, 18: 17, 19: 18, 20: 19, 21: 20,       # S types
        22: 21,                                         # P.3
        23: 22, 24: 23, 25: 24, 26: 25,               # halogens
        27: 32, 28: 33, 29: 34, 30: 35,               # SE, MG, SR, CU
        31: 36, 32: 37, 33: 38, 34: 39,               # MN, HG, CD, NI
        35: 28, 36: 29, 37: 30,                        # ZN, CA, FE
        38: 40, 39: 41, 40: 31,                        # CO.OH, DUMMY, SOL
    }
    return _sybyl_to_base_table.get(sybyl_type, 41)  # Dummy fallback


# ── Boltzmann constant (kcal/mol/K) ─────────────────────────────────────────

_kB_kcal: float = 0.0019872041


# ── contact table ───────────────────────────────────────────────────────────

@dataclass
class ContactTable:
    """Observed atom-type contact counts from a set of protein-ligand complexes.

    Used as input for knowledge-based potential derivation (Sippl reference
    ratio method).
    """

    ntypes: int
    counts: "np.ndarray"
    type_totals: "np.ndarray"
    n_structures: int
    distance_cutoff: float

    def save(self, path: str) -> None:
        """Serialize to JSON."""
        data = {
            "ntypes": self.ntypes,
            "counts": self.counts.tolist(),
            "type_totals": self.type_totals.tolist(),
            "n_structures": self.n_structures,
            "distance_cutoff": self.distance_cutoff,
        }
        with open(path, "w") as fh:
            json.dump(data, fh)

    @classmethod
    def load(cls, path: str) -> "ContactTable":
        """Deserialize from JSON."""
        with open(path) as fh:
            data = json.load(fh)
        return cls(
            ntypes=data["ntypes"],
            counts=np.array(data["counts"]),
            type_totals=np.array(data["type_totals"]),
            n_structures=data["n_structures"],
            distance_cutoff=data["distance_cutoff"],
        )


# ── knowledge-based trainer ─────────────────────────────────────────────────

class KnowledgeBasedTrainer:
    """Derive a knowledge-based potential from observed contact statistics.

    Implements the Sippl reference-ratio method: for each atom-type pair
    (i, j), the potential is ``-kBT * ln(N_obs / N_ref)`` where *N_ref* is
    the expected count under a uniform mixing assumption.
    """

    def __init__(self, ntypes: int, temperature: float = 300.0,
                 pseudocount: float = 1.0):
        self.ntypes = ntypes
        self.temperature = temperature
        self.pseudocount = pseudocount
        self._counts: Optional["np.ndarray"] = None
        self._type_totals: Optional["np.ndarray"] = None

    def add_contact_table(self, table: ContactTable) -> None:
        """Accumulate contacts from a :class:`ContactTable`."""
        if self._counts is None:
            self._counts = table.counts.astype(np.float64).copy()
            self._type_totals = table.type_totals.astype(np.float64).copy()
        else:
            self._counts += table.counts
            self._type_totals += table.type_totals

    def derive_potential(self) -> "EnergyMatrix":
        """Return an :class:`EnergyMatrix` with the derived potential."""
        if self._counts is None:
            raise RuntimeError("No contact tables added")
        total = self._type_totals.sum()
        freq = ((self._type_totals + self.pseudocount)
                / (total + self.pseudocount * self.ntypes))
        ref = np.outer(freq, freq)
        ref *= self._counts.sum()
        obs = self._counts + self.pseudocount
        ref_safe = np.where(ref > 0, ref, 1.0)
        ratio = obs / ref_safe
        kBT = _kB_kcal * self.temperature
        matrix = -kBT * np.log(ratio)
        matrix = (matrix + matrix.T) / 2.0
        return EnergyMatrix(self.ntypes, matrix)


# ── density function entry ───────────────────────────────────────────────────

@dataclass
class DensityPoint:
    """Single (x, y) point in a piecewise-linear density function."""
    x: float
    y: float


@dataclass
class MatrixEntry:
    """Energy matrix entry for a pair of atom types.

    For scalar entries (``is_scalar=True``), ``scalar_value`` holds the
    single energy value.  For density-function entries, ``density_points``
    stores the piecewise-linear (x, y) curve.
    """
    type1: int
    type2: int
    is_scalar: bool = True
    scalar_value: float = 0.0
    density_points: List[DensityPoint] = field(default_factory=list)

    def evaluate(self, relative_area: float) -> float:
        """Evaluate the energy for a given normalised contact area.

        Mirrors ``get_yval()`` in ``vcfunction.cpp``.
        """
        if self.is_scalar:
            return self.scalar_value

        if not self.density_points:
            return 0.0

        pts = self.density_points
        # Find bracketing interval
        if relative_area < pts[0].x:
            return 0.0  # no left-bound data
        for k in range(len(pts) - 1):
            if relative_area <= pts[k + 1].x:
                # linear interpolation
                dx = pts[k + 1].x - pts[k].x
                if abs(dx) < 1e-12:
                    return pts[k].y
                frac = (relative_area - pts[k].x) / dx
                return pts[k].y + frac * (pts[k + 1].y - pts[k].y)
        # beyond last point
        return pts[-1].y


# ── EnergyMatrix class ───────────────────────────────────────────────────────

class EnergyMatrix:
    """Read, write, and manipulate FlexAID energy matrices.

    Supports both legacy N×N scalar matrices and 256×256 binary blobs.

    Attributes:
        ntypes:  Number of atom types (e.g. 10, 40, or 256).
        matrix:  NumPy float64 array of shape ``(ntypes, ntypes)``.
                 Symmetric: ``matrix[i,j] == matrix[j,i]``.
        entries: Optional dict of ``(i,j) → MatrixEntry`` for density-function
                 entries.  Only populated when loaded from a ``.dat`` file that
                 contains density functions.
    """

    def __init__(self, ntypes: int,
                 matrix: "np.ndarray" = None,
                 entries: Optional[Dict[Tuple[int, int], MatrixEntry]] = None):
        if np is None:
            raise RuntimeError("NumPy is required for EnergyMatrix")
        self.ntypes = int(ntypes)
        if matrix is None:
            self.matrix = np.zeros((self.ntypes, self.ntypes), dtype=np.float64)
        else:
            self.matrix = np.asarray(matrix, dtype=np.float64)
        self.entries = entries or {}

    # ── legacy .dat I/O ──────────────────────────────────────────────────

    @classmethod
    def from_dat_file(cls, path: str) -> "EnergyMatrix":
        """Parse a legacy FlexAID ``.dat`` energy matrix file.

        Format: upper-triangle, N*(N+1)/2 lines.  Each line is either a
        single scalar or space-separated x,y pairs (density function).
        Lines may have a ``TYPE-TYPE = value`` prefix.

        Matches the C++ parser in ``read_emat.cpp``.
        """
        path = str(path)
        with open(path) as fh:
            raw_lines = fh.readlines()

        # Combine continuation lines (lines that don't contain a newline-
        # terminated complete entry are appended to the previous line).
        lines: List[str] = []
        current = ""
        for raw in raw_lines:
            current += raw
            if "\n" in current:
                lines.append(current)
                current = ""
        if current.strip():
            lines.append(current)

        nlines = len(lines)
        # Solve N*(N+1)/2 = nlines  →  N = (-1 + sqrt(1 + 8*nlines)) / 2
        discriminant = 1.0 + 8.0 * nlines
        ntypes = int((-1.0 + math.sqrt(discriminant)) / 2.0 + 0.001)

        if ntypes * (ntypes + 1) // 2 != nlines:
            raise ValueError(
                f"Invalid line count {nlines} in energy matrix file {path}; "
                f"expected N*(N+1)/2 lines for some integer N"
            )

        matrix = np.zeros((ntypes, ntypes), dtype=np.float64)
        entries: Dict[Tuple[int, int], MatrixEntry] = {}
        line_idx = 0

        for i in range(ntypes):
            for j in range(i, ntypes):
                line = lines[line_idx]
                line_idx += 1

                # Strip TYPE-TYPE prefix if present
                if "=" in line:
                    line = line[line.index("=") + 1:]
                tokens = line.split()

                entry = MatrixEntry(type1=i, type2=j)

                if len(tokens) == 1:
                    # Scalar entry
                    val = float(tokens[0])
                    entry.is_scalar = True
                    entry.scalar_value = val
                    matrix[i, j] = val
                    matrix[j, i] = val
                elif len(tokens) % 2 == 0:
                    # Density function (x,y pairs)
                    entry.is_scalar = False
                    for k in range(0, len(tokens), 2):
                        entry.density_points.append(
                            DensityPoint(x=float(tokens[k]),
                                         y=float(tokens[k + 1])))
                    # Store the mean y-value as the matrix scalar approximation
                    mean_y = sum(p.y for p in entry.density_points) / len(entry.density_points)
                    matrix[i, j] = mean_y
                    matrix[j, i] = mean_y
                else:
                    raise ValueError(
                        f"Invalid token count {len(tokens)} for pair ({i},{j}) "
                        f"in {path}"
                    )

                entries[(i, j)] = entry
                if i != j:
                    entries[(j, i)] = MatrixEntry(
                        type1=j, type2=i,
                        is_scalar=entry.is_scalar,
                        scalar_value=entry.scalar_value,
                        density_points=entry.density_points,
                    )

        return cls(ntypes, matrix, entries)

    def to_dat_file(self, path: str, labels: Optional[List[str]] = None) -> None:
        """Write legacy ``.dat`` format consumable by C++ ``read_emat.cpp``.

        Args:
            path:   Output file path.
            labels: Optional list of type labels (e.g. Roman numerals).
                    If None, Roman numerals are used for ntypes <= 20,
                    otherwise numeric indices.
        """
        if labels is None:
            if self.ntypes <= len(_ROMAN):
                labels = _ROMAN[:self.ntypes]
            else:
                labels = [str(i + 1) for i in range(self.ntypes)]

        with open(path, "w") as fh:
            for i in range(self.ntypes):
                for j in range(i, self.ntypes):
                    key = (i, j)
                    label = f"{labels[i]}-{labels[j]}"
                    # Right-align label to 10 chars
                    label_padded = label.rjust(10)

                    if key in self.entries and not self.entries[key].is_scalar:
                        # Density function
                        pts = self.entries[key].density_points
                        vals = " ".join(
                            f"{p.x:.4f} {p.y:.4f}" for p in pts)
                        fh.write(f"{label_padded} = {vals}\n")
                    else:
                        val = self.matrix[i, j]
                        fh.write(f"{label_padded} = {val:10.4f}\n")

    # ── 256×256 binary I/O ───────────────────────────────────────────────

    @classmethod
    def from_binary(cls, path: str) -> "EnergyMatrix":
        """Load a 256×256 matrix from binary blob (SHNN magic header).

        Binary format:
          - 4 bytes: ``SHNN`` magic
          - 4 bytes: uint32 version
          - 4 bytes: uint32 matrix dimension (256)
          - 256*256*4 bytes: float32 matrix data (row-major)
        """
        path = str(path)
        with open(path, "rb") as fh:
            magic = fh.read(4)
            if magic != SHNN_MAGIC:
                raise ValueError(
                    f"Invalid magic header: expected {SHNN_MAGIC!r}, got {magic!r}")
            version, dim = struct.unpack("<II", fh.read(8))
            if dim != MATRIX_256_SIZE:
                raise ValueError(
                    f"Expected {MATRIX_256_SIZE}×{MATRIX_256_SIZE} matrix, "
                    f"got {dim}×{dim}")
            data = fh.read(dim * dim * 4)
            if len(data) != dim * dim * 4:
                raise ValueError("Truncated binary matrix data")

        flat = np.frombuffer(data, dtype=np.float32).astype(np.float64)
        matrix = flat.reshape((dim, dim))
        return cls(dim, matrix)

    def to_binary(self, path: str) -> None:
        """Write 256×256 binary blob with SHNN magic header."""
        if self.ntypes != MATRIX_256_SIZE:
            raise ValueError(
                f"Binary format requires {MATRIX_256_SIZE}×{MATRIX_256_SIZE} "
                f"matrix, got {self.ntypes}×{self.ntypes}")
        with open(path, "wb") as fh:
            fh.write(SHNN_MAGIC)
            fh.write(struct.pack("<II", SHNN_VERSION, self.ntypes))
            fh.write(self.matrix.astype(np.float32).tobytes())

    # ── projection ───────────────────────────────────────────────────────

    def project_to_40(self) -> "EnergyMatrix":
        """Project a 256×256 matrix to 40×40 SYBYL via base_to_sybyl_parent().

        For each SYBYL type pair (s1, s2), averages over all 256-type codes
        whose base type maps to s1 and s2 respectively.

        Returns:
            New ``EnergyMatrix`` with ``ntypes=40``.
        """
        if self.ntypes != MATRIX_256_SIZE:
            raise ValueError("project_to_40 requires a 256×256 matrix")

        out = np.zeros((40, 40), dtype=np.float64)
        counts = np.zeros((40, 40), dtype=np.float64)

        for code_i in range(MATRIX_256_SIZE):
            base_i, _, _ = decode_256_type(code_i)
            sybyl_i = base_to_sybyl(base_i) - 1  # 0-indexed

            for code_j in range(MATRIX_256_SIZE):
                base_j, _, _ = decode_256_type(code_j)
                sybyl_j = base_to_sybyl(base_j) - 1

                if 0 <= sybyl_i < 40 and 0 <= sybyl_j < 40:
                    out[sybyl_i, sybyl_j] += self.matrix[code_i, code_j]
                    counts[sybyl_i, sybyl_j] += 1.0

        # Average where counts > 0
        mask = counts > 0
        out[mask] /= counts[mask]

        return EnergyMatrix(40, out)

    # ── lookup ───────────────────────────────────────────────────────────

    def lookup(self, type_i: int, type_j: int) -> float:
        """O(1) matrix lookup.

        Args:
            type_i: First atom type index (0-indexed).
            type_j: Second atom type index (0-indexed).

        Returns:
            Interaction energy value.
        """
        return float(self.matrix[type_i, type_j])

    def evaluate(self, type_i: int, type_j: int,
                 relative_area: float = 0.0) -> float:
        """Evaluate interaction energy, using density function if available.

        For scalar entries, returns the scalar value.
        For density-function entries, interpolates at ``relative_area``.
        """
        key = (type_i, type_j)
        if key in self.entries:
            return self.entries[key].evaluate(relative_area)
        return self.lookup(type_i, type_j)

    # ── utilities ────────────────────────────────────────────────────────

    @property
    def is_symmetric(self) -> bool:
        """Check if matrix is symmetric within floating-point tolerance."""
        return bool(np.allclose(self.matrix, self.matrix.T, atol=1e-6))

    def symmetrise(self) -> None:
        """Force symmetry: M = (M + M^T) / 2."""
        self.matrix = (self.matrix + self.matrix.T) / 2.0

    def __repr__(self) -> str:
        sym = "symmetric" if self.is_symmetric else "asymmetric"
        return f"<EnergyMatrix ntypes={self.ntypes} {sym}>"


# ── convenience functions ────────────────────────────────────────────────────

def parse_dat_file(path: str) -> Tuple[int, "np.ndarray"]:
    """Parse legacy ``.dat`` file and return ``(ntypes, matrix)``."""
    mat = EnergyMatrix.from_dat_file(path)
    return mat.ntypes, mat.matrix


def write_dat_file(path: str, ntypes: int, matrix: "np.ndarray",
                   labels: Optional[List[str]] = None) -> None:
    """Write a NumPy matrix as a legacy ``.dat`` file."""
    mat = EnergyMatrix(ntypes, matrix)
    mat.to_dat_file(path, labels=labels)
