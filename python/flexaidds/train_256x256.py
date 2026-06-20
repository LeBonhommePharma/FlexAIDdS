"""Offline training pipeline for the 256×256 soft contact matrix.

Parses PDBbind protein-ligand complexes (PDB + MOL2), assigns 256-types,
enumerates inter-molecular contacts via KD-tree, fits per-cell interaction
energies by ridge regression, then refines the full matrix via L-BFGS against
CASF-2016 Pearson r.  Validates via 256→40 projection against FlexAID's
existing SYBYL-based energy matrix.

Usage:
    python -m flexaidds.train_256x256 \\
        --pdbbind-dir /data/PDBbind/refined-set \\
        --output matrix_256x256.bin \\
        --validate-casf /data/CASF-2016 \\
        --project-40 legacy_40x40.dat

Dependencies: numpy, scipy (BSD-licensed).  No GPL dependencies.
"""

import argparse
import logging
import math
import os
import re
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from scipy.spatial import cKDTree
    from scipy.optimize import minimize
    from scipy.stats import pearsonr
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from .energy_matrix import (
    EnergyMatrix,
    SHNN_MAGIC,
    SHNN_VERSION,
    MATRIX_256_SIZE,
    encode_256_type,
    decode_256_type,
    base_to_sybyl,
    sybyl_to_base,
    SYBYL_TYPE_NAMES,
)

logger = logging.getLogger(__name__)

# ─── constants ───────────────────────────────────────────────────────────────

CONTACT_CUTOFF = 4.5    # Å — max distance for a contact pair
kB_kcal = 0.001987206   # kcal mol⁻¹ K⁻¹
TEMPERATURE = 298.15     # K
RIDGE_ALPHA = 1.0        # L2 regularisation strength


# ─── SYBYL MOL2 type → base type mapping ────────────────────────────────────

_MOL2_TYPE_MAP: Dict[str, int] = {
    "C.1": 0, "C.2": 1, "C.3": 2, "C.ar": 3, "C.cat": 4,
    "N.1": 5, "N.2": 6, "N.3": 7, "N.4": 8, "N.ar": 9, "N.am": 10, "N.pl3": 11,
    "O.2": 12, "O.3": 13, "O.co2": 14,
    "S.2": 16, "S.3": 17, "S.O": 18, "S.O2": 19,
    "P.3": 21,
    "F": 22, "Cl": 23, "Br": 24, "I": 25, "Se": 31,
    "Zn": 28, "Ca": 29, "Fe": 30,
    "Du": 31, "LP": 31,
}

_PDB_ELEMENT_MAP: Dict[str, int] = {
    "C": 2, "N": 7, "O": 13, "S": 17, "P": 21,
    "F": 22, "CL": 23, "BR": 24, "I": 25,
    "ZN": 28, "CA": 29, "FE": 30, "MG": 31,
}


# ─── data structures ────────────────────────────────────────────────────────

@dataclass
class Atom:
    """Parsed atom with 256-type assignment."""
    index: int
    name: str
    element: str
    x: float
    y: float
    z: float
    charge: float = 0.0
    base_type: int = 31
    type_256: int = 0

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclass
class ContactPair:
    """Observed contact between protein and ligand atoms."""
    type_a: int     # 256-type of protein atom
    type_b: int     # 256-type of ligand atom
    distance: float
    area: float = 1.0  # approximate; set to 1.0 for count-based training


@dataclass
class Complex:
    """Protein-ligand complex with experimental binding affinity."""
    pdb_code: str
    protein_atoms: List[Atom]
    ligand_atoms: List[Atom]
    contacts: List[ContactPair] = field(default_factory=list)
    pKd: float = 0.0       # -log10(Kd) experimental affinity
    deltaG: float = 0.0    # ΔG = RT ln Kd (kcal/mol)


# ─── parsers ─────────────────────────────────────────────────────────────────

def parse_pdb_atoms(pdb_path: str) -> List[Atom]:
    """Parse ATOM/HETATM records from a PDB file."""
    atoms = []
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                idx = int(line[6:11].strip())
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip()
                if not element:
                    element = name[0]
            except (ValueError, IndexError):
                continue

            base = _PDB_ELEMENT_MAP.get(element.upper(), 31)
            charge = 0.0
            donor = _is_hbond_donor_base(base, charge)
            acceptor = _is_hbond_acceptor_base(base, charge)
            t256 = encode_256_type(base, donor, acceptor)
            atoms.append(Atom(idx, name, element, x, y, z, charge, base, t256))
    return atoms


def parse_mol2_atoms(mol2_path: str) -> List[Atom]:
    """Parse @<TRIPOS>ATOM section from a MOL2 file."""
    atoms = []
    in_atom_section = False
    with open(mol2_path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom_section:
                break
            if not in_atom_section:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                idx = int(parts[0])
                name = parts[1]
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                sybyl_type = parts[5]
                charge = float(parts[8]) if len(parts) > 8 else 0.0
            except (ValueError, IndexError):
                continue

            base = _MOL2_TYPE_MAP.get(sybyl_type, 31)
            donor = _is_hbond_donor_base(base, charge)
            acceptor = _is_hbond_acceptor_base(base, charge)
            t256 = encode_256_type(base, donor, acceptor)
            element = sybyl_type.split(".")[0] if "." in sybyl_type else sybyl_type
            atoms.append(Atom(idx, name, element, x, y, z, charge, base, t256))
    return atoms


def _is_hbond_donor_base(base: int, charge: float) -> bool:
    """Mirror the v56 C++ heavy-atom donor role classifier."""
    return base in {7, 8, 10, 11, 13} or (base == 17 and charge > 0.3)


def _is_hbond_acceptor_base(base: int, charge: float) -> bool:
    """Mirror the v56 C++ heavy-atom acceptor role classifier."""
    if base in {5, 6, 7, 9, 11, 12, 13, 14, 15, 18, 19, 22}:
        return True
    if base in {16, 17, 20}:
        return charge < -0.3
    return False


def parse_pdbbind_index(index_path: str) -> Dict[str, float]:
    """Parse PDBbind INDEX/INDEX_general_PL_data.* for pdb_code → pKd."""
    affinity = {}
    with open(index_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            pdb_code = parts[0].strip().lower()
            try:
                pkd = float(parts[3])
            except ValueError:
                continue
            affinity[pdb_code] = pkd
    return affinity


# ─── contact enumeration ────────────────────────────────────────────────────

def enumerate_contacts(protein_atoms: List[Atom],
                       ligand_atoms: List[Atom],
                       cutoff: float = CONTACT_CUTOFF) -> List[ContactPair]:
    """Find protein-ligand contacts within cutoff using KD-tree."""
    if not HAS_SCIPY:
        return _enumerate_contacts_brute(protein_atoms, ligand_atoms, cutoff)

    prot_coords = np.array([[a.x, a.y, a.z] for a in protein_atoms])
    lig_coords = np.array([[a.x, a.y, a.z] for a in ligand_atoms])

    if len(prot_coords) == 0 or len(lig_coords) == 0:
        return []

    tree = cKDTree(prot_coords)
    contacts = []

    for j, lc in enumerate(lig_coords):
        neighbors = tree.query_ball_point(lc, cutoff)
        for i in neighbors:
            dist = np.linalg.norm(prot_coords[i] - lc)
            contacts.append(ContactPair(
                type_a=protein_atoms[i].type_256,
                type_b=ligand_atoms[j].type_256,
                distance=dist,
            ))
    return contacts


def _enumerate_contacts_brute(protein_atoms: List[Atom],
                               ligand_atoms: List[Atom],
                               cutoff: float) -> List[ContactPair]:
    """Brute-force O(N*M) fallback when scipy is unavailable."""
    contacts = []
    cutoff2 = cutoff * cutoff
    for pa in protein_atoms:
        for la in ligand_atoms:
            dx = pa.x - la.x
            dy = pa.y - la.y
            dz = pa.z - la.z
            d2 = dx*dx + dy*dy + dz*dz
            if d2 <= cutoff2:
                contacts.append(ContactPair(
                    type_a=pa.type_256,
                    type_b=la.type_256,
                    distance=math.sqrt(d2),
                ))
    return contacts


# ─── contact frequency matrix ───────────────────────────────────────────────

def build_contact_matrix(complexes: List[Complex]) -> np.ndarray:
    """Build 256×256 observed contact frequency matrix."""
    freq = np.zeros((256, 256), dtype=np.float64)
    for cpx in complexes:
        for c in cpx.contacts:
            freq[c.type_a, c.type_b] += 1.0
            if c.type_a != c.type_b:
                freq[c.type_b, c.type_a] += 1.0
    return freq


def build_reference_matrix(complexes: List[Complex]) -> np.ndarray:
    """Build volume-corrected random reference distribution.

    Reference state: product of marginal type frequencies × total contacts.
    This is the standard Sippl-style reference for knowledge-based potentials.
    """
    # Count type occurrences across all complexes
    type_count = np.zeros(256, dtype=np.float64)
    total_contacts = 0.0
    for cpx in complexes:
        for c in cpx.contacts:
            type_count[c.type_a] += 1.0
            type_count[c.type_b] += 1.0
            total_contacts += 1.0

    if total_contacts < 1.0:
        return np.ones((256, 256), dtype=np.float64)

    # Marginal probabilities
    p = type_count / type_count.sum()

    # Reference: f_ref[i][j] = total_contacts * p[i] * p[j]
    ref = total_contacts * np.outer(p, p)
    # Ensure no zeros (Laplace smoothing)
    ref = np.maximum(ref, 0.5)
    return ref


# ─── inverse Boltzmann (Sippl potential) ─────────────────────────────────────

def inverse_boltzmann(freq: np.ndarray, ref: np.ndarray,
                      temperature: float = TEMPERATURE) -> np.ndarray:
    """Compute knowledge-based potential via inverse Boltzmann.

    E(i,j) = -kT * ln(f_obs(i,j) / f_ref(i,j))
    """
    kT = kB_kcal * temperature
    ratio = np.maximum(freq, 0.5) / np.maximum(ref, 0.5)
    return -kT * np.log(ratio)


# ─── ridge regression per cell ───────────────────────────────────────────────

def ridge_fit(complexes: List[Complex],
              alpha: float = RIDGE_ALPHA) -> np.ndarray:
    """Fit 256×256 matrix values by per-cell ridge regression.

    For each cell (i,j), collects contact counts across all complexes
    and regresses against experimental ΔG.
    """
    n = len(complexes)
    if n == 0:
        return np.zeros((256, 256), dtype=np.float64)

    # Feature matrix: X[k, i*256+j] = contact count for pair (i,j) in complex k
    # This is very sparse, so we accumulate per-cell sums
    contact_sums = np.zeros((256, 256), dtype=np.float64)
    contact_sq = np.zeros((256, 256), dtype=np.float64)
    target_sums = np.zeros((256, 256), dtype=np.float64)
    cell_counts = np.zeros((256, 256), dtype=np.float64)

    for cpx in complexes:
        # Per-cell contact count for this complex
        cell_n = {}
        for c in cpx.contacts:
            key = (c.type_a, c.type_b)
            cell_n[key] = cell_n.get(key, 0) + 1

        for (i, j), count in cell_n.items():
            fc = float(count)
            contact_sums[i, j] += fc
            contact_sq[i, j] += fc * fc
            target_sums[i, j] += fc * cpx.deltaG
            cell_counts[i, j] += 1.0

    # Ridge solution per cell: w = (X^T X + α I)^-1 X^T y
    # For 1D per-cell regression: w = Σ(x_k * y_k) / (Σ(x_k²) + α)
    matrix = np.zeros((256, 256), dtype=np.float64)
    mask = cell_counts > 0
    denom = contact_sq[mask] + alpha
    matrix[mask] = target_sums[mask] / denom

    # Symmetrise
    matrix = (matrix + matrix.T) / 2.0
    return matrix


# ─── L-BFGS refinement against CASF Pearson r ───────────────────────────────

def lbfgs_refine(matrix: np.ndarray,
                 complexes: List[Complex],
                 max_iter: int = 200) -> np.ndarray:
    """Refine matrix via L-BFGS to maximise Pearson r against experimental ΔG.

    Objective: minimise -r(predicted_scores, experimental_deltaG)
    where predicted_score = Σ_contacts matrix[type_a, type_b]
    """
    if not HAS_SCIPY:
        logger.warning("scipy not available; skipping L-BFGS refinement")
        return matrix

    n_cpx = len(complexes)
    if n_cpx < 3:
        return matrix

    # Precompute contact feature vectors (sparse)
    exp_dg = np.array([cpx.deltaG for cpx in complexes])

    # Contact indices per complex
    contact_indices = []
    for cpx in complexes:
        indices = {}
        for c in cpx.contacts:
            key = c.type_a * 256 + c.type_b
            indices[key] = indices.get(key, 0) + 1
        contact_indices.append(indices)

    # Extract upper triangle for optimisation (symmetric matrix)
    tri_indices = []
    for i in range(256):
        for j in range(i, 256):
            tri_indices.append((i, j))
    n_params = len(tri_indices)

    # Initial params from upper triangle
    x0 = np.array([matrix[i, j] for i, j in tri_indices], dtype=np.float64)

    def predict(params):
        scores = np.zeros(n_cpx)
        for k, ci in enumerate(contact_indices):
            score = 0.0
            for idx, count in ci.items():
                i, j = divmod(idx, 256)
                # Map to upper triangle index
                ii, jj = min(i, j), max(i, j)
                tri_idx = ii * 256 - ii * (ii + 1) // 2 + jj
                if tri_idx < n_params:
                    score += params[tri_idx] * count
            scores[k] = score
        return scores

    def objective(params):
        scores = predict(params)
        if np.std(scores) < 1e-12:
            return 0.0
        r, _ = pearsonr(scores, exp_dg)
        return -r  # minimise negative correlation

    logger.info("L-BFGS refinement: %d parameters, %d complexes", n_params, n_cpx)
    result = minimize(objective, x0, method='L-BFGS-B',
                      options={'maxiter': max_iter})

    # Reconstruct full matrix from optimised upper triangle
    opt_matrix = np.zeros((256, 256), dtype=np.float64)
    for idx, (i, j) in enumerate(tri_indices):
        opt_matrix[i, j] = result.x[idx]
        opt_matrix[j, i] = result.x[idx]

    if result.success:
        logger.info("L-BFGS converged: final r = %.4f", -result.fun)
    else:
        logger.warning("L-BFGS did not converge: %s", result.message)

    return opt_matrix


# ─── CASF-2016 validation ───────────────────────────────────────────────────

def validate_casf(matrix: np.ndarray,
                  complexes: List[Complex]) -> Dict[str, float]:
    """Compute CASF-2016-style metrics: Pearson r, RMSE, Kendall tau."""
    if not complexes:
        return {"pearson_r": 0.0, "rmse": 0.0}

    predicted = []
    experimental = []
    for cpx in complexes:
        score = 0.0
        for c in cpx.contacts:
            score += matrix[c.type_a, c.type_b]
        predicted.append(score)
        experimental.append(cpx.deltaG)

    predicted = np.array(predicted)
    experimental = np.array(experimental)

    if len(predicted) < 3 or np.std(predicted) < 1e-12:
        return {"pearson_r": 0.0, "rmse": float(np.std(experimental)),
                "n_complexes": len(complexes)}

    r, p_value = pearsonr(predicted, experimental) if HAS_SCIPY else (0.0, 1.0)
    rmse = float(np.sqrt(np.mean((predicted - experimental) ** 2)))

    return {
        "pearson_r": float(r),
        "p_value": float(p_value),
        "rmse": rmse,
        "n_complexes": len(complexes),
    }


# ─── 256→40 projection validation ───────────────────────────────────────────

def validate_projection(matrix_256: np.ndarray,
                        reference_dat_path: str) -> Dict[str, float]:
    """Compare 256→40 projection against a reference FlexAID .dat matrix."""
    em = EnergyMatrix(256, matrix_256)
    proj = em.project_to_40()
    ref = EnergyMatrix.from_dat_file(reference_dat_path)

    # Truncate to min dimension
    n = min(proj.ntypes, ref.ntypes)
    proj_sub = proj.matrix[:n, :n]
    ref_sub = ref.matrix[:n, :n]

    # Correlation between upper triangles
    tri = np.triu_indices(n)
    p_vals = proj_sub[tri]
    r_vals = ref_sub[tri]

    if np.std(p_vals) < 1e-12 or np.std(r_vals) < 1e-12:
        return {"projection_r": 0.0, "projection_ntypes": n,
                "projection_rmse": float(np.sqrt(np.mean((p_vals - r_vals) ** 2)))}

    if HAS_SCIPY:
        r, _ = pearsonr(p_vals, r_vals)
    else:
        # Numpy-only correlation
        r = float(np.corrcoef(p_vals, r_vals)[0, 1])

    return {
        "projection_r": float(r),
        "projection_ntypes": n,
        "projection_rmse": float(np.sqrt(np.mean((p_vals - r_vals) ** 2))),
    }


# ─── full training pipeline ─────────────────────────────────────────────────

@dataclass
class TrainingConfig:
    """Configuration for the training pipeline."""
    pdbbind_dir: str = ""
    output_path: str = "matrix_256x256.bin"
    casf_dir: str = ""
    reference_dat: str = ""
    contact_cutoff: float = CONTACT_CUTOFF
    ridge_alpha: float = RIDGE_ALPHA
    lbfgs_maxiter: int = 200
    temperature: float = TEMPERATURE
    seed: int = 42


def load_pdbbind_complexes(pdbbind_dir: str,
                            cutoff: float = CONTACT_CUTOFF) -> List[Complex]:
    """Load all complexes from a PDBbind refined-set directory.

    Expected structure:
        pdbbind_dir/
            INDEX/INDEX_general_PL_data.2020
            1a0q/1a0q_protein.pdb
            1a0q/1a0q_ligand.mol2
            ...
    """
    pdbbind = Path(pdbbind_dir)
    if not pdbbind.is_dir():
        raise FileNotFoundError(f"PDBbind directory not found: {pdbbind_dir}")

    # Find index file
    index_file = None
    for pattern in ["INDEX/*PL_data*", "INDEX/*general*", "index/*"]:
        candidates = list(pdbbind.glob(pattern))
        if candidates:
            index_file = str(candidates[0])
            break

    affinities = {}
    if index_file:
        affinities = parse_pdbbind_index(index_file)
        logger.info("Loaded %d affinities from %s", len(affinities), index_file)

    complexes = []
    for subdir in sorted(pdbbind.iterdir()):
        if not subdir.is_dir():
            continue
        code = subdir.name.lower()

        pdb_files = list(subdir.glob("*_protein.pdb")) + list(subdir.glob("*_pocket.pdb"))
        mol2_files = list(subdir.glob("*_ligand.mol2"))

        if not pdb_files or not mol2_files:
            continue

        try:
            prot_atoms = parse_pdb_atoms(str(pdb_files[0]))
            lig_atoms = parse_mol2_atoms(str(mol2_files[0]))
            contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff)
        except Exception as e:
            logger.debug("Skipping %s: %s", code, e)
            continue

        if not contacts:
            continue

        pkd = affinities.get(code, 0.0)
        dg = -kB_kcal * TEMPERATURE * math.log(10) * pkd if pkd != 0.0 else 0.0

        cpx = Complex(
            pdb_code=code,
            protein_atoms=prot_atoms,
            ligand_atoms=lig_atoms,
            contacts=contacts,
            pKd=pkd,
            deltaG=dg,
        )
        complexes.append(cpx)

    logger.info("Loaded %d complexes from %s", len(complexes), pdbbind_dir)
    return complexes


def train_matrix(config: TrainingConfig) -> EnergyMatrix:
    """Run the full training pipeline.

    Steps:
        1. Load PDBbind complexes
        2. Enumerate contacts
        3. Build observed and reference frequency matrices
        4. Compute inverse Boltzmann (Sippl potential)
        5. Ridge regression per cell
        6. L-BFGS refinement against experimental ΔG
        7. Validate (CASF-2016 if available, 256→40 projection if reference)
        8. Save binary blob

    Returns:
        Trained EnergyMatrix (256×256).
    """
    np.random.seed(config.seed)

    # Step 1: Load data
    complexes = load_pdbbind_complexes(config.pdbbind_dir,
                                       config.contact_cutoff)
    if not complexes:
        raise RuntimeError("No complexes loaded from PDBbind directory")

    # Step 2–3: Build frequency matrices
    freq = build_contact_matrix(complexes)
    ref = build_reference_matrix(complexes)

    # Step 4: Inverse Boltzmann
    sippl = inverse_boltzmann(freq, ref, config.temperature)
    logger.info("Inverse Boltzmann matrix computed (non-zero cells: %d)",
                np.count_nonzero(sippl))

    # Step 5: Ridge regression
    ridge_matrix = ridge_fit(complexes, config.ridge_alpha)
    logger.info("Ridge regression fit completed")

    # Combine: 70% Sippl + 30% ridge
    combined = 0.7 * sippl + 0.3 * ridge_matrix
    combined = (combined + combined.T) / 2.0  # ensure symmetry

    # Step 6: L-BFGS refinement
    refined = lbfgs_refine(combined, complexes, config.lbfgs_maxiter)

    # Step 7: Validate
    metrics = validate_casf(refined, complexes)
    logger.info("Training set metrics: Pearson r=%.4f, RMSE=%.4f",
                metrics.get("pearson_r", 0), metrics.get("rmse", 0))

    if config.casf_dir and os.path.isdir(config.casf_dir):
        casf_complexes = load_pdbbind_complexes(config.casf_dir,
                                                 config.contact_cutoff)
        if casf_complexes:
            casf_metrics = validate_casf(refined, casf_complexes)
            logger.info("CASF-2016 metrics: Pearson r=%.4f, RMSE=%.4f",
                        casf_metrics.get("pearson_r", 0),
                        casf_metrics.get("rmse", 0))

    if config.reference_dat and os.path.isfile(config.reference_dat):
        proj_metrics = validate_projection(refined, config.reference_dat)
        logger.info("256→40 projection: r=%.4f, RMSE=%.4f",
                    proj_metrics.get("projection_r", 0),
                    proj_metrics.get("projection_rmse", 0))

    # Step 8: Save
    em = EnergyMatrix(256, refined)
    em.to_binary(config.output_path)
    logger.info("Saved trained matrix to %s", config.output_path)

    # Also save legacy .dat projection if requested
    if config.reference_dat:
        proj = em.project_to_40()
        proj_path = config.output_path.replace(".bin", "_40x40.dat")
        proj.to_dat_file(proj_path)
        logger.info("Saved 40×40 projection to %s", proj_path)

    return em


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Train 256×256 soft contact matrix from PDBbind data")
    parser.add_argument("--pdbbind-dir", required=True,
                        help="Path to PDBbind refined-set directory")
    parser.add_argument("--output", default="matrix_256x256.bin",
                        help="Output binary matrix path")
    parser.add_argument("--validate-casf", default="",
                        help="Path to CASF-2016 directory for validation")
    parser.add_argument("--project-40", default="",
                        help="Path to reference 40-type .dat for projection validation")
    parser.add_argument("--cutoff", type=float, default=CONTACT_CUTOFF,
                        help="Contact distance cutoff (Å)")
    parser.add_argument("--ridge-alpha", type=float, default=RIDGE_ALPHA,
                        help="Ridge regression regularisation")
    parser.add_argument("--lbfgs-maxiter", type=int, default=200,
                        help="L-BFGS maximum iterations")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    config = TrainingConfig(
        pdbbind_dir=args.pdbbind_dir,
        output_path=args.output,
        casf_dir=args.validate_casf,
        reference_dat=args.project_40,
        contact_cutoff=args.cutoff,
        ridge_alpha=args.ridge_alpha,
        lbfgs_maxiter=args.lbfgs_maxiter,
    )

    train_matrix(config)


if __name__ == "__main__":
    main()
