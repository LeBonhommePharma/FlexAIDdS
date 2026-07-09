"""High-level docking interface for FlexAID∆S.

Provides Pythonic API for molecular docking workflows.
"""

import subprocess
import shutil
import re
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Any, Sequence
from dataclasses import dataclass, replace

from .thermodynamics import Thermodynamics, StatMechEngine, kB_kcal
from .dift import (
    DiFTEngine,
    RotatableBondTorsion,
    TorsionalScore,
    make_bond_torsion,
    score_torsional,
)

try:
    from . import _core
except ImportError:
    _core = None


@dataclass
class Pose:
    """Single docked pose within a binding mode.
    
    Attributes:
        index: Pose index in GA population
        energy: Binding energy (CF score) in kcal/mol
        rmsd: RMSD to reference structure (if available)
        coordinates: Atomic coordinates (Nx3 array)
        boltzmann_weight: Statistical weight in ensemble
    """
    index: int
    energy: float
    rmsd: Optional[float] = None
    coordinates: Optional[np.ndarray] = None
    boltzmann_weight: float = 0.0
    
    def __repr__(self) -> str:
        rmsd_str = f" RMSD={self.rmsd:.2f}" if self.rmsd is not None else ""
        return (
            f"<Pose {self.index} E={self.energy:.3f}{rmsd_str} "
            f"w={self.boltzmann_weight:.4g}>"
        )

    def to_dict(self) -> dict:
        return {
            'index': self.index,
            'energy_kcal_mol': self.energy,
            'rmsd_angstrom': self.rmsd,
            'boltzmann_weight': self.boltzmann_weight,
        }


class BindingMode:
    """Binding mode: cluster of docked poses with thermodynamic scoring.

    A binding mode represents a distinct local minimum on the binding energy
    landscape, characterized by an ensemble of similar poses.

    The configurational thermodynamics come from the pose-energy ensemble. In
    addition, a binding mode can carry the DiFT (Discrete Fourier Transform)
    torsional potentials of the ligand's rotatable bonds; when attached, the
    per-bond torsional free-energy contribution ``V_tors(φ) − T·S_tors`` folds
    into ``free_energy``, ``enthalpy``, and ``entropy``. This replaces the
    classical crude per-rotatable-bond count penalty with a first-principles
    statistical-mechanical ΔS term derived from the same Fourier spectrum that
    parametrizes the torsional energy. See :mod:`flexaidds.dift`.

    Example:
        >>> mode = results.binding_modes[0]  # top-ranked mode
        >>> thermo = mode.get_thermodynamics()
        >>> print(f"ΔG = {thermo.free_energy:.2f} kcal/mol")
        >>> print(f"ΔH = {thermo.mean_energy:.2f}, TΔS = {thermo.entropy_term:.2f}")

        >>> # Fold in ligand torsional entropy from QM/CG dihedral profiles:
        >>> mode.set_torsional_profiles([qm_scan_bond0, qm_scan_bond1],
        ...                             dihedral_angles_rad=[phi0, phi1])
        >>> mode.torsional_free_energy  # kcal/mol (energy + confinement −TΔS)
    """

    def __init__(self, cpp_binding_mode=None, temperature: float = 300.0):
        """Initialize binding mode.

        Args:
            cpp_binding_mode: C++ BindingMode object (internal use, optional).
            temperature: Temperature in Kelvin for thermodynamic calculations.
        """
        self._cpp_mode = cpp_binding_mode
        self._poses: List[Pose] = []
        self._temperature: float = temperature
        self._cached_thermo: Optional[Thermodynamics] = None
        # Receptor-bound ions and cofactors present in the complex.
        # Each entry is a string "RESNAME:CHAIN:RESNUM", e.g. "MG:A:101".
        self.receptor_cofactors: List[str] = []
        # DiFT torsional potentials of the ligand's rotatable bonds (optional).
        # When empty, torsional contributions are exactly zero and the binding
        # mode behaves as a pure configurational ensemble.
        self._torsional_bonds: List[RotatableBondTorsion] = []
        # Representative dihedral state (radians) used for the torsional energy
        # term. May be None: the confinement −TΔS penalty is state-independent
        # and is still applied, while the energy term contributes zero.
        self._torsional_dihedrals: Optional[List[float]] = None
        self._cached_torsional: Optional[TorsionalScore] = None

    def _invalidate_cache(self) -> None:
        self._cached_thermo = None
        self._cached_torsional = None

    def _compute_python_thermo(self) -> Thermodynamics:
        """Compute thermodynamics from pose energies using StatMechEngine."""
        if self._cached_thermo is not None:
            return self._cached_thermo
        if not self._poses:
            return Thermodynamics(
                temperature=self._temperature, log_Z=0.0,
                free_energy=float('inf'), mean_energy=float('inf'),
                mean_energy_sq=float('inf'), heat_capacity=0.0,
                entropy=0.0, std_energy=0.0,
            )
        engine = StatMechEngine(self._temperature)
        for pose in self._poses:
            engine.add_sample(pose.energy)
        self._cached_thermo = engine.compute()
        return self._cached_thermo

    # ── DiFT torsional contribution ─────────────────────────────────────────

    def set_torsional_potentials(
            self,
            bonds: Sequence[RotatableBondTorsion],
            dihedral_angles_rad: Optional[Sequence[float]] = None) -> None:
        """Attach DiFT torsional potentials for the ligand's rotatable bonds.

        Each bond carries a DiFT-parametrized :class:`~flexaidds.dift.\
TorsionalPotential` (from a QM scan or a Boltzmann-inverted coarse-grained
        dihedral histogram). Once attached, this binding mode's ``free_energy``,
        ``enthalpy``, and ``entropy`` include the torsional contribution.

        Args:
            bonds: Per-rotatable-bond DiFT potentials.
            dihedral_angles_rad: Representative dihedral value (radians) of each
                bond, in the same order as *bonds* — typically the dihedral
                state of the mode's best pose. Drives the torsional *energy*
                term; if *None*, only the state-independent confinement −TΔS
                penalty is applied (the energy term is zero). Must match the
                length of *bonds* to be used.
        """
        self._torsional_bonds = list(bonds)
        self._torsional_dihedrals = (
            list(dihedral_angles_rad) if dihedral_angles_rad is not None else None)
        self._cached_torsional = None

    def set_torsional_profiles(
            self,
            profiles: Sequence[Sequence[float]],
            dihedral_angles_rad: Optional[Sequence[float]] = None,
            temperature_K: Optional[float] = None,
            max_multiplicity: int = 6) -> None:
        """Attach torsional potentials by DiFT-parametrizing raw profiles.

        Convenience wrapper over :func:`~flexaidds.dift.make_bond_torsion`:
        each entry of *profiles* is an M-point torsional energy profile over
        [0, 2π) (a QM scan or a Boltzmann-inverted CG histogram) which is
        transformed and Shannon-collapse truncated on the spot.

        Args:
            profiles: One torsional profile per rotatable bond.
            dihedral_angles_rad: Representative dihedral of each bond (radians);
                see :meth:`set_torsional_potentials`.
            temperature_K: Temperature for the DiFT fit; defaults to this
                mode's temperature.
            max_multiplicity: Anti-overfit cap on Fourier multiplicity.
        """
        temp = self._temperature if temperature_K is None else temperature_K
        bonds = [
            make_bond_torsion(profile, gene_index=i, temperature_K=temp,
                              max_multiplicity=max_multiplicity)
            for i, profile in enumerate(profiles)
        ]
        self.set_torsional_potentials(bonds, dihedral_angles_rad)

    @property
    def torsional_score(self) -> TorsionalScore:
        """Decomposed torsional contribution (kcal/mol) of the attached bonds.

        ``energy`` is Σ V_tors,b(φ_b) relative to each well minimum (zero unless
        representative dihedral angles were supplied); ``minus_TS`` is the
        confinement −TΔS penalty, well-defined from the potentials alone.
        Returns an all-zero score when no potentials are attached.
        """
        if self._cached_torsional is not None:
            return self._cached_torsional

        bonds = self._torsional_bonds
        angles = self._torsional_dihedrals
        if not bonds:
            score = TorsionalScore()
        elif angles is not None and len(angles) == len(bonds):
            # Full energy + entropy from the representative dihedral state.
            score = score_torsional(bonds, angles, self._temperature)
        else:
            # No (or mismatched) dihedral state: apply only the state-independent
            # confinement −TΔS penalty; leave the energy term at zero.
            engine = DiFTEngine(self._temperature)
            score = TorsionalScore()
            for bond in bonds:
                score.minus_TS += engine.thermodynamics(bond.potential).minus_TS
                score.n_bonds += 1
        self._cached_torsional = score
        return score

    @property
    def torsional_free_energy(self) -> float:
        """Torsional free-energy contribution ΔG_tors (kcal/mol).

        ``Σ_b [V_tors,b(φ_b) − T·S_tors,b]`` — energy at the representative
        dihedral state plus the confinement entropy penalty. Zero when no
        DiFT potentials are attached.
        """
        return self.torsional_score.total()

    @property
    def torsional_entropy(self) -> float:
        """Torsional entropy S_tors (kcal mol⁻¹ K⁻¹), ≤ 0 (a confinement loss).

        Derived from the same Fourier spectra as the torsional energy; zero
        when no DiFT potentials are attached.
        """
        # minus_TS = −T·S_tors  ⇒  S_tors = −minus_TS / T.
        return -self.torsional_score.minus_TS / self._temperature

    def get_thermodynamics(self) -> Thermodynamics:
        """Get full thermodynamic properties of this binding mode.

        When DiFT torsional potentials are attached, the returned free energy,
        enthalpy (``mean_energy``), and entropy include the torsional
        contribution on top of the configurational ensemble.

        Returns:
            Thermodynamics object with F, S, H, Cv, etc.
        """
        if self._cpp_mode is not None:
            thermo_cpp = self._cpp_mode.get_thermodynamics()
            base = Thermodynamics(
                temperature=thermo_cpp.temperature,
                log_Z=thermo_cpp.log_Z,
                free_energy=thermo_cpp.free_energy,
                mean_energy=thermo_cpp.mean_energy,
                mean_energy_sq=thermo_cpp.mean_energy_sq,
                heat_capacity=thermo_cpp.heat_capacity,
                entropy=thermo_cpp.entropy,
                std_energy=thermo_cpp.std_energy,
            )
        else:
            base = self._compute_python_thermo()
        return self._with_torsional(base)

    def _with_torsional(self, base: Thermodynamics) -> Thermodynamics:
        """Fold the torsional contribution into a configurational Thermodynamics.

        Returns *base* unchanged when no torsional potentials are attached, so
        the pure-configurational path is untouched.
        """
        score = self.torsional_score
        if score.n_bonds == 0:
            return base
        return replace(
            base,
            free_energy=base.free_energy + score.total(),
            mean_energy=base.mean_energy + score.energy,
            entropy=base.entropy + self.torsional_entropy,
        )

    @property
    def free_energy(self) -> float:
        """Helmholtz free energy F = -kT ln Z (kcal/mol), incl. torsional term."""
        if self._cpp_mode:
            base = self._cpp_mode.get_free_energy()
        else:
            base = self._compute_python_thermo().free_energy
        return base + self.torsional_free_energy

    @property
    def enthalpy(self) -> float:
        """Boltzmann-weighted average energy ⟨E⟩ (kcal/mol), incl. torsional."""
        if self._cpp_mode:
            base = self._cpp_mode.compute_enthalpy()
        else:
            base = self._compute_python_thermo().mean_energy
        return base + self.torsional_score.energy

    @property
    def entropy(self) -> float:
        """Entropy S (kcal mol⁻¹ K⁻¹) — configurational plus torsional."""
        if self._cpp_mode:
            base = self._cpp_mode.compute_entropy()
        else:
            base = self._compute_python_thermo().entropy
        return base + self.torsional_entropy

    @property
    def n_poses(self) -> int:
        """Number of poses in this binding mode."""
        if self._cpp_mode:
            return self._cpp_mode.get_BindingMode_size()
        return len(self._poses)

    def __len__(self) -> int:
        return self.n_poses

    def __repr__(self) -> str:
        return (f"<BindingMode n_poses={self.n_poses} "
                f"F={self.free_energy:.2f} H={self.enthalpy:.2f} "
                f"S={self.entropy:.5f}>")


class BindingPopulation:
    """Collection of binding modes from a docking run.
    
    Provides ensemble-level analysis and ranking of binding modes.
    """
    
    def __init__(self, modes: Optional[List[BindingMode]] = None,
                 temperature: float = 300.0):
        self._modes: List[BindingMode] = list(modes) if modes else []
        self._temperature: float = temperature
    
    def add_mode(self, mode: BindingMode) -> None:
        """Add a binding mode to the population."""
        self._modes.append(mode)
    
    def rank_by_free_energy(self) -> List[BindingMode]:
        """Return binding modes sorted by free energy (best first)."""
        return sorted(self._modes, key=lambda m: m.free_energy)
    
    def compute_global_thermodynamics(self) -> Thermodynamics:
        """Compute thermodynamics over all binding modes.

        Aggregates all pose energies from all modes into a single
        canonical ensemble.

        Returns:
            Global ensemble thermodynamics
        """
        engine = StatMechEngine(self._temperature)
        for mode in self._modes:
            for pose in mode._poses:
                engine.add_sample(pose.energy)
        return engine.compute()

    def compute_super_cluster_thermodynamics(self) -> Thermodynamics:
        """Compute thermodynamics using only the super-cluster subset.

        Extracts the dominant energy basin via SuperCluster, then
        computes canonical ensemble thermodynamics on the filtered set.

        Returns:
            Thermodynamics for the super-cluster subset.
        """
        from .supercluster import SuperCluster
        all_energies = [p.energy for m in self._modes for p in m._poses]
        if not all_energies:
            return Thermodynamics(
                temperature=self._temperature, log_Z=0.0,
                free_energy=float('inf'), mean_energy=float('inf'),
                mean_energy_sq=float('inf'), heat_capacity=0.0,
                entropy=0.0, std_energy=0.0,
            )
        sc = SuperCluster(all_energies)
        filtered = sc.filter_energies()
        engine = StatMechEngine(self._temperature)
        for e in filtered:
            engine.add_sample(e)
        return engine.compute()
    
    def get_shannon_entropy(self) -> float:
        """Population-level Shannon configurational entropy.

        S = -kB * sum(p_i * ln(p_i)) over all poses across all modes,
        where p_i are Boltzmann probabilities.

        Returns:
            Shannon entropy in kcal/mol/K
        """
        import math

        all_energies = []
        for mode in self._modes:
            for pose in mode._poses:
                all_energies.append(pose.energy)

        if not all_energies:
            return 0.0

        beta = 1.0 / (kB_kcal * self._temperature)
        # Log-sum-exp for numerical stability
        neg_beta_e = [-beta * e for e in all_energies]
        max_val = max(neg_beta_e)
        log_Z = max_val + math.log(sum(math.exp(v - max_val) for v in neg_beta_e))

        shannon_S = 0.0
        for e in all_energies:
            log_p = -beta * e - log_Z
            p = math.exp(log_p)
            if p > 1e-30:
                shannon_S -= p * log_p
        shannon_S *= kB_kcal

        return shannon_S

    def get_deltaG_matrix(self) -> List[List[float]]:
        """ΔG matrix between all pairs of binding modes.

        matrix[i][j] = F_i - F_j. Anti-symmetric: matrix[i][j] = -matrix[j][i].

        Returns:
            n x n matrix of pairwise ΔG values (kcal/mol)
        """
        n = len(self._modes)
        energies = [m.free_energy for m in self._modes]
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                dg = energies[i] - energies[j]
                matrix[i][j] = dg
                matrix[j][i] = -dg
        return matrix

    @property
    def n_modes(self) -> int:
        """Number of binding modes."""
        return len(self._modes)
    
    def __len__(self) -> int:
        return self.n_modes
    
    def __getitem__(self, index: int) -> BindingMode:
        return self._modes[index]
    
    def __iter__(self):
        return iter(self._modes)
    
    def __repr__(self) -> str:
        return f"<BindingPopulation n_modes={self.n_modes} T={self._temperature}K>"


class Docking:
    """High-level interface for FlexAID∆S molecular docking.
    
    Example:
        >>> docking = Docking("config.inp")
        >>> results = docking.run()
        >>> top_mode = results.binding_modes[0]
        >>> print(f"Best ΔG: {top_mode.free_energy:.2f} kcal/mol")
    """
    
    def __init__(self, config_file: str):
        """Initialize docking from configuration file.
        
        Args:
            config_file: Path to FlexAID .inp config file
        """
        self.config_file = Path(config_file)
        if not self.config_file.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")
        
        self._config: Dict[str, Any] = {}
        self._parse_config()
    
    def _parse_config(self) -> None:
        """Parse FlexAID config file.

        The format is fixed-width: the first 6 characters are the keyword,
        character 7 is a space delimiter, and the remainder of the line is the
        value.  Lines that start with '#' or are blank are ignored.

        Keywords that may appear multiple times (OPTIMZ, FLEXSC) are collected
        into lists.  All other keywords map to a single string value (or a
        boolean ``True`` for flag-only keywords such as EXCHET, ROTOBS, etc.).

        After parsing, ``self._config`` is populated with keys matching the
        6-character keyword names used by the C++ FlexAID engine.
        """
        # Keywords whose value is the rest of the line (path / string).
        _string_keys = {
            "PDBNAM", "INPLIG", "RNGOPT", "METOPT", "BPKENM", "COMPLF",
            "VCTSCO", "IMATRX", "DEFTYP", "CONSTR", "NMAAMP", "NMAEIG",
            "RMSDST", "DEPSPA", "STATEP", "TEMPOP", "CLUSTA",
        }
        # Keywords whose value is a float.
        _float_keys = {
            "ACSWEI", "CLRMSD", "PERMEA", "INTRAF", "VARDIS", "VARANG",
            "VARDIH", "VARFLX", "SLVPEN", "DEECLA", "ROTPER", "SPACER",
        }
        # Keywords whose value is an integer.
        _int_keys = {
            "NMAMOD", "MAXRES", "TEMPER", "NRGOUT",
        }
        # Keywords that appear multiple times; values are collected into a list.
        _list_keys = {"OPTIMZ", "FLEXSC"}
        # Flag-only keywords: presence means True, no value expected.
        _flag_keys = {
            "DEEFLX", "ROTOBS", "NORMAR", "USEACS", "EXCHET", "INCHOH",
            "NOINTR", "OMITBU", "VINDEX", "HTPMOD", "OUTRNG", "USEDEE",
            "NRGSUI", "SCOLIG", "SCOOUT", "ROTOUT", "SUPCLU",
        }

        # Initialise list accumulators so callers can always iterate them.
        for key in _list_keys:
            self._config[key] = []

        with open(self.config_file) as fh:
            for raw_line in fh:
                line = raw_line.rstrip("\n").rstrip("\r")

                # Skip blank lines and comments.
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                # The keyword occupies exactly the first 6 characters.
                if len(line) < 6:
                    continue
                keyword = line[:6].strip()
                value_str = line[7:].strip() if len(line) > 7 else ""

                if keyword in _flag_keys:
                    self._config[keyword] = True
                elif keyword in _list_keys:
                    self._config[keyword].append(value_str)
                elif keyword in _float_keys:
                    try:
                        self._config[keyword] = float(value_str.split()[0])
                    except (ValueError, IndexError):
                        self._config[keyword] = value_str
                elif keyword in _int_keys:
                    try:
                        self._config[keyword] = int(value_str.split()[0])
                    except (ValueError, IndexError):
                        self._config[keyword] = value_str
                elif keyword in _string_keys:
                    self._config[keyword] = value_str
                else:
                    # Unknown keyword: store raw value string.
                    self._config[keyword] = value_str

    @property
    def receptor(self) -> Optional[str]:
        """Path to receptor PDB file (PDBNAM keyword)."""
        return self._config.get("PDBNAM")

    @property
    def ligand(self) -> Optional[str]:
        """Path to ligand input file (INPLIG keyword)."""
        return self._config.get("INPLIG")

    @property
    def temperature(self) -> Optional[int]:
        """Simulation temperature in Kelvin (TEMPER keyword)."""
        return self._config.get("TEMPER")

    @property
    def optimization_method(self) -> Optional[str]:
        """Optimization method, e.g. 'GA' (METOPT keyword)."""
        return self._config.get("METOPT")
    
    def run(self, binary: Optional[str] = None,
            timeout: int = 3600, **kwargs) -> BindingPopulation:
        """Execute docking via the FlexAID C++ binary and parse results.

        Locates the ``FlexAID`` binary (in PATH, project build/, or explicit
        *binary* argument), invokes it with this config file, waits for
        completion, then parses all ``*_N_M.pdb`` output files written by
        ``output_Population()`` to reconstruct a ``BindingPopulation``.

        Args:
            binary:  Path to FlexAID executable.  If *None*, searches PATH and
                     common build locations (``build/FlexAID``,
                     ``../build/FlexAID``).
            timeout: Wall-clock timeout in seconds (default 3600).
            **kwargs: Ignored; reserved for future keyword overrides.

        Returns:
            BindingPopulation populated from the PDB REMARK lines written by
            ``output_BindingMode()`` / ``output_Population()``.

        Raises:
            FileNotFoundError: binary not found.
            RuntimeError:      FlexAID exited non-zero or produced no output.
        """
        # ── 1. Locate binary ─────────────────────────────────────────────────
        exe = self._find_binary(binary)

        # ── 2. Invoke FlexAID ────────────────────────────────────────────────
        cmd = [str(exe), str(self.config_file)]
        try:
            result = subprocess.run(
                cmd,
                cwd=self.config_file.parent,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"FlexAID timed out after {timeout}s"
            ) from exc

        if result.returncode != 0:
            raise RuntimeError(
                f"FlexAID exited with code {result.returncode}.\n"
                f"stdout: {result.stdout[-2000:]}\n"
                f"stderr: {result.stderr[-2000:]}"
            )

        # ── 3. Discover output PDBs ──────────────────────────────────────────
        # output_BindingMode writes files named <prefix>_<minPoints>_<mode>.pdb
        # Collect all candidate PDB files in the working directory.
        work_dir = self.config_file.parent
        pdb_files = sorted(work_dir.glob("*_*.pdb"),
                           key=lambda p: p.stat().st_mtime)

        if not pdb_files:
            raise RuntimeError(
                "FlexAID completed but no PDB output files were found in "
                f"{work_dir}. Check the config file NRGOUT / output settings."
            )

        # ── 4. Parse PDB REMARK lines into BindingModes ──────────────────────
        temperature = self._config.get("TEMPER", 300) or 300
        modes: List[BindingMode] = []
        seen_modes: Dict[int, BindingMode] = {}

        for pdb_path in pdb_files:
            mode_info = self._parse_remark_pdb(pdb_path, temperature)
            if mode_info is None:
                continue
            mode_idx, pose = mode_info
            if mode_idx not in seen_modes:
                seen_modes[mode_idx] = BindingMode(temperature=float(temperature))
            seen_modes[mode_idx]._poses.append(pose)

        # Sort modes by free energy (ascending → most favourable first)
        modes = sorted(seen_modes.values(),
                       key=lambda m: m.free_energy)

        return BindingPopulation(modes, temperature=float(temperature))

    # ── helpers ───────────────────────────────────────────────────────────────

    def _find_binary(self, binary: Optional[str]) -> Path:
        """Locate the FlexAID executable."""
        if binary is not None:
            p = Path(binary)
            if not p.is_file():
                raise FileNotFoundError(f"Specified FlexAID binary not found: {binary}")
            return p

        # Search order: PATH → project-relative build dirs
        in_path = shutil.which("FlexAID")
        if in_path:
            return Path(in_path)

        candidates = [
            self.config_file.parent / "FlexAID",
            self.config_file.parent / "build" / "FlexAID",
            Path(__file__).parents[3] / "build" / "FlexAID",
        ]
        for c in candidates:
            if c.is_file():
                return c

        raise FileNotFoundError(
            "FlexAID binary not found in PATH or build/. "
            "Build with 'cmake --build build' or pass binary= argument."
        )

    @staticmethod
    def _parse_remark_pdb(
            pdb_path: Path, temperature: float) -> Optional[tuple]:
        """Parse a single output PDB written by output_BindingMode().

        Extracts mode index, CF, RMSD, and per-pose energy from REMARK lines.
        Returns (mode_index, Pose) or None if the file lacks FlexAID remarks.
        """
        mode_idx   = None
        cf_val     = None
        rmsd_val   = None
        freq       = 1

        try:
            text = pdb_path.read_text(errors="replace")
        except OSError:
            return None

        for line in text.splitlines():
            if not line.startswith("REMARK"):
                continue
            # "REMARK Binding Mode:N Best CF in Binding Mode:X …"
            m = re.search(
                r"Binding Mode:(\d+).*?Best CF in Binding Mode:\s*([-\d.]+)"
                r".*?Binding Mode Frequency:(\d+)",
                line)
            if m:
                mode_idx = int(m.group(1))
                cf_val   = float(m.group(2))
                freq     = int(m.group(3))
            # "REMARK 0.12345 RMSD to ref. structure …"
            m2 = re.search(r"REMARK\s+([\d.]+)\s+RMSD to ref\.", line)
            if m2 and rmsd_val is None:
                rmsd_val = float(m2.group(1))

        if mode_idx is None or cf_val is None:
            return None

        import math
        beta = 1.0 / (0.001987206 * float(temperature))
        bw   = math.exp(-beta * cf_val)

        pose = Pose(
            index=mode_idx,
            energy=cf_val,
            rmsd=rmsd_val,
            boltzmann_weight=bw,
        )
        return mode_idx, pose
    
    def __repr__(self) -> str:
        return f"<Docking config={self.config_file.name}>"
