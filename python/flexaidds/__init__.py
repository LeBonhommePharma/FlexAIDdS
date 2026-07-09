"""flexaidds: Python bindings and analysis tools for FlexAID∆S."""

from .models import BindingModeResult, DockingResult, PoseResult
from .results import load_results
from .io import is_ion, _ION_RESNAMES, Atom, PDBStructure, read_pdb, write_pdb
from .docking import Docking, BindingMode, BindingPopulation, Pose
from .encom import ENCoMEngine, FrequencyCalibration, NormalMode, VibrationalEntropy
from .tencm import (
    TorsionalENM, TorsionalNormalMode, Conformer, FullThermoResult,
    compute_shannon_entropy, compute_torsional_vibrational_entropy,
    run_shannon_thermo_stack,
)
from .dift import (
    DiFTEngine, FourierTerm, TorsionalPotential, TorsionalThermo,
    RotatableBondTorsion, TorsionalScore, score_torsional, make_bond_torsion,
    spectral_entropy,
)
from .__version__ import __version__ as __version__
from .updater import check_for_updates, UpdateInfo
from .boltz2 import (
    Boltz2Client,
    Boltz2PredictionResult,
    Boltz2AffinityResult,
    Boltz2Polymer,
    Boltz2Ligand,
    PocketConstraint,
    PocketContact,
    Boltz2Error,
)
from .benchmark import (
    BenchmarkSystem,
    MethodResult,
    SystemBenchmarkResult,
    BenchmarkResult,
    BenchmarkSummary,
    run_benchmark,
    load_benchmark_dataset,
    save_benchmark_dataset,
    auto_workers,
)

# Pure-Python thermodynamics (always available)
from .thermodynamics import (
    StatMechEngine,
    Thermodynamics,
    ThermodynamicBreakdown,
    kB_kcal,
    kB_SI,
    deltaG_standard_to_Kd_M,
    Kd_M_to_deltaG_standard,
    StabilityCurve,
    gibbs_helmholtz_dG,
    kirchhoff_dH,
    kirchhoff_dS,
)

# Pure-Python μVT competitive binding (mirrors GrandPartitionFunction; no C++ required)
from .grand_canonical import (
    CompetitiveSite,
    OccupancyPoint,
    set_concentration,
    plot_occupancy_curve,
)

# entropy.help audit schema (A1.1) — pure Python, always available
from .schemas.thermo_audit import (
    ThermodynamicOutput,
    ThermodynamicOutputDC,
    TotalSampledPartitionFunction,
    TotalSampledPartitionFunctionDC,
    Provenance,
    ProvenanceDC,
    make_total_sampled_output,
)

# C++ extension — optional: pure-Python helpers work without it
try:
    from ._core import (
        BoltzmannLUT,
        Replica,
        State,
        TIPoint,
        WHAMBin,
        TemperatureScanPoint,
        DeltaCpFit,
        kB_kcal,  # noqa: F811  (more precise C++ value)
        kB_SI,    # noqa: F811
    )
    # Override with C++ StatMechEngine when available
    from ._core import StatMechEngine as _CppStatMechEngine  # noqa: F811
    from ._core import Thermodynamics as _CppThermodynamics  # noqa: F811
    from ._core import ThermodynamicBreakdown as _CppThermodynamicBreakdown  # noqa: F811
    from ._core import ENCoMEngine as _CppENCoMEngine        # noqa: F811
    from ._core import NormalMode as _CppNormalMode           # noqa: F811
    from ._core import VibrationalEntropy as _CppVibrationalEntropy  # noqa: F811
    HAS_CORE_BINDINGS = True
except ImportError:
    # Fallback when C++ extension is not built
    from ._fallback_types import (
        BoltzmannLUT, Replica, State, TIPoint, WHAMBin,
        TemperatureScanPoint, DeltaCpFit,
    )
    kB_kcal = 0.001987206   # kcal mol⁻¹ K⁻¹
    kB_SI = 1.380649e-23    # J K⁻¹
    HAS_CORE_BINDINGS = False

# GA hyperparameter optimizer
from .optimize import GAOptimizer, OptimizationResult

# ML rescoring bridge
from .ml_rescore import (
    VoronoiGraphExtractor,
    ShannonProfileExtractor,
    FeatureBuilder,
    ThermoFeatures,
    MLRescorer,
)

# Dataset runner (requires numpy + pyyaml optional deps)
try:
    from .dataset_runner import (
        DatasetRunner,
        DatasetConfig,
        DatasetResult,
        BenchmarkReport,
    )
except ImportError:
    pass  # numpy/pyyaml not installed; dataset_runner not available

from .supercluster import SuperCluster
from .tencom_results import FlexModeResult, FlexPopulationResult, parse_tencom_pdb, parse_tencom_json
from .energy_matrix import (
    EnergyMatrix,
    MatrixEntry,
    DensityPoint,
    encode_256_type,
    decode_256_type,
    base_to_sybyl,
    sybyl_to_base,
    parse_dat_file,
    write_dat_file,
    SYBYL_TYPE_NAMES,
    SYBYL_RADII,
)
from .reporting import (
    generate_pymol_script,
    generate_markdown_report,
    generate_temperature_scan_plot,
    write_all_reports,
)

# Imagine / publication figure + animation generation (Grok Imagine integration, post-hoc only)
from .figures import (
    prepare_publication_figures,
    build_imagine_cover_prompt,
    build_imagine_animation_prompt,
    check_gate6_passed,
    extract_best_mode_summary,
    # NRDD / FlexAID∆S journal cover generation (dramatic entropy-enthalpy, E-E index, full branding)
    NRDDCoverParams,
    build_nrdd_cover_prompt,
    generate_flexaids_nrdd_cover,
)


def dock(
    receptor: str,
    ligand: str,
    *,
    binding_site: str = "auto",
    compute_entropy: bool = True,
    temperature: float = 300.0,
    config: str = "",
    binary: str = None,
    timeout: int = 3600,
) -> BindingPopulation:
    """High-level docking interface.

    Generates a FlexAID config, runs the C++ engine, and returns results
    as a :class:`BindingPopulation`.

    Args:
        receptor:        Path to receptor PDB file.
        ligand:          Path to ligand MOL2 file.
        binding_site:    Binding site specification (default ``"auto"``).
        compute_entropy: Enable thermodynamic entropy calculation.
        temperature:     Simulation temperature in Kelvin.
        config:          Optional path to a FlexAID ``.inp`` config file.
                         When provided, ``receptor`` and ``ligand`` are ignored
                         and the config is used directly.
        binary:          Path to the FlexAID executable (auto-detected if None).
        timeout:         Wall-clock timeout in seconds.

    Returns:
        :class:`BindingPopulation` with ranked binding modes.

    Example:
        >>> results = flexaidds.dock('receptor.pdb', 'ligand.mol2')
        >>> for mode in results.rank_by_free_energy():
        ...     print(f"Mode: ΔG={mode.free_energy:.2f} kcal/mol")
    """
    import tempfile
    from pathlib import Path

    if config:
        docking = Docking(config)
        return docking.run(binary=binary, timeout=timeout)

    # Generate a temporary .inp config from receptor + ligand
    receptor_path = Path(receptor).resolve()
    ligand_path = Path(ligand).resolve()
    if not receptor_path.is_file():
        raise FileNotFoundError(f"Receptor not found: {receptor}")
    if not ligand_path.is_file():
        raise FileNotFoundError(f"Ligand not found: {ligand}")

    lines = [
        f"PDBNAM {receptor_path}",
        f"INPLIG {ligand_path}",
        f"TEMPER {int(temperature)}",
        "METOPT GA",
        "COMPLF VCT",
    ]
    if not compute_entropy:
        lines.append("TEMPER 0")
    if binding_site != "auto":
        lines.append(f"RNGOPT {binding_site}")

    tmp_dir = Path(tempfile.mkdtemp(prefix="flexaidds_"))
    cfg_path = tmp_dir / "dock.inp"
    cfg_path.write_text("\n".join(lines) + "\n")

    docking = Docking(str(cfg_path))
    return docking.run(binary=binary, timeout=timeout)


__all__ = [
    # Version
    "__version__",
    # Data models
    "PoseResult",
    "BindingModeResult",
    "DockingResult",
    # Result loading
    "load_results",
    # I/O utilities
    "is_ion",
    "_ION_RESNAMES",
    "Atom",
    "PDBStructure",
    "read_pdb",
    "write_pdb",
    # Docking
    "Docking",
    "BindingMode",
    "BindingPopulation",
    "Pose",
    "dock",
    # ENCoM
    "ENCoMEngine",
    "FrequencyCalibration",
    "NormalMode",
    "VibrationalEntropy",
    # Torsional ENCoM
    "TorsionalENM",
    "TorsionalNormalMode",
    "Conformer",
    "FullThermoResult",
    "compute_shannon_entropy",
    "compute_torsional_vibrational_entropy",
    "run_shannon_thermo_stack",
    # Thermodynamics (pure-Python or C++ override)
    "StatMechEngine",
    "Thermodynamics",
    "kB_kcal",
    "kB_SI",
    # Grand-canonical competitive binding (pure-Python μVT helpers)
    "CompetitiveSite",
    "OccupancyPoint",
    "set_concentration",
    "plot_occupancy_curve",
    # Availability flag
    "HAS_CORE_BINDINGS",
    # Core types (C++ when available, pure-Python fallback otherwise)
    "State",
    "BoltzmannLUT",
    "Replica",
    "WHAMBin",
    "TIPoint",
    "TemperatureScanPoint",
    "DeltaCpFit",
    # Affinity converters (Task 6/7, safe, with validation)
    "deltaG_standard_to_Kd_M",
    "Kd_M_to_deltaG_standard",
    "StabilityCurve",
    "gibbs_helmholtz_dG",
    "kirchhoff_dH",
    "kirchhoff_dS",
    # Reporting (Task 9 - JSON driven only)
    "generate_pymol_script",
    "generate_markdown_report",
    "generate_temperature_scan_plot",
    "write_all_reports",
    # Updater
    "check_for_updates",
    "UpdateInfo",
    # Boltz2
    "Boltz2Client",
    "Boltz2PredictionResult",
    "Boltz2AffinityResult",
    "Boltz2Polymer",
    "Boltz2Ligand",
    "PocketConstraint",
    "PocketContact",
    "Boltz2Error",
    # Benchmark
    "BenchmarkSystem",
    "MethodResult",
    "SystemBenchmarkResult",
    "BenchmarkResult",
    "BenchmarkSummary",
    "run_benchmark",
    "load_benchmark_dataset",
    "save_benchmark_dataset",
    "auto_workers",
    # SuperCluster
    "SuperCluster",
    # tENCoM results
    "FlexModeResult",
    "FlexPopulationResult",
    "parse_tencom_pdb",
    "parse_tencom_json",
    # Energy matrix I/O and 256-type encoding (always available — pure Python)
    "EnergyMatrix",
    "MatrixEntry",
    "DensityPoint",
    "encode_256_type",
    "decode_256_type",
    "base_to_sybyl",
    "sybyl_to_base",
    "parse_dat_file",
    "write_dat_file",
    "SYBYL_TYPE_NAMES",
    "SYBYL_RADII",
    # Dataset runner
    "DatasetRunner",
    "DatasetConfig",
    "DatasetResult",
    "BenchmarkReport",
    # GA hyperparameter optimizer
    "GAOptimizer",
    "OptimizationResult",
    # ML rescoring bridge
    "VoronoiGraphExtractor",
    "ShannonProfileExtractor",
    "FeatureBuilder",
    "ThermoFeatures",
    "MLRescorer",
]
