// core_bindings.cpp — pybind11 bindings for FlexAID∆S C++ core
//
// Exposes:
//   - statmech::StatMechEngine, Thermodynamics, BoltzmannLUT
//   - BindingMode / BindingPopulation (Phase 1/2 StatMech API)
//   - encom::ENCoMEngine, VibrationalEntropy, NormalMode (Phase 3)
//
// Build: See python/setup.py and CMakeLists.txt with -DBUILD_PYTHON_BINDINGS=ON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../../LIB/statmech.h"
#ifdef FLEXAIDS_FULL_GA
#include "../../LIB/BindingMode.h"
#endif
#include "../../LIB/encom.h"
#include "../../LIB/fast_optics.hpp"
#include "../../LIB/Spectrophore.h"
#include "../../LIB/BindingResidues.h"

namespace py = pybind11;
using namespace statmech;

// ──────────────────────────────────────────────────────────────────────────────
// Helper: Convert C++ std::vector to NumPy array (zero-copy view when possible)
// ──────────────────────────────────────────────────────────────────────────────

template <typename T>
py::array_t<T> to_numpy(const std::vector<T>& vec) {
    return py::array_t<T>(vec.size(), vec.data());
}

// ──────────────────────────────────────────────────────────────────────────────
// Module definition
// ──────────────────────────────────────────────────────────────────────────────

PYBIND11_MODULE(_core, m) {
    m.doc() = "FlexAID∆S C++ core bindings: statistical mechanics and docking engine";
    
    // Physical constants
    m.attr("kB_kcal") = kB_kcal;
    m.attr("kB_SI")   = kB_SI;
    
    // ═══════════════════════════════════════════════════════════════════════
    // Thermodynamics data structures
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<State>(m, "State", "Micr ostate with energy and degeneracy")
        .def(py::init<>())
        .def_readwrite("energy", &State::energy, "Energy in kcal/mol")
        .def_readwrite("count",  &State::count,  "Degeneracy/multiplicity")
        .def("__repr__", [](const State& s) {
            return "<State energy=" + std::to_string(s.energy) + 
                   " count=" + std::to_string(s.count) + ">";
        });
    
    py::class_<Thermodynamics>(m, "Thermodynamics", 
        "Complete thermodynamic properties of an ensemble")
        .def(py::init<>())
        .def_readwrite("temperature",    &Thermodynamics::temperature,    "K")
        .def_readwrite("log_Z",          &Thermodynamics::log_Z,          "ln(partition function)")
        .def_readwrite("free_energy",    &Thermodynamics::free_energy,    "Helmholtz F (kcal/mol)")
        .def_readwrite("mean_energy",    &Thermodynamics::mean_energy,    "⟨E⟩ (kcal/mol)")
        .def_readwrite("mean_energy_sq", &Thermodynamics::mean_energy_sq, "⟨E²⟩")
        .def_readwrite("heat_capacity",  &Thermodynamics::heat_capacity,  "Cv (kcal mol⁻¹ K⁻²)")
        .def_readwrite("entropy",        &Thermodynamics::entropy,        "S (kcal mol⁻¹ K⁻¹)")
        .def_readwrite("std_energy",     &Thermodynamics::std_energy,     "σ_E (kcal/mol)")
        .def("__repr__", [](const Thermodynamics& t) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<Thermodynamics T=%.1fK F=%.3f H=%.3f S=%.5f Cv=%.3f>",
                t.temperature, t.free_energy, t.mean_energy, 
                t.entropy, t.heat_capacity);
            return std::string(buf);
        });
    
    py::class_<Replica>(m, "Replica", "Parallel tempering replica")
        .def(py::init<>())
        .def_readwrite("id",              &Replica::id)
        .def_readwrite("temperature",     &Replica::temperature)
        .def_readwrite("beta",            &Replica::beta)
        .def_readwrite("current_energy",  &Replica::current_energy);
    
    py::class_<WHAMBin>(m, "WHAMBin", "WHAM histogram bin with free energy")
        .def(py::init<>())
        .def_readwrite("coord_center",  &WHAMBin::coord_center)
        .def_readwrite("count",         &WHAMBin::count)
        .def_readwrite("free_energy",   &WHAMBin::free_energy);
    
    py::class_<TIPoint>(m, "TIPoint", "Thermodynamic integration data point")
        .def(py::init<>())
        .def_readwrite("lambda",       &TIPoint::lambda)
        .def_readwrite("dV_dlambda",   &TIPoint::dV_dlambda);
    
    // ═══════════════════════════════════════════════════════════════════════
    // StatMechEngine: core thermodynamics calculator
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<StatMechEngine>(m, "StatMechEngine", 
        "Statistical mechanics engine for conformational ensembles")
        .def(py::init<double>(), 
            py::arg("temperature_K") = 300.0,
            "Initialize engine at given temperature (default 300K)")
        
        // ─── Ensemble construction ───
        .def("add_sample", &StatMechEngine::add_sample,
            py::arg("energy"), py::arg("multiplicity") = 1,
            "Add a sampled configuration with energy (kcal/mol) and multiplicity")
        .def("clear", &StatMechEngine::clear, "Remove all samples")
        
        // ─── Thermodynamic analysis ───
        .def("compute", &StatMechEngine::compute,
            "Compute full thermodynamics (F, S, H, Cv, etc.)")
        .def("boltzmann_weights", &StatMechEngine::boltzmann_weights,
            "Return Boltzmann weights for all samples (same order as insertion)")
        .def("delta_G", &StatMechEngine::delta_G,
            py::arg("reference"),
            "Compute ΔG relative to another ensemble")
        
        // ─── Advanced methods ───
        .def_static("init_replicas", &StatMechEngine::init_replicas,
            py::arg("temperatures"),
            "Initialize parallel tempering replicas at given temperatures")
        .def_static("attempt_swap", 
            [](Replica& a, Replica& b) {
                // Python RNG not compatible with std::mt19937, use C++ RNG
                thread_local std::mt19937 rng{std::random_device{}()};
                return StatMechEngine::attempt_swap(a, b, rng);
            },
            py::arg("replica_a"), py::arg("replica_b"),
            "Attempt Metropolis swap between two replicas (returns True if accepted)")
        .def_static("boltzmann_pmf", &StatMechEngine::boltzmann_pmf,
            py::arg("energies"), py::arg("coordinates"),
            py::arg("temperature"), py::arg("n_bins"),
            py::arg("max_iter") = 1000, py::arg("tolerance") = 1e-6,
            "Boltzmann-reweighted PMF along a reaction coordinate")
        .def_static("thermodynamic_integration", 
            &StatMechEngine::thermodynamic_integration,
            py::arg("points"),
            "Compute ΔG via thermodynamic integration (trapezoidal rule)")
        .def_static("helmholtz", &StatMechEngine::helmholtz,
            py::arg("energies"), py::arg("temperature"),
            "Compute Helmholtz free energy from raw energy vector")
        
        // ─── Properties ───
        .def_property_readonly("temperature", &StatMechEngine::temperature)
        .def_property_readonly("beta", &StatMechEngine::beta)
        .def_property_readonly("size", &StatMechEngine::size)
        .def("__len__", &StatMechEngine::size)
        .def("__repr__", [](const StatMechEngine& e) {
            return "<StatMechEngine T=" + std::to_string(e.temperature()) + 
                   "K n_samples=" + std::to_string(e.size()) + ">";
        });
    
    // ═══════════════════════════════════════════════════════════════════════
    // BoltzmannLUT: fast lookup table
    // ═══════════════════════════════════════════════════════════════════════
    
    py::class_<BoltzmannLUT>(m, "BoltzmannLUT",
        "Pre-tabulated Boltzmann factors for O(1) inner-loop evaluation")
        .def(py::init<double, double, double, int>(),
            py::arg("beta"), py::arg("e_min"), py::arg("e_max"),
            py::arg("n_bins") = 10000,
            "Initialize LUT for energy range [e_min, e_max]")
        .def("__call__", &BoltzmannLUT::operator(),
            py::arg("energy"),
            "Look up exp(-β E) for given energy");
    
#ifdef FLEXAIDS_FULL_GA
    // ═══════════════════════════════════════════════════════════════════════
    // BindingMode: pose cluster with thermodynamic scoring
    // (requires full GA infrastructure — only available in CMake build)
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<BindingMode>(m, "BindingMode",
        "Binding mode: cluster of docked poses with thermodynamic analysis")
        // Legacy interface (backward compatibility)
        .def("compute_energy", &BindingMode::compute_energy,
            "Helmholtz free energy F = H - TS (kcal/mol)")
        .def("compute_entropy", &BindingMode::compute_entropy,
            "Configurational entropy S (kcal mol⁻¹ K⁻¹)")
        .def("compute_enthalpy", &BindingMode::compute_enthalpy,
            "Boltzmann-weighted average energy ⟨E⟩ (kcal/mol)")

        // New StatMech API
        .def("get_thermodynamics", &BindingMode::get_thermodynamics,
            "Full thermodynamic properties (F, S, H, Cv, σ_E)")
        .def("get_free_energy", &BindingMode::get_free_energy,
            "Alias for compute_energy()")
        .def("get_heat_capacity", &BindingMode::get_heat_capacity,
            "Heat capacity Cv (kcal mol⁻¹ K⁻²)")
        .def("get_boltzmann_weights", &BindingMode::get_boltzmann_weights,
            "Boltzmann weights for all poses in this mode")
        .def("get_BindingMode_size", &BindingMode::get_BindingMode_size,
            "Number of poses in this binding mode")
        .def("get_pose", &BindingMode::get_pose,
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Access pose by index (bounds-checked)")
        .def("__len__", &BindingMode::get_BindingMode_size)
        .def("__repr__", [](const BindingMode& m) {
            auto thermo = m.get_thermodynamics();
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<BindingMode n_poses=%d F=%.3f H=%.3f S=%.5f>",
                m.get_BindingMode_size(),
                thermo.free_energy, thermo.mean_energy, thermo.entropy);
            return std::string(buf);
        });

    // ═══════════════════════════════════════════════════════════════════════
    // BindingPopulation: global ensemble thermodynamics
    // ═══════════════════════════════════════════════════════════════════════
    py::class_<BindingPopulation>(m, "BindingPopulation",
        "Collection of binding modes from a docking run, with global ensemble analysis")
        .def("get_population_size", &BindingPopulation::get_Population_size,
            "Number of distinct binding modes")
        .def("get_binding_mode",
            static_cast<const BindingMode& (BindingPopulation::*)(int) const>(
                &BindingPopulation::get_binding_mode),
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Access binding mode by index (bounds-checked)")
        .def("compute_delta_G",
            [](const BindingPopulation& pop, const BindingMode& m1, const BindingMode& m2) {
                return pop.compute_delta_G(m1, m2);
            },
            py::arg("mode1"), py::arg("mode2"),
            "ΔG between two binding modes (kcal/mol); positive = mode1 less favoured")
        .def("get_global_ensemble", &BindingPopulation::get_global_ensemble,
            "StatMechEngine aggregating all poses across all binding modes")
        .def("get_super_cluster_ensemble", &BindingPopulation::get_super_cluster_ensemble,
            "StatMechEngine with super-cluster filtered poses (dominant energy basin only)")
        .def("get_shannon_entropy", &BindingPopulation::get_shannon_entropy,
            "Population-level Shannon configurational entropy S = -kB * sum(p_i * ln(p_i)) (kcal/mol/K)")
        .def("get_deltaG_matrix", &BindingPopulation::get_deltaG_matrix,
            "Full ΔG matrix between all pairs of binding modes (kcal/mol); matrix[i][j] = F_i - F_j")
        .def("__len__", &BindingPopulation::get_Population_size)
        .def("__repr__", [](const BindingPopulation& p) {
            return "<BindingPopulation n_modes=" +
                   std::to_string(const_cast<BindingPopulation&>(p).get_Population_size()) + ">";
        });
#endif  // FLEXAIDS_FULL_GA

    // ═══════════════════════════════════════════════════════════════════════
    // ENCoM: normal mode analysis + vibrational entropy (Phase 3)
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<encom::NormalMode>(m, "NormalMode",
        "Normal mode from ENCoM elastic network calculation")
        .def(py::init([](int index, double eigenvalue, double frequency,
                         std::vector<double> eigenvector) {
            encom::NormalMode nm;
            nm.index = index;
            nm.eigenvalue = eigenvalue;
            nm.frequency = frequency;
            nm.eigenvector = std::move(eigenvector);
            return nm;
        }),
            py::arg("index") = 0,
            py::arg("eigenvalue") = 0.0,
            py::arg("frequency") = 0.0,
            py::arg("eigenvector") = std::vector<double>(),
            "Create a normal mode with optional initial values")
        .def_readwrite("index",      &encom::NormalMode::index,      "1-based mode index")
        .def_readwrite("eigenvalue", &encom::NormalMode::eigenvalue,  "λ (ENCoM arbitrary units)")
        .def_readwrite("frequency",  &encom::NormalMode::frequency,
            "sqrt(lambda) in model scale unless calibrated")
        .def_readwrite("eigenvector",&encom::NormalMode::eigenvector, "Displacement vector (3N)")
        .def("__repr__", [](const encom::NormalMode& nm) {
            return "<NormalMode " + std::to_string(nm.index) +
                   " λ=" + std::to_string(nm.eigenvalue) + ">";
        });

    py::class_<encom::VibrationalEntropy>(m, "VibrationalEntropy",
        "Quasi-harmonic vibrational entropy from ENCoM normal modes")
        .def(py::init([](double S_vib_kcal_mol_K, double S_vib_J_mol_K,
                         double omega_eff, int n_modes, double temperature,
                         double dG_vib_kcal_mol, double eigenvalue_to_omega,
                         bool calibrated, const std::string& calibration_label,
                         const std::string& calibration_provenance) {
            encom::VibrationalEntropy vs;
            vs.S_vib_kcal_mol_K = S_vib_kcal_mol_K;
            vs.S_vib_J_mol_K = S_vib_J_mol_K;
            vs.omega_eff = omega_eff;
            vs.n_modes = n_modes;
            vs.temperature = temperature;
            vs.eigenvalue_to_omega = eigenvalue_to_omega;
            vs.calibrated = calibrated;
            vs.calibration_label = calibration_label;
            vs.calibration_provenance = calibration_provenance;
            // dG_vib is stored implicitly as -T*S; if explicitly given, back-compute
            // (otherwise default to -T*S_vib)
            (void)dG_vib_kcal_mol; // used via property below
            return vs;
        }),
            py::arg("S_vib_kcal_mol_K") = 0.0,
            py::arg("S_vib_J_mol_K") = 0.0,
            py::arg("omega_eff") = 0.0,
            py::arg("n_modes") = 0,
            py::arg("temperature") = 300.0,
            py::arg("dG_vib_kcal_mol") = 0.0,
            py::arg("eigenvalue_to_omega") = 1.0,
            py::arg("calibrated") = false,
            py::arg("calibration_label") = "model-scale",
            py::arg("calibration_provenance") =
                "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied",
            "Create a VibrationalEntropy result")
        .def_readwrite("S_vib_kcal_mol_K", &encom::VibrationalEntropy::S_vib_kcal_mol_K,
            "S_vib in kcal mol⁻¹ K⁻¹")
        .def_readwrite("S_vib_J_mol_K",    &encom::VibrationalEntropy::S_vib_J_mol_K,
            "S_vib in J mol⁻¹ K⁻¹")
        .def_readwrite("omega_eff",        &encom::VibrationalEntropy::omega_eff,
            "Effective frequency ω_eff (rad/s)")
        .def_readwrite("n_modes",          &encom::VibrationalEntropy::n_modes,
            "Number of non-trivial normal modes (3N − 6)")
        .def_readwrite("temperature",      &encom::VibrationalEntropy::temperature, "K")
        .def_readwrite("eigenvalue_to_omega", &encom::VibrationalEntropy::eigenvalue_to_omega,
            "rad/s per sqrt(model eigenvalue)")
        .def_readwrite("calibrated", &encom::VibrationalEntropy::calibrated,
            "True when a physical/empirical frequency calibration was supplied")
        .def_readwrite("calibration_label", &encom::VibrationalEntropy::calibration_label)
        .def_readwrite("calibration_provenance", &encom::VibrationalEntropy::calibration_provenance)
        .def_property_readonly("absolute_claim_allowed",
            [](const encom::VibrationalEntropy& vs) {
                return vs.absolute_claim_allowed();
            })
        .def_property_readonly("dG_vib_kcal_mol", [](const encom::VibrationalEntropy& vs) {
            return -vs.temperature * vs.S_vib_kcal_mol_K;
        }, "−T·S_vib vibrational free energy correction (kcal/mol)")
        .def_property_readonly("TS_vib_kcal_mol", [](const encom::VibrationalEntropy& vs) {
            return vs.temperature * vs.S_vib_kcal_mol_K;
        }, "T·S_vib (kcal/mol)")
        .def("free_energy_correction", [](const encom::VibrationalEntropy& vs) {
            return -vs.temperature * vs.S_vib_kcal_mol_K;
        }, "−T·S_vib vibrational free energy correction (kcal/mol)")
        .def("__repr__", [](const encom::VibrationalEntropy& vs) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<VibrationalEntropy n_modes=%d S_vib=%.6f kcal/(mol·K) T=%.1fK>",
                vs.n_modes, vs.S_vib_kcal_mol_K, vs.temperature);
            return std::string(buf);
        });

    // ENCoMEngine wrapper that stores eigenvalue_cutoff for instance usage
    struct ENCoMEngineWrapper {
        double eigenvalue_cutoff;
        double eigenvalue_to_omega;
        bool calibrated;
        std::string calibration_label;
        std::string calibration_provenance;
        ENCoMEngineWrapper(
            double cutoff = 1e-6,
            double scale = 1.0,
            bool is_calibrated = false,
            const std::string& label = "model-scale",
            const std::string& provenance =
                "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied")
            : eigenvalue_cutoff(cutoff),
              eigenvalue_to_omega(scale),
              calibrated(is_calibrated),
              calibration_label(label),
              calibration_provenance(provenance) {}
    };

    py::class_<ENCoMEngineWrapper>(m, "ENCoMEngine",
        "ENCoM quasi-harmonic entropy calculator")
        .def(py::init<double, double, bool, const std::string&, const std::string&>(),
            py::arg("eigenvalue_cutoff") = 1e-6,
            py::arg("eigenvalue_to_omega") = 1.0,
            py::arg("calibrated") = false,
            py::arg("calibration_label") = "model-scale",
            py::arg("calibration_provenance") =
                "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied",
            "Initialize ENCoM engine with optional eigenvalue cutoff")
        .def_readwrite("eigenvalue_cutoff", &ENCoMEngineWrapper::eigenvalue_cutoff)
        .def_readwrite("eigenvalue_to_omega", &ENCoMEngineWrapper::eigenvalue_to_omega)
        .def_readwrite("calibrated", &ENCoMEngineWrapper::calibrated)
        .def_readwrite("calibration_label", &ENCoMEngineWrapper::calibration_label)
        .def_readwrite("calibration_provenance", &ENCoMEngineWrapper::calibration_provenance)
        .def("compute_vibrational_entropy",
            [](const ENCoMEngineWrapper& self,
               const std::vector<encom::NormalMode>& modes,
               double temperature) {
                encom::FrequencyCalibration calibration;
                calibration.eigenvalue_to_omega = self.eigenvalue_to_omega;
                calibration.calibrated = self.calibrated;
                calibration.label = self.calibration_label;
                calibration.provenance = self.calibration_provenance;
                return encom::ENCoMEngine::compute_vibrational_entropy(
                    modes, temperature, calibration, self.eigenvalue_cutoff);
            },
            py::arg("modes"),
            py::arg("temperature") = 300.0,
            "Compute quasi-harmonic S_vib from a list of NormalMode objects")
        .def_static("load_modes",
            &encom::ENCoMEngine::load_modes,
            py::arg("eigenvalue_file"), py::arg("eigenvector_file"),
            "Load normal modes from ENCoM output files")
        .def_static("total_entropy",
            &encom::ENCoMEngine::total_entropy,
            py::arg("S_conf_kcal_mol_K"), py::arg("S_vib_kcal_mol_K"),
            "S_total = S_conf + S_vib  (kcal mol⁻¹ K⁻¹)")
        .def_static("free_energy_with_vibrations",
            &encom::ENCoMEngine::free_energy_with_vibrations,
            py::arg("F_electronic"), py::arg("S_vib_kcal_mol_K"), py::arg("temperature_K"),
            "F_total = F_elec − T·S_vib  (kcal/mol)");

    // ═══════════════════════════════════════════════════════════════════════
    // FastOPTICS: lightweight super-cluster extraction
    // ═══════════════════════════════════════════════════════════════════════

    py::enum_<fast_optics::ClusterMode>(m, "ClusterMode",
        "Clustering mode for FastOPTICS extraction")
        .value("FULL_OPTICS", fast_optics::ClusterMode::FULL_OPTICS, "Full OPTICS hierarchy walk")
        .value("SUPER_CLUSTER_ONLY", fast_optics::ClusterMode::SUPER_CLUSTER_ONLY,
               "Fast super-cluster extraction (~40%% faster)")
        .export_values();

    py::class_<fast_optics::Point>(m, "Point", "N-dimensional point for OPTICS clustering")
        .def(py::init<>())
        .def(py::init([](std::vector<double> c) {
            fast_optics::Point p; p.coords = std::move(c); return p;
        }), py::arg("coords"))
        .def_readwrite("coords", &fast_optics::Point::coords);

    py::class_<fast_optics::Reachability>(m, "Reachability",
        "OPTICS ordering entry with point index and reachability distance")
        .def_readonly("index", &fast_optics::Reachability::index)
        .def_readonly("reach", &fast_optics::Reachability::reach);

    py::class_<fast_optics::FastOPTICS>(m, "FastOPTICSLight",
        "Lightweight FastOPTICS clustering with super-cluster extraction")
        .def(py::init<const std::vector<fast_optics::Point>&, int, double>(),
            py::arg("points"), py::arg("min_pts") = 4, py::arg("eps") = 0.0,
            "Build OPTICS ordering from points")
        .def("get_ordering", &fast_optics::FastOPTICS::getOrdering,
            py::return_value_policy::reference_internal,
            "Return the OPTICS reachability ordering")
        .def("extract_super_cluster", &fast_optics::FastOPTICS::extractSuperCluster,
            py::arg("mode") = fast_optics::ClusterMode::FULL_OPTICS,
            "Extract cluster indices (SUPER_CLUSTER_ONLY for fast mode)");

    // ──────────────────────────────────────────────────────────────────────
    // Spectrophore — 144-float 3D fingerprint
    // ──────────────────────────────────────────────────────────────────────

    m.attr("SPECTROPHORE_DESCRIPTOR_SIZE") = spectrophore::DESCRIPTOR_SIZE;
    m.attr("SPECTROPHORE_N_ANGULAR") = spectrophore::N_ANGULAR;
    m.attr("SPECTROPHORE_N_RADIAL") = spectrophore::N_RADIAL;
    m.attr("SPECTROPHORE_N_PROPERTIES") = spectrophore::N_PROPERTIES;

    py::class_<spectrophore::Spectrophore>(m, "Spectrophore",
        "144-float 3D molecular fingerprint (4 properties × 12 angular × 3 radial)")
        .def(py::init<>())
        .def("tanimoto", &spectrophore::Spectrophore::tanimoto,
            py::arg("other"),
            "Tanimoto similarity to another Spectrophore (0–1)")
        .def("euclidean", &spectrophore::Spectrophore::euclidean,
            py::arg("other"),
            "Euclidean distance to another Spectrophore")
        .def_property_readonly("values", [](const spectrophore::Spectrophore& s) {
            return py::array_t<float>({spectrophore::DESCRIPTOR_SIZE},
                                      {sizeof(float)},
                                      s.values, py::cast(s));
        }, "144-element descriptor array (read-only view)")
        .def_property_readonly("centroid", [](const spectrophore::Spectrophore& s) {
            return py::array_t<float>({3}, {sizeof(float)}, s.centroid, py::cast(s));
        }, "3D centroid used for computation");

    py::class_<spectrophore::SimpleAtom>(m, "SpectrophoreAtom",
        "Lightweight atom for spectrophore computation")
        .def(py::init<>())
        .def(py::init([](float x, float y, float z, float radius, float charge) {
            return spectrophore::SimpleAtom{x, y, z, radius, charge};
        }), py::arg("x"), py::arg("y"), py::arg("z"),
            py::arg("radius") = 1.7f, py::arg("charge") = 0.0f)
        .def_readwrite("x", &spectrophore::SimpleAtom::x)
        .def_readwrite("y", &spectrophore::SimpleAtom::y)
        .def_readwrite("z", &spectrophore::SimpleAtom::z)
        .def_readwrite("radius", &spectrophore::SimpleAtom::radius)
        .def_readwrite("charge", &spectrophore::SimpleAtom::charge);

    m.def("compute_spectrophore_from_atoms",
        [](py::array_t<float> coords, py::array_t<float> radii,
           py::array_t<float> charges, py::array_t<float> center) {
            auto c = coords.unchecked<2>();
            auto r = radii.unchecked<1>();
            auto q = charges.unchecked<1>();
            auto ctr = center.unchecked<1>();
            int n = static_cast<int>(c.shape(0));
            std::vector<spectrophore::SimpleAtom> atoms(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                atoms[static_cast<size_t>(i)] = {c(i, 0), c(i, 1), c(i, 2),
                                                  r(i), q(i)};
            }
            float cen[3] = {ctr(0), ctr(1), ctr(2)};
            return spectrophore::compute_from_atoms(atoms.data(), n, cen);
        },
        py::arg("coords"), py::arg("radii"), py::arg("charges"),
        py::arg("center"),
        "Compute spectrophore from atom coordinates, radii, and charges");

    // ──────────────────────────────────────────────────────────────────────
    // BindingResidues — key residue identification from MIF scores
    // ──────────────────────────────────────────────────────────────────────

    py::class_<binding_residues::ResidueContribution>(m, "ResidueContribution",
        "A protein residue's contribution to binding (from MIF analysis)")
        .def_readonly("res_index", &binding_residues::ResidueContribution::res_index)
        .def_property_readonly("name", [](const binding_residues::ResidueContribution& r) {
            return std::string(r.name);
        })
        .def_readonly("number", &binding_residues::ResidueContribution::number)
        .def_property_readonly("chain", [](const binding_residues::ResidueContribution& r) {
            return std::string(1, r.chain);
        })
        .def_readonly("mif_score", &binding_residues::ResidueContribution::mif_score)
        .def_readonly("contact_count", &binding_residues::ResidueContribution::contact_count)
        .def_readonly("min_distance", &binding_residues::ResidueContribution::min_distance)
        .def("__repr__", [](const binding_residues::ResidueContribution& r) {
            return std::string(r.name) + " " + std::to_string(r.number) +
                   ":" + r.chain + " (MIF=" +
                   std::to_string(r.mif_score) + ")";
        });
}
