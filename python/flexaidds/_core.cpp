// _core.cpp — pybind11 bindings for FlexAID∆S C++ core
//
// Exposes all standalone thermodynamic modules:
//   - statmech::StatMechEngine, Thermodynamics, BoltzmannLUT
//   - statmech data structures: State, Replica, WHAMBin, TIPoint
//   - encom::ENCoMEngine, VibrationalEntropy, NormalMode
//   - Physical constants (kB_kcal, kB_SI)
//
// Note: BindingMode / BindingPopulation require the full GA infrastructure
// (gaboom.h, fileio.h) and are not exposed here.  Use the CMake build with
// -DBUILD_PYTHON_BINDINGS=ON for the extended module that includes them.
//
// Build: pip install -e python/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "statmech.h"
#include "encom.h"
#include "tENCoM/tencm.h"
#include "ShannonThermoStack/ShannonThermoStack.h"

namespace py = pybind11;
using namespace statmech;

PYBIND11_MODULE(_core, m) {
    m.doc() = "FlexAID\u0394S C++ core: statistical mechanics, "
              "Boltzmann lookup tables, and ENCoM vibrational entropy";

    // Physical constants
    m.attr("kB_kcal") = kB_kcal;
    m.attr("kB_SI")   = kB_SI;

    // ═══════════════════════════════════════════════════════════════════════
    // Thermodynamics data structures
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<State>(m, "State", "Microstate with energy and degeneracy")
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
        .def_readwrite("mean_energy",    &Thermodynamics::mean_energy,    "<E> (kcal/mol)")
        .def_readwrite("mean_energy_sq", &Thermodynamics::mean_energy_sq, "<E^2>")
        .def_readwrite("heat_capacity",  &Thermodynamics::heat_capacity,  "Cv (kcal/mol/K^2)")
        .def_readwrite("entropy",        &Thermodynamics::entropy,        "S (kcal/mol/K)")
        .def_readwrite("std_energy",     &Thermodynamics::std_energy,     "sigma_E (kcal/mol)")
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

    // Note: TIPoint.lambda is accessed as "lambda_val" to avoid Python keyword conflict
    py::class_<TIPoint>(m, "TIPoint", "Thermodynamic integration data point")
        .def(py::init<>())
        .def_readwrite("lambda_val",  &TIPoint::lambda)
        .def_readwrite("dV_dlambda",  &TIPoint::dV_dlambda);

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
            "Compute full thermodynamics (F, S, H, Cv, etc.)",
            py::call_guard<py::gil_scoped_release>())
        .def("boltzmann_weights", &StatMechEngine::boltzmann_weights,
            "Return Boltzmann weights for all samples (same order as insertion)",
            py::call_guard<py::gil_scoped_release>())
        .def("delta_G", &StatMechEngine::delta_G,
            py::arg("reference"),
            "Compute delta-G relative to another ensemble",
            py::call_guard<py::gil_scoped_release>())

        // ─── Advanced methods ───
        // Note: std::span parameters require lambda wrappers for pybind11
        .def_static("init_replicas",
            [](const std::vector<double>& temps) {
                return StatMechEngine::init_replicas(temps);
            },
            py::arg("temperatures"),
            "Initialize parallel tempering replicas at given temperatures")
        .def_static("attempt_swap",
            [](Replica& a, Replica& b) {
                thread_local std::mt19937 rng{std::random_device{}()};
                return StatMechEngine::attempt_swap(a, b, rng);
            },
            py::arg("replica_a"), py::arg("replica_b"),
            "Attempt Metropolis swap between two replicas (returns True if accepted)")
        .def_static("boltzmann_pmf",
            [](const std::vector<double>& energies,
               const std::vector<double>& coordinates,
               double temperature, int n_bins, int max_iter, double tolerance) {
                py::gil_scoped_release release;
                return StatMechEngine::boltzmann_pmf(energies, coordinates,
                    temperature, n_bins, max_iter, tolerance);
            },
            py::arg("energies"), py::arg("coordinates"),
            py::arg("temperature"), py::arg("n_bins"),
            py::arg("max_iter") = 1000, py::arg("tolerance") = 1e-6,
            "Boltzmann-reweighted PMF along a reaction coordinate")
        .def_static("thermodynamic_integration",
            [](const std::vector<TIPoint>& points) {
                py::gil_scoped_release release;
                return StatMechEngine::thermodynamic_integration(points);
            },
            py::arg("points"),
            "Compute delta-G via thermodynamic integration (trapezoidal rule)")
        .def_static("helmholtz",
            [](const std::vector<double>& energies, double temperature) {
                py::gil_scoped_release release;
                return StatMechEngine::helmholtz(energies, temperature);
            },
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
            "Look up exp(-beta * E) for given energy");

    // ═══════════════════════════════════════════════════════════════════════
    // ENCoM: normal mode analysis + vibrational entropy
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
        .def_readwrite("index",       &encom::NormalMode::index,       "1-based mode index")
        .def_readwrite("eigenvalue",  &encom::NormalMode::eigenvalue,  "eigenvalue (ENCoM units)")
        .def_readwrite("frequency",   &encom::NormalMode::frequency,   "omega = sqrt(eigenvalue)")
        .def_readwrite("eigenvector", &encom::NormalMode::eigenvector, "Displacement vector (3N)")
        .def("__repr__", [](const encom::NormalMode& nm) {
            return "<NormalMode " + std::to_string(nm.index) +
                   " eigenvalue=" + std::to_string(nm.eigenvalue) + ">";
        });

    py::class_<encom::VibrationalEntropy>(m, "VibrationalEntropy",
        "Quasi-harmonic vibrational entropy from ENCoM normal modes")
        .def(py::init([](double S_vib_kcal_mol_K, double S_vib_J_mol_K,
                         double omega_eff, int n_modes, double temperature,
                         double dG_vib_kcal_mol) {
            encom::VibrationalEntropy vs;
            vs.S_vib_kcal_mol_K = S_vib_kcal_mol_K;
            vs.S_vib_J_mol_K = S_vib_J_mol_K;
            vs.omega_eff = omega_eff;
            vs.n_modes = n_modes;
            vs.temperature = temperature;
            (void)dG_vib_kcal_mol;
            return vs;
        }),
            py::arg("S_vib_kcal_mol_K") = 0.0,
            py::arg("S_vib_J_mol_K") = 0.0,
            py::arg("omega_eff") = 0.0,
            py::arg("n_modes") = 0,
            py::arg("temperature") = 300.0,
            py::arg("dG_vib_kcal_mol") = 0.0,
            "Create a VibrationalEntropy result")
        .def_readwrite("S_vib_kcal_mol_K", &encom::VibrationalEntropy::S_vib_kcal_mol_K,
            "S_vib in kcal/(mol*K)")
        .def_readwrite("S_vib_J_mol_K",    &encom::VibrationalEntropy::S_vib_J_mol_K,
            "S_vib in J/(mol*K)")
        .def_readwrite("omega_eff",        &encom::VibrationalEntropy::omega_eff,
            "Effective frequency omega_eff (rad/s)")
        .def_readwrite("n_modes",          &encom::VibrationalEntropy::n_modes,
            "Number of non-trivial normal modes (3N - 6)")
        .def_readwrite("temperature",      &encom::VibrationalEntropy::temperature, "K")
        .def_property_readonly("dG_vib_kcal_mol", [](const encom::VibrationalEntropy& vs) {
            return -vs.temperature * vs.S_vib_kcal_mol_K;
        }, "-T*S_vib vibrational free energy correction (kcal/mol)")
        .def_property_readonly("TS_vib_kcal_mol", [](const encom::VibrationalEntropy& vs) {
            return vs.temperature * vs.S_vib_kcal_mol_K;
        }, "T*S_vib (kcal/mol)")
        .def("free_energy_correction", [](const encom::VibrationalEntropy& vs) {
            return -vs.temperature * vs.S_vib_kcal_mol_K;
        }, "-T*S_vib vibrational free energy correction (kcal/mol)")
        .def("__repr__", [](const encom::VibrationalEntropy& vs) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<VibrationalEntropy n_modes=%d S_vib=%.6f kcal/(mol*K) T=%.1fK>",
                vs.n_modes, vs.S_vib_kcal_mol_K, vs.temperature);
            return std::string(buf);
        });

    // ENCoMEngine wrapper that stores eigenvalue_cutoff for instance usage
    struct ENCoMEngineWrapper {
        double eigenvalue_cutoff;
        ENCoMEngineWrapper(double cutoff = 1e-6) : eigenvalue_cutoff(cutoff) {}
    };

    py::class_<ENCoMEngineWrapper>(m, "ENCoMEngine",
        "ENCoM quasi-harmonic entropy calculator")
        .def(py::init<double>(),
            py::arg("eigenvalue_cutoff") = 1e-6,
            "Initialize ENCoM engine with optional eigenvalue cutoff")
        .def_readwrite("eigenvalue_cutoff", &ENCoMEngineWrapper::eigenvalue_cutoff)
        .def("compute_vibrational_entropy",
            [](const ENCoMEngineWrapper& self,
               const std::vector<encom::NormalMode>& modes,
               double temperature) {
                py::gil_scoped_release release;
                return encom::ENCoMEngine::compute_vibrational_entropy(
                    modes, temperature, self.eigenvalue_cutoff);
            },
            py::arg("modes"),
            py::arg("temperature") = 300.0,
            "Compute quasi-harmonic S_vib from a list of NormalMode objects")
        .def_static("load_modes",
            &encom::ENCoMEngine::load_modes,
            py::arg("eigenvalue_file"), py::arg("eigenvector_file"),
            "Load normal modes from ENCoM output files",
            py::call_guard<py::gil_scoped_release>())
        .def_static("total_entropy",
            &encom::ENCoMEngine::total_entropy,
            py::arg("S_conf_kcal_mol_K"), py::arg("S_vib_kcal_mol_K"),
            "S_total = S_conf + S_vib  (kcal/(mol*K))")
        .def_static("free_energy_with_vibrations",
            &encom::ENCoMEngine::free_energy_with_vibrations,
            py::arg("F_electronic"), py::arg("S_vib_kcal_mol_K"), py::arg("temperature_K"),
            "F_total = F_elec - T*S_vib  (kcal/mol)");

    // ═══════════════════════════════════════════════════════════════════════
    // TorsionalENM: backbone flexibility via torsional normal modes
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<tencm::NormalMode>(m, "TorsionalNormalMode",
        "Torsional normal mode from TorsionalENM")
        .def(py::init<>())
        .def_readwrite("eigenvalue",  &tencm::NormalMode::eigenvalue,
            "Stiffness in kcal/(mol*rad^2)")
        .def_readwrite("eigenvector", &tencm::NormalMode::eigenvector,
            "Displacement over torsion DOFs");

    py::class_<tencm::Conformer>(m, "Conformer",
        "Perturbed backbone conformation from TorsionalENM sampling")
        .def_readonly("delta_theta",    &tencm::Conformer::delta_theta)
        .def_readonly("strain_energy",  &tencm::Conformer::strain_energy);

    py::class_<tencm::TorsionalENM>(m, "TorsionalENM",
        "Torsional elastic network model for backbone flexibility")
        .def(py::init<>())
        .def("build_from_pdb",
            [](tencm::TorsionalENM& self, const std::string& pdb_path,
               float cutoff, float k0) {
                py::gil_scoped_release release;
                // Parse PDB for Cα coordinates, build model
                // Read PDB file, extract atoms/residues, then build
                std::vector<std::array<float,3>> ca_coords;
                std::ifstream ifs(pdb_path);
                if (!ifs.is_open())
                    throw std::runtime_error("Cannot open PDB: " + pdb_path);

                std::string line;
                while (std::getline(ifs, line)) {
                    if (line.size() >= 54 && line.substr(0,4) == "ATOM") {
                        std::string aname = line.substr(12, 4);
                        // Strip spaces
                        aname.erase(0, aname.find_first_not_of(' '));
                        aname.erase(aname.find_last_not_of(' ') + 1);
                        if (aname == "CA") {
                            float x = std::stof(line.substr(30, 8));
                            float y = std::stof(line.substr(38, 8));
                            float z = std::stof(line.substr(46, 8));
                            ca_coords.push_back({x, y, z});
                        }
                    }
                }
                if (ca_coords.size() < 4)
                    throw std::runtime_error("Need at least 4 Cα atoms, found " +
                                             std::to_string(ca_coords.size()));

                // Build directly from Cα coordinates using internal methods
                // We use a simplified path: create minimal atom/residue stubs
                int n = static_cast<int>(ca_coords.size());
                std::vector<atom> atoms(n);
                std::vector<resid> residues(n + 1); // 1-based indexing

                for (int i = 0; i < n; ++i) {
                    atoms[i].coor[0] = ca_coords[i][0];
                    atoms[i].coor[1] = ca_coords[i][1];
                    atoms[i].coor[2] = ca_coords[i][2];
                    strncpy(atoms[i].name, " CA ", 5);
                    atoms[i].number = i + 1;

                    residues[i+1].fatm = (int*)malloc(sizeof(int));
                    residues[i+1].latm = (int*)malloc(sizeof(int));
                    residues[i+1].fatm[0] = i;
                    residues[i+1].latm[0] = i;
                    strncpy(residues[i+1].name, "ALA", 4);
                }

                self.build(atoms.data(), residues.data() + 1, n, cutoff, k0);

                // Clean up malloc'd fatm/latm
                for (int i = 0; i < n; ++i) {
                    free(residues[i+1].fatm);
                    free(residues[i+1].latm);
                }
            },
            py::arg("pdb_path"),
            py::arg("cutoff") = 9.0f,
            py::arg("k0") = 1.0f,
            "Build torsional ENM from a PDB file (reads Cα atoms)")
        .def("is_built", &tencm::TorsionalENM::is_built)
        .def("n_residues", &tencm::TorsionalENM::n_residues)
        .def("n_bonds", &tencm::TorsionalENM::n_bonds)
        .def("modes", &tencm::TorsionalENM::modes,
            py::return_value_policy::reference_internal)
        .def("bfactors", &tencm::TorsionalENM::bfactors,
            py::arg("temperature") = 300.0f,
            "Predicted Cα B-factors at given temperature (Å²)",
            py::call_guard<py::gil_scoped_release>())
        .def("__repr__", [](const tencm::TorsionalENM& t) {
            return "<TorsionalENM n_res=" + std::to_string(t.n_residues()) +
                   " n_bonds=" + std::to_string(t.n_bonds()) +
                   " built=" + (t.is_built() ? "true" : "false") + ">";
        });

    // ═══════════════════════════════════════════════════════════════════════
    // ShannonThermoStack: combined entropy pipeline
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<shannon_thermo::FullThermoResult>(m, "FullThermoResult",
        "Combined thermodynamic result from ShannonThermoStack")
        .def_readonly("deltaG",              &shannon_thermo::FullThermoResult::deltaG)
        .def_readonly("shannonEntropy",      &shannon_thermo::FullThermoResult::shannonEntropy)
        .def_readonly("torsionalVibEntropy", &shannon_thermo::FullThermoResult::torsionalVibEntropy)
        .def_readonly("entropyContribution", &shannon_thermo::FullThermoResult::entropyContribution)
        .def_readonly("report",              &shannon_thermo::FullThermoResult::report)
        .def("__repr__", [](const shannon_thermo::FullThermoResult& r) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "<FullThermoResult dG=%.3f shannon=%.3f vib=%.6f>",
                r.deltaG, r.shannonEntropy, r.torsionalVibEntropy);
            return std::string(buf);
        });

    m.def("compute_shannon_entropy",
        &shannon_thermo::compute_shannon_entropy,
        py::arg("values"), py::arg("num_bins") = 20,
        "Compute Shannon entropy of continuous values (bits)",
        py::call_guard<py::gil_scoped_release>());

    m.def("compute_torsional_vibrational_entropy",
        &shannon_thermo::compute_torsional_vibrational_entropy,
        py::arg("modes"), py::arg("temperature_K") = 298.15,
        "Compute torsional vibrational entropy from normal modes",
        py::call_guard<py::gil_scoped_release>());

    m.def("run_shannon_thermo_stack",
        &shannon_thermo::run_shannon_thermo_stack,
        py::arg("stat_engine"), py::arg("tencm_model"),
        py::arg("base_deltaG"), py::arg("temperature_K") = 298.15,
        "Run full ShannonThermoStack pipeline",
        py::call_guard<py::gil_scoped_release>());
}
