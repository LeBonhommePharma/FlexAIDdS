// bindings_matrix.cpp — pybind11 bindings for the 256×256 soft contact matrix
//
// Exposes SoftContactMatrix, atom_typing_256, and ShannonMatrixScorer to Python
// with zero-copy numpy views for matrix data.
//
// Build: Compiled as part of _core extension when FLEXAIDS_USE_256_MATRIX=ON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../../LIB/atom_typing_256.h"
#include "../../LIB/soft_contact_matrix.h"
#include "../../LIB/shannon_matrix_scorer.h"

namespace py = pybind11;

// Forward declaration: this function is called from core_bindings.cpp
// to register the matrix submodule.  When FLEXAIDS_USE_256_MATRIX is defined,
// the linker includes this TU.

void register_matrix_bindings(py::module_& m) {
    // ═══════════════════════════════════════════════════════════════════════
    // atom_typing_256 — 8-bit encoding
    // ═══════════════════════════════════════════════════════════════════════

    auto m_at = m.def_submodule("atom256", "256-type atom encoding");

    // atom256::encode is overloaded (plus deleted guard overloads for stale
    // v55 call sites), so &atom256::encode is ambiguous for pybind11's
    // template deduction. Bind via a lambda pinned to the real signature.
    m_at.def("encode",
        [](uint8_t base_type, bool donor, bool acceptor) {
            return atom256::encode(base_type, donor, acceptor);
        },
        py::arg("base_type"), py::arg("donor"), py::arg("acceptor"),
        "Encode 8-bit atom type from (base_type, donor, acceptor)");

    m_at.def("get_base", &atom256::get_base,
        py::arg("code"), "Extract base type (bits 0-5)");
    m_at.def("get_donor", &atom256::get_hbond_donor,
        py::arg("code"), "Extract H-bond donor flag (bit 7)");
    m_at.def("get_acceptor", &atom256::get_hbond_acceptor,
        py::arg("code"), "Extract H-bond acceptor flag (bit 6)");
    m_at.def("get_hbond", &atom256::get_hbond,
        py::arg("code"), "True if donor or acceptor (bits 6-7)");

    m_at.def("sybyl_to_base", &atom256::sybyl_to_base,
        py::arg("sybyl_type"),
        "Map SYBYL type (1-40) to base type (0-63)");
    m_at.def("base_to_sybyl_parent", &atom256::base_to_sybyl_parent,
        py::arg("base_type"),
        "Map base type (0-63) to SYBYL parent (1-40)");
    m_at.def("base_type_name", &atom256::base_type_name,
        py::arg("base"), "Human-readable name for base type");
    m_at.def("encode_from_sybyl", &atom256::encode_from_sybyl,
        py::arg("sybyl_type"), py::arg("partial_charge"),
        py::arg("n_hydrogens"),
        py::arg("has_heteroatom_neighbor") = false,
        py::arg("is_bridgehead") = false,
        "Full encoding from SYBYL type + charge + structural context");

    m_at.attr("BASE_TYPE_COUNT") = py::int_(static_cast<int>(atom256::BASE_TYPE_COUNT));

    // ═══════════════════════════════════════════════════════════════════════
    // SoftContactMatrix — 256×256 interaction matrix
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<scm::SoftContactMatrix>(m, "SoftContactMatrix",
        "256×256 cache-aligned soft contact interaction matrix (256 KB)")
        .def(py::init([]() {
            auto mat = new scm::SoftContactMatrix();
            mat->zero();
            return mat;
        }), "Create zero-initialised matrix")

        .def("lookup", &scm::SoftContactMatrix::lookup,
            py::arg("type_i"), py::arg("type_j"),
            "O(1) pairwise interaction energy lookup")

        .def("set", &scm::SoftContactMatrix::set,
            py::arg("type_i"), py::arg("type_j"), py::arg("value"),
            "Set interaction energy for a pair")

        // Zero-copy numpy view
        .def("as_numpy", [](scm::SoftContactMatrix& self) {
            return py::array_t<float>(
                {scm::MATRIX_DIM, scm::MATRIX_DIM},
                {scm::MATRIX_DIM * sizeof(float), sizeof(float)},
                self.data,
                py::cast(self)  // prevent dealloc while view is alive
            );
        }, py::return_value_policy::reference_internal,
        "Return matrix as writable numpy float32 array (256×256, zero-copy)")

        // Load from numpy
        .def("from_numpy", [](scm::SoftContactMatrix& self,
                               py::array_t<float, py::array::c_style |
                                                   py::array::forcecast> arr) {
            auto buf = arr.request();
            if (buf.size != scm::MATRIX_SIZE)
                throw std::runtime_error("Expected 65536 elements");
            std::memcpy(self.data, buf.ptr, scm::MATRIX_BYTES);
        }, py::arg("array"), "Load matrix from numpy float32 array")

        .def("symmetrise", &scm::SoftContactMatrix::symmetrise,
            "Force symmetry: M = (M + M^T) / 2")

        // Batch scoring
        .def("score_contacts",
            [](const scm::SoftContactMatrix& self,
               py::array_t<uint8_t> type_a,
               py::array_t<uint8_t> type_b,
               py::array_t<float> areas) {
                auto ba = type_a.request();
                auto bb = type_b.request();
                auto ba2 = areas.request();
                int n = std::min({(int)ba.size, (int)bb.size, (int)ba2.size});
                return self.score_contacts(
                    static_cast<const uint8_t*>(ba.ptr),
                    static_cast<const uint8_t*>(bb.ptr),
                    static_cast<const float*>(ba2.ptr), n);
            },
            py::arg("type_a"), py::arg("type_b"), py::arg("areas"),
            "Score N contacts: sum of matrix[a[k]][b[k]] * area[k] (AVX2 when available)")

        // Binary I/O
        .def("save", &scm::SoftContactMatrix::save,
            py::arg("path"), "Save to binary file (SHNN header)")
        .def("load", &scm::SoftContactMatrix::load,
            py::arg("path"), "Load from binary file (SHNN header)")

        // 256→40 projection
        .def("project_to_40x40", [](const scm::SoftContactMatrix& self) {
            auto proj = self.project_to_40x40();
            return py::array_t<float>(
                {40, 40},
                {40 * sizeof(float), sizeof(float)},
                proj.data()
            );
        }, "Project 256×256 → 40×40 SYBYL matrix (returns numpy float32)")

        .def("__repr__", [](const scm::SoftContactMatrix&) {
            return "<SoftContactMatrix 256×256 (256 KB aligned)>";
        });

    // ═══════════════════════════════════════════════════════════════════════
    // FastOPTICS super-cluster detection
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<scm::FOPTICSResult>(m, "FOPTICSResult",
        "FastOPTICS clustering result on matrix row vectors")
        .def_readonly("order", &scm::FOPTICSResult::order)
        .def_readonly("reachability", &scm::FOPTICSResult::reachability)
        .def_readonly("cluster_labels", &scm::FOPTICSResult::cluster_labels)
        .def_readonly("n_clusters", &scm::FOPTICSResult::n_clusters)
        .def("__repr__", [](const scm::FOPTICSResult& r) {
            return "<FOPTICSResult n_clusters=" + std::to_string(r.n_clusters) + ">";
        });

    m.def("find_super_clusters",
        [](const scm::SoftContactMatrix& mat, int min_pts,
           int n_projections, uint32_t seed) {
            return scm::find_super_clusters(mat, min_pts, n_projections, seed);
        },
        py::arg("matrix"),
        py::arg("min_pts") = 5,
        py::arg("n_projections") = 20,
        py::arg("seed") = 42,
        "Run FastOPTICS on matrix row vectors to detect super-clusters");

    m.def("apply_supercluster_bias",
        [](scm::SoftContactMatrix& mat, const scm::FOPTICSResult& clusters,
           float alpha, float sigma) {
            scm::apply_supercluster_bias(mat, clusters, alpha, sigma);
        },
        py::arg("matrix"), py::arg("clusters"),
        py::arg("alpha") = 0.3f, py::arg("sigma") = 2.0f,
        "Apply Gaussian supercluster bias modulation to matrix");

    // ═══════════════════════════════════════════════════════════════════════
    // ShannonMatrixScorer — two-term scoring engine
    // ═══════════════════════════════════════════════════════════════════════

    py::class_<scorer::Contact>(m, "Contact",
        "Single atom-atom contact for scoring")
        .def(py::init([](uint8_t ta, uint8_t tb, float area, float dist,
                         float ra, float rb, float qa, float qb) {
            return scorer::Contact{ta, tb, area, dist, ra, rb, qa, qb};
        }),
            py::arg("type_a"), py::arg("type_b"),
            py::arg("area") = 1.0f, py::arg("distance") = 3.0f,
            py::arg("radius_a") = 1.7f, py::arg("radius_b") = 1.7f,
            py::arg("charge_a") = 0.0f, py::arg("charge_b") = 0.0f)
        .def_readwrite("type_a", &scorer::Contact::type_a)
        .def_readwrite("type_b", &scorer::Contact::type_b)
        .def_readwrite("area", &scorer::Contact::area)
        .def_readwrite("distance", &scorer::Contact::distance)
        .def_readwrite("charge_a", &scorer::Contact::charge_a)
        .def_readwrite("charge_b", &scorer::Contact::charge_b);

    py::class_<scorer::PoseScore>(m, "PoseScore",
        "Result of scoring a single pose")
        .def_readonly("matrix_score", &scorer::PoseScore::matrix_score)
        .def_readonly("lj_score", &scorer::PoseScore::lj_score)
        .def_readonly("elec_score", &scorer::PoseScore::elec_score)
        .def_readonly("total_score", &scorer::PoseScore::total_score)
        .def_readonly("survived_filter", &scorer::PoseScore::survived_filter);

    py::class_<scorer::EnsembleResult>(m, "EnsembleResult",
        "Result of scoring a pose ensemble")
        .def_readonly("deltaG", &scorer::EnsembleResult::deltaG)
        .def_readonly("shannonEntropy", &scorer::EnsembleResult::shannonEntropy)
        .def_readonly("meanScore", &scorer::EnsembleResult::meanScore)
        .def_readonly("n_survivors", &scorer::EnsembleResult::n_survivors)
        .def_readonly("n_superclusters", &scorer::EnsembleResult::n_superclusters)
        .def_readonly("scores", &scorer::EnsembleResult::scores);

    py::class_<scorer::ShannonMatrixScorer>(m, "ShannonMatrixScorer",
        "Two-term scoring: matrix pre-filter → LJ+Coulomb → Shannon entropy → ΔG")
        .def(py::init<const scm::SoftContactMatrix&, double, float>(),
            py::arg("matrix"),
            py::arg("temperature_K") = 298.15,
            py::arg("filter_threshold") = 0.0f)
        .def("set_temperature", &scorer::ShannonMatrixScorer::set_temperature)
        .def("set_filter_threshold", &scorer::ShannonMatrixScorer::set_filter_threshold)
        .def("set_dielectric", &scorer::ShannonMatrixScorer::set_dielectric)
        .def("matrix_score", &scorer::ShannonMatrixScorer::matrix_score,
            py::arg("contacts"), "Fast matrix-only score")
        .def("analytic_score", &scorer::ShannonMatrixScorer::analytic_score,
            py::arg("contacts"), "LJ + Coulomb analytic score")
        .def("score_pose", &scorer::ShannonMatrixScorer::score_pose,
            py::arg("contacts"), "Full two-stage pose scoring")
        .def("score_ensemble", &scorer::ShannonMatrixScorer::score_ensemble,
            py::arg("pose_contacts"), "Score ensemble with Shannon entropy")
        .def("__repr__", [](const scorer::ShannonMatrixScorer&) {
            return "<ShannonMatrixScorer two-term 256×256>";
        });

    // Constants
    m.attr("MATRIX_DIM") = py::int_(scm::MATRIX_DIM);
    m.attr("MATRIX_SIZE") = py::int_(scm::MATRIX_SIZE);
}
