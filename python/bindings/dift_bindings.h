// dift_bindings.h — pybind11 bindings for the DiFT torsional engine
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Shared by both _core builds:
//   - python/flexaidds/_core.cpp        (setup.py / pip install -e)
//   - python/bindings/core_bindings.cpp (CMake -DBUILD_PYTHON_BINDINGS=ON)
//
// Exposes dift::DiFTEngine and friends. The pure-Python flexaidds.dift module
// mirrors this API and is the always-available fallback; these bindings are a
// speed path that must produce numerically identical results.
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdio>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "../../LIB/DiFT/DiFT.h"

// Registers the DiFT classes/functions on an existing pybind11 module.
inline void register_dift_bindings(pybind11::module_& m) {
    namespace py = pybind11;
    using namespace dift;

    m.attr("dift_kB_kcal") = kB_kcal;

    // ── FourierTerm ──────────────────────────────────────────────────────────
    py::class_<FourierTerm>(m, "FourierTerm",
        "One cosine term: amplitude * cos(multiplicity*phi - phase)")
        .def(py::init<>())
        .def_readwrite("multiplicity", &FourierTerm::multiplicity, "n - integer frequency")
        .def_readwrite("amplitude",    &FourierTerm::amplitude,    "A_n (kcal/mol)")
        .def_readwrite("phase",        &FourierTerm::phase,        "omega_n (radians)")
        .def_readwrite("power",        &FourierTerm::power,        "A_n^2 spectral power")
        .def("__repr__", [](const FourierTerm& t) {
            char buf[96];
            std::snprintf(buf, sizeof(buf), "<FourierTerm n=%d A=%.4g w=%.4g>",
                          t.multiplicity, t.amplitude, t.phase);
            return std::string(buf);
        });

    // ── TorsionalPotential ───────────────────────────────────────────────────
    py::class_<TorsionalPotential>(m, "TorsionalPotential",
        "Analytical torsional potential V(phi) = mean + sum A_n cos(n phi - w_n)")
        .def(py::init<>())
        .def_readwrite("terms",            &TorsionalPotential::terms)
        .def_readwrite("mean",             &TorsionalPotential::mean)
        .def_readwrite("v_min",            &TorsionalPotential::v_min)
        .def_readwrite("n_samples",        &TorsionalPotential::n_samples)
        .def_readwrite("r_squared",        &TorsionalPotential::r_squared)
        .def_readwrite("spectral_entropy", &TorsionalPotential::spectral_entropy)
        .def_readwrite("effective_modes",  &TorsionalPotential::effective_modes)
        .def_readwrite("refinement_iters", &TorsionalPotential::refinement_iters)
        .def("evaluate", &TorsionalPotential::evaluate, py::arg("phi"),
             "Absolute potential at angle phi (radians)")
        .def("relative", &TorsionalPotential::relative, py::arg("phi"),
             "Potential relative to its global minimum")
        .def("sample", &TorsionalPotential::sample, py::arg("n"),
             "Sample V(phi) on n uniform points over [0, 2pi)")
        .def_property_readonly("n_terms", &TorsionalPotential::n_terms)
        .def("__repr__", [](const TorsionalPotential& p) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "<TorsionalPotential terms=%zu R2=%.4f N_eff=%.3f>",
                p.n_terms(), p.r_squared, p.effective_modes);
            return std::string(buf);
        });

    // ── TorsionalThermo ──────────────────────────────────────────────────────
    py::class_<TorsionalThermo>(m, "TorsionalThermo",
        "Per-bond torsional thermodynamics (excess, vs a free rotor)")
        .def(py::init<>())
        .def_readwrite("temperature",        &TorsionalThermo::temperature)
        .def_readwrite("partition_function", &TorsionalThermo::partition_function)
        .def_readwrite("free_energy",        &TorsionalThermo::free_energy)
        .def_readwrite("mean_energy",        &TorsionalThermo::mean_energy)
        .def_readwrite("entropy",            &TorsionalThermo::entropy)
        .def_readwrite("minus_TS",           &TorsionalThermo::minus_TS)
        .def("__repr__", [](const TorsionalThermo& t) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<TorsionalThermo T=%.1fK z=%.4f F=%.4f S=%.6f -TS=%.4f>",
                t.temperature, t.partition_function, t.free_energy,
                t.entropy, t.minus_TS);
            return std::string(buf);
        });

    // ── DiFTEngine ───────────────────────────────────────────────────────────
    py::class_<DiFTEngine>(m, "DiFTEngine",
        "Discrete Fourier Transform torsional parametrization engine")
        .def(py::init<double>(), py::arg("temperature_K") = 300.0)
        .def("transform",
            [](const DiFTEngine& self, const std::vector<double>& profile) {
                double mean = 0.0;
                auto spec = self.transform(std::span<const double>(profile), mean);
                return std::make_tuple(spec, mean);
            },
            py::arg("profile"),
            "Forward DiFT -> (spectrum, mean)")
        .def("parametrize",
            [](const DiFTEngine& self, const std::vector<double>& profile,
               int max_multiplicity) {
                return self.parametrize(std::span<const double>(profile),
                                        max_multiplicity);
            },
            py::arg("profile"), py::arg("max_multiplicity") = 0,
            "Forward transform + Shannon-collapse spectral truncation")
        .def("refine",
            [](const DiFTEngine& self, const std::vector<double>& qm,
               const std::vector<double>& mm_initial, double lambda,
               double r2_target, int max_iter, int max_multiplicity) {
                return self.refine(std::span<const double>(qm),
                                   std::span<const double>(mm_initial),
                                   lambda, r2_target, max_iter, max_multiplicity);
            },
            py::arg("qm"), py::arg("mm_initial"), py::arg("lambda_") = 0.5,
            py::arg("r2_target") = 0.98, py::arg("max_iter") = 50,
            py::arg("max_multiplicity") = 6,
            "Iterative QM-MM refinement (V_{i+1} = V_i + lambda*D_i)")
        .def("boltzmann_invert",
            [](const DiFTEngine& self, const std::vector<double>& histogram) {
                return self.boltzmann_invert(std::span<const double>(histogram));
            },
            py::arg("histogram"),
            "Boltzmann-invert a CG dihedral histogram into an energy profile")
        .def("thermodynamics", &DiFTEngine::thermodynamics, py::arg("potential"),
            "Per-bond torsional thermodynamics from a parametrized potential")
        .def_static("circular_mean",
            [](const std::vector<double>& angles) {
                return DiFTEngine::circular_mean(std::span<const double>(angles));
            },
            py::arg("angles"), "Directional mean of phase offsets")
        .def_static("r_squared",
            [](const std::vector<double>& observed,
               const std::vector<double>& model) {
                return DiFTEngine::r_squared(std::span<const double>(observed),
                                             std::span<const double>(model));
            },
            py::arg("observed"), py::arg("model"),
            "Coefficient of determination R^2")
        .def_property_readonly("temperature", &DiFTEngine::temperature)
        .def("set_temperature", &DiFTEngine::set_temperature, py::arg("T_K"))
        .def("__repr__", [](const DiFTEngine& e) {
            return "<DiFTEngine T=" + std::to_string(e.temperature()) + "K>";
        });

    // ── free function ────────────────────────────────────────────────────────
    m.def("dift_spectral_entropy",
        [](const std::vector<FourierTerm>& spectrum) {
            return spectral_entropy(std::span<const FourierTerm>(spectrum));
        },
        py::arg("spectrum"),
        "Spectral Shannon entropy H_spec (nats) of a power spectrum");
}
