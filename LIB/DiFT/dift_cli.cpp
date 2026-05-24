// dift_cli.cpp - command-line DiFT torsional parametrization
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "DiFT.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;
constexpr double kPi = 3.1415926535897932384626433832795;

struct Options {
    std::string input;
    std::string json_output;
    std::string sidecar_output;
    int max_multiplicity = 6;
    double temperature_K = 300.0;
    bool histogram = false;
};

struct ScanData {
    std::vector<double> angles_deg;
    std::vector<double> values;
    bool has_angles = false;
};

[[noreturn]] void usage(int exit_code) {
    std::ostream& out = exit_code == 0 ? std::cout : std::cerr;
    out << "Usage: dift_torsion_fit --input scan.csv --json fit.json [options]\n"
        << "\n"
        << "Options:\n"
        << "  --sidecar FILE            Write FlexAIDdS .dift sidecar terms\n"
        << "  --max-multiplicity N      Highest Fourier multiplicity to retain (default: 6)\n"
        << "  --temperature K           Temperature for thermodynamics (default: 300)\n"
        << "  --histogram               Treat values as counts and Boltzmann-invert first\n"
        << "  --help                    Show this help\n"
        << "\n"
        << "Input may be one numeric column (values on an implicit uniform grid) or\n"
        << "two numeric columns: angle_degrees,value. Angle grids must be uniform over\n"
        << "[0, 360) and must not repeat the 360-degree endpoint.\n";
    std::exit(exit_code);
}

Options parse_args(int argc, char** argv) {
    Options opts;
    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);
        auto need_value = [&](std::string_view flag) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(flag) + " needs a value");
            return argv[++i];
        };
        if (arg == "--help" || arg == "-h") usage(0);
        else if (arg == "--input") opts.input = need_value(arg);
        else if (arg == "--json") opts.json_output = need_value(arg);
        else if (arg == "--sidecar") opts.sidecar_output = need_value(arg);
        else if (arg == "--max-multiplicity") opts.max_multiplicity = std::stoi(need_value(arg));
        else if (arg == "--temperature") opts.temperature_K = std::stod(need_value(arg));
        else if (arg == "--histogram") opts.histogram = true;
        else throw std::invalid_argument("unknown option: " + std::string(arg));
    }
    if (opts.input.empty()) throw std::invalid_argument("--input is required");
    if (opts.json_output.empty()) throw std::invalid_argument("--json is required");
    if (opts.max_multiplicity < 0) throw std::invalid_argument("--max-multiplicity must be >= 0");
    if (!(opts.temperature_K > 0.0) || !std::isfinite(opts.temperature_K)) {
        throw std::invalid_argument("--temperature must be finite and > 0");
    }
    return opts;
}

std::string trim(std::string s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::vector<std::string> split_numeric_fields(std::string line) {
    for (char& c : line) {
        if (c == ',' || c == ';' || c == '\t') c = ' ';
    }
    std::istringstream in(line);
    std::vector<std::string> fields;
    std::string token;
    while (in >> token) fields.push_back(token);
    return fields;
}

bool parse_double(const std::string& s, double& out) {
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && std::isfinite(out);
}

ScanData read_scan(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("cannot open input: " + filename);

    ScanData data;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto fields = split_numeric_fields(line);
        if (fields.empty()) continue;

        double a = 0.0, v = 0.0;
        if (fields.size() == 1) {
            if (!parse_double(fields[0], v)) {
                if (data.values.empty()) continue; // tolerate one header before data
                throw std::runtime_error("non-numeric value on line " + std::to_string(line_no));
            }
            if (data.has_angles) {
                throw std::runtime_error("mixed one-column and two-column data");
            }
            data.values.push_back(v);
        } else {
            if (!parse_double(fields[0], a) || !parse_double(fields[1], v)) {
                if (data.values.empty()) continue; // tolerate one header before data
                throw std::runtime_error("non-numeric angle/value on line " + std::to_string(line_no));
            }
            if (!data.values.empty() && !data.has_angles) {
                throw std::runtime_error("mixed one-column and two-column data");
            }
            data.has_angles = true;
            data.angles_deg.push_back(a);
            data.values.push_back(v);
        }
    }

    if (data.values.size() < 3) throw std::runtime_error("need at least 3 torsion samples");
    return data;
}

void validate_angles(const std::vector<double>& angles) {
    if (angles.empty()) return;
    const std::size_t n = angles.size();
    if (n < 3) throw std::runtime_error("need at least 3 torsion angle samples");
    for (double a : angles) {
        if (!std::isfinite(a)) throw std::runtime_error("angle grid contains non-finite values");
        if (a < -1e-9 || a >= 360.0 - 1e-9) {
            throw std::runtime_error("angle grid must cover [0, 360) without a duplicate 360 endpoint");
        }
    }
    std::vector<double> sorted = angles;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (std::abs(sorted[i] - sorted[i - 1]) < 1e-9) {
            throw std::runtime_error("angle grid contains duplicate angles");
        }
    }
    const double expected = 360.0 / static_cast<double>(n);
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        const double step = sorted[i] - sorted[i - 1];
        if (std::abs(step - expected) > 1e-6) {
            throw std::runtime_error("angle grid must be uniformly spaced over [0, 360)");
        }
    }
    const double wrap_step = sorted.front() + 360.0 - sorted.back();
    if (std::abs(wrap_step - expected) > 1e-6) {
        throw std::runtime_error("angle grid must cover a complete [0, 360) period");
    }
}

std::vector<double> order_by_angle(const ScanData& data) {
    if (!data.has_angles) return data.values;
    validate_angles(data.angles_deg);
    std::vector<std::size_t> order(data.values.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return data.angles_deg[a] < data.angles_deg[b];
    });
    std::vector<double> values;
    values.reserve(order.size());
    for (std::size_t idx : order) values.push_back(data.values[idx]);
    return values;
}

std::vector<double> full_model(const dift::TorsionalPotential& pot, std::size_t n) {
    return pot.sample(static_cast<int>(n));
}

double rmse(const std::vector<double>& obs, const std::vector<double>& model) {
    if (obs.size() != model.size() || obs.empty()) return 0.0;
    double ss = 0.0;
    for (std::size_t i = 0; i < obs.size(); ++i) {
        const double d = obs[i] - model[i];
        ss += d * d;
    }
    return std::sqrt(ss / static_cast<double>(obs.size()));
}

double max_abs_error(const std::vector<double>& obs, const std::vector<double>& model) {
    double out = 0.0;
    for (std::size_t i = 0; i < obs.size() && i < model.size(); ++i) {
        out = std::max(out, std::abs(obs[i] - model[i]));
    }
    return out;
}

void write_json(const Options& opts,
                const std::vector<double>& profile,
                const dift::TorsionalPotential& pot,
                const dift::TorsionalThermo& thermo) {
    auto model = full_model(pot, profile.size());
    std::ofstream out(opts.json_output);
    if (!out) throw std::runtime_error("cannot write JSON: " + opts.json_output);

    out << std::setprecision(12);
    out << "{\n";
    out << "  \"method\": \"DiFT\",\n";
    out << "  \"input\": \"" << opts.input << "\",\n";
    out << "  \"n_samples\": " << profile.size() << ",\n";
    out << "  \"temperature_K\": " << thermo.temperature << ",\n";
    out << "  \"max_multiplicity\": " << opts.max_multiplicity << ",\n";
    out << "  \"profile_type\": \"" << (opts.histogram ? "boltzmann_inverted_histogram" : "energy_scan") << "\",\n";
    out << "  \"mean_kcal_mol\": " << pot.mean << ",\n";
    out << "  \"v_min_kcal_mol\": " << pot.v_min << ",\n";
    out << "  \"r_squared\": " << pot.r_squared << ",\n";
    out << "  \"rmse_kcal_mol\": " << rmse(profile, model) << ",\n";
    out << "  \"max_abs_error_kcal_mol\": " << max_abs_error(profile, model) << ",\n";
    out << "  \"spectral_entropy_nats\": " << pot.spectral_entropy << ",\n";
    out << "  \"effective_modes\": " << pot.effective_modes << ",\n";
    out << "  \"thermodynamics\": {\n";
    out << "    \"partition_function\": " << thermo.partition_function << ",\n";
    out << "    \"free_energy_kcal_mol\": " << thermo.free_energy << ",\n";
    out << "    \"mean_energy_kcal_mol\": " << thermo.mean_energy << ",\n";
    out << "    \"entropy_kcal_mol_K\": " << thermo.entropy << ",\n";
    out << "    \"minus_TS_kcal_mol\": " << thermo.minus_TS << "\n";
    out << "  },\n";
    out << "  \"terms\": [\n";
    for (std::size_t i = 0; i < pot.terms.size(); ++i) {
        const auto& t = pot.terms[i];
        out << "    {\"multiplicity\": " << t.multiplicity
            << ", \"amplitude_kcal_mol\": " << t.amplitude
            << ", \"phase_rad\": " << t.phase
            << ", \"phase_deg\": " << (t.phase * 180.0 / kPi)
            << ", \"power\": " << t.power << "}";
        out << (i + 1 == pot.terms.size() ? "\n" : ",\n");
    }
    out << "  ]\n";
    out << "}\n";
}

void write_sidecar(const Options& opts,
                   const dift::TorsionalPotential& pot,
                   const dift::TorsionalThermo& thermo) {
    if (opts.sidecar_output.empty()) return;
    std::ofstream out(opts.sidecar_output);
    if (!out) throw std::runtime_error("cannot write sidecar: " + opts.sidecar_output);
    out << std::setprecision(12);
    out << "# FlexAIDdS DiFT torsion prior v1\n";
    out << "source " << opts.input << "\n";
    out << "samples " << pot.n_samples << "\n";
    out << "temperature_K " << thermo.temperature << "\n";
    out << "mean_kcal_mol " << pot.mean << "\n";
    out << "v_min_kcal_mol " << pot.v_min << "\n";
    out << "r_squared " << pot.r_squared << "\n";
    out << "minus_TS_kcal_mol " << thermo.minus_TS << "\n";
    out << "# term multiplicity amplitude_kcal_mol phase_rad\n";
    for (const auto& t : pot.terms) {
        out << "term " << t.multiplicity << ' ' << t.amplitude << ' ' << t.phase << "\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parse_args(argc, argv);
        const ScanData scan = read_scan(opts.input);
        std::vector<double> profile = order_by_angle(scan);
        for (double v : profile) {
            if (!std::isfinite(v)) throw std::runtime_error("profile contains non-finite values");
        }

        dift::DiFTEngine engine(opts.temperature_K);
        if (opts.histogram) profile = engine.boltzmann_invert(profile);
        auto pot = engine.parametrize(profile, opts.max_multiplicity);
        auto thermo = engine.thermodynamics(pot);

        write_json(opts, profile, pot, thermo);
        write_sidecar(opts, pot, thermo);
        std::cout << "DiFT fit: " << pot.n_terms() << " terms, R2=" << pot.r_squared
                  << ", JSON=" << opts.json_output;
        if (!opts.sidecar_output.empty()) std::cout << ", sidecar=" << opts.sidecar_output;
        std::cout << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "dift_torsion_fit: " << e.what() << "\n";
        return 2;
    }
}
