// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
#ifndef FLEXAIDDS_DATASET_THERMO_LOG_H
#define FLEXAIDDS_DATASET_THERMO_LOG_H

#include <string>

namespace dataset {

// Parses the stable thermodynamic-summary fields consumed by DatasetRunner.
// Both the historical physical-looking labels and the schema-v2 proxy labels
// are accepted so that relabeling output cannot silently change CSV numerics.
//
// Selectivity matters as much as recognition. Two producers print an ensemble
// summary to the same stdout stream:
//
//   gaboom's post-GA block   "F-like proxy          F~  = ..."   (was
//                            "Helmholtz free energy  F  = ...")
//   top.cpp's Post-GA block  "F-like proxy   = ..."              (was
//                            "Free energy F  = ...")
//
// Before the relabel, only gaboom's spelling matched the free-energy and
// entropy patterns; top.cpp's "Free energy F" / "Entropy S" did not. Matching
// the shared "F-like proxy" / "S-like" stems would newly admit top.cpp's block
// and change result.csv, so the free-energy and entropy matchers additionally
// require gaboom's distinctive `F~` marker and its "diagnostic" spelling.
//
// Mean energy is deliberately NOT narrowed: at HEAD both producers matched the
// "Mean energy" substring, so both must keep matching to preserve the existing
// last-writer-wins outcome.
struct DatasetThermoLog {
    float free_energy = 0.0f;
    float mean_energy = 0.0f;
    float configurational_entropy = 0.0f;
    bool have_free_energy = false;
    bool have_mean_energy = false;
    bool have_configurational_entropy = false;

    void consume(const std::string& line) {
        if (is_ensemble_free_energy(line)) {
            assign_if_parsed(line, free_energy, have_free_energy);
        } else if (contains_any(line, "Mean CF proxy", "Mean CF") ||
                   line.find("Mean energy") != std::string::npos) {
            assign_if_parsed(line, mean_energy, have_mean_energy);
        } else if (contains_any(line, "S-like diagnostic", "Entropy (conf)")) {
            assign_if_parsed(line, configurational_entropy,
                             have_configurational_entropy);
        }
    }

private:
    // gaboom's ensemble line, in either vocabulary. The `F~` marker and the
    // "Helmholtz" spelling are both absent from top.cpp's block, which is what
    // keeps that block out of the CSV exactly as it was before the relabel.
    static bool is_ensemble_free_energy(const std::string& line) {
        const bool proxy_form = line.find("F-like proxy") != std::string::npos &&
                                line.find("F~") != std::string::npos;
        return proxy_form || line.find("Helmholtz free energy") != std::string::npos;
    }

    // A later line that matches a label but carries an unparsable value must
    // not revoke an earlier successful parse. The inline parsers this helper
    // replaced only ever set their flag inside a successful try.
    static void assign_if_parsed(const std::string& line,
                                 float& destination,
                                 bool& flag) {
        float parsed = 0.0f;
        if (parse_value(line, parsed)) {
            destination = parsed;
            flag = true;
        }
    }

    static bool contains_any(const std::string& line,
                             const char* first,
                             const char* second) {
        return line.find(first) != std::string::npos ||
               line.find(second) != std::string::npos;
    }

    static bool parse_value(const std::string& line, float& destination) {
        const auto equals = line.find('=');
        if (equals == std::string::npos) return false;
        try {
            destination = std::stof(line.substr(equals + 1));
            return true;
        } catch (...) {
            return false;
        }
    }
};

}  // namespace dataset

#endif  // FLEXAIDDS_DATASET_THERMO_LOG_H
