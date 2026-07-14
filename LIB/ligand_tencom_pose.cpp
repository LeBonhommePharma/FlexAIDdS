#include "PoseBust/Loaders.h"
#include "VibEntropy.h"
#include "flexaid.h"
#include "tENCoM/tencm.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::string json_escape(const std::string& value) {
    std::string out;
    out.reserve(value.size());
    for (const char ch : value) {
        if (ch == '\\' || ch == '"') out.push_back('\\');
        out.push_back(ch);
    }
    return out;
}

void usage(const char* program) {
    std::cerr << "Usage: " << program
              << " --pose ligand.sdf [--output metrics.json]\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string pose_path;
    std::string output_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--pose" && i + 1 < argc) {
            pose_path = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    if (pose_path.empty()) {
        usage(argv[0]);
        return 2;
    }

    flexaids::posebust::Molecule molecule;
    std::string error;
    if (!flexaids::posebust::load_sdf(pose_path, molecule, &error)) {
        std::cerr << "ligand_tencom_pose: " << error << "\n";
        return 3;
    }

    std::vector<atom> atoms(molecule.atoms.size());
    for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
        atoms[i] = atom{};
        atoms[i].coor[0] = molecule.atoms[i].x;
        atoms[i].coor[1] = molecule.atoms[i].y;
        atoms[i].coor[2] = molecule.atoms[i].z;
        std::strncpy(atoms[i].element, molecule.atoms[i].element.c_str(),
                     sizeof(atoms[i].element) - 1);
    }

    tencm::TorsionalENM model;
    model.build_from_ligand(atoms.data(), 0, static_cast<int>(atoms.size()));
    if (!model.is_built()) {
        std::cerr << "ligand_tencom_pose: ligand Cartesian ANM did not build\n";
        return 4;
    }

    double max_eigenvalue = 0.0;
    for (const auto& mode : model.modes()) {
        if (std::isfinite(mode.eigenvalue))
            max_eigenvalue = std::max(max_eigenvalue, mode.eigenvalue);
    }
    const double cutoff = std::max(1e-10, max_eigenvalue * 1e-8);
    std::vector<double> eigenvalues;
    eigenvalues.reserve(model.modes().size());
    for (const auto& mode : model.modes()) {
        if (std::isfinite(mode.eigenvalue) && mode.eigenvalue > cutoff)
            eigenvalues.push_back(mode.eigenvalue);
    }
    if (eigenvalues.empty()) {
        std::cerr << "ligand_tencom_pose: Eigen returned no positive modes\n";
        return 5;
    }

    const auto entropy = vibentropy::compute_vib_entropy_collapse({eigenvalues});
    std::ostringstream json;
    json << std::setprecision(10)
         << "{\n"
         << "  \"pose\": \"" << json_escape(pose_path) << "\",\n"
         << "  \"engine\": \"tENCoM/Eigen\",\n"
         << "  \"model\": \"ligand_cartesian_anm\",\n"
         << "  \"quantity\": \"eigenvalue_spectrum_shannon_entropy\",\n"
         << "  \"units\": \"nats\",\n"
         << "  \"n_atoms\": " << molecule.atoms.size() << ",\n"
         << "  \"n_heavy_atoms\": " << molecule.n_heavy() << ",\n"
         << "  \"n_positive_modes\": " << eigenvalues.size() << ",\n"
         << "  \"H_vib_nats\": " << entropy.H_pop << ",\n"
         << "  \"eigen_status\": \"ok\",\n"
         << "  \"absolute_free_energy_calibrated\": false\n"
         << "}\n";

    if (output_path.empty()) {
        std::cout << json.str();
    } else {
        std::ofstream output(output_path);
        if (!output) {
            std::cerr << "ligand_tencom_pose: cannot write " << output_path << "\n";
            return 6;
        }
        output << json.str();
    }
    return 0;
}
