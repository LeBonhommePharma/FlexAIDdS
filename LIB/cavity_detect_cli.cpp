// cavity_detect_cli.cpp -- native FlexAIDdS CavityDetector sphere exporter
//
// Generates FlexAID-compatible cleft sphere PDB files from the in-tree
// CavityDetector. With --ligand, the selected cleft is the detected cleft whose
// probe sphere is nearest the ligand centroid, matching the "occupied cavity"
// benchmark protocol without calling external Get_Cleft.

#include "CavityDetect/CavityDetect.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Options {
    std::string receptor;
    std::string ligand;
    std::string out;
    std::string out_dir;
    std::string prefix = "cleft";
    int top = 1;
    float min_radius = 1.4f;
    float max_radius = 4.0f;
    float site_cutoff = 15.0f;
};

void usage(const char* argv0) {
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " --receptor receptor.pdb --out cleft.pdb [--ligand ligand.sdf]\n"
        << "  " << argv0 << " --receptor receptor.pdb --out-dir cavities --prefix 1G9V --top 3\n"
        << "\nOptions:\n"
        << "  --receptor <pdb>       Receptor PDB used for native cavity detection\n"
        << "  --ligand <sdf|pdb>     Optional ligand used to select occupied cleft\n"
        << "  --out <pdb>            Write one selected sphere PDB\n"
        << "  --out-dir <dir>        Write top-N ranked sphere PDBs\n"
        << "  --prefix <name>        Prefix for --out-dir files (default: cleft)\n"
        << "  --top <N>              Number of ranked clefts for --out-dir (default: 1)\n"
        << "  --min-radius <A>       SURFNET probe lower bound (default: 1.4)\n"
        << "  --max-radius <A>       SURFNET probe upper bound (default: 4.0)\n"
        << "  --site-cutoff <A>      Receptor atom cutoff around ligand centroid (default: 15; 0 disables)\n";
}

bool starts_with(const std::string& s, const char* prefix) {
    return s.rfind(prefix, 0) == 0;
}

std::string lower_ext(const std::string& path) {
    std::string ext = fs::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext;
}

std::optional<std::array<double, 3>> read_pdb_centroid(const std::string& path) {
    std::ifstream in(path);
    if (!in) return std::nullopt;

    std::array<double, 3> sum{0.0, 0.0, 0.0};
    std::size_t n = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!starts_with(line, "ATOM  ") && !starts_with(line, "HETATM")) continue;
        if (line.size() < 54) continue;
        try {
            sum[0] += std::stod(line.substr(30, 8));
            sum[1] += std::stod(line.substr(38, 8));
            sum[2] += std::stod(line.substr(46, 8));
            ++n;
        } catch (...) {
            continue;
        }
    }
    if (n == 0) return std::nullopt;
    const double inv = 1.0 / static_cast<double>(n);
    return std::array<double, 3>{sum[0] * inv, sum[1] * inv, sum[2] * inv};
}

std::optional<std::array<double, 3>> read_sdf_centroid(const std::string& path) {
    std::ifstream in(path);
    if (!in) return std::nullopt;

    std::string line;
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(in, line)) return std::nullopt;
    }
    if (!std::getline(in, line) || line.size() < 3) return std::nullopt;

    int atom_count = 0;
    try {
        atom_count = std::stoi(line.substr(0, 3));
    } catch (...) {
        return std::nullopt;
    }
    if (atom_count <= 0) return std::nullopt;

    std::array<double, 3> sum{0.0, 0.0, 0.0};
    int n = 0;
    for (int i = 0; i < atom_count && std::getline(in, line); ++i) {
        std::istringstream iss(line);
        double x = 0.0, y = 0.0, z = 0.0;
        if (iss >> x >> y >> z) {
            sum[0] += x;
            sum[1] += y;
            sum[2] += z;
            ++n;
        }
    }
    if (n == 0) return std::nullopt;
    const double inv = 1.0 / static_cast<double>(n);
    return std::array<double, 3>{sum[0] * inv, sum[1] * inv, sum[2] * inv};
}

std::optional<std::array<double, 3>> read_ligand_centroid(const std::string& path) {
    const std::string ext = lower_ext(path);
    if (ext == ".sdf" || ext == ".mol") return read_sdf_centroid(path);
    return read_pdb_centroid(path);
}

std::optional<std::array<double, 3>> pdb_line_xyz(const std::string& line) {
    if (line.size() < 54) return std::nullopt;
    try {
        return std::array<double, 3>{
            std::stod(line.substr(30, 8)),
            std::stod(line.substr(38, 8)),
            std::stod(line.substr(46, 8)),
        };
    } catch (...) {
        return std::nullopt;
    }
}

double dist2(const float center[3], const std::array<double, 3>& point) {
    const double dx = static_cast<double>(center[0]) - point[0];
    const double dy = static_cast<double>(center[1]) - point[1];
    const double dz = static_cast<double>(center[2]) - point[2];
    return dx * dx + dy * dy + dz * dz;
}

double dist2(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

std::size_t write_site_filtered_receptor(
    const std::string& receptor,
    const std::array<double, 3>& centroid,
    double cutoff,
    const fs::path& out_path)
{
    std::ifstream in(receptor);
    if (!in) throw std::runtime_error("cannot open receptor: " + receptor);
    if (!out_path.parent_path().empty()) {
        fs::create_directories(out_path.parent_path());
    }
    std::ofstream out(out_path);
    if (!out) throw std::runtime_error("cannot write site receptor: " + out_path.string());

    const double cutoff2 = cutoff * cutoff;
    std::size_t kept = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!starts_with(line, "ATOM  ") && !starts_with(line, "HETATM")) continue;
        auto xyz = pdb_line_xyz(line);
        if (!xyz || dist2(*xyz, centroid) > cutoff2) continue;
        out << line << "\n";
        ++kept;
    }
    out << "END\n";
    return kept;
}

fs::path site_receptor_path(const Options& opt) {
    if (!opt.out.empty()) {
        return fs::path(opt.out).replace_extension(".site.pdb");
    }
    return fs::path(opt.out_dir) / (opt.prefix + "_site.pdb");
}

int select_occupied_cleft(
    const std::vector<cavity_detect::DetectedCleft>& clefts,
    const std::array<double, 3>& ligand_centroid)
{
    int best_id = clefts.front().id;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& cleft : clefts) {
        for (const auto& sphere : cleft.spheres) {
            const double d2 = dist2(sphere.center, ligand_centroid);
            if (d2 < best) {
                best = d2;
                best_id = cleft.id;
            }
        }
    }
    return best_id;
}

Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* label) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("missing value for ") + label);
            }
            return argv[++i];
        };

        if (arg == "--receptor") opt.receptor = need_value("--receptor");
        else if (arg == "--ligand") opt.ligand = need_value("--ligand");
        else if (arg == "--out") opt.out = need_value("--out");
        else if (arg == "--out-dir") opt.out_dir = need_value("--out-dir");
        else if (arg == "--prefix") opt.prefix = need_value("--prefix");
        else if (arg == "--top") opt.top = std::max(1, std::stoi(need_value("--top")));
        else if (arg == "--min-radius") opt.min_radius = std::stof(need_value("--min-radius"));
        else if (arg == "--max-radius") opt.max_radius = std::stof(need_value("--max-radius"));
        else if (arg == "--site-cutoff") opt.site_cutoff = std::stof(need_value("--site-cutoff"));
        else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    if (opt.receptor.empty()) throw std::runtime_error("--receptor is required");
    if (opt.out.empty() && opt.out_dir.empty()) throw std::runtime_error("--out or --out-dir is required");
    if (opt.min_radius <= 0.0f || opt.max_radius < opt.min_radius) {
        throw std::runtime_error("invalid radius bounds");
    }
    return opt;
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    try {
        opt = parse_args(argc, argv);
    } catch (const std::exception& exc) {
        std::cerr << "cavity_detect_cli: " << exc.what() << "\n";
        usage(argv[0]);
        return 2;
    }

    try {
        std::optional<std::array<double, 3>> ligand_centroid;
        if (!opt.ligand.empty()) {
            ligand_centroid = read_ligand_centroid(opt.ligand);
            if (!ligand_centroid) {
                std::cerr << "cavity_detect_cli: could not read ligand centroid: "
                          << opt.ligand << "\n";
                return 4;
            }
        }

        std::string receptor_for_detection = opt.receptor;
        if (ligand_centroid && opt.site_cutoff > 0.0f) {
            const fs::path site_path = site_receptor_path(opt);
            const std::size_t kept = write_site_filtered_receptor(
                opt.receptor, *ligand_centroid, opt.site_cutoff, site_path);
            if (kept < 2) {
                std::cerr << "cavity_detect_cli: ligand-centered receptor filter kept "
                          << kept << " atoms; increase --site-cutoff\n";
                return 5;
            }
            receptor_for_detection = site_path.string();
            std::cerr << "cavity_detect_cli: kept " << kept
                      << " receptor atoms within " << opt.site_cutoff
                      << " A at " << site_path << "\n";
        }

        cavity_detect::CavityDetector detector;
        detector.load_from_pdb(receptor_for_detection);
        detector.detect(opt.min_radius, opt.max_radius);
        const auto& clefts = detector.clefts();
        if (clefts.empty()) {
            std::cerr << "cavity_detect_cli: no clefts detected in "
                      << receptor_for_detection << "\n";
            return 3;
        }

        if (!opt.out.empty()) {
            int selected_id = clefts.front().id;
            if (ligand_centroid) {
                selected_id = select_occupied_cleft(clefts, *ligand_centroid);
            }
            const fs::path out_path(opt.out);
            if (!out_path.parent_path().empty()) {
                fs::create_directories(out_path.parent_path());
            }
            detector.write_sphere_pdb(opt.out, selected_id);
            std::cerr << "cavity_detect_cli: wrote cleft " << selected_id
                      << " to " << opt.out << "\n";
        }

        if (!opt.out_dir.empty()) {
            fs::create_directories(opt.out_dir);
            const int n = std::min<int>(opt.top, static_cast<int>(clefts.size()));
            for (int i = 0; i < n; ++i) {
                const auto& cleft = clefts[static_cast<std::size_t>(i)];
                fs::path out = fs::path(opt.out_dir) /
                    (opt.prefix + "_sph_" + std::to_string(i + 1) + ".pdb");
                detector.write_sphere_pdb(out.string(), cleft.id);
                std::cerr << "cavity_detect_cli: wrote ranked cleft " << (i + 1)
                          << " id " << cleft.id << " to " << out << "\n";
            }
        }
    } catch (const std::exception& exc) {
        std::cerr << "cavity_detect_cli: " << exc.what() << "\n";
        return 1;
    }
    return 0;
}
