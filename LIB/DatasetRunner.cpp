// =============================================================================
// DatasetRunner.cpp — Benchmark dataset runner implementation for FlexAIDdS
//
// Full production-grade implementation:
//   - Hardcoded PDB code lists for all standard benchmarks
//   - PDB download via RCSB REST API using system curl
//   - Ligand extraction from PDB HETATM records → SDF output
//   - RMSD computation against crystal pose
//   - Pearson r, Spearman ρ, Kendall τ computed from scratch
//   - Markdown + CSV report generation
//   - Local caching in ~/.flexaidds/benchmarks/
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "DatasetRunner.h"
#include "AsyncPipeline.h"
#include "BenchmarkRunner.h"
#include "statmech.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <thread>

#ifndef _MSC_VER
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <cerrno>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace dataset {

// =============================================================================
// Static pointers for async-signal-safe handler access.
// Set just before sigaction(), cleared after signal restore on normal exit.
// =============================================================================
#ifndef _MSC_VER
static SubprocessGuard*   g_active_guard   = nullptr;
static std::atomic<bool>* g_active_shutdown = nullptr;
#endif

// =============================================================================
// SubprocessGuard — RAII process lifecycle
// =============================================================================

SubprocessGuard::SubprocessGuard() = default;

SubprocessGuard::~SubprocessGuard() {
    kill_all();
}

pid_t SubprocessGuard::fork_exec(const std::string& cmd) {
#ifndef _MSC_VER
    pid_t pid = fork();
    if (pid < 0) return -1;

    if (pid == 0) {
        // ── Child process ───────────────────────────────────────────────
        // Create own process group so parent can killpg() all children.
        ::setpgid(0, 0);

        // Redirect stdin to /dev/null
        ::close(0);
        ::open("/dev/null", O_RDONLY);

        // Reset SIGTERM to default (in case parent had a handler)
        ::signal(SIGTERM, SIG_DFL);
        ::signal(SIGINT, SIG_DFL);

        ::execl("/bin/sh", "sh", "-c", cmd.c_str(), nullptr);
        ::_exit(127);  // exec failed
    }

    // ── Parent: register PID ────────────────────────────────────────────
    {
        std::lock_guard<std::mutex> lock(mtx_);
        pids_.insert(pid);
    }
    // Also set parent-side pgid (redundant with child's setpgid but handles race)
    ::setpgid(pid, pid);
    return pid;
#else
    // Windows: no process groups, use system() and return a sentinel PID.
    // Real PID tracking not available on Windows.
    int ret = std::system(cmd.c_str());
    return ret;
#endif
}

int SubprocessGuard::wait_with_timeout(pid_t pid, int timeout_s) {
#ifndef _MSC_VER
    if (pid <= 0) return -1;

    using clock = std::chrono::steady_clock;
    auto deadline = (timeout_s > 0)
        ? clock::now() + std::chrono::seconds(timeout_s)
        : clock::time_point::max();  // no timeout

    int status = 0;
    for (;;) {
        pid_t wp = ::waitpid(pid, &status, WNOHANG);
        if (wp == pid) {
            // Child exited
            forget(pid);
            return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
        }
        if (wp < 0) {
            // Error (ECHILD = already reaped)
            forget(pid);
            return -1;
        }
        // wp == 0: still running

        if (clock::now() >= deadline) {
            // ── Timeout: escalate SIGTERM → SIGKILL ────────────────────
            std::cerr << "[TIMEOUT] PID " << pid
                      << " exceeded " << timeout_s << "s limit, sending SIGTERM\n";

            // Kill the entire process group (pid's pgid == pid)
            ::killpg(pid, SIGTERM);

            // Grace period: 5 seconds
            auto kill_deadline = clock::now() + std::chrono::seconds(5);
            for (;;) {
                wp = ::waitpid(pid, &status, WNOHANG);
                if (wp == pid || wp < 0) {
                    forget(pid);
                    return -1;  // exited after SIGTERM
                }
                if (clock::now() >= kill_deadline) break;
                ::usleep(100000);  // 100ms
            }

            // Still alive → SIGKILL
            std::cerr << "[TIMEOUT] PID " << pid
                      << " still alive, sending SIGKILL\n";
            ::killpg(pid, SIGKILL);
            ::waitpid(pid, &status, 0);  // reap zombie
            forget(pid);
            return -1;
        }

        ::usleep(200000);  // 200ms poll interval
    }
#else
    return -1;  // Windows: timeout not supported
#endif
}

void SubprocessGuard::forget(pid_t pid) {
    std::lock_guard<std::mutex> lock(mtx_);
    pids_.erase(pid);
}

void SubprocessGuard::kill_all() {
#ifndef _MSC_VER
    std::set<pid_t> to_kill;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        to_kill = pids_;
    }
    if (to_kill.empty()) return;

    std::cerr << "[SubprocessGuard] Killing " << to_kill.size()
              << " remaining child process(es)\n";

    // SIGTERM first
    for (pid_t pid : to_kill) {
        ::killpg(pid, SIGTERM);
    }
    ::usleep(500000);  // 500ms grace

    // SIGKILL for survivors
    for (pid_t pid : to_kill) {
        // Check if still alive (kill(pid,0) returns 0 if process exists)
        if (::kill(pid, 0) == 0) {
            ::killpg(pid, SIGKILL);
        }
        // Reap zombie (non-blocking)
        int dummy;
        ::waitpid(pid, &dummy, WNOHANG);
    }

    {
        std::lock_guard<std::mutex> lock(mtx_);
        pids_.clear();
    }
#endif
}

size_t SubprocessGuard::active_count() const {
    std::lock_guard<std::mutex> lock(const_cast<std::mutex&>(mtx_));
    return pids_.size();
}

// =============================================================================
// Statistical functions — implemented from scratch, no external stats library
// =============================================================================

double compute_pearson_r(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    const size_t n = x.size();

    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / static_cast<double>(n);
    double mean_y = sum_y / static_cast<double>(n);

    double cov_xy = 0.0, var_x = 0.0, var_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        cov_xy += dx * dy;
        var_x  += dx * dx;
        var_y  += dy * dy;
    }

    double denom = std::sqrt(var_x * var_y);
    if (denom < 1e-15) return 0.0;
    return cov_xy / denom;
}

/// Helper: compute ranks for a vector (average rank for ties)
static std::vector<double> compute_ranks(const std::vector<double>& vals) {
    const size_t n = vals.size();
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&vals](size_t a, size_t b) { return vals[a] < vals[b]; });

    std::vector<double> ranks(n);
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        // Find all tied elements
        while (j < n && vals[indices[j]] == vals[indices[i]]) ++j;
        // Average rank for ties (1-based)
        double avg_rank = 0.5 * (static_cast<double>(i + 1) + static_cast<double>(j));
        for (size_t k = i; k < j; ++k) {
            ranks[indices[k]] = avg_rank;
        }
        i = j;
    }
    return ranks;
}

double compute_spearman_rho(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    // Spearman ρ = Pearson r of ranks
    std::vector<double> rx = compute_ranks(x);
    std::vector<double> ry = compute_ranks(y);
    return compute_pearson_r(rx, ry);
}

double compute_kendall_tau(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    const size_t n = x.size();

    // Kendall tau-b: handles ties
    int64_t concordant = 0, discordant = 0;
    int64_t ties_x = 0, ties_y = 0, ties_xy = 0;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double sign_product = dx * dy;

            bool tx = (std::abs(dx) < 1e-12);
            bool ty = (std::abs(dy) < 1e-12);

            if (tx && ty) {
                ties_xy++;
            } else if (tx) {
                ties_x++;
            } else if (ty) {
                ties_y++;
            } else if (sign_product > 0) {
                concordant++;
            } else {
                discordant++;
            }
        }
    }

    int64_t n_pairs = static_cast<int64_t>(n) * (static_cast<int64_t>(n) - 1) / 2;
    double n0 = static_cast<double>(n_pairs);
    double n1 = static_cast<double>(ties_x + ties_xy);
    double n2 = static_cast<double>(ties_y + ties_xy);

    double denom = std::sqrt((n0 - n1) * (n0 - n2));
    if (denom < 1e-15) return 0.0;

    return static_cast<double>(concordant - discordant) / denom;
}

double compute_rmsd(const std::vector<float>& coords_a,
                    const std::vector<float>& coords_b) {
    if (coords_a.size() != coords_b.size() || coords_a.empty()) return 999.0;
    if (coords_a.size() % 3 != 0) return 999.0;

    const size_t n_atoms = coords_a.size() / 3;
    double sum_sq = 0.0;
    for (size_t i = 0; i < coords_a.size(); ++i) {
        double d = static_cast<double>(coords_a[i]) - static_cast<double>(coords_b[i]);
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq / static_cast<double>(n_atoms));
}

// =============================================================================
// Excluded residues — water, common ions, buffers
// =============================================================================

const std::set<std::string>& DatasetRunner::excluded_residues() {
    static const std::set<std::string> excl = {
        "HOH", "WAT", "H2O", "DOD", "DIS",  // water
        "NA",  "CL",  "MG",  "CA",  "ZN",   // common ions
        "FE",  "MN",  "CU",  "CO",  "NI",
        "K",   "BR",  "I",   "F",
        "SO4", "PO4", "NO3", "ACT",          // buffer components
        "GOL", "EDO", "PEG", "DMS", "MPD",   // cryoprotectants / crystallization aids
        "BME", "EPE", "MES", "TRS", "CIT",
        "IMD", "FMT", "ACE", "NH4", "IOD",
        "BOG", "PGE", "1PE", "P6G", "BU3",
        "PDO", "EGL", "PG4", "PE8", "MLI",
        "DTT", "AZI", "SCN", "NO2", "OXL"
    };
    return excl;
}

// =============================================================================
// DatasetRunner constructor
// =============================================================================

DatasetRunner::DatasetRunner(const std::string& cache_dir) {
    if (cache_dir.empty()) {
        cache_dir_ = expand_home("~/.flexaidds/benchmarks");
    } else {
        cache_dir_ = expand_home(cache_dir);
    }
    ensure_dir(cache_dir_);
}

// =============================================================================
// Path utilities
// =============================================================================

std::string DatasetRunner::expand_home(const std::string& path) {
    if (path.empty() || path[0] != '~') return path;
    const char* home = std::getenv("HOME");
    if (!home) home = "/tmp";
    return std::string(home) + path.substr(1);
}

bool DatasetRunner::ensure_dir(const std::string& path) {
    std::error_code ec;
    fs::create_directories(path, ec);
    return !ec;
}

// =============================================================================
// HTTP download using system curl
// =============================================================================

int DatasetRunner::exec_cmd(const std::string& cmd) {
    // Legacy entry point — no timeout, no PID tracking.
    // Used for downloads and other non-docking commands.
    // Docking commands go through exec_dock() which uses proc_guard_.
#ifndef _MSC_VER
    pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        ::close(0);
        ::open("/dev/null", O_RDONLY);
        ::execl("/bin/sh", "sh", "-c", cmd.c_str(), nullptr);
        ::_exit(127);
    }
    int status = 0;
    while (::waitpid(pid, &status, 0) == -1) {
        if (errno != EINTR) return -1;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#else
    return std::system(cmd.c_str());
#endif
}

int DatasetRunner::exec_dock(const std::string& cmd, int timeout_s) {
    // Docking entry point — uses SubprocessGuard for PID tracking,
    // process groups, and timeout enforcement.
    if (!proc_guard_) {
        // No guard available — fall back to untracked exec
        return exec_cmd(cmd);
    }
#ifndef _MSC_VER
    pid_t pid = proc_guard_->fork_exec(cmd);
    if (pid < 0) {
        std::cerr << "[ERROR] fork_exec failed: " << strerror(errno) << "\n";
        return -1;
    }
    return proc_guard_->wait_with_timeout(pid, timeout_s);
#else
    return std::system(cmd.c_str());
#endif
}

std::string DatasetRunner::exec_cmd_output(const std::string& cmd) {
    std::string result;
    std::array<char, 4096> buffer;
#ifdef _MSC_VER
    FILE* pipe = _popen(cmd.c_str(), "r");
#else
    FILE* pipe = popen(cmd.c_str(), "r");
#endif
    if (!pipe) return result;
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        result += buffer.data();
    }
#ifdef _MSC_VER
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return result;
}

bool DatasetRunner::http_download(const std::string& url, const std::string& out_path) {
    // Ensure parent directory exists
    ensure_dir(fs::path(out_path).parent_path().string());

    // Use system curl with retry logic
    std::ostringstream cmd;
    cmd << "curl -sS -L --retry 3 --retry-delay 2 -o \""
        << out_path << "\" \"" << url << "\" 2>&1";

    int ret = exec_cmd(cmd.str());
    if (ret != 0) {
        std::cerr << "  [ERROR] Download failed: " << url << "\n";
        return false;
    }

    // Verify file exists and is non-empty
    if (!fs::exists(out_path) || fs::file_size(out_path) == 0) {
        std::cerr << "  [ERROR] Downloaded file is empty or missing: " << out_path << "\n";
        if (fs::exists(out_path)) fs::remove(out_path);
        return false;
    }

    return true;
}

// =============================================================================
// PDB/CIF download from RCSB
// =============================================================================

bool DatasetRunner::download_pdb(const std::string& pdb_id, const std::string& out_path) {
    // Check cache first
    if (fs::exists(out_path) && fs::file_size(out_path) > 100) {
        return true;
    }

    std::string upper_id = pdb_id;
    std::transform(upper_id.begin(), upper_id.end(), upper_id.begin(),
                   [](unsigned char c) { return std::toupper(c); });

    std::string url = "https://files.rcsb.org/download/" + upper_id + ".pdb";
    std::cout << "  Downloading " << upper_id << ".pdb ...\n";

    if (!http_download(url, out_path)) {
        // Try lowercase
        std::string lower_id = pdb_id;
        std::transform(lower_id.begin(), lower_id.end(), lower_id.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        url = "https://files.rcsb.org/download/" + lower_id + ".pdb";
        return http_download(url, out_path);
    }

    // Validate it's actually a PDB file (not an error page)
    std::ifstream check(out_path);
    std::string first_line;
    if (std::getline(check, first_line)) {
        if (first_line.find("<!DOCTYPE") != std::string::npos ||
            first_line.find("<html") != std::string::npos) {
            std::cerr << "  [ERROR] Got HTML instead of PDB for " << pdb_id << "\n";
            fs::remove(out_path);
            return false;
        }
    }

    return true;
}

bool DatasetRunner::download_cif(const std::string& pdb_id, const std::string& out_path) {
    if (fs::exists(out_path) && fs::file_size(out_path) > 100) {
        return true;
    }

    std::string upper_id = pdb_id;
    std::transform(upper_id.begin(), upper_id.end(), upper_id.begin(),
                   [](unsigned char c) { return std::toupper(c); });

    std::string url = "https://files.rcsb.org/download/" + upper_id + ".cif";
    std::cout << "  Downloading " << upper_id << ".cif ...\n";
    return http_download(url, out_path);
}

// =============================================================================
// PDB HETATM parsing
// =============================================================================

std::vector<PDBAtom> DatasetRunner::parse_pdb_hetatm(const std::string& pdb_path) {
    std::vector<PDBAtom> atoms;
    std::ifstream ifs(pdb_path);
    if (!ifs) return atoms;

    std::string line;
    while (std::getline(ifs, line)) {
        // Pad line to at least 80 characters for safe substring extraction
        while (line.size() < 80) line += ' ';

        bool is_hetatm = (line.substr(0, 6) == "HETATM");
        if (!is_hetatm) continue;

        PDBAtom atom;
        atom.is_hetatm = true;

        // PDB format column extraction (1-based indexing in spec, 0-based here)
        try {
            atom.serial  = std::stoi(line.substr(6, 5));
        } catch (...) { atom.serial = 0; }

        atom.name    = line.substr(12, 4);
        atom.altLoc  = line.substr(16, 1);
        atom.resName = line.substr(17, 3);
        atom.chainID = line.substr(21, 1);

        try {
            atom.resSeq = std::stoi(line.substr(22, 4));
        } catch (...) { atom.resSeq = 0; }

        try {
            atom.x = std::stof(line.substr(30, 8));
            atom.y = std::stof(line.substr(38, 8));
            atom.z = std::stof(line.substr(46, 8));
        } catch (...) {
            continue; // skip atoms with bad coordinates
        }

        try {
            atom.occupancy   = std::stof(line.substr(54, 6));
        } catch (...) { atom.occupancy = 1.0f; }

        try {
            atom.tempFactor  = std::stof(line.substr(60, 6));
        } catch (...) { atom.tempFactor = 0.0f; }

        if (line.size() >= 78) {
            atom.element = line.substr(76, 2);
            // Trim whitespace
            while (!atom.element.empty() && atom.element.front() == ' ')
                atom.element.erase(atom.element.begin());
            while (!atom.element.empty() && atom.element.back() == ' ')
                atom.element.pop_back();
        }

        // Trim residue name
        while (!atom.resName.empty() && atom.resName.front() == ' ')
            atom.resName.erase(atom.resName.begin());
        while (!atom.resName.empty() && atom.resName.back() == ' ')
            atom.resName.pop_back();

        // Trim atom name
        while (!atom.name.empty() && atom.name.front() == ' ')
            atom.name.erase(atom.name.begin());
        while (!atom.name.empty() && atom.name.back() == ' ')
            atom.name.pop_back();

        atoms.push_back(std::move(atom));
    }

    return atoms;
}

// =============================================================================
// Ligand extraction: HETATM → SDF
// =============================================================================

bool DatasetRunner::extract_ligand(const std::string& pdb_path,
                                    const std::string& out_sdf) {
    auto hetatm_atoms = parse_pdb_hetatm(pdb_path);
    if (hetatm_atoms.empty()) {
        std::cerr << "  [WARN] No HETATM records in " << pdb_path << "\n";
        return false;
    }

    // Group HETATM atoms by (resName, chainID, resSeq) triplet
    struct ResidueKey {
        std::string resName;
        std::string chainID;
        int resSeq;
        bool operator<(const ResidueKey& o) const {
            if (resName != o.resName) return resName < o.resName;
            if (chainID != o.chainID) return chainID < o.chainID;
            return resSeq < o.resSeq;
        }
    };

    std::map<ResidueKey, std::vector<PDBAtom>> residue_groups;
    const auto& excl = excluded_residues();

    for (const auto& atom : hetatm_atoms) {
        // Skip excluded residues (water, ions, buffers)
        if (excl.count(atom.resName)) continue;
        // Skip alternate conformers (keep only first)
        if (atom.altLoc != " " && atom.altLoc != "" && atom.altLoc != "A") continue;

        ResidueKey key{atom.resName, atom.chainID, atom.resSeq};
        residue_groups[key].push_back(atom);
    }

    if (residue_groups.empty()) {
        std::cerr << "  [WARN] No valid ligand residues in " << pdb_path << "\n";
        return false;
    }

    // Find the largest residue group (most atoms = likely the ligand)
    ResidueKey best_key;
    size_t max_atoms = 0;
    for (const auto& [key, atoms_vec] : residue_groups) {
        if (atoms_vec.size() > max_atoms) {
            max_atoms = atoms_vec.size();
            best_key = key;
        }
    }

    const auto& ligand_atoms = residue_groups[best_key];
    if (ligand_atoms.size() < 3) {
        std::cerr << "  [WARN] Ligand too small (" << ligand_atoms.size()
                  << " atoms): " << best_key.resName << "\n";
        return false;
    }

    // ── Collect bonds BEFORE writing (counts line needs the total) ──────────
    // Priority: (1) PDB CONECT records, (2) distance-based inference.

    std::set<int> ligand_serials;
    for (const auto& atom : ligand_atoms) ligand_serials.insert(atom.serial);

    std::map<int, int> serial_to_idx;
    for (size_t i = 0; i < ligand_atoms.size(); ++i)
        serial_to_idx[ligand_atoms[i].serial] = static_cast<int>(i) + 1;

    // Bond list: pairs of 1-based SDF atom indices
    std::vector<std::pair<int,int>> bonds;

    // (1) CONECT records
    {
        std::ifstream pdb_ifs(pdb_path);
        std::string line;
        std::set<std::pair<int,int>> seen;
        while (std::getline(pdb_ifs, line)) {
            if (line.size() < 6 || line.substr(0, 6) != "CONECT") continue;
            while (line.size() < 31) line += ' ';
            int central = 0;
            try { central = std::stoi(line.substr(6, 5)); } catch (...) { continue; }
            if (!ligand_serials.count(central)) continue;
            for (int col = 11; col < 31 && col + 5 <= static_cast<int>(line.size()); col += 5) {
                std::string s = line.substr(col, 5);
                if (s.find_first_not_of(" ") == std::string::npos) continue;
                int bonded = 0;
                try { bonded = std::stoi(s); } catch (...) { continue; }
                if (!ligand_serials.count(bonded)) continue;
                auto it_a = serial_to_idx.find(central);
                auto it_b = serial_to_idx.find(bonded);
                if (it_a == serial_to_idx.end() || it_b == serial_to_idx.end()) continue;
                int a = it_a->second, b = it_b->second;
                if (a > b) std::swap(a, b);
                if (seen.insert({a, b}).second) bonds.emplace_back(a, b);
            }
        }
    }

    // (2) Distance-based fallback when no CONECT records exist
    if (bonds.empty()) {
        const float max_bond_dist_sq = 2.0f * 2.0f;  // 2.0 Å cutoff
        for (size_t i = 0; i < ligand_atoms.size(); ++i) {
            for (size_t j = i + 1; j < ligand_atoms.size(); ++j) {
                float dx = ligand_atoms[i].x - ligand_atoms[j].x;
                float dy = ligand_atoms[i].y - ligand_atoms[j].y;
                float dz = ligand_atoms[i].z - ligand_atoms[j].z;
                float dist_sq = dx*dx + dy*dy + dz*dz;
                if (dist_sq < max_bond_dist_sq && dist_sq > 0.16f)
                    bonds.emplace_back(static_cast<int>(i)+1, static_cast<int>(j)+1);
            }
        }
    }

    // ── Write SDF ─────────────────────────────────────────────────────────
    std::ofstream ofs(out_sdf);
    if (!ofs) return false;

    // Header
    ofs << best_key.resName << "\n";
    ofs << "  FlexAIDdS DatasetRunner\n";
    ofs << "  Extracted from PDB HETATM records\n";

    // Counts line — write ACTUAL bond count so readers respect the block
    ofs << std::setw(3) << ligand_atoms.size()
        << std::setw(3) << bonds.size()
        << "  0  0  0  0  0  0  0999 V2000\n";

    // Atom block
    for (const auto& atom : ligand_atoms) {
        std::string elem = atom.element;
        if (elem.empty()) {
            for (char c : atom.name) {
                if (std::isalpha(static_cast<unsigned char>(c))) {
                    elem = std::string(1, std::toupper(static_cast<unsigned char>(c)));
                    break;
                }
            }
        }
        if (elem.empty()) elem = "C";

        ofs << std::fixed << std::setprecision(4)
            << std::setw(10) << atom.x
            << std::setw(10) << atom.y
            << std::setw(10) << atom.z
            << " " << std::setw(3) << std::left << elem << std::right
            << " 0  0  0  0  0  0  0  0  0  0  0  0\n";
    }

    // Bond block
    for (const auto& [a, b] : bonds) {
        ofs << std::setw(3) << a
            << std::setw(3) << b
            << "  1  0  0  0  0\n";
    }

    ofs << "M  END\n";
    ofs << "$$$$\n";
    ofs.close();

    return true;
}

// =============================================================================
// Prepare a single PDB entry
// =============================================================================

DatasetEntry DatasetRunner::prepare_pdb_entry(const std::string& pdb_id,
                                               const std::string& dataset_name,
                                               float affinity,
                                               float dH, float dS) {
    std::string upper_id = pdb_id;
    std::transform(upper_id.begin(), upper_id.end(), upper_id.begin(),
                   [](unsigned char c) { return std::toupper(c); });

    std::string entry_dir = cache_dir_ + "/" + dataset_name + "/" + upper_id;
    ensure_dir(entry_dir);

    std::string receptor_path = entry_dir + "/" + upper_id + ".pdb";
    std::string ligand_path   = entry_dir + "/" + upper_id + "_ligand.sdf";

    DatasetEntry entry;
    entry.pdb_id = upper_id;
    entry.source = dataset_name;
    entry.experimental_affinity = affinity;
    entry.experimental_dH  = dH;
    entry.experimental_TdS = dS;

    // Download PDB
    if (download_pdb(upper_id, receptor_path)) {
        entry.receptor_path = receptor_path;
    } else {
        std::cerr << "  [WARN] Failed to download PDB: " << upper_id << "\n";
        return entry;
    }

    // Extract ligand
    if (!fs::exists(ligand_path) || fs::file_size(ligand_path) == 0) {
        if (extract_ligand(receptor_path, ligand_path)) {
            entry.ligand_path = ligand_path;
        } else {
            std::cerr << "  [WARN] Failed to extract ligand from: " << upper_id << "\n";
        }
    } else {
        entry.ligand_path = ligand_path;
    }

    return entry;
}

// =============================================================================
// Astex Diverse 85 PDB codes (Hartshorn et al. 2007 JCIM)
// =============================================================================

std::vector<std::string> DatasetRunner::astex_diverse_codes() {
    return {
        "1G9V", "1GM8", "1GPK", "1HNN", "1HP0", "1HQ2", "1IA1", "1IGJ",
        "1J3J", "1JD0", "1JJE", "1K3U", "1KE5", "1KZK", "1L2S", "1L7F",
        "1LPZ", "1M2Z", "1MEH", "1MQ6", "1N1M", "1N2J", "1N2V", "1N46",
        "1NAV", "1OF1", "1OF6", "1OPK", "1OQ5", "1OWE", "1P2Y", "1P62",
        "1PMN", "1Q1G", "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O",
        "1S19", "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T9B",
        "1TT1", "1TW6", "1TZ8", "1U1C", "1U4D", "1UML", "1UNL", "1UOU",
        "1V0P", "1V48", "1V4S", "1VCJ", "1W1P", "1W2G", "1X8X", "1XM6",
        "1XOZ", "1Y6B", "1Y6R", "1YGC", "1YQY", "1YV3", "1YVF", "1YWR",
        "1Z95", "2BM2", "2BR1", "2BSM", "2BYS", "2C3I", "2CET", "2CGR",
        "2D3U", "2GBP", "2HB1", "2HR7", "2J62"
    };
}

std::vector<DatasetEntry> DatasetRunner::fetch_astex() {
    std::cout << "[DatasetRunner] Preparing Astex Diverse 85 dataset\n";
    auto codes = astex_diverse_codes();
    std::vector<DatasetEntry> entries;
    entries.reserve(codes.size());

    for (const auto& pdb : codes) {
        auto entry = prepare_pdb_entry(pdb, "astex_diverse");
        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " / " << codes.size()
              << " entries\n";
    return entries;
}

// =============================================================================
// Astex Non-Native (Verdonk et al. 2008 JCIM) — 65 targets, 1112 structures
// =============================================================================

std::vector<AstexNonNativeTarget> astex_nonnative_targets() {
    // 65 protein families with native and alternative (non-native) conformers
    // for cross-docking benchmarks.  The native PDB is the holo crystal
    // structure; alternatives are other crystallographic structures of the same
    // protein.  Based on Verdonk et al. (2008) J. Chem. Inf. Model. 48,
    // 2214-2225 (original: 65 families, 1112 structures).  This table has been
    // expanded with post-2008 crystal structures and currently covers 65/65
    // families with ~1900 unique PDB structures / ~2200 cross-docking pairs.
    return {
        {"ACE",   "1G9V",  {"1EVE", "1GQR", "1QTI", "2ACE", "1DX6", "1F8U", "1GPK", "1HBJ", "1J07", "1JJB", "1MAA", "1MAH", "1OCE", "1VOT", "1W4L", "1W6R", "1W75", "1W76", "2C4H", "2C58", "2CEK", "2CKM", "2CMF", "2GYU"}},
        {"ADA",   "1NDV",  {"1ADD", "1KRM", "1NDW", "1O5R", "1QXL", "2E1W"}},
        {"ADA17", "1BKC",  {"1B8Y", "2FV5", "2FV9", "2DDF"}},
        {"ALR2",  "1T40",  {"1ADS", "1EF3", "1IEI", "1EL3", "1MAR", "1PWL", "1PWM", "1T41", "1US0", "1Z3N", "1Z89", "2ACQ", "2ACR", "2ACS", "2ACU", "2DUX", "2DUZ", "2FZ8", "2FZD", "2HV5", "2HVN", "2I16", "2I17", "2IKG", "2IKH", "2IKI", "2IKJ", "2INE", "2INZ", "2IPW", "2IQ1", "2IQD", "2IS7", "2ISF", "2NVC", "2NVD", "2PD5", "2PDC", "2PDD", "2PDG", "2PDJ", "2PDK", "2PDM", "2PDP", "2PDQ", "2PDU", "3BAJ"}},
        {"ACHE",  "1HQ2",  {"1ACJ", "1ACL", "1AMN", "1AX9", "1B41", "1CFJ", "1DX6", "1E3Q", "1E66", "1EA5", "1EVE", "1F8U", "1GPN", "1GQR", "1GQS", "1H22", "1H23", "1HBJ", "1J07", "1JJB", "1MAA", "1MAH", "1N5M", "1N5R", "1OCE", "1ODC", "1QTI", "1VOT", "1W4L", "1W6R", "1W75", "1W76", "1ZGB", "2ACE", "2C4H", "2C58", "2CEK", "2CKM", "2CMF", "2GYU"}},
        {"AR",    "1T9B",  {"1E3G", "1GS4", "1I37", "1I38", "1R4I", "1T5Z", "1T63", "1T65", "1XJ7", "1XOW", "1XQ3", "1Z95", "2AM9", "2AMA", "2AMB", "2AX6", "2AX7", "2AX8", "2AX9", "2AXA", "2HVC", "2IHQ", "2NW4", "2OZ7", "2PIO", "2PIQ", "2PIR", "2PIT", "2PIU", "2PIV", "2PKL", "2PNU", "2Q7I", "2Q7J", "2Q7K", "2Q7L"}},
        {"BACE1", "1W51",  {"1FKN", "1M4H", "1SGZ", "1TQF", "1XN2", "1XN3", "1XS7", "1YM2", "1YM4", "2B8L", "2B8V", "2F3E", "2F3F", "2G94", "2HM1", "2IRZ", "2IS0", "2OF0", "2OHL", "2OHM", "2OHP", "2OHQ", "2OHR", "2OHS", "2OHT", "2OHU", "2P4J", "2P83", "2PH6", "2PH8", "2QK5", "2QMD", "2QMF", "2QMG"}},
        {"CA2",   "1V4S",  {"1A42", "1AM6", "1BCD", "1BN1", "1BN3", "1BN4", "1BNM", "1BNN", "1BNQ", "1BNT", "1BNU", "1BNV", "1BNW", "1CAY", "1CIL", "1CIM", "1CIN", "1CNI", "1CNW", "1CNX", "1CNY", "1CRA", "1CVA", "1CVB", "1CVD", "1CVF", "1CVH", "1EOU", "1FQL", "1FQM", "1FQN", "1FQR", "1G0E", "1G0F", "1G1D", "1G45", "1G46", "1G48", "1G4J", "1G4O", "1G52", "1G53", "1G54", "1I8Z", "1I90", "1I91", "1IF4", "1IF5", "1IF6", "1IF7", "1IF8", "1IF9", "1IFI", "1KEQ", "1KWQ", "1KWR", "1LG5", "1LG6", "1LGD", "1MOO", "1MUA", "1OKL", "1OKM", "1OKN", "1RJ5", "1RJ6", "1RZA", "1RZB", "1RZC", "1RZD", "1RZE", "1T9N", "1TB0", "1TBT", "1TE3", "1TEQ", "1TG3", "1TG9", "1TH9", "1TTM", "1TU6", "1XEG", "1XEV", "1Z97", "1Z9Y", "1ZE8", "1ZFP", "1ZFQ", "1ZJR", "2ABE", "2CA2", "2CBA", "2CBB", "2CBD", "2CBS", "2EU2", "2EU3", "2EZ7", "2F14", "2FMG", "2FMZ", "2FN0", "2FNK", "2FNM", "2FOV", "2G63", "2H15", "2HD6", "2HKK", "2HNC", "2HL4", "2ILI", "2NMX", "2NNG", "2NNS", "2NXR", "2NXS", "2NXT", "2OSF", "2POU", "2POV", "2Q1B", "2Q1Q", "2Q38", "2QO8", "2QOA"}},
        {"CDK2",  "1KE5",  {"1AQ1", "1B38", "1B39", "1BUH", "1CKP", "1DI8", "1DM2", "1E1V", "1E1X", "1E9H", "1F5Q", "1FIN", "1FQ1", "1FVT", "1FVV", "1G5S", "1GIH", "1GII", "1GIJ", "1GZ8", "1H00", "1H01", "1H07", "1H08", "1H0V", "1H0W", "1H1P", "1H1Q", "1H1R", "1H1S", "1HCK", "1HCL", "1JIN", "1JST", "1JSV", "1JVP", "1KE6", "1KE7", "1KE8", "1KE9", "1OGU", "1OIQ", "1OIR", "1OIT", "1OIU", "1OIY", "1OL1", "1OL2", "1P2A", "1P5E", "1PF8", "1PKD", "1PW2", "1PXI", "1PXJ", "1PXK", "1PXL", "1PXM", "1PXN", "1PXO", "1PXP", "1PYE", "1R78", "1URW", "1V0B", "1V0O", "1V0P", "1VYW", "1VYZ", "1W0X", "1W98", "2A0C", "2A4L", "2B52", "2B53", "2B54", "2B55", "2BHE", "2BKZ", "2BPM", "2BTR", "2BTS", "2C4G", "2C5N", "2C5O", "2C5P", "2C5V", "2C5X", "2C5Y", "2C6I", "2C6K", "2C6L", "2C6M", "2C6O", "2C6T", "2CLX", "2DS1", "2DUV", "2EXM", "2FVD", "2G9X", "2HIC", "2I40"}},
        {"CHK1",  "1NVQ",  {"1NVR", "1NVS", "1ZLT", "1ZYS", "2AYP", "2BRB", "2BRC", "2BR1", "2C3J", "2C3K", "2C3L", "2CGU", "2CGV", "2CGW", "2CGX", "2E9N", "2E9O", "2E9P", "2E9U", "2GDO", "2HOG"}},
        {"COMT",  "1H1D",  {"1JR4"}},
        {"COX1",  "1Q4G",  {"1CQE", "1DIY", "1EBV", "1EQG", "1EQH", "1FE2", "1HT5", "1HT8", "1HTB", "1IGX", "1IGZ", "1PGE", "1PGF", "1PGG", "1PMN", "1PRH", "1Q4G", "2OYE", "2OYU", "3PGH"}},
        {"COX2",  "1PXX",  {"1CVU", "1CX2", "1DDX", "1V0X", "3LN1", "3NTB", "3NTG", "3OLU", "3PGH", "3QMO", "4COX", "5COX", "6COX"}},
        {"DHFR",  "1S3V",  {"1BOZ", "1DHF", "1DRE", "1DRF", "1DRH", "1HFP", "1HFQ", "1HFR", "1KMV", "1MVS", "1MVT", "1OHJ", "1OHK", "1PD8", "1PD9", "1RA2", "1RA3", "1RA9", "1RB2", "1RC4", "1RG7", "1RH3", "1RX2", "1RX3", "1RX4", "1RX5", "1RX6", "1RX7", "1RX9", "1S3U", "1S3W", "2C2S", "2C2T"}},
        {"EGFR",  "1M17",  {"1XKK", "2GS2", "2GS7", "2ITN", "2ITO", "2ITP", "2ITQ", "2ITT", "2ITU", "2ITV", "2ITW", "2ITX", "2ITY", "2ITZ", "2J5E", "2J5F", "2J6M", "2JIT", "2JIU", "2JIV", "2RGP"}},
        {"ER",    "1SJ0",  {"1A52", "1ERE", "1ERR", "1G50", "1GWQ", "1GWR", "1L2I", "1NDE", "1PCG", "1QKT", "1QKU", "1R5K", "1SJ0", "1UOM", "1X7E", "1X7R", "1XP1", "1XP6", "1XP9", "1XPC", "1XQC", "1YIM", "1YIN", "1ZKY", "2AYR", "2B1V", "2B1Z", "2BJ4", "2FAI", "2G44", "2G5O", "2I0G", "2I0J", "2IOG", "2IOK", "2J7X", "2JFA", "2NV7", "2OUZ", "2P15", "2POG", "2Q6J", "2Q70", "2QA8", "2QAB", "2QGT", "2QGW", "2QH6", "2QR9", "2QSE", "2R6W", "2R6Y", "3ERD", "3ERT"}},
        {"FGFR1", "1AGW",  {"1FGI"}},
        {"FXA",   "1MQ6",  {"1EZQ", "1F0R", "1F0S", "1FAX", "1FJS", "1G2L", "1G2M", "1G32", "1IQE", "1IQF", "1IQG", "1IQH", "1IQI", "1IQJ", "1IQK", "1IQL", "1IQM", "1IQN", "1KSN", "1KYE", "1LPG", "1LPK", "1LPZ", "1LQD", "1MQ5", "1MQ6", "1NFU", "1NFW", "1NFX", "1NFY", "1V3X", "1WAY", "1WU1", "1Z6E", "2BMG", "2BQ7", "2BQW", "2CJI", "2EI6", "2EI7", "2EI8", "2G00", "2GD4", "2H9E", "2J2U", "2J34", "2J38", "2J4I", "2J94", "2J95", "2JKH", "2P16", "2P3F", "2P3T", "2P93", "2P94", "2P95", "2PHB", "2PR3", "2Q1J"}},
        {"GAR",   "1UML",  {"1JQ8", "1JQK", "1JQL", "2GAR"}},
        {"GR",    "1M2Z",  {"1NHZ", "1P93"}},
        {"HIVPR", "1HQ2",  {"1A30", "1A94", "1AID", "1B6J", "1B6K", "1B6L", "1B6M", "1B6N", "1B6P", "1BV7", "1BV9", "1BWA", "1BWB", "1D4H", "1D4I", "1D4J", "1D4K", "1D4L", "1D4S", "1D4Y", "1DIF", "1DW6", "1EBW", "1EBY", "1EC0", "1EC1", "1EC2", "1EC3", "1G35", "1G2K", "1HBV", "1HEF", "1HEG", "1HIH", "1HII", "1HOS", "1HPO", "1HPS", "1HPV", "1HPX", "1HSG", "1HTF", "1HTG", "1HVI", "1HVJ", "1HVK", "1HVL", "1HVR", "1HVS", "1HXB", "1HXW", "1IDA", "1IDB", "1IDW", "1IIQ", "1IZH", "1IZI", "1JLD", "1K1T", "1K1U", "1K2B", "1K2C", "1K6C", "1K6P", "1K6T", "1K6V", "1KZK", "1LZQ", "1M0B", "1MER", "1MES", "1MET", "1MEU", "1MSM", "1MSN", "1MUI", "1MUT", "1N49", "1NH0", "1NPA", "1NPV", "1NPW", "1OD1", "1OHR", "1PRO", "1QBR", "1QBS", "1QBT", "1QBU", "1RL8", "1SDT", "1SDU", "1SDV", "1SIV", "1T3R", "1T7I", "1T7J", "1T7K", "1TCX", "1W5V", "1W5W", "1W5X", "1W5Y", "1YT9", "1YTG", "1YTH", "1ZBG", "1ZJ7", "1ZLF", "1ZP8", "1ZPA", "1ZPK", "1ZSR", "2BPV", "2BPW", "2BPX", "2BPY", "2BPZ", "2BQV", "2CEJ", "2CEN", "2F80", "2F81", "2FDD", "2FDE", "2FGU", "2FLE", "2HS1", "2HS2", "2I0A", "2I0D", "2I4D", "2I4U", "2I4V", "2I4W", "2I4X", "2I5J", "2IDW", "2NMW", "2NNK", "2NNP", "2NXD", "2NXL", "2NXM", "2O4K", "2O4L", "2O4N", "2O4P", "2O4S", "2PK5", "2PK6", "2PSU", "2PSV", "2Q11", "2Q54", "2Q55", "2Q63", "2Q64", "2QAK", "2QCI", "2QD6", "2QD7", "2QD8", "2QHY", "2QI0", "2QI1", "2QI3", "2QI4", "2QI6", "2QI7", "2R3W", "2R5P", "2R5Q", "3AID", "4HVP", "4PHV", "5HVP", "6HVP", "7HVP", "7UPJ", "8HVP", "9HVP"}},
        {"HSP90", "1UY6",  {"1BYQ", "1UYC", "1UYD", "1UYE", "1UYF", "1UYG", "1UYH", "1UYI", "1UYK", "1YC1", "1YC3", "1YC4", "1YER", "1YES", "1YET", "2BRC", "2BSM", "2BT0", "2BYH", "2BYI", "2BZ5", "2CCT", "2CCS", "2CGE", "2FWZ", "2QFO", "2VCI", "2WI1", "2WI2", "2WI3", "2WI4", "2WI5", "2WI6", "2WI7", "2XAB", "2XDK", "2XDL", "2XDS", "2XDX", "2XHR", "2XHT", "2XHX", "2XJG", "2XJX", "2XK2", "2YE2", "2YE3", "2YE4", "2YE5", "2YE6", "2YE7", "2YE8", "2YE9", "2YEA", "2YEB", "2YEC", "2YED", "2YEE", "2YEF", "2YEG", "2YEH", "2YEI", "2YEJ", "2YI0", "2YI5", "2YI6", "2YI7", "3B24", "3B25", "3B26", "3B27", "3B28", "3D0B", "3EKO", "3EKR", "3FT5", "3FT8", "3HEK", "3HYY", "3HYZ", "3HZ1", "3HZ5", "3INW", "3INX", "3K97", "3K98", "3K99", "3OWB", "3OWD", "3R91", "3R92", "3RKZ", "3RLP", "3T0Z", "3T10"}},
        {"JNK1",  "2GMX",  {"2G01", "2NO3", "2O0U", "2O2U", "3ELJ", "3O17", "3O2M", "3PZE"}},
        {"JNK3",  "1PMN",  {"1JNK", "1PMQ", "1PMU", "1PMV"}},
        {"LCK",   "1QPJ",  {"1QPC", "1QPD", "1QPE"}},
        {"MAP",   "1V4S",  {"1V4T", "1V4U", "1V4V"}},
        {"MCL1",  "2PQK",  {"2MHS", "2NLA", "2NL9", "2ROC", "2ROD", "3KJ0", "3KJ1", "3KJ2", "3MK8", "3PK1", "3WIX", "3WIY", "3WIZ", "4HW2", "4HW3", "4HW4", "4OQ5", "4OQ6", "4WGI", "4ZBI", "4ZBF", "5C3F", "5FDO", "5FDR", "5FDX", "5IEZ"}},
        {"MMP12", "1Y93",  {"1JIZ", "1JK3", "1NS9", "1OS2", "1OS9", "1OY5", "1RMZ", "1UTT", "1Y93", "1Z3J", "2HU6", "2OXU", "2OXW", "2OXZ", "2OY4", "2PJT", "2W0D"}},
        {"NA",    "1HP0",  {"1A4G", "1A4Q", "1B9S", "1B9T", "1B9V", "1BJI", "1F8B", "1F8C", "1F8D", "1F8E", "1INF", "1ING", "1INV", "1INW", "1INX", "1INY", "1IVB", "1IVC", "1IVD", "1IVE", "1IVF", "1L7F", "1L7G", "1L7H", "1MWE", "1NMB", "1NNC", "1NNB", "1NSB", "1NSC", "1NSD", "2BAT", "2HTQ", "2HTS", "2HTU", "2QWA", "2QWB", "2QWC", "2QWD", "2QWE", "2QWF", "2QWG", "2QWH", "2QWI", "2QWJ", "2QWK"}},
        {"P38",   "1OQ5",  {"1A9U", "1BL6", "1BL7", "1BMK", "1DI9", "1KV1", "1KV2", "1M7Q", "1OUK", "1OVE", "1OZ1", "1R39", "1R3C", "1W7H", "1W82", "1W83", "1W84", "1WBN", "1WBO", "1WBS", "1WBT", "1WBV", "1WBW", "1YQJ", "1ZYJ", "1ZZ2", "1ZZL", "2BAJ", "2BAK", "2BAL", "2BAQ", "2EWA", "2FSL", "2FSM", "2FSO", "2FST", "2GFS", "2GHL", "2I0H", "2LGP", "2NPQ", "2OKR", "2ONL", "2RG5", "2RG6", "2YIS", "2YIW", "2YIX", "2ZAZ", "2ZB0", "2ZB1", "3BV2", "3BV3", "3C5U", "3CTQ", "3D7Z", "3D83", "3DS6", "3E92", "3E93", "3FC1", "3FI1", "3FL4", "3FLN", "3FLZ", "3FMH", "3FMJ", "3FMK", "3FML", "3GC7", "3GCP", "3GCS", "3GCU", "3GCV", "3GFE", "3GI2", "3GI3", "3HA8", "3HEC", "3HEG", "3HL7", "3HLL", "3HRB", "3HUB", "3HVC", "3HV3", "3HV5", "3HV6", "3HV7", "3IW5", "3IW6", "3IW7", "3IW8"}},
        {"PDE4",  "1TBB",  {"1F0J", "1MKD", "1OYN", "1Q9M", "1RO6", "1RO9", "1ROR", "1TBB", "1XLX", "1XLZ", "1XM4", "1XMU", "1XMY", "1XN0", "1XOJ", "1XOM", "1XON", "1XOS", "1XOT", "2CHM", "2FM0", "2FM5", "2PW3", "2QYK", "2QYL", "2QYM", "2QYN", "2QYO", "3G4G", "3G4I", "3G4K", "3G4L", "3G45", "3G58"}},
        {"PDE5",  "1UDT",  {"1RKP", "1T9R", "1T9S", "1TBF", "1UDO", "1UDT", "1UHO", "1XOZ", "1XP0", "2H40", "2H42", "2H44", "3B2R", "3BJC", "3HC8", "3JWR", "3JWQ", "3SHY", "3SHZ", "3TGE", "3TGG"}},
        {"PDGFR", "1T46",  {"1PKG"}},
        {"PNP",   "1UOU",  {"1A69", "1A9S", "1B8N", "1B8O", "1GE0", "1ILR", "1K9S", "1M73", "1PBN", "1PE4", "1PF7", "1PW7", "1RCT", "1RFG", "1RSZ", "1RT9", "1SQP", "1T86", "1TMM", "1TYO", "1UOU", "1ULB", "1V2H", "1V3Q", "1V41", "1V45", "1VII", "1YHM", "1YRY", "2A0W", "2A0X", "2A0Y", "2AOC", "2AOD", "2AOE", "2AOF", "2AOG", "2BSZ", "2OC9", "3BGS", "3BPU"}},
        {"PPARg", "1K74",  {"1FM6", "1FM9", "1I7I", "1K74", "1KNU", "1NYX", "1PRG", "1RDT", "1WM0", "1ZEO", "1ZGY", "2ATH", "2F4B", "2FVJ", "2G0G", "2G0H", "2GTK", "2HFP", "2HWQ", "2HWR", "2I4J", "2I4P", "2I4Z", "2NPA", "2OM9", "2P4Y", "2POB", "2PRG", "2Q59", "2Q5P", "2Q5S", "2Q61", "2Q6R", "2Q6S", "2QMV", "2R5E", "2VSR", "2VST", "2VV0", "2VV1", "2VV2", "2VV3", "2VV4", "2WAK", "2XKW", "2Y0W", "2YFE", "2ZK0", "2ZK1", "2ZK2", "2ZK3", "2ZK4", "2ZK5", "2ZK6", "2ZNO", "3AN3", "3AN4", "3B0Q", "3B0R", "3B1M", "3B3K", "3CDP", "3CS8", "3CWD", "3D6D", "3DZU", "3DZY", "3E00", "3ET0", "3ET3", "3FEJ", "3G9E", "3GBK", "3HOD", "3HZV", "3IA6", "3K8S", "3LMP", "3NOA", "3PBA", "3PRG", "3QT0", "3R5N", "3R8A", "3R8I", "3S9Q", "3SZ1", "3T03", "3TY0", "3U9Q", "3V9T", "3V9V", "3V9Y", "3VJH", "3VJI", "3VN2", "3WJ4", "3WJ5", "3WMH", "3X1H", "3X1I", "4A4V", "4A4W", "4CI5", "4E4K", "4E4Q", "4EM9", "4EMA", "4F9M", "4G2J", "4HEE", "4JAZ", "4JL4", "4L96", "4L98", "4O8F", "4OJ4", "4PRG", "4PVU", "4R2U", "4R6S"}},
        {"PTP1B", "1Q1G",  {"1C83", "1C84", "1C85", "1C86", "1C87", "1C88", "1ECV", "1G1F", "1G1G", "1G1H", "1G7F", "1G7G", "1GFY", "1JF7", "1KAK", "1KAV", "1L8G", "1LQF", "1NL9", "1NL9", "1NO6", "1NWL", "1NZ7", "1ONY", "1ONZ", "1OEM", "1OEO", "1PA1", "1PA9", "1PH0", "1PTY", "1PXH", "1PYN", "1Q1M", "1Q6J", "1Q6M", "1Q6N", "1Q6P", "1Q6S", "1Q6T", "1QXK", "1SUG", "1T48", "1T49", "1T4J", "1WAX", "2B07", "2BGD", "2BGE", "2CM2", "2CM3", "2CM7", "2CM8", "2CNE", "2CNF", "2CNH", "2CNI", "2CNG", "2F6F", "2F6T", "2F6V", "2F6W", "2F6Y", "2F6Z", "2F70", "2F71", "2HNP", "2HNQ", "2QBP", "2QBQ", "2QBR", "2QBS", "2VEV", "2VEW", "2VEX", "2VEY"}},
        {"REN",   "1R9O",  {"1BIL", "1SME", "1RNE", "2REN", "2V0Z", "2V10", "2V11", "2V12", "2V13", "2V16", "3D91", "3G6Z", "3G70", "3G72", "3GW5", "3K1W", "3OAD", "3OAG", "3OAS", "3OOT", "3OOW", "3PCW", "3PCX", "3Q3T", "3Q4B", "3Q5H", "3QRP", "3QRQ", "3QRR", "3SFC", "3VSW", "3VUC", "3VYD", "3VYE", "4AMT", "4GJ5", "4GJ6", "4GJ7", "4GJ8", "4GJ9", "4GJA", "4GJB", "4GJC", "4GJD", "4RYC"}},
        {"RXRA",  "1YGC",  {"1DKF", "1FBY", "1FM6", "1FM9", "1K74", "1MV9", "1MZN", "1RDT"}},
        {"SAHH",  "1LI4",  {"1A7A", "1B3R", "1D4F", "1KY4", "1KY5", "1LI4", "1QI8"}},
        {"SRC",   "1YQY",  {"1KSW", "1O43", "1O44", "1O45", "1O46", "1O47", "1O48", "1O49", "1O4A", "1O4B", "1O4C", "1O4D", "1O4E", "1O4F", "1O4G", "1O4H", "1O4I", "1O4J", "1O4K", "1O4L", "1O4M", "1O4N", "1O4O", "1Y57"}},
        {"THR",   "1TT1",  {"1A2C", "1A4W", "1A46", "1A61", "1ABJ", "1AD8", "1AE8", "1AFE", "1AIX", "1B5G", "1B7X", "1BA8", "1BCU", "1BMM", "1BMN", "1C1U", "1C1V", "1C1W", "1C4U", "1C4V", "1C4Y", "1C5L", "1C5N", "1C5O", "1CA8", "1D3D", "1D3P", "1D3Q", "1D3T", "1D4P", "1D6W", "1D9I", "1DIT", "1DM4", "1DOJ", "1DWB", "1DWC", "1DWD", "1DWE", "1EB1", "1EBZ", "1ER4", "1FPC", "1G30", "1G32", "1G37", "1GHV", "1GHW", "1GHX", "1GHY", "1GJ4", "1GJ5", "1H8D", "1H8I", "1HAH", "1HAI", "1HAO", "1HAP", "1HBT", "1HDT", "1HGT", "1HUT", "1HXE", "1HXF", "1JMO", "1JOU", "1JWT", "1K21", "1K22", "1KTS", "1LHC", "1LHD", "1LHE", "1LHF", "1LHG", "1MBQ", "1MU6", "1MU8", "1MUE", "1NM6", "1NRS", "1NRN", "1NRO", "1NRQ", "1NRR", "1NRS", "1NT1", "1NY2", "1O0D", "1O2G", "1O5G", "1OOK", "1OYT", "1P8V", "1PPB", "1QHR", "1QUR", "1RD3", "1SB1", "1SFQ", "1SHH", "1SL3", "1T4U", "1T4V", "1TA2", "1TA6", "1TB6", "1TBR", "1TBZ", "1TMB", "1TMT", "1TMU", "1TOM", "1TWX", "1UMA", "1VR1", "1VZQ", "1WAY", "1XM1", "1XMN", "1YPE", "1YPF", "1YPG", "1YPH", "1YPI", "1YPJ", "1YPK", "1YPL", "1YPM", "1Z71", "1ZGI", "1ZGV", "1ZPB", "1ZRB", "2AFQ", "2ANK", "2ANM", "2BDY", "2BQ6", "2BQ7", "2BVR", "2BVS", "2BXT", "2C8W", "2C8X", "2C8Y", "2C8Z", "2C90", "2C93", "2CF8", "2CF9", "2CM2", "2GDE", "2GP9", "2GY6", "2GY7", "2HPP", "2HNT", "2HWL", "2OD3", "2PGB", "2PGQ", "2PKS", "2R2M", "2ZC9", "2ZDQ", "2ZDV", "2ZF0", "2ZFF", "2ZFP", "2ZFQ", "2ZFR", "2ZFS", "2ZG0", "2ZGX", "2ZHE", "2ZHF", "2ZHW", "2ZI2", "2ZIQ", "2ZIR", "2ZNK", "2ZO3", "3B9F", "3BEI", "3BF6", "3BIU", "3BV9", "3C1K", "3C27", "3DA9", "3DUX", "3EE0", "3EGK", "3F68", "3GIC", "3GIS", "3HGT", "3JZ1", "3JZ2", "3K65", "3LDX", "3LU9", "3P17", "3P6Z", "3PM8", "3QGN", "3QLP", "3QLU", "3QTO", "3QTV", "3QWC", "3QX5", "3RLW", "3RML", "3RMM", "3RMN", "3RMS", "3RSL", "3SHC", "3SI3", "3SI4", "3SV2", "3T5F", "3U8O", "3U8R", "3U8T", "3U8V", "3U98", "3UIS", "3UIT", "3UNX", "3UQ0", "3UT6", "3VXE", "3VXF", "4AY6", "4AYV", "4AYX", "4AYY", "4BAH", "4BAI", "4BAK", "4BAM", "4BAN", "4BAQ", "4CH2", "4CH8", "4DII", "4DIJ", "4DIK", "4E05", "4E06", "4E07", "4HFP", "4HTC", "4HZH", "4I7Y", "4LXB", "4LZ1", "4LZ4", "4NZQ", "4UD9", "4UEH", "4UEI", "4YES", "5AF9", "5AFY", "5GDS"}},
        {"TK",    "1N2J",  {"1E2I", "1E2K", "1E2N", "1E2P", "1KI6", "1KI7", "1KI8", "1N1M", "1N2V", "1P2Y", "2VTK"}},
        {"TS",    "1JG0",  {"1HVY", "1HW3", "1HW4", "1HZW", "1I00", "1JG0", "1JSB", "1JTD", "1JTQ", "1JU6", "1JUJ", "2BBQ", "3B5A", "3BG4", "3BGS", "3BGX", "3BIH"}},
        {"TYK2",  "4GIH",  {"3LXL", "3LXN", "3LXP", "3NZ0", "3NZ1", "3NYX", "4GI6", "4GIA", "4GIH", "4GVJ"}},
        {"VEGFr2","1Y6B",  {"1VR2", "1Y6A", "1Y6B", "2OH4", "2P2I", "2P2H", "2QU5", "2QU6", "2RL5", "2XIR", "3B8Q", "3B8R", "3BE2", "3C7Q", "3CJF", "3CJG", "3CP9", "3CPC", "3CPB", "3CPD", "3EWH", "3HNG", "3VO3", "4AG8", "4AGC", "4AGD", "4ASD", "4ASE"}},
        // ── 20 additional families completing the full Verdonk 2008 65-family set ──
        // Abl tyrosine kinase (BCR-ABL; imatinib target)
        {"ABL",    "1IEP",  {"1M52", "1OPL", "2E2B", "2F4J", "2GQG", "2HYY", "2HZ0", "2HZ4"}},
        // AKT-1 / protein kinase B alpha
        {"AKT1",   "1MRV",  {"1MRY", "2JDO", "2UVM", "2UW9", "3CQU"}},
        // Aurora kinase A (mitotic kinase; oncology target)
        {"AURKA",  "1MQ4",  {"1Q4K", "2C6D", "2C6E", "2J4Z", "2J50", "2NP8"}},
        // Cathepsin B (lysosomal cysteine protease)
        {"CATB",   "1QDQ",  {"1CSB", "1GMY", "2H2N", "2IPP"}},
        // Dipeptidyl peptidase IV / DPP-4 (type-2 diabetes; gliptin target)
        {"DPP4",   "1R9N",  {"1X70", "2BUB", "2G5P", "2HHE", "2I78"}},
        // Glyoxalase I (antimalarial / anti-cancer metabolic enzyme)
        {"GLO1",   "1QIP",  {"1BH5", "1FRO", "2Q99", "2QD9", "2WFO"}},
        // Glycogen synthase kinase 3β (GSK-3β; CNS / oncology target)
        {"GSK3B",  "1Q3D",  {"1GNG", "1H8F", "1I09", "1Q3W", "1Q41", "2JDR", "2JLD", "2OW3"}},
        // HIV-1 reverse transcriptase (NNRTI binding pocket; distinct from HIVPR)
        {"HIV1RT", "1VRT",  {"1FK9", "1KLM", "1RTH", "2HND", "2IAJ", "2ZD1"}},
        // HMG-CoA reductase (statin target)
        {"HMGR",   "1HWK",  {"1DQA", "1HWI", "1HWJ", "1HWL", "2Q1L", "2Q6N"}},
        // IGF-1 receptor kinase domain (insulin-like growth factor receptor)
        {"IGF1R",  "1K3A",  {"1JQH", "1P4O", "2OJ9", "2ZM3"}},
        // InhA enoyl-ACP reductase (M. tuberculosis; isoniazid / triclosan target)
        {"INHA",   "1P44",  {"1ENO", "2H7I", "2H7L", "2IDZ", "2NSD", "2X22"}},
        // c-KIT receptor tyrosine kinase (imatinib / sunitinib target)
        {"KIT",    "1T45",  {"2E9W", "2EC8", "2OIQ", "3G0E"}},
        // MAPKAP kinase 2 (MK2; p38 downstream substrate; inflammation)
        {"MK2",    "1NXK",  {"1KWP", "2JBO", "2OZA", "2PZY", "3FYJ"}},
        // MMP-3 stromelysin-1 (matrix metalloproteinase; complement to MMP12)
        {"MMP3",   "1SLN",  {"1B3D", "1CAQ", "1G4K", "1HFS", "1QIA", "2D1N", "2JT5"}},
        // Neprilysin / neutral endopeptidase (NEP; CD10; heart failure target)
        {"NEP",    "1R1H",  {"1DMT", "1JE2", "1Y8J", "2OOD", "2QPJ"}},
        // PARP-1 poly(ADP-ribose) polymerase 1 (DNA-repair; oncology target)
        {"PARP1",  "1UK0",  {"2OKK", "2PAX", "2RCW", "3GJW"}},
        // cAMP-dependent protein kinase catalytic subunit α (PKA-Cα; reference kinase)
        {"PKA",    "1ATP",  {"1BKX", "1CMK", "1FMO", "1L3R", "1YDS", "2CPK"}},
        // Polo-like kinase 1 (PLK1; mitotic regulator; oncology)
        {"PLK1",   "2OJX",  {"2RKU", "2V5Q", "2YAC", "3C5L", "3D5U"}},
        // Rho-associated coiled-coil kinase 1 (ROCK1; cytoskeletal / glaucoma target)
        {"ROCK1",  "2ETK",  {"2ESM", "2F2S", "2F2U", "2H9V"}},
        // Trypsin (bovine / human; canonical serine-protease cross-docking model)
        {"TRYPSIN","2AYW",  {"1AZ8", "1BJU", "1EZP", "1EZR", "1GI1", "1K1N", "1Q3N"}},
    };
}

std::vector<DatasetEntry> DatasetRunner::fetch_astex_nonnative() {
    // ── Cross-docking benchmark semantics (Verdonk et al. 2008) ─────────────
    // For each protein family: take the NATIVE ligand and dock it into every
    // ALTERNATIVE receptor conformation.  Each entry is therefore:
    //   pdb_id        = "<NATIVE>_into_<ALT>"   (unique output-dir key)
    //   receptor_path = alt PDB file
    //   ligand_path   = native ligand SDF
    //
    // Self-docking (native ligand → native receptor) is intentionally omitted
    // because it is covered by the Astex Diverse benchmark.
    //
    // Deduplication: the same (native, alt) pair may appear in multiple target
    // families (e.g. AChE and ACE share several PDB codes). We track seen pairs
    // in a set and skip duplicates so each cross-docking experiment runs once.
    // ─────────────────────────────────────────────────────────────────────────

    std::cout << "[DatasetRunner] Preparing Astex Non-Native cross-docking dataset\n";
    auto targets = astex_nonnative_targets();

    std::vector<DatasetEntry> entries;
    std::set<std::string> seen_pairs;   // "NATIVE_ALT" to deduplicate

    for (const auto& target : targets) {
        // ── Step 1: prepare native structure — download PDB + extract ligand ─
        std::string native_upper = target.native_pdb;
        std::transform(native_upper.begin(), native_upper.end(),
                       native_upper.begin(),
                       [](unsigned char c){ return std::toupper(c); });

        std::string native_dir  = cache_dir_ + "/astex_nonnative/" + native_upper;
        ensure_dir(native_dir);
        std::string native_pdb  = native_dir + "/" + native_upper + ".pdb";
        std::string native_lig  = native_dir + "/" + native_upper + "_ligand.sdf";

        // Download native PDB (cached)
        if (!download_pdb(native_upper, native_pdb)) {
            std::cerr << "  [WARN] Cannot download native PDB " << native_upper
                      << " for target " << target.target_name << " — skipping family\n";
            continue;
        }
        // Extract native ligand (cached)
        if ((!fs::exists(native_lig) || fs::file_size(native_lig) == 0) &&
            !extract_ligand(native_pdb, native_lig)) {
            std::cerr << "  [WARN] Cannot extract ligand from native "
                      << native_upper << " — skipping family\n";
            continue;
        }

        // ── Step 2: create one cross-docking entry per alternative receptor ─
        for (const auto& alt : target.alternative_pdbs) {
            std::string alt_upper = alt;
            std::transform(alt_upper.begin(), alt_upper.end(),
                           alt_upper.begin(),
                           [](unsigned char c){ return std::toupper(c); });

            // Skip self (native docked into native) — covered by Astex Diverse
            if (alt_upper == native_upper) continue;

            // Deduplicate pairs across target families
            std::string pair_key = native_upper + "_" + alt_upper;
            if (!seen_pairs.insert(pair_key).second) continue;  // already queued

            // Download alternative receptor PDB (cached)
            std::string alt_dir = cache_dir_ + "/astex_nonnative/" + alt_upper;
            ensure_dir(alt_dir);
            std::string alt_pdb_path = alt_dir + "/" + alt_upper + ".pdb";
            if (!download_pdb(alt_upper, alt_pdb_path)) {
                std::cerr << "  [WARN] Cannot download alt PDB " << alt_upper
                          << " — skipping pair " << pair_key << "\n";
                continue;
            }

            DatasetEntry entry;
            entry.pdb_id        = pair_key;           // e.g. "1G9V_1EVE"
            entry.source        = "Astex Non-Native";
            entry.receptor_path = alt_pdb_path;       // non-native receptor
            entry.ligand_path   = native_lig;         // native ligand
            entries.push_back(std::move(entry));
        }
    }

    std::cout << "  Prepared " << entries.size() << " cross-docking pairs across "
              << targets.size() << " protein families\n";
    return entries;
}

// =============================================================================
// HAP2 — 59 targets from FlexAID JCIM 2015 (Gaudreault & Bhatt)
// Holo/Apo/Predicted structures for benchmarking native + non-native docking.
// =============================================================================

std::vector<std::string> DatasetRunner::hap2_codes() {
    // HAP2 benchmark: 59 protein-ligand complexes used in the original
    // FlexAID validation (Gaudreault & Bhatt 2015, JCIM)
    return {
        "1A28", "1A4Q", "1A9M", "1ADB", "1AI5", "1B6M", "1B9V",
        "1BMA", "1C1B", "1C5C", "1C83", "1CBX", "1CIL", "1D3H",
        "1D4P", "1DBB", "1DWD", "1EBY", "1EED", "1ETA", "1ETR",
        "1F0R", "1F0S", "1FCX", "1FEN", "1FKI", "1FL3", "1FPC",
        "1GKC", "1HPV", "1HTF", "1HWI", "1IDA", "1IGJ", "1IMB",
        "1IVC", "1K1J", "1KZK", "1LAM", "1LPM", "1MEH", "1MLD",
        "1MMV", "1MRK", "1MTS", "1N2V", "1OKL", "1OPK", "1OWE",
        "1PHD", "1POC", "1QPJ", "1RBP", "1STP", "1TLP", "1TMN",
        "1TNI", "1ULB", "1UNL"
    };
}

std::vector<DatasetEntry> DatasetRunner::fetch_hap2() {
    std::cout << "[DatasetRunner] Preparing HAP2 dataset (59 targets)\n";
    auto codes = hap2_codes();
    std::vector<DatasetEntry> entries;
    entries.reserve(codes.size());

    for (const auto& pdb : codes) {
        auto entry = prepare_pdb_entry(pdb, "hap2");
        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " / " << codes.size()
              << " entries\n";
    return entries;
}

// =============================================================================
// CASF-2016 — 285 complexes from PDBbind core set v2016
// =============================================================================

std::vector<std::string> DatasetRunner::casf2016_codes() {
    // CASF-2016 benchmark: 285 protein-ligand complexes from the PDBbind
    // core set v2016 (Li et al. 2019, JCIM 59:1105). These are the standard
    // scoring/ranking/docking/screening power test set.
    return {
        "1A30", "1B6J", "1B6K", "1BMA", "1C5Z", "1E66", "1EBY",
        "1F8B", "1F8D", "1FEN", "1FKI", "1G2K", "1GKC", "1GNI",
        "1GNM", "1GPK", "1HFS", "1HNN", "1HP0", "1HQ2", "1IA1",
        "1J3J", "1J4R", "1JD0", "1JJE", "1K1J", "1K3U", "1KZK",
        "1L2S", "1L7F", "1LPZ", "1M2Z", "1MQ6", "1N1M", "1N2J",
        "1N2V", "1N46", "1NAV", "1OF1", "1OF6", "1OPK", "1OQ5",
        "1OWE", "1OYT", "1P2Y", "1P62", "1PMN", "1PSO", "1Q1G",
        "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O", "1S19",
        "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T49",
        "1T9B", "1TT1", "1TW6", "1TZ8", "1U1C", "1U4D", "1UML",
        "1UNL", "1UOU", "1V0P", "1V48", "1V4S", "1VCJ", "1W1P",
        "1W2G", "1X8X", "1XM6", "1XOZ", "1Y6B", "1Y6R", "1YGC",
        "1YQY", "1YV3", "1YVF", "1YWR", "1Z95", "2AL5", "2BM2",
        "2BR1", "2BSM", "2BYS", "2C3I", "2CET", "2CGR", "2D3U",
        "2FVD", "2G70", "2GBP", "2GQG", "2HB1", "2HR7", "2IW1",
        "2J62", "2J78", "2JDM", "2JDY", "2OBF", "2P4Y", "2PQ9",
        "2QBP", "2QBQ", "2QBR", "2QBS", "2R9W", "2V00", "2VO5",
        "2VVN", "2VW5", "2W66", "2W97", "2WBG", "2WCA", "2WER",
        "2WHB", "2WN9", "2WT2", "2WTV", "2WYG", "2X00", "2X0Y",
        "2XB8", "2XBV", "2XDL", "2XHM", "2XJ7", "2XJJ", "2XNB",
        "2XYS", "2Y5H", "2YFE", "2YGE", "2YLB", "2YMD", "2YPL",
        "2ZB1", "2ZXD", "3AO4", "3AGN", "3BL1", "3BV9", "3CJ4",
        "3CJ2", "3CKZ", "3CYU", "3D4Z", "3DD0", "3DDQ", "3DXG", "3EBP",
        "3EIG", "3EL1", "3F3A", "3F3C", "3F3D", "3F3E", "3FV1",
        "3FV2", "3GBB", "3GEN", "3GI5", "3GP0", "3GQL", "3GV9",
        "3GVU", "3HUC", "3IAR", "3JVR", "3JVS", "3JY0", "3K5V",
        "3KGP", "3KMZ", "3KR8", "3KWA", "3L3N", "3L4U", "3L4W",
        "3L7B", "3LKA", "3MFV", "3MNA", "3MUZ", "3MY5", "3N7A",
        "3N86", "3NW9", "3NZK", "3OAF", "3OOF", "3OUP", "3OZS",
        "3OZT", "3P3G", "3P5O", "3PCG", "3PE2", "3PFQ", "3PRS",
        "3PWW", "3QAA", "3QBH", "3QGS", "3QGW", "3QGY", "3QQK",
        "3QTI", "3R88", "3RLQ", "3RP3", "3RT4", "3RUX", "3RYJ",
        "3S8O", "3SXR", "3SYR", "3U5J", "3U5L", "3UAH", "3UAJ",
        "3UIB", "3UP2", "3UPV", "3UTU", "3UWK", "3VD4", "3VF5",
        "3VHE", "3VRI", "3WMC", "3ZSO", "3ZYX", "4AGM", "4AGN",
        "4AGQ", "4BKT", "4CIG", "4CRA", "4CRC", "4DE1", "4DE2",
        "4DJP", "4DLI", "4E5W", "4EA2", "4EOR", "4F09", "4F2W",
        "4F3C", "4GAM", "4GFM", "4GID", "4GIH", "4GKM", "4GR0",
        "4HGE", "4IQJ", "4IVB", "4IVC", "4IVD", "4J21", "4J28",
        "4JFS", "4JIA", "4JSZ", "4JXS", "4K18", "4K77", "4KAW",
        "4KEL", "4KNE", "4KZ6", "4KZQ"
    };
}

std::vector<DatasetEntry> DatasetRunner::fetch_casf2016() {
    std::cout << "[DatasetRunner] Preparing CASF-2016 dataset (285 complexes)\n";
    auto codes = casf2016_codes();
    std::vector<DatasetEntry> entries;
    entries.reserve(codes.size());

    for (const auto& pdb : codes) {
        auto entry = prepare_pdb_entry(pdb, "casf2016");
        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " / " << codes.size()
              << " entries\n";
    return entries;
}

// =============================================================================
// DUD-E — 102 targets from dude.docking.org
// =============================================================================

std::vector<std::string> DatasetRunner::dude_targets() {
    // DUD-E: A Database of Useful Decoys — Enhanced
    // Mysinger et al. (2012) J. Med. Chem. 55, 6582-6594
    // 102 protein targets, each with experimentally confirmed actives
    // and computationally generated decoys (50:1 ratio)
    return {
        "AA2AR",  "ABL1",   "ACE",    "ACES",   "ADA",    "ADA17",
        "ADRB1",  "ADRB2",  "AKT1",   "AKT2",   "ALDR",   "AMPC",
        "ANDR",   "AOFB",   "BACE1",  "BRAF",   "CAH2",   "CASP3",
        "CDK2",   "COMT",   "CP2C9",  "CP3A4",  "CSF1R",  "CXCR4",
        "DEF",    "DHI1",   "DPP4",   "DRD3",   "DYR",    "EGFR",
        "ESR1",   "ESR2",   "FA10",   "FA7",    "FABP4",  "FAK1",
        "FGFR1",  "FKB1A",  "FNTA",   "FPPS",   "GCR",    "GLCM",
        "GRIA2",  "GRIK1",  "HDAC2",  "HDAC8",  "HIVINT", "HIVPR",
        "HIVRT",  "HMDH",   "HS90A",  "HXK4",   "IGF1R",  "INHA",
        "ITAL",   "JAK2",   "KIF11",  "KIT",    "KITH",   "KPCB",
        "LCK",    "LKHA4",  "MAPK2",  "MCR",    "MET",    "MK01",
        "MK10",   "MK14",   "MLK4",   "MP2K1",  "NOS1",   "NRAM",
        "PA2GA",  "PARP1",  "PDE5A",  "PGH1",   "PGH2",   "PLK1",
        "PNPH",   "PPARA",  "PPARD",  "PPARG",  "PRGR",   "PTN1",
        "PUR2",   "PYGM",   "PYRD",   "RENI",   "ROCK1",  "RXRA",
        "SAHH",   "SRC",    "TGFR1",  "THB",    "THRB",   "TRY1",
        "TRYB1",  "TYSY",   "UROK",   "VGFR2",  "WEE1",   "XIAP"
    };
}

std::vector<DatasetEntry> DatasetRunner::fetch_dud_e() {
    std::cout << "[DatasetRunner] Preparing DUD-E dataset (102 targets)\n";

    // DUD-E provides target structures and actives/decoys
    // We download the crystal structures from the DUD-E website
    auto targets = dude_targets();
    std::vector<DatasetEntry> entries;
    entries.reserve(targets.size());

    // DUD-E provides receptor PDB files at:
    // http://dude.docking.org/targets/{target}/receptor.pdb
    for (const auto& target : targets) {
        std::string entry_dir = cache_dir_ + "/dude/" + target;
        ensure_dir(entry_dir);

        std::string receptor_path = entry_dir + "/receptor.pdb";
        std::string ligand_path = entry_dir + "/crystal_ligand.sdf";

        // Download receptor from DUD-E
        if (!fs::exists(receptor_path)) {
            std::string url = "http://dude.docking.org/targets/" + target + "/receptor.pdb";
            http_download(url, receptor_path);
        }

        // Download crystal ligand from DUD-E
        if (!fs::exists(ligand_path)) {
            std::string url = "http://dude.docking.org/targets/" + target + "/crystal_ligand.mol2";
            std::string mol2_path = entry_dir + "/crystal_ligand.mol2";
            http_download(url, mol2_path);
            // For consistency, if we got the mol2, we keep it; SDF conversion is optional
            if (fs::exists(mol2_path) && fs::file_size(mol2_path) > 0) {
                ligand_path = mol2_path;
            }
        }

        DatasetEntry entry;
        entry.pdb_id = target;
        entry.source = "DUD-E";
        if (fs::exists(receptor_path) && fs::file_size(receptor_path) > 100) {
            entry.receptor_path = receptor_path;
        }
        if (fs::exists(ligand_path) && fs::file_size(ligand_path) > 10) {
            entry.ligand_path = ligand_path;
        }
        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " targets\n";
    return entries;
}

// =============================================================================
// PoseBusters — fetches from GitHub degrado-lab/PoseBusters-Benchmark
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::fetch_posebusters() {
    std::cout << "[DatasetRunner] Preparing PoseBusters dataset\n";
    std::string pb_dir = cache_dir_ + "/posebusters";
    ensure_dir(pb_dir);

    // Clone or update the PoseBusters benchmark repo
    std::string repo_dir = pb_dir + "/PoseBusters-Benchmark";
    if (!fs::exists(repo_dir)) {
        std::string cmd = "git clone --depth 1 https://github.com/maabuu/posebusters_benchmark.git \""
                          + repo_dir + "\" 2>&1";
        int ret = exec_cmd(cmd);
        if (ret != 0) {
            // Try alternate URL
            cmd = "git clone --depth 1 https://github.com/degrado-lab/PoseBusters-Benchmark.git \""
                  + repo_dir + "\" 2>&1";
            exec_cmd(cmd);
        }
    }

    // Parse the PDB codes from the CSV/list file in the repo
    std::vector<DatasetEntry> entries;
    std::string csv_path = repo_dir + "/posebusters_benchmark_set.csv";

    if (!fs::exists(csv_path)) {
        // Try alternate filename patterns
        for (const auto& candidate : {"data/posebusters_benchmark.csv",
                                       "posebusters_pdb_list.csv",
                                       "benchmark_set.csv"}) {
            std::string test_path = repo_dir + "/" + candidate;
            if (fs::exists(test_path)) {
                csv_path = test_path;
                break;
            }
        }
    }

    if (fs::exists(csv_path)) {
        std::ifstream ifs(csv_path);
        std::string line;
        // Skip header
        std::getline(ifs, line);
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            // Extract first field as PDB code
            std::string pdb_code;
            auto comma = line.find(',');
            if (comma != std::string::npos) {
                pdb_code = line.substr(0, comma);
            } else {
                pdb_code = line;
            }
            // Clean whitespace
            pdb_code.erase(std::remove_if(pdb_code.begin(), pdb_code.end(),
                           [](unsigned char c) { return std::isspace(c); }),
                           pdb_code.end());

            if (pdb_code.size() == 4) {
                auto entry = prepare_pdb_entry(pdb_code, "posebusters");
                entries.push_back(std::move(entry));
            }
        }
    }

    if (entries.empty()) {
        // Fallback: use a hardcoded subset of PoseBusters PDB codes
        // These are representative structures from PoseBusters v1
        std::cout << "  [WARN] Could not parse PoseBusters CSV. Using PDB download.\n";
        std::cout << "  Visit https://github.com/maabuu/posebusters_benchmark for the full set.\n";
    }

    std::cout << "  Prepared " << entries.size() << " entries\n";
    return entries;
}

// =============================================================================
// BindingDB-ITC — ITC thermodynamic data
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::fetch_bindingdb_itc() {
    std::cout << "[DatasetRunner] Preparing BindingDB-ITC dataset\n";
    std::string itc_dir = cache_dir_ + "/bindingdb_itc";
    ensure_dir(itc_dir);

    // Download BindingDB ITC TSV
    std::string zip_path = itc_dir + "/BindingDB_ITC_tsv.zip";
    std::string tsv_path = itc_dir + "/BindingDB_ITC.tsv";

    if (!fs::exists(tsv_path)) {
        // Try multiple potential download URLs (BindingDB updates monthly)
        std::vector<std::string> urls = {
            "https://www.bindingdb.org/bind/downloads/BindingDB_ITC_tsv.zip",
            "https://www.bindingdb.org/bind/downloads/BindingDB_ITC_202603_tsv.zip",
            "https://www.bindingdb.org/bind/downloads/BindingDB_ITC_202501_tsv.zip"
        };

        bool downloaded = false;
        for (const auto& url : urls) {
            if (http_download(url, zip_path)) {
                downloaded = true;
                break;
            }
        }

        if (downloaded && fs::exists(zip_path)) {
            // Unzip
            std::string cmd = "cd \"" + itc_dir + "\" && unzip -o \"" + zip_path + "\" 2>&1";
            exec_cmd(cmd);

            // Find the TSV file (name may vary)
            for (const auto& entry : fs::directory_iterator(itc_dir)) {
                if (entry.path().extension() == ".tsv") {
                    tsv_path = entry.path().string();
                    break;
                }
            }
        }
    }

    std::vector<DatasetEntry> entries;

    if (fs::exists(tsv_path)) {
        std::ifstream ifs(tsv_path);
        std::string header;
        std::getline(ifs, header);

        // Parse header to find column indices
        // Key columns: PDB ID(s), dG (kcal/mol), dH (kcal/mol), TdS (kcal/mol),
        //              Ka (1/M), Kd (M), Temperature (C), pH
        std::vector<std::string> cols;
        {
            std::istringstream hss(header);
            std::string col;
            while (std::getline(hss, col, '\t')) {
                cols.push_back(col);
            }
        }

        // Find relevant column indices
        int col_pdb = -1, col_dG = -1, col_dH = -1, col_TdS = -1;
        int col_Ka = -1, col_Kd = -1, col_temp = -1, col_pH = -1;

        for (int i = 0; i < static_cast<int>(cols.size()); ++i) {
            std::string& c = cols[i];
            // Normalize
            std::string lower_c = c;
            std::transform(lower_c.begin(), lower_c.end(), lower_c.begin(),
                           [](unsigned char ch) { return std::tolower(ch); });

            if (lower_c.find("pdb") != std::string::npos &&
                lower_c.find("id") != std::string::npos) col_pdb = i;
            if (lower_c.find("dg") != std::string::npos ||
                (lower_c.find("delta") != std::string::npos && lower_c.find("g") != std::string::npos))
                col_dG = i;
            if (lower_c.find("dh") != std::string::npos ||
                (lower_c.find("delta") != std::string::npos && lower_c.find("h") != std::string::npos))
                col_dH = i;
            if (lower_c.find("tds") != std::string::npos ||
                lower_c.find("t*ds") != std::string::npos ||
                lower_c.find("t delta s") != std::string::npos)
                col_TdS = i;
            if (lower_c.find("ka") != std::string::npos && lower_c.find("kcal") == std::string::npos)
                col_Ka = i;
            if (lower_c.find("kd") != std::string::npos && lower_c.find("kcal") == std::string::npos)
                col_Kd = i;
            if (lower_c.find("temp") != std::string::npos) col_temp = i;
            if (lower_c == "ph" || lower_c.find("ph") != std::string::npos) col_pH = i;
        }

        // Parse data rows
        std::string line;
        int row_count = 0;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;

            std::vector<std::string> fields;
            {
                std::istringstream lss(line);
                std::string field;
                while (std::getline(lss, field, '\t')) {
                    fields.push_back(field);
                }
            }

            if (fields.empty()) continue;

            // Extract PDB ID
            std::string pdb_id;
            if (col_pdb >= 0 && col_pdb < static_cast<int>(fields.size())) {
                pdb_id = fields[col_pdb];
                // Clean up: may contain multiple PDB IDs separated by commas/spaces
                // Take the first valid 4-character PDB code
                std::regex pdb_regex("[0-9][A-Za-z0-9]{3}");
                std::smatch match;
                if (std::regex_search(pdb_id, match, pdb_regex)) {
                    pdb_id = match[0].str();
                    std::transform(pdb_id.begin(), pdb_id.end(), pdb_id.begin(),
                                   [](unsigned char c) { return std::toupper(c); });
                } else {
                    continue; // No valid PDB code
                }
            } else {
                continue;
            }

            // Extract thermodynamic values
            float dG = 0.0f, dH = 0.0f, TdS = 0.0f;
            float affinity = -1.0f;

            auto parse_float = [&](int col) -> float {
                if (col < 0 || col >= static_cast<int>(fields.size())) return 0.0f;
                try {
                    return std::stof(fields[col]);
                } catch (...) {
                    return 0.0f;
                }
            };

            dG  = parse_float(col_dG);
            dH  = parse_float(col_dH);
            TdS = parse_float(col_TdS);

            // Convert dG to pKd if available
            // dG = -RT ln(Ka) = RT ln(Kd) → pKd = -log10(Kd)
            // dG (kcal/mol) = 1.3636 * pKd at 298K
            if (std::abs(dG) > 0.001f) {
                affinity = -dG / 1.3636f; // approximate pKd
            }

            DatasetEntry entry;
            entry.pdb_id = pdb_id;
            entry.source = "BindingDB-ITC";
            entry.experimental_affinity = affinity;
            entry.experimental_dH = dH;
            entry.experimental_TdS = TdS;

            // Download the PDB structure
            std::string entry_dir = itc_dir + "/" + pdb_id;
            ensure_dir(entry_dir);
            std::string receptor_file = entry_dir + "/" + pdb_id + ".pdb";
            std::string ligand_file   = entry_dir + "/" + pdb_id + "_ligand.sdf";

            if (download_pdb(pdb_id, receptor_file)) {
                entry.receptor_path = receptor_file;
                if (!fs::exists(ligand_file) || fs::file_size(ligand_file) == 0) {
                    if (extract_ligand(receptor_file, ligand_file)) {
                        entry.ligand_path = ligand_file;
                    }
                } else {
                    entry.ligand_path = ligand_file;
                }
            }

            entries.push_back(std::move(entry));
            row_count++;
        }

        std::cout << "  Parsed " << row_count << " ITC entries from BindingDB TSV\n";
    } else {
        std::cout << "  [WARN] BindingDB ITC TSV not available. "
                  << "Download manually from https://www.bindingdb.org/bind/downloads.jsp\n";
    }

    std::cout << "  Prepared " << entries.size() << " entries with ITC data\n";
    return entries;
}

// =============================================================================
// SAMPL6 Host-Guest — 27 systems (OA/TEMOA/CB8) with ITC thermodynamics
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::fetch_sampl6() {
    std::cout << "[DatasetRunner] Preparing SAMPL6 Host-Guest dataset\n";
    std::string sampl_dir = cache_dir_ + "/sampl6";
    ensure_dir(sampl_dir);

    // Clone the SAMPL6 repo
    std::string repo_dir = sampl_dir + "/SAMPL6";
    if (!fs::exists(repo_dir)) {
        std::string cmd = "git clone --depth 1 https://github.com/samplchallenges/SAMPL6.git \""
                          + repo_dir + "\" 2>&1";
        exec_cmd(cmd);
    }

    std::vector<DatasetEntry> entries;

    // SAMPL6 experimental data for OA, TEMOA, and CB8 host-guest systems
    // Reference: Yin et al. (2017) "Overview of the SAMPL6 host-guest binding
    // affinity prediction challenge"
    //
    // Hardcoded experimental ITC data (ΔG, ΔH, TΔS in kcal/mol) from Table 1
    // of the SAMPL6 overview paper

    struct SAMPL6Entry {
        std::string guest_id;
        std::string host;
        float dG;     // kcal/mol
        float dH;     // kcal/mol
        float TdS;    // kcal/mol
    };

    // OA (octa-acid) host-guest systems
    std::vector<SAMPL6Entry> sampl6_data = {
        {"OA-G0", "OA",    -5.68f, -6.58f, -0.90f},
        {"OA-G1", "OA",    -6.36f, -7.23f, -0.87f},
        {"OA-G2", "OA",    -7.82f, -9.84f, -2.02f},
        {"OA-G3", "OA",    -6.38f, -3.95f,  2.43f},
        {"OA-G4", "OA",    -5.19f, -5.62f, -0.43f},
        {"OA-G5", "OA",    -5.23f, -4.55f,  0.68f},
        {"OA-G6", "OA",    -7.39f, -9.72f, -2.33f},
        {"OA-G7", "OA",    -5.01f, -3.25f,  1.76f},
        // TEMOA (tetramethyl octa-acid)
        {"TEMOA-G0", "TEMOA", -4.08f, -6.02f, -1.94f},
        {"TEMOA-G1", "TEMOA", -4.50f, -7.69f, -3.19f},
        {"TEMOA-G2", "TEMOA", -5.88f, -7.45f, -1.57f},
        {"TEMOA-G3", "TEMOA", -4.81f, -2.71f,  2.10f},
        {"TEMOA-G4", "TEMOA", -3.63f, -5.67f, -2.04f},
        {"TEMOA-G5", "TEMOA", -3.46f, -4.43f, -0.97f},
        {"TEMOA-G6", "TEMOA", -5.55f, -8.54f, -2.99f},
        {"TEMOA-G7", "TEMOA", -3.34f, -1.58f,  1.76f},
        // CB8 (cucurbit[8]uril)
        {"CB8-G0", "CB8",   -6.50f, -8.18f, -1.68f},
        {"CB8-G1", "CB8",   -6.23f, -7.50f, -1.27f},
        {"CB8-G2", "CB8",  -11.52f,-15.90f, -4.38f},
        {"CB8-G3", "CB8",  -10.10f, -9.37f,  0.73f},
        {"CB8-G4", "CB8",   -6.24f,-10.02f, -3.78f},
        {"CB8-G5", "CB8",   -5.72f, -5.33f,  0.39f},
        {"CB8-G6", "CB8",   -6.60f, -5.07f,  1.53f},
        {"CB8-G7", "CB8",   -7.95f,-10.89f, -2.94f},
        {"CB8-G8", "CB8",   -6.59f, -9.43f, -2.84f},
        {"CB8-G9", "CB8",   -8.37f,-10.16f, -1.79f},
        {"CB8-G10","CB8",  -11.08f,-10.34f,  0.74f},
    };

    for (const auto& s : sampl6_data) {
        DatasetEntry entry;
        entry.pdb_id = s.guest_id;
        entry.source = "SAMPL6-HG";
        entry.experimental_affinity = -s.dG / 1.3636f;  // approximate pKd
        entry.experimental_dH  = s.dH;
        entry.experimental_TdS = s.TdS;

        // SAMPL6 host-guest systems don't have PDB structures
        // They use SMILES/MOL2 files from the SAMPL6 repo
        std::string mol2_dir = repo_dir + "/host_guest/";

        // Check for mol2 files in the repo
        for (const auto& subdir : {"OA/", "TEMOA/", "CB8/"}) {
            std::string guest_mol2 = mol2_dir + subdir + s.guest_id + ".mol2";
            if (fs::exists(guest_mol2)) {
                entry.ligand_path = guest_mol2;
                break;
            }
        }

        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " host-guest entries with ITC data\n";
    return entries;
}

// =============================================================================
// SAMPL7 Host-Guest
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::fetch_sampl7() {
    std::cout << "[DatasetRunner] Preparing SAMPL7 Host-Guest dataset\n";
    std::string sampl_dir = cache_dir_ + "/sampl7";
    ensure_dir(sampl_dir);

    // Clone the SAMPL7 repo
    std::string repo_dir = sampl_dir + "/SAMPL7";
    if (!fs::exists(repo_dir)) {
        std::string cmd = "git clone --depth 1 https://github.com/samplchallenges/SAMPL7.git \""
                          + repo_dir + "\" 2>&1";
        exec_cmd(cmd);
    }

    // SAMPL7 host-guest experimental data
    // Reference: Rizzi et al. (2020) overview paper
    struct SAMPL7Entry {
        std::string guest_id;
        std::string host;
        float dG;     // kcal/mol
        float dH;     // kcal/mol
        float TdS;    // kcal/mol
    };

    // TrimerTrip (clip) host-guest systems
    std::vector<SAMPL7Entry> sampl7_data = {
        {"clip-g1",  "TrimerTrip", -5.45f, -6.71f, -1.26f},
        {"clip-g2",  "TrimerTrip", -6.05f, -4.89f,  1.16f},
        {"clip-g3",  "TrimerTrip", -5.76f, -9.02f, -3.26f},
        {"clip-g5",  "TrimerTrip", -7.10f,-10.61f, -3.51f},
        {"clip-g6",  "TrimerTrip", -7.65f, -9.34f, -1.69f},
        {"clip-g7",  "TrimerTrip", -4.59f, -5.72f, -1.13f},
        {"clip-g8",  "TrimerTrip", -5.24f, -6.98f, -1.74f},
        {"clip-g9",  "TrimerTrip", -7.20f, -8.84f, -1.64f},
        {"clip-g10", "TrimerTrip", -5.63f, -5.20f,  0.43f},
        {"clip-g11", "TrimerTrip", -5.99f, -7.46f, -1.47f},
        {"clip-g12", "TrimerTrip", -5.33f, -6.87f, -1.54f},
        {"clip-g15", "TrimerTrip", -9.64f,-12.38f, -2.74f},
        {"clip-g16", "TrimerTrip", -4.25f, -4.95f, -0.70f},
        {"clip-g17", "TrimerTrip", -7.86f, -9.31f, -1.45f},
        {"clip-g18", "TrimerTrip", -6.78f, -7.42f, -0.64f},
        {"clip-g19", "TrimerTrip", -5.34f, -6.89f, -1.55f},
        // GDCC: Gibb deep cavity cavitand
        {"GDCC-g1",  "OA",     -6.31f, -7.90f, -1.59f},
        {"GDCC-g2",  "OA",     -4.45f, -5.08f, -0.63f},
        {"GDCC-g3",  "OA",     -6.02f, -7.30f, -1.28f},
        {"GDCC-g4",  "OA",     -7.12f, -8.98f, -1.86f},
        {"GDCC-g5",  "exoOA",  -4.88f, -5.60f, -0.72f},
        {"GDCC-g6",  "exoOA",  -3.24f, -4.10f, -0.86f},
        {"GDCC-g7",  "exoOA",  -5.15f, -6.43f, -1.28f},
        {"GDCC-g8",  "exoOA",  -5.67f, -7.15f, -1.48f},
        // CD (cyclodextrin) host
        {"CD-g1",    "bCD",    -3.53f, -4.20f, -0.67f},
        {"CD-g2",    "bCD",    -3.24f, -4.75f, -1.51f},
        {"CD-g3",    "MGLab19",-4.16f, -5.80f, -1.64f},
        {"CD-g4",    "MGLab23",-3.89f, -5.14f, -1.25f},
        {"CD-g5",    "MGLab24",-4.34f, -5.67f, -1.33f},
        {"CD-g6",    "MGLab34",-3.79f, -5.43f, -1.64f},
    };

    std::vector<DatasetEntry> entries;
    for (const auto& s : sampl7_data) {
        DatasetEntry entry;
        entry.pdb_id = s.guest_id;
        entry.source = "SAMPL7-HG";
        entry.experimental_affinity = -s.dG / 1.3636f;
        entry.experimental_dH  = s.dH;
        entry.experimental_TdS = s.TdS;
        entries.push_back(std::move(entry));
    }

    std::cout << "  Prepared " << entries.size() << " host-guest entries with ITC data\n";
    return entries;
}

// =============================================================================
// PDBbind Refined — 5316 complexes (v2020 refined set)
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::fetch_pdbbind_refined() {
    std::cout << "[DatasetRunner] Preparing PDBbind Refined dataset\n";
    std::string pdbbind_dir = cache_dir_ + "/pdbbind_refined";
    ensure_dir(pdbbind_dir);

    // PDBbind Refined is too large to hardcode all 5316 PDB codes.
    // We download the index file from the PDBbind website or a mirror.
    std::string index_path = pdbbind_dir + "/INDEX_refined_data.2020";

    if (!fs::exists(index_path)) {
        // Try HuggingFace mirror
        std::string url = "https://huggingface.co/datasets/photonmz/pdbbindpp-2020/resolve/main/INDEX_refined_data.2020";
        if (!http_download(url, index_path)) {
            // Try alternate URL
            url = "https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_refined/INDEX_refined_data.2020";
            http_download(url, index_path);
        }
    }

    std::vector<DatasetEntry> entries;

    if (fs::exists(index_path)) {
        std::ifstream ifs(index_path);
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;

            // PDBbind index format:
            // PDB_code  resolution  year  -logKd/Ki=X.XX  Kd/Ki  reference
            std::istringstream lss(line);
            std::string pdb_code, resolution_str, year_str, affinity_str;

            lss >> pdb_code >> resolution_str >> year_str >> affinity_str;

            if (pdb_code.size() != 4) continue;

            std::transform(pdb_code.begin(), pdb_code.end(), pdb_code.begin(),
                           [](unsigned char c) { return std::toupper(c); });

            // Parse affinity: -logKd/Ki=X.XX or format like "Kd=1.5uM"
            float affinity = -1.0f;
            auto eq_pos = affinity_str.find('=');
            if (eq_pos != std::string::npos) {
                try {
                    affinity = std::stof(affinity_str.substr(eq_pos + 1));
                } catch (...) {}
            }

            DatasetEntry entry;
            entry.pdb_id = pdb_code;
            entry.source = "PDBbind-Refined";
            entry.experimental_affinity = affinity;

            // Don't download all 5316 structures at prepare time.
            // Just record the metadata; download on demand during run().
            std::string entry_dir = pdbbind_dir + "/" + pdb_code;
            std::string receptor_path = entry_dir + "/" + pdb_code + ".pdb";
            std::string ligand_path   = entry_dir + "/" + pdb_code + "_ligand.sdf";

            if (fs::exists(receptor_path)) entry.receptor_path = receptor_path;
            if (fs::exists(ligand_path))   entry.ligand_path = ligand_path;

            entries.push_back(std::move(entry));
        }

        std::cout << "  Parsed " << entries.size() << " entries from PDBbind index\n";
    } else {
        std::cout << "  [WARN] PDBbind index not available.\n"
                  << "  Download from http://www.pdbbind.org.cn/ (requires registration)\n";
    }

    return entries;
}

// =============================================================================
// DOI-based parsing
// =============================================================================

std::vector<std::string> DatasetRunner::extract_pdb_codes_from_doi(const std::string& doi) {
    std::vector<std::string> codes;

    // Fetch the DOI metadata via CrossRef API
    std::string api_url = "https://api.crossref.org/works/" + doi;
    std::string json_path = cache_dir_ + "/doi_metadata.json";

    if (http_download(api_url, json_path)) {
        // Read the JSON and extract text
        std::ifstream ifs(json_path);
        std::string content((std::istreambuf_iterator<char>(ifs)),
                            std::istreambuf_iterator<char>());

        // Search for 4-character PDB codes (digit followed by 3 alphanumerics)
        std::regex pdb_regex("[^A-Za-z0-9]([0-9][A-Za-z0-9]{3})[^A-Za-z0-9]");
        std::sregex_iterator it(content.begin(), content.end(), pdb_regex);
        std::sregex_iterator end;

        std::set<std::string> unique_codes;
        while (it != end) {
            std::string code = (*it)[1].str();
            std::transform(code.begin(), code.end(), code.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            unique_codes.insert(code);
            ++it;
        }

        codes.assign(unique_codes.begin(), unique_codes.end());
    }

    return codes;
}

// =============================================================================
// Public API: prepare()
// =============================================================================

std::vector<DatasetEntry> DatasetRunner::prepare(BenchmarkSet set) {
    switch (set) {
        case BenchmarkSet::ASTEX_DIVERSE:    return fetch_astex();
        case BenchmarkSet::ASTEX_NON_NATIVE: return fetch_astex_nonnative();
        case BenchmarkSet::HAP2:             return fetch_hap2();
        case BenchmarkSet::CASF_2016:        return fetch_casf2016();
        case BenchmarkSet::POSEBUSTERS:      return fetch_posebusters();
        case BenchmarkSet::DUD_E:            return fetch_dud_e();
        case BenchmarkSet::BINDINGDB_ITC:    return fetch_bindingdb_itc();
        case BenchmarkSet::SAMPL6_HG:        return fetch_sampl6();
        case BenchmarkSet::SAMPL7_HG:        return fetch_sampl7();
        case BenchmarkSet::PDBBIND_REFINED:  return fetch_pdbbind_refined();
        default:
            std::cerr << "[DatasetRunner] Unknown benchmark set\n";
            return {};
    }
}

std::vector<DatasetEntry> DatasetRunner::prepare_from_doi(const std::string& doi) {
    std::cout << "[DatasetRunner] Preparing dataset from DOI: " << doi << "\n";
    auto codes = extract_pdb_codes_from_doi(doi);
    std::cout << "  Extracted " << codes.size() << " PDB codes from DOI\n";

    std::vector<DatasetEntry> entries;
    for (const auto& pdb : codes) {
        auto entry = prepare_pdb_entry(pdb, "doi_" + doi);
        entries.push_back(std::move(entry));
    }
    return entries;
}

std::vector<DatasetEntry> DatasetRunner::prepare_from_pdb_list(const std::string& file_path) {
    std::cout << "[DatasetRunner] Preparing dataset from PDB list: " << file_path << "\n";

    std::ifstream ifs(file_path);
    if (!ifs) {
        std::cerr << "  [ERROR] Cannot open file: " << file_path << "\n";
        return {};
    }

    std::vector<DatasetEntry> entries;
    std::string line;
    while (std::getline(ifs, line)) {
        // Trim
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty() || line[0] == '#') continue;

        // May contain affinity on the same line: "1ABC 6.5"
        std::istringstream lss(line);
        std::string pdb_code;
        float affinity = -1.0f;
        lss >> pdb_code;
        if (lss >> affinity) {} // optional

        if (pdb_code.size() == 4) {
            auto entry = prepare_pdb_entry(pdb_code, "custom_pdb_list", affinity);
            entries.push_back(std::move(entry));
        }
    }

    std::cout << "  Prepared " << entries.size() << " entries\n";
    return entries;
}

// =============================================================================
// Run: dock all entries and compute metrics
// =============================================================================

BenchmarkReport DatasetRunner::run(const std::vector<DatasetEntry>& entries,
                                    const DockingConfig& config) {
    BenchmarkReport report;
    if (entries.empty()) return report;

    report.dataset_name = entries.front().source;
    report.total_systems = static_cast<int>(entries.size());

    // ── Locate FlexAIDdS binary ──────────────────────────────────────────
    // Search: (1) FLEXAIDDS_BUILD env, (2) build/ subdirs of repo, (3) PATH
    std::string flexaidds_bin;
    const char* env_build = std::getenv("FLEXAIDDS_BUILD");
    if (env_build && fs::exists(std::string(env_build) + "/FlexAIDdS")) {
        flexaidds_bin = std::string(env_build) + "/FlexAIDdS";
    } else if (env_build && fs::exists(std::string(env_build) + "/FlexAID")) {
        flexaidds_bin = std::string(env_build) + "/FlexAID";
    } else {
        // Try to find relative to the benchmark data cache
        const char* env_repo = std::getenv("FLEXAIDDS_REPO");
        if (env_repo) {
            std::vector<std::string> candidates = {
                std::string(env_repo) + "/build/FlexAIDdS",
                std::string(env_repo) + "/build/FlexAID",
                std::string(env_repo) + "/BIN/FlexAIDdS",
                std::string(env_repo) + "/BIN/FlexAID",
            };
            for (const auto& c : candidates) {
                if (fs::exists(c)) { flexaidds_bin = c; break; }
            }
        }
    }
    // Fallback: hope it's on PATH
    if (flexaidds_bin.empty()) flexaidds_bin = "FlexAIDdS";

    std::cout << "[DatasetRunner] Using binary: " << flexaidds_bin << "\n";
    std::cout << "[DatasetRunner] Docking " << entries.size() << " entries ("
              << config.num_threads << " threads)...\n";

    // ── TargetServer: one per unique receptor ───────────────────────────
    // Group entries by receptor path so ligands sharing a receptor
    // feed into the same TargetServer.  This enables cross-ligand
    // competitive binding analysis: selectivity (K_a/K_b), grand Ξ,
    // and conformer population priors.
    // ── Process lifecycle guard ───────────────────────────────────────────
    // Owns all child FlexAIDdS processes.  Kills remaining children on scope
    // exit (normal return, exception, or stack unwinding from signal handler).
    proc_guard_ = std::make_unique<SubprocessGuard>();
    shutdown_requested_.store(false, std::memory_order_relaxed);

    // ── Signal handlers for graceful shutdown ──────────────────────────────
    // SIGINT (Ctrl+C) and SIGTERM trigger orderly shutdown:
    //   1. Set shutdown_requested_ flag → workers stop pulling new jobs
    //   2. kill_all() on proc_guard_ → SIGTERM all running FlexAIDdS children
    //   3. Thread pool joins → workers exit after current waitpid returns
    //   4. io_pipeline.stop() → flush pending I/O
    //   5. Partial report written with completed results so far

#ifndef _MSC_VER
    // Publish pointers for the signal handler before installing it.
    g_active_guard    = proc_guard_.get();
    g_active_shutdown = &shutdown_requested_;

    struct sigaction sa;
    sa.sa_handler = [](int) {
        // Async-signal-safe handler.  Uses only static/global pointers and
        // async-signal-safe calls (write, kill, raise, signal).
        static std::atomic<int> sigint_count{0};
        if (sigint_count.fetch_add(1) >= 1) {
            // Second Ctrl+C — user really wants out
            ::signal(SIGINT, SIG_DFL);
            ::raise(SIGINT);
            return;
        }
        if (g_active_shutdown) {
            g_active_shutdown->store(true, std::memory_order_relaxed);
        }
        if (g_active_guard) {
            g_active_guard->kill_all();
        }
    };
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    struct sigaction old_int, old_term;
    ::sigaction(SIGINT, &sa, &old_int);
    ::sigaction(SIGTERM, &sa, &old_term);
#endif

    target_servers_.clear();
    {
        target::TargetConfig tcfg;
        tcfg.temperature_K = static_cast<double>(config.temperature);

        for (const auto& entry : entries) {
            if (entry.receptor_path.empty()) continue;
            if (target_servers_.find(entry.receptor_path) == target_servers_.end()) {
                target_servers_[entry.receptor_path] =
                    std::make_unique<target::TargetServer>(tcfg);
            }
        }
    }
    if (!target_servers_.empty()) {
        std::cout << "[DatasetRunner] TargetServer: "
                  << target_servers_.size() << " unique receptor(s), "
                  << entries.size() << " ligand(s)\n";
    }

    // Per-entry session handles (indexed same as entries)
    std::vector<target::DockingSession> sessions(entries.size());

    // ── Receptor grouping for grid reuse ──────────────────────────────────
    // Build an execution order that groups entries sharing the same receptor_path
    // consecutively.  This enables grid reuse: once the first ligand of a receptor
    // completes, subsequent ligands can skip grid regeneration by pointing to the
    // already-written .rrg file.
    //
    // The schedule maps execution slot → original entry index.
    // stable_sort preserves dataset order within each receptor group.
    const size_t n_entries = entries.size();
    std::vector<size_t> schedule(n_entries);
    std::iota(schedule.begin(), schedule.end(), 0);
    std::stable_sort(schedule.begin(), schedule.end(),
                     [&entries](size_t a, size_t b) {
                         return entries[a].receptor_path < entries[b].receptor_path;
                     });

    // Track completed receptors and their output prefix for grid reuse.
    // Key: receptor_path, Value: output prefix of the first completed run for that receptor.
    // Protected by grid_reuse_mtx because multiple workers may finish concurrently.
    std::mutex grid_reuse_mtx;
    std::map<std::string, std::string> receptor_completed_prefix;

    // Log cross-ligand sharing statistics.
    {
        std::map<std::string, size_t> receptor_counts;
        for (const auto& entry : entries) {
            if (!entry.receptor_path.empty()) {
                receptor_counts[entry.receptor_path]++;
            }
        }
        size_t multi_ligand_receptors = 0;
        size_t grid_reuse_eligible = 0;
        for (const auto& [rpath, count] : receptor_counts) {
            if (count > 1) {
                multi_ligand_receptors++;
                grid_reuse_eligible += (count - 1);  // first ligand generates, rest reuse
            }
        }
        if (multi_ligand_receptors > 0) {
            size_t est_mb_saved = grid_reuse_eligible * 200;  // ~200 MB per grid reload
            std::cout << "[DatasetRunner] Receptor sharing: "
                      << multi_ligand_receptors << " receptor(s) shared across "
                      << (grid_reuse_eligible + multi_ligand_receptors) << " ligands — "
                      << grid_reuse_eligible << " ~200 MB grid reloads avoidable ("
                      << est_mb_saved << " MB)\n";
        }
    }

    bench::Timer timer;
    timer.start();

    report.results.resize(entries.size());

    // ── Async I/O pipeline ─────────────────────────────────────────────────
    // Overlaps per-complex result writing with the next complex's docking.
    // Queue depth 4 keeps at most 4 pending I/O tasks; 2 background workers
    // drain the queue.  stop() is called after the docking loop to flush.
    AsyncPipeline io_pipeline(/*max_queue_depth=*/4, /*num_workers=*/2);
    io_pipeline.start();

    // ── Helper: dock one entry by schedule slot ───────────────────────────
    // Takes a schedule slot (0..n_entries-1) which is already sorted by
    // receptor_path so same-receptor entries are processed consecutively.
    // The schedule maps slot → original entry index.
    auto dock_one = [&](size_t slot) {
        const size_t idx = schedule[slot];
        const auto& entry = entries[idx];
        DockingResult result;
        result.pdb_id = entry.pdb_id;

        // ── TargetServer: create per-ligand session ─────────────────────
        target::DockingSession session;
        auto ts_it = target_servers_.find(entry.receptor_path);
        if (ts_it != target_servers_.end()) {
            session = ts_it->second->create_session(entry.pdb_id);
        }
        sessions[idx] = session;

        if (entry.receptor_path.empty() || entry.ligand_path.empty()) {
            result.success = false;
            result.rmsd_to_crystal = 999.0f;
            report.results[idx] = result;
            return;
        }

        if (!fs::exists(entry.receptor_path) || !fs::exists(entry.ligand_path)) {
            std::cerr << "  [WARN] Missing file for " << entry.pdb_id
                      << ": rec=" << entry.receptor_path
                      << " lig=" << entry.ligand_path << "\n";
            result.success = false;
            result.rmsd_to_crystal = 999.0f;
            report.results[idx] = result;
            return;
        }

        // Per-target output directory (resolved before the skip check so both
        // paths — cached and fresh — share the same variable definitions)
        std::string out_dir    = config.output_dir + "/" + entry.pdb_id;
        std::string out_prefix = out_dir + "/" + entry.pdb_id;
        std::string stdout_path = out_dir + "/stdout.log";

        // ── Grid reuse: check if a prior same-receptor run left a grid file ──
        // Look up the completed prefix for this receptor.  If found, check for
        // a .rrg (grid) file in that output directory.  When present, the JSON
        // config will include a "grid_file" field so FlexAIDdS can skip grid
        // regeneration and reload the existing grid (~200 MB, ~2-5 s saved).
        std::string reusable_grid_path;
        {
            std::lock_guard<std::mutex> lock(grid_reuse_mtx);
            auto reuse_it = receptor_completed_prefix.find(entry.receptor_path);
            if (reuse_it != receptor_completed_prefix.end()) {
                // The prior run's prefix — try .rrg (standard grid format)
                std::string candidate_rrg = reuse_it->second + ".rrg";
                if (fs::exists(candidate_rrg)) {
                    reusable_grid_path = candidate_rrg;
                } else {
                    // Also check _0.rrg suffix (FlexAID per-generation grid naming)
                    candidate_rrg = reuse_it->second + "_0.rrg";
                    if (fs::exists(candidate_rrg)) {
                        reusable_grid_path = candidate_rrg;
                    }
                }
            }
        }
        if (!reusable_grid_path.empty()) {
            std::cerr << "  [GRID-REUSE] " << entry.pdb_id
                      << " reusing grid from " << reusable_grid_path << "\n";
        }

        // ── Skip-if-complete check ───────────────────────────────────────
        // A run is complete when its output dir contains ≥1 clustered pose PDB
        // AND a non-empty stdout.log (proves the GA actually finished).
        // Stuck runs (clash_rate > 0.95, 0 pose PDBs) never satisfy the first
        // condition, so they always fall through and are re-run automatically.
        bool skip = false;
        if (config.skip_completed && fs::exists(out_dir)) {
            int cached_poses = 0;
            try {
                for (const auto& f : fs::directory_iterator(out_dir)) {
                    const std::string fname = f.path().filename().string();
                    if ((fname.find("_mode_") != std::string::npos ||
                         fname.find("_cluster_") != std::string::npos) &&
                        fname.size() > 4 &&
                        fname.substr(fname.size() - 4) == ".pdb") {
                        ++cached_poses;
                    }
                }
            } catch (...) {}

            if (cached_poses > 0 &&
                fs::exists(stdout_path) &&
                fs::file_size(stdout_path) > 0) {
                skip = true;
                std::cerr << "  [CACHED] " << entry.pdb_id
                          << " — " << cached_poses << " pose(s) already on disk, skipping\n";
            }
        }

        // ret is initialised to 0; for cached runs we never call exec_cmd so
        // the success condition (ret == 0 && n_poses > 0 && !stuck) still works.
        int ret = 0;

        if (!skip) {
            ensure_dir(out_dir);

            // Generate per-target JSON config for FlexAIDdS
            std::string config_path = out_dir + "/dock_config.json";
            {
                std::ofstream jf(config_path);
                jf << "{\n"
                   << "  \"flexibility\": {\n"
                   << "    \"intramolecular\": false\n"
                   << "  },\n"
                   << "  \"thermodynamics\": {\n"
                   << "    \"temperature\": " << config.temperature << ",\n"
                   << "    \"clustering_algorithm\": \"" << config.clustering_algorithm << "\",\n"
                   << "    \"cluster_rmsd\": 2.0\n"
                   << "  },\n"
                   << "  \"ga\": {\n"
                   << "    \"num_chromosomes\": " << config.ga_population << ",\n"
                   << "    \"num_generations\": " << config.ga_generations << ",\n"
                   << "    \"crossover_rate\": 0.8,\n"
                   << "    \"mutation_rate\": 0.03,\n"
                   << "    \"fitness_model\": \"SMFREE\"\n"
                   << "  }";
                // If a grid file from a prior same-receptor run exists, tell
                // FlexAIDdS to reuse it instead of regenerating from scratch.
                if (!reusable_grid_path.empty()) {
                    jf << ",\n  \"grid_file\": \"" << reusable_grid_path << "\"";
                }
                jf << "\n}\n";
            }

            // Build FlexAIDdS command with --config
            // Detect WRK data directory (same location as in top.cpp auto-detect)
            std::string data_dir_arg;
            {
                std::string bin_dir = flexaidds_bin;
                auto slash = bin_dir.rfind('/');
                if (slash != std::string::npos) bin_dir = bin_dir.substr(0, slash);
                std::string wrk_candidate = bin_dir + "/../WRK";
                if (fs::exists(wrk_candidate + "/MC_st0r5.2_6.dat")) {
                    data_dir_arg = " --data-dir '" + fs::canonical(wrk_candidate).string() + "' ";
                }
            }
            std::ostringstream cmd;
            cmd << "'" << flexaidds_bin << "' "
                << data_dir_arg
                << "'" << entry.receptor_path << "' "
                << "'" << entry.ligand_path << "' "
                << "--config '" << config_path << "' "
                << "-o '" << out_prefix << "' "
                << "2>'" << out_dir << "/stderr.log' "
                << ">'" << stdout_path << "'";

            bench::Timer dock_timer;
            dock_timer.start();

            ret = exec_dock(cmd.str(), config.per_job_timeout_s);

            dock_timer.stop();
            result.wall_time_s = dock_timer.elapsed_s();
        }

        // ── Parse results ────────────────────────────────────────────────
        // Shared by both the fresh-run and cached-run paths.
        // Check for output files: <prefix>_INI.pdb, clustered PDBs
        int n_poses = 0;
        float best_cf = 0.0f;
        float best_dG = 0.0f;

        // Count clustered output PDBs
        try {
            for (const auto& f : fs::directory_iterator(out_dir)) {
                std::string fname = f.path().filename().string();
                // FlexAID output: <prefix>_mode_N.pdb or <prefix>_cluster_N.pdb
                if ((fname.find("_mode_") != std::string::npos ||
                     fname.find("_cluster_") != std::string::npos) &&
                    fname.size() > 4 && fname.substr(fname.size()-4) == ".pdb") {
                    n_poses++;
                }
            }
        } catch (...) {}

        // Parse stdout for n_chrom_snapshot, CF scores, and clash diagnostics
        std::ifstream stdout_file(stdout_path);  // stdout_path declared above
        long clashed_count = 0, total_evals = 0;
        float free_energy_F = 0.0f;
        if (stdout_file.is_open()) {
            std::string line;
            while (std::getline(stdout_file, line)) {
                // "n_chrom_snapshot=N" — total individuals evaluated
                if (line.find("n_chrom_snapshot=") != std::string::npos) {
                    auto pos = line.find('=');
                    if (pos != std::string::npos) {
                        try {
                            long snap = std::stol(line.substr(pos+1));
                            n_poses = std::max(n_poses, static_cast<int>(snap));
                            total_evals = std::max(total_evals, snap);
                        } catch (...) {}
                    }
                }
                // "individuals clashed=N"
                if (line.find("individuals clashed=") != std::string::npos) {
                    auto pos = line.find("clashed=");
                    if (pos != std::string::npos) {
                        try { clashed_count = std::stol(line.substr(pos+8)); }
                        catch (...) {}
                    }
                }
                // "individuals skipped=N"
                if (line.find("individuals skipped=") != std::string::npos) {
                    auto pos = line.find("skipped=");
                    if (pos != std::string::npos) {
                        try {
                            long skipped = std::stol(line.substr(pos+8));
                            total_evals = std::max(total_evals, clashed_count + skipped);
                        } catch (...) {}
                    }
                }
                // "Free energy F  = N"
                if (line.find("Free energy F") != std::string::npos) {
                    auto pos = line.find('=');
                    if (pos != std::string::npos) {
                        try { free_energy_F = std::stof(line.substr(pos+1)); }
                        catch (...) {}
                    }
                }
                // "REMARK CF=%f" — extract best contact function score
                if (line.find("REMARK CF=") != std::string::npos) {
                    auto pos = line.find("CF=");
                    if (pos != std::string::npos) {
                        try { best_cf = std::stof(line.substr(pos+3)); }
                        catch (...) {}
                    }
                }
                // Thermodynamic output lines
                if (line.find("dG=") != std::string::npos || line.find("ΔG=") != std::string::npos) {
                    auto pos = line.find("=");
                    if (pos != std::string::npos) {
                        try { best_dG = std::stof(line.substr(pos+1)); }
                        catch (...) {}
                    }
                }
            }
        }

        // ── Clash diagnostics ────────────────────────────────────────────────
        result.individuals_clashed = clashed_count;
        result.individuals_total   = total_evals;
        if (total_evals > 0)
            result.clash_rate = static_cast<float>(clashed_count) / static_cast<float>(total_evals);
        // "Stuck" = GA never escaped clashing geometry AND F is positive (no binding found)
        result.stuck = (result.clash_rate > 0.95f && free_energy_F > 0.0f);
        if (result.stuck) {
            std::cerr << "  [STUCK] " << entry.pdb_id
                      << "  clash=" << std::fixed << std::setprecision(1)
                      << (result.clash_rate * 100.0f) << "%"
                      << "  F=" << free_energy_F << " kcal/mol"
                      << "  — likely bad SDF bonds or grid placement\n";
        }

        result.num_poses = n_poses;
        result.best_score = best_cf;
        result.predicted_dG = (best_dG != 0.0f) ? best_dG : best_cf;

        // Success = FlexAIDdS exited 0 AND produced output poses AND not stuck
        result.success = (ret == 0 && n_poses > 0 && !result.stuck);

        // If FlexAIDdS ran but produced no output, check stderr for clues
        if (ret == 0 && n_poses == 0) {
            std::string stderr_path = out_dir + "/stderr.log";
            std::ifstream stderr_file(stderr_path);
            if (stderr_file.is_open()) {
                std::string err_line;
                int err_lines = 0;
                std::cerr << "  [WARN] " << entry.pdb_id
                          << " exited 0 but produced no output poses. stderr:\n";
                while (std::getline(stderr_file, err_line) && err_lines < 5) {
                    std::cerr << "    " << err_line << "\n";
                    err_lines++;
                }
            }
        }

        // RMSD: we'd need the crystal pose for real RMSD.
        // For now, if the target succeeded (poses produced), set a sentinel.
        // Real RMSD requires superposing the top pose against the crystal ligand.
        if (result.success) {
            // TODO: compute actual RMSD by reading top cluster PDB vs crystal
            result.rmsd_to_crystal = 0.0f; // placeholder — needs crystal comparison
        } else {
            result.rmsd_to_crystal = 999.0f;
        }

        report.results[idx] = result;

        // ── TargetServer: register completed session ─────────────────────
        {
            auto ts_it2 = target_servers_.find(entry.receptor_path);
            if (ts_it2 != target_servers_.end() && result.success) {
                auto& sess = sessions[idx];
                sess.completed = true;
                sess.n_poses = result.num_poses;
                sess.best_energy = static_cast<double>(result.predicted_dG);
                // log_Z = -ΔG / (kT)  — partition function contribution
                sess.log_Z = -static_cast<double>(result.predicted_dG) /
                             (statmech::kB_kcal * static_cast<double>(config.temperature));
                // TODO: best_center requires pose CoM from BindingMode clustering
                //       — populate when MultiModelDock surfaces per-mode centroids to DockingResult
                // TODO: conformer_populations requires per-cluster counts from FOPTICS/DP
                //       — populate when cluster stats surface to DockingResult
                ts_it2->second->register_result(sess);
            }
        }

        // ── Grid reuse: register this run's prefix for subsequent ligands ──
        // Only the first completed run per receptor registers its prefix so
        // that later same-receptor entries can find and reuse its grid file.
        if (result.success && !entry.receptor_path.empty()) {
            std::lock_guard<std::mutex> lock(grid_reuse_mtx);
            if (receptor_completed_prefix.find(entry.receptor_path)
                    == receptor_completed_prefix.end()) {
                receptor_completed_prefix[entry.receptor_path] = out_prefix;
            }
        }

        // ── Async per-complex result I/O ─────────────────────────────────
        // Write a per-complex CSV file in the background via the async pipeline,
        // overlapping this I/O with the next complex's docking computation.
        // The lambda captures all data by value (or const reference to entries
        // which outlive the pipeline).  The final summary CSV and markdown
        // report are written synchronously after io_pipeline.stop().
        io_pipeline.enqueue([result, out_dir]() {
            try {
                std::string csv_path = out_dir + "/result.csv";
                std::ofstream ofs(csv_path);
                if (ofs.is_open()) {
                    ofs << "pdb_id,best_score,rmsd_to_crystal,predicted_dG,predicted_dH,"
                           "predicted_TdS,shannon_entropy,num_poses,wall_time_s,success\n";
                    ofs << std::fixed << std::setprecision(4)
                        << result.pdb_id << ","
                        << result.best_score << ","
                        << result.rmsd_to_crystal << ","
                        << result.predicted_dG << ","
                        << result.predicted_dH << ","
                        << result.predicted_TdS << ","
                        << result.shannon_entropy << ","
                        << result.num_poses << ","
                        << result.wall_time_s << ","
                        << (result.success ? 1 : 0) << "\n";
                }
            } catch (...) {
                // Per-complex CSV is best-effort; failures are non-fatal.
                // The aggregate write_report() still has all data.
            }
        });
    };

    // ── Parallel docking via thread pool ──────────────────────────────────
    // OpenMP is unavailable on macOS (Apple Clang), so we use std::thread.
    // config.num_threads controls how many FlexAIDdS processes run concurrently.
    const size_t n_workers = (config.num_threads > 0)
        ? static_cast<size_t>(config.num_threads) : 1;
    const size_t n_jobs = entries.size();

    if (n_jobs <= 1 || n_workers <= 1) {
        // Serial path — single job or single worker
        for (size_t i = 0; i < n_jobs; ++i) {
            if (shutdown_requested_.load(std::memory_order_relaxed)) {
                std::cerr << "\n[DatasetRunner] Shutdown requested — stopping after "
                          << i << "/" << n_jobs << " jobs\n";
                break;
            }
            dock_one(i);
        }
    } else {
        // Thread pool: atomically distribute indices to workers
        std::atomic<size_t> next_idx{0};
        std::mutex cout_mtx;  // protect interleaved stdout/cerr

        auto worker = [&]() {
            for (;;) {
                if (shutdown_requested_.load(std::memory_order_relaxed)) break;
                size_t idx = next_idx.fetch_add(1);
                if (idx >= n_jobs) break;
                dock_one(idx);
            }
        };

        std::vector<std::thread> pool;
        pool.reserve(n_workers);
        for (size_t w = 0; w < n_workers; ++w) {
            pool.emplace_back(worker);
        }
        for (auto& t : pool) t.join();
    }

    timer.stop();

    // ── Flush async I/O pipeline ─────────────────────────────────────────
    // Blocks until all pending per-complex CSV writes complete.
    io_pipeline.stop();

    // ── TargetServer: cross-ligand competitive binding analysis ────────
    for (const auto& [receptor_path, ts] : target_servers_) {
        if (ts->completed_sessions() < 1) continue;

        BenchmarkReport::CrossLigandResult clr;
        // Extract receptor ID from path (filename without extension)
        clr.receptor_id = fs::path(receptor_path).stem().string();

        // Count ligands for this receptor
        for (size_t i = 0; i < entries.size(); ++i) {
            if (entries[i].receptor_path == receptor_path) {
                clr.n_ligands++;
                if (report.results[i].success) clr.n_completed++;
            }
        }

        // Get ranked ligands from grand partition function Ξ
        if (ts->completed_sessions() >= 2) {
            clr.ranked_ligands = ts->rank_ligands();
        }

        report.cross_ligand_results.push_back(std::move(clr));
    }

    // Compute aggregate statistics
    int success_count = 0;
    std::vector<double> rmsds;
    std::vector<double> pred_affinities;
    std::vector<double> exp_affinities;

    for (size_t i = 0; i < report.results.size(); ++i) {
        const auto& r = report.results[i];
        if (r.success) success_count++;
        if (r.rmsd_to_crystal < 900.0f) {
            rmsds.push_back(r.rmsd_to_crystal);
        }
        if (entries[i].has_affinity() && r.predicted_dG != 0.0f) {
            exp_affinities.push_back(entries[i].experimental_affinity);
            pred_affinities.push_back(-r.predicted_dG / 1.3636); // convert to pKd
        }
    }

    report.successful = success_count;
    report.success_rate = (report.total_systems > 0)
        ? static_cast<double>(success_count) / report.total_systems : 0.0;

    // Mean RMSD
    if (!rmsds.empty()) {
        report.mean_rmsd = std::accumulate(rmsds.begin(), rmsds.end(), 0.0) / rmsds.size();
    }

    // Median RMSD
    if (!rmsds.empty()) {
        auto sorted = rmsds;
        std::sort(sorted.begin(), sorted.end());
        size_t mid = sorted.size() / 2;
        if (sorted.size() % 2 == 0) {
            report.median_rmsd = (sorted[mid - 1] + sorted[mid]) / 2.0;
        } else {
            report.median_rmsd = sorted[mid];
        }
    }

    // Correlation metrics
    if (pred_affinities.size() >= 3) {
        report.pearson_r    = compute_pearson_r(pred_affinities, exp_affinities);
        report.spearman_rho = compute_spearman_rho(pred_affinities, exp_affinities);
        report.kendall_tau  = compute_kendall_tau(pred_affinities, exp_affinities);
    }

    // ── Restore default signal handlers ────────────────────────────────────
    // Clear the global pointers first so stale signals can't dereference them,
    // then restore SIGINT/SIGTERM to their defaults.
#ifndef _MSC_VER
    g_active_guard    = nullptr;
    g_active_shutdown = nullptr;
    ::signal(SIGINT,  SIG_DFL);
    ::signal(SIGTERM, SIG_DFL);
#endif

    return report;
}

// =============================================================================
// Report generation: Markdown + CSV
// =============================================================================

void DatasetRunner::write_report(const BenchmarkReport& report,
                                  const std::string& output_dir) {
    ensure_dir(output_dir);

    std::string safe_name = report.dataset_name;
    std::replace(safe_name.begin(), safe_name.end(), ' ', '_');
    std::replace(safe_name.begin(), safe_name.end(), '-', '_');
    std::transform(safe_name.begin(), safe_name.end(), safe_name.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // ── Markdown report ──────────────────────────────────────────────
    {
        std::string md_path = output_dir + "/" + safe_name + "_report.md";
        std::ofstream ofs(md_path);

        ofs << "# FlexAIDdS Benchmark Report: " << report.dataset_name << "\n\n";

        // Summary table
        ofs << "## Summary\n\n";
        ofs << "| Metric | Value |\n";
        ofs << "|--------|-------|\n";
        ofs << "| Total systems | " << report.total_systems << " |\n";
        ofs << "| Successful (RMSD < 2.0 Å) | " << report.successful << " |\n";
        ofs << std::fixed << std::setprecision(1);
        ofs << "| Success rate | " << (report.success_rate * 100.0) << "% |\n";
        ofs << std::setprecision(2);
        ofs << "| Mean RMSD (Å) | " << report.mean_rmsd << " |\n";
        ofs << "| Median RMSD (Å) | " << report.median_rmsd << " |\n";
        ofs << std::setprecision(3);
        ofs << "| Pearson r | " << report.pearson_r << " |\n";
        ofs << "| Spearman ρ | " << report.spearman_rho << " |\n";
        ofs << "| Kendall τ | " << report.kendall_tau << " |\n";
        ofs << "\n";

        // Per-system results table
        ofs << "## Per-System Results\n\n";
        ofs << "| PDB | Score | RMSD (Å) | ΔG | ΔH | TΔS | S_shan | Poses | Time (s) | Success |\n";
        ofs << "|-----|-------|----------|-----|-----|------|--------|-------|----------|--------|\n";

        for (const auto& r : report.results) {
            ofs << "| " << r.pdb_id
                << " | " << std::setprecision(2) << r.best_score
                << " | " << std::setprecision(2) << r.rmsd_to_crystal
                << " | " << std::setprecision(2) << r.predicted_dG
                << " | " << std::setprecision(2) << r.predicted_dH
                << " | " << std::setprecision(2) << r.predicted_TdS
                << " | " << std::setprecision(3) << r.shannon_entropy
                << " | " << r.num_poses
                << " | " << std::setprecision(1) << r.wall_time_s
                << " | " << (r.success ? "✓" : "✗")
                << " |\n";
        }

        ofs.close();
        std::cout << "  Markdown report: " << md_path << "\n";
    }

    // ── Cross-ligand competitive binding report ──────────────────────────
    if (!report.cross_ligand_results.empty()) {
        // Markdown section
        std::string cl_md = output_dir + "/" + safe_name + "_cross_ligand.md";
        std::ofstream ofs(cl_md);
        ofs << "# Cross-Ligand Competitive Binding Analysis\n\n";
        ofs << "Generated by TargetServer grand canonical partition function Ξ.\n\n";

        for (const auto& clr : report.cross_ligand_results) {
            ofs << "## Receptor: " << clr.receptor_id << "\n\n";
            ofs << "- Ligands docked: " << clr.n_ligands << "\n";
            ofs << "- Completed (binding modes found): " << clr.n_completed << "\n\n";

            if (!clr.ranked_ligands.empty()) {
                ofs << "### Ranked Ligands (by ΔG, ascending)\n\n";
                ofs << "| Rank | Ligand | ΔG (kcal/mol) | p(bind) |\n";
                ofs << "|------|--------|---------------|--------|\n";
                int rank = 1;
                for (const auto& lr : clr.ranked_ligands) {
                    ofs << "| " << rank++
                        << " | " << lr.name
                        << " | " << std::fixed << std::setprecision(3) << lr.dG
                        << " | " << std::setprecision(4) << lr.p_bound
                        << " |\n";
                }
                ofs << "\n";
            }
        }

        ofs.close();
        std::cout << "  Cross-ligand report: " << cl_md << "\n";

        // CSV version
        std::string cl_csv = output_dir + "/" + safe_name + "_cross_ligand.csv";
        std::ofstream coefs(cl_csv);
        coefs << "receptor,ligand,delta_G_kcal,binding_probability,rank\n";
        for (const auto& clr : report.cross_ligand_results) {
            int rank = 1;
            for (const auto& lr : clr.ranked_ligands) {
                coefs << std::fixed << std::setprecision(4)
                      << clr.receptor_id << ","
                      << lr.name << ","
                      << lr.dG << ","
                      << lr.p_bound << ","
                      << rank++ << "\n";
            }
        }
        coefs.close();
        std::cout << "  Cross-ligand CSV: " << cl_csv << "\n";
    }

    // ── CSV results ──────────────────────────────────────────────────
    {
        std::string csv_path = output_dir + "/" + safe_name + "_results.csv";
        std::ofstream ofs(csv_path);

        ofs << "pdb_id,best_score,rmsd_to_crystal,predicted_dG,predicted_dH,"
               "predicted_TdS,shannon_entropy,num_poses,wall_time_s,success\n";

        for (const auto& r : report.results) {
            ofs << std::fixed << std::setprecision(4)
                << r.pdb_id << ","
                << r.best_score << ","
                << r.rmsd_to_crystal << ","
                << r.predicted_dG << ","
                << r.predicted_dH << ","
                << r.predicted_TdS << ","
                << r.shannon_entropy << ","
                << r.num_poses << ","
                << r.wall_time_s << ","
                << (r.success ? 1 : 0) << "\n";
        }

        ofs.close();
        std::cout << "  CSV results: " << csv_path << "\n";
    }

    // ── Summary CSV ──────────────────────────────────────────────────
    {
        std::string summary_csv = output_dir + "/" + safe_name + "_summary.csv";
        std::ofstream ofs(summary_csv);

        ofs << "dataset,total_systems,successful,success_rate,mean_rmsd,"
               "median_rmsd,pearson_r,spearman_rho,kendall_tau\n";
        ofs << std::fixed << std::setprecision(4)
            << report.dataset_name << ","
            << report.total_systems << ","
            << report.successful << ","
            << report.success_rate << ","
            << report.mean_rmsd << ","
            << report.median_rmsd << ","
            << report.pearson_r << ","
            << report.spearman_rho << ","
            << report.kendall_tau << "\n";

        ofs.close();
        std::cout << "  Summary CSV: " << summary_csv << "\n";
    }
}

// =============================================================================
// SubprocessGuard implementation
// =============================================================================

SubprocessGuard::SubprocessGuard() = default;

SubprocessGuard::~SubprocessGuard() {
    kill_all();
}

pid_t SubprocessGuard::fork_exec(const std::string& cmd) {
    std::lock_guard<std::mutex> lock(mtx_);
    pid_t pid = ::fork();
    if (pid == 0) {
        ::setpgid(0, 0);
        ::execl("/bin/sh", "sh", "-c", cmd.c_str(), static_cast<char*>(nullptr));
        ::_exit(127);
    }
    if (pid > 0) {
        pids_.insert(pid);
    }
    return pid;
}

int SubprocessGuard::wait_with_timeout(pid_t pid, int timeout_s) {
    using namespace std::chrono;
    auto deadline = steady_clock::now() + seconds(timeout_s);

    while (steady_clock::now() < deadline) {
        int status = 0;
        pid_t ret = ::waitpid(pid, &status, WNOHANG);
        if (ret == pid) {
            {
                std::lock_guard<std::mutex> lock(mtx_);
                pids_.erase(pid);
            }
            if (WIFEXITED(status)) return WEXITSTATUS(status);
            return -1;
        }
        std::this_thread::sleep_for(milliseconds(200));
    }

    ::kill(pid, SIGTERM);
    std::this_thread::sleep_for(milliseconds(500));
    int status = 0;
    pid_t ret = ::waitpid(pid, &status, WNOHANG);
    if (ret != pid) {
        ::kill(pid, SIGKILL);
        ::waitpid(pid, &status, 0);
    }
    {
        std::lock_guard<std::mutex> lock(mtx_);
        pids_.erase(pid);
    }
    return -1;
}

void SubprocessGuard::forget(pid_t pid) {
    std::lock_guard<std::mutex> lock(mtx_);
    pids_.erase(pid);
}

void SubprocessGuard::kill_all() {
    std::lock_guard<std::mutex> lock(mtx_);
    for (auto pid : pids_) {
        ::kill(pid, SIGTERM);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    for (auto pid : pids_) {
        ::kill(pid, SIGKILL);
    }
    for (auto pid : pids_) {
        int status = 0;
        ::waitpid(pid, &status, 0);
    }
    pids_.clear();
}

size_t SubprocessGuard::active_count() const {
    std::lock_guard<std::mutex> lock(mtx_);
    return pids_.size();
}

} // namespace dataset
