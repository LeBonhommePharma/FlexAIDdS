// =============================================================================
// probe_cf.cpp — frozen-pose CF scoring diagnostic (no GA, no search)
//
// Scores a receptor + a FROZEN ligand pose through the exact production
// vcfunction()/ic2cf() path and prints a decomposed CF breakdown as JSON.
//
// Design
// ------
// The CF decomposition must be BYTE-FAITHFUL to what the docking engine
// computes during a real run, otherwise it is useless as a regression
// instrument.  Reproducing that in-process would require duplicating ~300
// lines of top.cpp init (ProcessLigand SYBYL typing enrichment, formal-charge
// assignment, type256 donor/acceptor population, apply_config, cleft/grid/
// optres setup) — every one a silent-divergence risk.
//
// Instead probe_cf DRIVES the engine's own already-wired frozen-pose scorer:
//   FLEXAIDDS_SCORE_NATIVE=1 + FLEXAIDDS_NATIVE_ONLY=1 make top.cpp run the
//   full production init and then call score_native_pose() (native_score.cpp),
//   which overrides the ligand atoms[].coor[] with the pose read from
//   FLEXAIDDS_RMSDST and calls vcfunction() DIRECTLY — no GA — then std::exit()s
//   before the search.  It prints:
//       [NATIVE_CF] cf=<total> breakdown=com:..,wal:..,sas:..,con:..,elec:..,hbond:..,gist_desolv:..,metal_coord:..,entropy:..,pb_clash:..
//   probe_cf captures that line and reformats it as JSON.
//
// This wraps the existing score_native_pose() diagnostic (per the task's
// "wrap, don't reimplement" guidance) and therefore reproduces the engine's
// numbers exactly (validated: 1G9V native = -49.612113, com:-122.9833).
//
// Modes
// -----
//   --mode direct     : engine loads --ligand (topology+typing), score_native
//                       overrides its coords with the pose (FLEXAIDDS_RMSDST).
//                       Scores the cartesian pose coordinates directly.
//   --mode roundtrip  : engine loads the POSE itself as the ligand, builds its
//                       internal coordinates, and score_native reconstructs via
//                       ic2cf(FA->opt_par) — exercising the encoder.  The
//                       [NATIVE-SEED-RMSD] round-trip error is reported too.
//
// CF term mapping (engine cf_str -> JSON)
// ---------------------------------------
//   cf_total  = get_cf_evalue()          (sum of the terms below)
//   cf_vct    = com   (VCT complementarity; same field as cf_com)
//   cf_com    = com
//   cf_con    = con   (constraint term)
//   cf_wal    = wal   (soft/hard wall)
//   cf_sas    = sas   (solvent-accessible-surface penalty)
//   cf_hbond  = hbond (angular H-bond energy)
//   cf_clash  = pb_clash (PoseBust-basis physical-realism clash penalty)
//
// Copyright 2026 Le Bonhomme Pharma.  SPDX-License-Identifier: Apache-2.0
// =============================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cctype>
#include <unistd.h>

namespace fs = std::filesystem;

// Escape a string for embedding in a JSON string literal: backslash, quote,
// the shorthand control escapes, and any remaining control char < 0x20 as \u00XX.
static std::string json_escape(const std::string& s) {
    std::string r;
    r.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
            case '\\': r += "\\\\"; break;
            case '\"': r += "\\\""; break;
            case '\b': r += "\\b";  break;
            case '\f': r += "\\f";  break;
            case '\n': r += "\\n";  break;
            case '\r': r += "\\r";  break;
            case '\t': r += "\\t";  break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    r += buf;
                } else {
                    r += static_cast<char>(c);
                }
        }
    }
    return r;
}

// ── CLI ──────────────────────────────────────────────────────────────────────
struct Args {
    std::string receptor;
    std::string ligand;      // topology SDF/MOL2 (defaults to pose if omitted)
    std::string pose;        // coords to score (SDF/PDB/MOL2); defaults to ligand
    std::string mode = "direct";
    std::string config;      // dock_config.json (optional, passed as -c)
    std::string binary;      // FlexAIDdS executable
    std::string data_dir;    // dir holding MC_st0r5.2_6.dat / AMINO.def
    std::string pdb = "?";   // label for JSON
    bool keep_tmp = false;
};

static void usage(const char* p) {
    std::fprintf(stderr,
        "probe_cf — frozen-pose CF scoring diagnostic (no GA)\n\n"
        "Usage: %s --receptor REC.pdb --pose POSE.{sdf|pdb} [options]\n\n"
        "  --receptor PATH   receptor PDB (required)\n"
        "  --pose PATH       ligand pose to score: SDF, or PDB (ligand extracted)\n"
        "  --ligand PATH     ligand topology SDF/MOL2 (default: --pose if it is SDF)\n"
        "  --mode MODE       direct | roundtrip                (default: direct)\n"
        "  --config PATH     dock_config.json (passed to engine as -c)\n"
        "  --binary PATH     FlexAIDdS executable (default: autodetect near cwd)\n"
        "  --data-dir PATH   dir with MC_st0r5.2_6.dat/AMINO.def (default: binary dir)\n"
        "  --pdb ID          label emitted in JSON\n"
        "  --keep-tmp        keep the generated pose SDF\n",
        p);
}

// ── pose coordinate parsing ──────────────────────────────────────────────────
struct XYZ { double x, y, z; };

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Read the V2000 atom block of the first molecule of an SDF/MOL file.
static bool parse_sdf_coords(const std::string& path, std::vector<XYZ>& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string l;
    std::vector<std::string> lines;
    while (std::getline(f, l)) lines.push_back(l);
    if (lines.size() < 4) return false;
    int na = std::atoi(lines[3].substr(0, 3).c_str());
    if (na <= 0 || (int)lines.size() < 4 + na) return false;
    for (int i = 0; i < na; ++i) {
        const std::string& a = lines[4 + i];
        if (a.size() < 30) return false;
        out.push_back({ std::atof(a.substr(0, 10).c_str()),
                        std::atof(a.substr(10, 10).c_str()),
                        std::atof(a.substr(20, 10).c_str()) });
    }
    return true;
}

// Extract ligand-atom coordinates from a PDB pose.  Prefers records whose
// residue name matches `resname` (the SDF title); falls back to the last
// `want` ATOM/HETATM records if that count does not match.
static bool parse_pdb_lig_coords(const std::string& path, const std::string& resname,
                                 int want, std::vector<XYZ>& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string l;
    std::vector<XYZ> named, all;
    while (std::getline(f, l)) {
        if (l.rfind("ATOM", 0) != 0 && l.rfind("HETATM", 0) != 0) continue;
        if (l.size() < 54) continue;
        XYZ p{ std::atof(l.substr(30, 8).c_str()),
               std::atof(l.substr(38, 8).c_str()),
               std::atof(l.substr(46, 8).c_str()) };
        all.push_back(p);
        std::string rn = trim(l.substr(17, 3));
        if (!resname.empty() && rn == resname) named.push_back(p);
    }
    if ((int)named.size() == want) { out = named; return true; }
    if ((int)all.size() >= want) {
        std::fprintf(stderr,
            "WARN [probe_cf]: ligand residue-name match failed; falling back to "
            "last %d ATOM/HETATM records — verify the pose PDB is ligand-only or "
            "pass an SDF pose\n", want);
        out.assign(all.end() - want, all.end());   // last `want` records
        return true;
    }
    return false;
}

// Read an SDF's title (line 0) and atom count (line 3).
static bool read_sdf_meta(const std::string& path, std::string& title, int& na,
                          std::vector<std::string>& lines) {
    std::ifstream f(path);
    if (!f) return false;
    std::string l;
    while (std::getline(f, l)) lines.push_back(l);
    if (lines.size() < 4) return false;
    title = trim(lines[0]);
    na = std::atoi(lines[3].substr(0, 3).c_str());
    return na > 0 && (int)lines.size() >= 4 + na;
}

// Write `tmpl` SDF with the atom-block XYZ replaced by `coords`.
static bool write_templated_sdf(const std::vector<std::string>& tmpl, int na,
                                const std::vector<XYZ>& coords,
                                const std::string& out_path) {
    if ((int)coords.size() != na) return false;
    std::ofstream o(out_path);
    if (!o) return false;
    for (int i = 0; i < 4; ++i) o << tmpl[i] << "\n";
    for (int i = 0; i < na; ++i) {
        const std::string& row = tmpl[4 + i];
        std::string rest = row.size() > 30 ? row.substr(30) : "";
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%10.4f%10.4f%10.4f",
                      coords[i].x, coords[i].y, coords[i].z);
        o << buf << rest << "\n";
    }
    for (size_t i = 4 + na; i < tmpl.size(); ++i) o << tmpl[i] << "\n";
    return true;
}

// ── engine invocation ────────────────────────────────────────────────────────
static std::string shq(const std::string& s) {
    std::string r = "'";
    for (char c : s) { if (c == '\'') r += "'\\''"; else r += c; }
    return r + "'";
}

int main(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (k == "--receptor") a.receptor = next();
        else if (k == "--ligand")   a.ligand   = next();
        else if (k == "--pose")     a.pose     = next();
        else if (k == "--mode")     a.mode     = next();
        else if (k == "--config" || k == "-c") a.config = next();
        else if (k == "--binary")   a.binary   = next();
        else if (k == "--data-dir") a.data_dir = next();
        else if (k == "--pdb")      a.pdb      = next();
        else if (k == "--keep-tmp") a.keep_tmp = true;
        else if (k == "-h" || k == "--help") { usage(argv[0]); return 0; }
        else { std::fprintf(stderr, "unknown arg: %s\n", k.c_str()); usage(argv[0]); return 2; }
    }

    if (a.receptor.empty()) { std::fprintf(stderr, "ERROR: --receptor required\n"); return 2; }
    if (a.pose.empty() && a.ligand.empty()) {
        std::fprintf(stderr, "ERROR: at least one of --pose / --ligand required\n"); return 2;
    }
    if (a.mode != "direct" && a.mode != "roundtrip") {
        std::fprintf(stderr, "ERROR: --mode must be direct|roundtrip\n"); return 2;
    }
    if (a.ligand.empty()) a.ligand = a.pose;   // pose doubles as topology
    if (a.pose.empty())   a.pose   = a.ligand; // score the ligand's own coords

    // Resolve binary / data-dir.
    if (a.binary.empty()) {
        for (const char* c : { "build/FlexAIDdS", "build_ltofix/FlexAIDdS",
                               "./FlexAIDdS", "FlexAIDdS" }) {
            if (fs::exists(c)) { a.binary = c; break; }
        }
    }
    if (a.binary.empty() || !fs::exists(a.binary)) {
        std::fprintf(stderr, "ERROR: FlexAIDdS binary not found (use --binary)\n"); return 2;
    }
    a.binary = fs::absolute(a.binary).string();
    if (a.data_dir.empty()) a.data_dir = fs::path(a.binary).parent_path().string();

    // Determine the ligand-topology SDF (must be SDF for pose-coord templating).
    // If --ligand is MOL2 there is no coord override — engine scores it as-is.
    bool ligand_is_sdf = false;
    {
        std::string e = fs::path(a.ligand).extension().string();
        for (auto& c : e) c = std::tolower(c);
        ligand_is_sdf = (e == ".sdf" || e == ".mol");
    }

    // Build the RMSDST pose SDF (coords to score).
    std::string rmsdst = a.pose;
    std::string tmp_pose;
    bool made_tmp = false;

    bool pose_same_as_ligand = (fs::weakly_canonical(a.pose) == fs::weakly_canonical(a.ligand));
    if (!pose_same_as_ligand) {
        if (!ligand_is_sdf) {
            std::fprintf(stderr,
                "ERROR: pose-coordinate override requires an SDF --ligand template "
                "(got '%s'). Provide the crystal/topology ligand as SDF.\n",
                a.ligand.c_str());
            return 2;
        }
        std::string title; int na = 0; std::vector<std::string> tmpl;
        if (!read_sdf_meta(a.ligand, title, na, tmpl)) {
            std::fprintf(stderr, "ERROR: cannot read ligand SDF: %s\n", a.ligand.c_str()); return 2;
        }
        std::vector<XYZ> coords;
        std::string pe = fs::path(a.pose).extension().string();
        for (auto& c : pe) c = std::tolower(c);
        bool ok;
        if (pe == ".pdb" || pe == ".ent")
            ok = parse_pdb_lig_coords(a.pose, title, na, coords);
        else
            ok = parse_sdf_coords(a.pose, coords);
        if (!ok || (int)coords.size() != na) {
            std::fprintf(stderr,
                "ERROR: could not extract %d pose atoms from %s (got %zu)\n",
                na, a.pose.c_str(), coords.size());
            return 2;
        }
        tmp_pose = (fs::temp_directory_path() /
                    ("probe_cf_pose_" + std::to_string(::getpid()) + ".sdf")).string();
        if (!write_templated_sdf(tmpl, na, coords, tmp_pose)) {
            std::fprintf(stderr, "ERROR: failed writing templated pose SDF\n"); return 2;
        }
        rmsdst = tmp_pose;
        made_tmp = true;
    }

    // In roundtrip mode the engine must LOAD the pose as the ligand so its IC is
    // built from the pose coordinates and ic2cf() reconstructs them.
    std::string engine_ligand = (a.mode == "roundtrip") ? rmsdst : a.ligand;

    // Compose the command (stderr+stdout merged so we can grep the diagnostic).
    std::ostringstream cmd;
    cmd << "FLEXAIDDS_SCORE_NATIVE=1 FLEXAIDDS_NATIVE_ONLY=1 "
        << "FLEXAIDDS_RMSDST=" << shq(rmsdst) << " "
        << shq(a.binary) << " " << shq(a.receptor) << " " << shq(engine_ligand);
    if (!a.config.empty()) cmd << " -c " << shq(a.config);
    cmd << " --data-dir " << shq(a.data_dir) << " 2>&1";

    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (!pipe) { std::fprintf(stderr, "ERROR: popen failed\n"); return 3; }

    std::string native_line, seed_line;
    char line[4096];
    while (std::fgets(line, sizeof(line), pipe)) {
        std::string s(line);
        if (s.find("[NATIVE_CF]") != std::string::npos && s.find("breakdown=") != std::string::npos)
            native_line = s;
        else if (s.find("[NATIVE-SEED-RMSD]") != std::string::npos)
            seed_line = s;
    }
    int rc = pclose(pipe);

    if (made_tmp && !a.keep_tmp) std::remove(tmp_pose.c_str());

    if (native_line.empty()) {
        std::fprintf(stderr,
            "ERROR: no [NATIVE_CF] line captured (engine rc=%d). "
            "Re-run with --keep-tmp and check the engine output.\n", rc);
        return 4;
    }

    // Parse: [NATIVE_CF] cf=<t> breakdown=com:<>,wal:<>,sas:<>,con:<>,elec:<>,hbond:<>,gist_desolv:<>,metal_coord:<>,entropy:<>,pb_clash:<>
    double cf_total = 0, com = 0, wal = 0, sas = 0, con = 0, elec = 0, hbond = 0;
    double gist_desolv = 0, metal_coord = 0, entropy = 0, pb_clash = 0;
    {
        size_t p = native_line.find("cf=");
        if (p != std::string::npos) cf_total = std::atof(native_line.c_str() + p + 3);
        auto grab = [&](const char* key, double& dst) {
            std::string k = std::string(key) + ":";
            size_t q = native_line.find(k);
            if (q != std::string::npos) dst = std::atof(native_line.c_str() + q + k.size());
        };
        grab("com", com); grab("wal", wal); grab("sas", sas);
        grab("con", con); grab("elec", elec); grab("hbond", hbond);
        grab("gist_desolv", gist_desolv); grab("metal_coord", metal_coord);
        grab("entropy", entropy); grab("pb_clash", pb_clash);
    }

    double roundtrip_rmsd = -1.0;
    if (!seed_line.empty()) {
        size_t p = seed_line.find("RMSD = ");
        if (p != std::string::npos) roundtrip_rmsd = std::atof(seed_line.c_str() + p + 7);
    }

    // Emit JSON.
    std::printf("{");
    std::printf("\"pdb\": \"%s\", ", json_escape(a.pdb).c_str());
    std::printf("\"mode\": \"%s\", ", a.mode.c_str());
    std::printf("\"cf_total\": %.6f, ", cf_total);
    std::printf("\"cf_vct\": %.6f, ", com);
    std::printf("\"cf_com\": %.6f, ", com);
    std::printf("\"cf_con\": %.6f, ", con);
    std::printf("\"cf_hbond\": %.6f, ", hbond);
    std::printf("\"cf_sas\": %.6f, ", sas);
    std::printf("\"cf_wal\": %.6f, ", wal);
    std::printf("\"cf_elec\": %.6f, ", elec);
    std::printf("\"cf_gist_desolv\": %.6f, ", gist_desolv);
    std::printf("\"cf_metal_coord\": %.6f, ", metal_coord);
    std::printf("\"cf_entropy\": %.6f, ", entropy);
    std::printf("\"cf_clash\": %.6f, ", pb_clash);
    if (a.mode == "roundtrip")
        std::printf("\"roundtrip_rmsd\": %.4f, ", roundtrip_rmsd);
    std::printf("\"receptor\": \"%s\", ", json_escape(a.receptor).c_str());
    std::printf("\"pose\": \"%s\"", json_escape(a.pose).c_str());
    std::printf("}\n");
    return 0;
}
