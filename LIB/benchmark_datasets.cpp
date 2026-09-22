// =============================================================================
// benchmark_datasets.cpp — Standalone benchmark dataset runner executable
//
// Usage:
//   benchmark_datasets --benchmark astex [--output results/] [--threads 8]
//   benchmark_datasets --benchmark casf2016
//   benchmark_datasets --benchmark all
//   benchmark_datasets --benchmark doi:10.1021/acs.jcim.3c00817
//   benchmark_datasets --benchmark pdb_list:my_targets.txt
//
// Supported benchmarks:
//   astex, astex_nonnative, hap2, casf2016, posebusters, dude,
//   bindingdb_itc, sampl6, sampl7, pdbbind, all
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

// PROCESS EXIT CODES — a shell driver must be able to branch WITHOUT parsing
// prose. Until 2026-09-20 every runtime failure left here as exit 2 with the
// single line "ERROR: Incomplete docking", so an environment fault that killed
// the engine at startup was indistinguishable from nine genuinely hard targets.
// It took six hours and a manual stderr.log read to tell them apart.
//
//     0   every target completed at runtime level (poses written, clean exit)
//     1   CLI / usage error                         [pre-existing, unchanged]
//     2   generic incomplete docking, or Fleet chunk with no target results
//                                                   [pre-existing, unchanged —
//         still the value for any class this file has not learned to name, so
//         a driver testing `rc -eq 2` keeps working]
//     3   Fleet result publication failure          [pre-existing, unchanged]
//
//    20   engine_startup_abort    engine terminated through its own top-level
//                                 handler with zero poses written
//    21   engine_crash_midrun     same handler fired, but poses already existed
//    22   timeout                 child killed after per_job_timeout_s
//    23   engine_killed_signal    child did not exit normally, well under timeout
//    24   engine_nonzero_exit     clean exit, nonzero status, no handler line
//    25   ga_stuck_clashes        GA drowned in clashes (DockingResult::stuck)
//    26   zero_poses_written      exit 0, ran to completion, wrote no poses
//    27   engine_not_launched     no child was ever started for this target
//    28   cached_row_incomplete   skip-cache replayed a row marked incomplete
//    29   no_target_results       report carries no systems at all
//    30   result_count_mismatch   results.size() != total_systems
//
// `rc != 0` remains the universal failure test; the granular codes are additive.
// The dominant class is the most diagnostic-urgent one present, not the most
// frequent: one engine_startup_abort in a chunk of eighty-five is an
// environment fault that invalidates every target launched after it, whereas
// eighty-four zero_poses_written rows are a docking result.
// =============================================================================

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// The diagnosis classifier below is deliberately free of DatasetRunner /
// BenchmarkRunner / FleetRunner so that tests/test_startup_diagnosis.cpp can
// #include this translation unit and link against libc++ alone. Defining
// FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY compiles the classifier and
// nothing else. Two guarded regions, both closed at the end of the file.
#ifndef FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY
#include "DatasetRunner.h"
#include "BenchmarkRunner.h"
#include "FleetRunner.h"
#include "fs_safe.h"
#endif

namespace fs = std::filesystem;

// =============================================================================
// startup_diag — name the failure class instead of the symptom
//
// WHY THIS EXISTS
// ---------------
// benchmark_runtime_exit_code() (LIB/DatasetRunner.h) collapses six distinct
// conditions into one nonzero return, and this file printed one sentence for
// all of them. On 2026-09-20 nine cells of an 85-target Astex campaign died
// because $TMPDIR named a directory an agent harness had swept; the engine
// threw out of read_input/direct_input, LIB/top.cpp's catch(std::exception)
// printed "Fatal error: filesystem error: in temp_directory_path: ..." to the
// child's stderr.log and returned 1, and the parent reported "Incomplete
// docking". Nothing was wrong with those nine targets.
//
// The repair is not a better sentence. It is a classification the parent can
// make from evidence it already holds, written to a machine-readable sidecar
// and mapped onto a distinct process exit code, so a driver branches on rc and
// a human reads the class name — neither has to grep a child log to find out
// that the failure was environmental.
//
// WHAT THE CLASSIFIER IS ALLOWED TO USE
// -------------------------------------
// Only two kinds of evidence, both already on hand:
//   (1) the per-target runtime fields of DockingResult — docking_exit_code,
//       num_poses, docking_completed, stuck, wall_time_s;
//   (2) two filesystem facts about <output_dir>/<pdb_id>/ (the layout fixed at
//       DatasetRunner.cpp:6330) — whether the directory and its stderr.log
//       exist, and whether that log's tail carries the engine's own top-level
//       handler prefix.
//
// The vocabulary is closed and every term is decidable from the above. Two
// terms that suggest themselves are deliberately ABSENT because this site
// cannot produce them:
//
//   * poses_written_but_all_filtered — a target with poses > 0, exit 0 and
//     !stuck has docking_completed == true and never trips
//     benchmark_runtime_exit_code at all. Pose-level rejection is recorded in
//     DockingResult::rmsd_fail_reason, which is a different axis. A class here
//     would be unreachable.
//   * input_missing — BenchmarkReport carries no input paths and main() has no
//     entries vector, so the cause cannot be read. The observable is that no
//     child process was ever started, which is named engine_not_launched.
//     Calling that input_missing would assert a cause from an absence.
// =============================================================================
namespace startup_diag {

enum class Class : int {
    None = 0,
    NoTargetResults,
    ResultCountMismatch,
    EngineNotLaunched,
    EngineStartupAbort,
    EngineCrashMidrun,
    Timeout,
    EngineKilledSignal,
    EngineNonzeroExit,
    GaStuckClashes,
    ZeroPosesWritten,
    CachedRowIncomplete,
};

// Stable strings. These land in a CSV column and a JSON field; downstream
// readers match on them, so they are renamed only with the readers.
inline const char* to_string(Class c) {
    switch (c) {
        case Class::None:                return "none";
        case Class::NoTargetResults:     return "no_target_results";
        case Class::ResultCountMismatch: return "result_count_mismatch";
        case Class::EngineNotLaunched:   return "engine_not_launched";
        case Class::EngineStartupAbort:  return "engine_startup_abort";
        case Class::EngineCrashMidrun:   return "engine_crash_midrun";
        case Class::Timeout:             return "timeout";
        case Class::EngineKilledSignal:  return "engine_killed_signal";
        case Class::EngineNonzeroExit:   return "engine_nonzero_exit";
        case Class::GaStuckClashes:      return "ga_stuck_clashes";
        case Class::ZeroPosesWritten:    return "zero_poses_written";
        case Class::CachedRowIncomplete: return "cached_row_incomplete";
    }
    return "unclassified";
}

// Exit-code map. The table lives in the banner at the top of this file; this
// function is the single implementation of it. 0/1/2/3 are pre-existing and
// are not reassigned.
inline int exit_code_for(Class c) {
    switch (c) {
        case Class::None:                return 0;
        case Class::EngineStartupAbort:  return 20;
        case Class::EngineCrashMidrun:   return 21;
        case Class::Timeout:             return 22;
        case Class::EngineKilledSignal:  return 23;
        case Class::EngineNonzeroExit:   return 24;
        case Class::GaStuckClashes:      return 25;
        case Class::ZeroPosesWritten:    return 26;
        case Class::EngineNotLaunched:   return 27;
        case Class::CachedRowIncomplete: return 28;
        case Class::NoTargetResults:     return 29;
        case Class::ResultCountMismatch: return 30;
    }
    return 2;  // unnamed class keeps the legacy generic code
}

// Diagnostic urgency, not frequency. An environment fault ranks above any
// docking outcome because it invalidates everything launched after it.
inline int severity_rank(Class c) {
    switch (c) {
        case Class::None:                return 0;
        case Class::CachedRowIncomplete: return 1;
        case Class::ZeroPosesWritten:    return 2;
        case Class::GaStuckClashes:      return 3;
        case Class::EngineNonzeroExit:   return 4;
        case Class::Timeout:             return 5;
        case Class::EngineKilledSignal:  return 6;
        case Class::EngineCrashMidrun:   return 7;
        case Class::ResultCountMismatch: return 8;
        case Class::NoTargetResults:     return 9;
        case Class::EngineNotLaunched:   return 10;
        case Class::EngineStartupAbort:  return 11;
    }
    return 0;
}

// The runtime facts the classifier reads, lifted out of DockingResult so the
// classifier stays independent of DatasetRunner.h and remains unit-testable.
struct TargetFacts {
    std::string pdb_id;
    int    docking_exit_code{-1};
    int    num_poses{0};
    bool   docking_completed{false};
    bool   stuck{false};
    double wall_time_s{0.0};
};

struct ChildEvidence {
    bool out_dir_exists{false};
    bool stderr_log_exists{false};
    bool engine_terminal_handler{false};   // engine died through its own catch
    std::string stderr_log_path;
    std::vector<std::string> stderr_tail;  // last few lines, sanitized
};

// Literals emitted by the engine's own top-level handlers in LIB/top.cpp:
//   3890  fprintf(stderr, "%sFlexAID Error:%s %s\n", ...)   FlexAIDException
//   3892  fprintf(stderr, "Fatal error: %s\n", e.what());   std::exception
// Both run immediately before main() returns, so their presence means the
// engine terminated through an exception rather than finishing. Matching an
// exact substring the engine itself prints is not prose-parsing: these are
// format strings in this repository, and a test pins them.
inline bool line_has_engine_terminal_handler(const std::string& line) {
    return line.find("Fatal error: ") != std::string::npos ||
           line.find("FlexAID Error:") != std::string::npos;
}

// Read the tail of the child's stderr.log. Never throws: every filesystem call
// uses the error_code overload, and an unreadable log degrades to "no
// evidence" rather than taking the parent down — which is the whole lesson of
// the incident this file is hardening against.
//
// The window is bounded at kTailBytes because a healthy run's stderr.log is
// megabytes. That bound is safe for the marker search: both handlers fire
// immediately before process exit, so a terminal marker is always within the
// last few lines. A "Fatal error: " further back would not be the terminal one.
inline ChildEvidence read_child_evidence(const std::string& out_dir,
                                         std::size_t tail_lines = 6) {
    constexpr std::streamoff kTailBytes = 8192;
    ChildEvidence ev;

    std::error_code ec;
    ev.out_dir_exists = fs::is_directory(fs::path(out_dir), ec) && !ec;

    const std::string log_path = out_dir + "/stderr.log";
    ev.stderr_log_path = log_path;

    ec.clear();
    const bool is_file = fs::is_regular_file(fs::path(log_path), ec) && !ec;
    if (!is_file) return ev;
    ev.stderr_log_exists = true;

    std::ifstream in(log_path, std::ios::binary | std::ios::ate);
    if (!in) return ev;
    const std::streamoff size = static_cast<std::streamoff>(in.tellg());
    const std::streamoff start = (size > kTailBytes) ? (size - kTailBytes) : 0;
    in.seekg(start, std::ios::beg);

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        lines.push_back(line);
    }
    // A truncated first line is an artifact of the byte window, not content.
    if (start > 0 && !lines.empty()) lines.erase(lines.begin());

    for (const auto& l : lines) {
        if (line_has_engine_terminal_handler(l)) { ev.engine_terminal_handler = true; break; }
    }

    const std::size_t keep = (lines.size() > tail_lines) ? tail_lines : lines.size();
    for (std::size_t i = lines.size() - keep; i < lines.size(); ++i) {
        std::string s = lines[i];
        if (s.size() > 300) s = s.substr(0, 300) + " …";
        // Strip control bytes (the engine colourises some handlers) so the tail
        // is safe to paste into a CSV cell or a JSON string.
        std::string clean;
        clean.reserve(s.size());
        for (unsigned char ch : s) {
            if (ch == '\t') { clean.push_back(' '); continue; }
            if (ch >= 0x20 || ch >= 0x80) clean.push_back(static_cast<char>(ch));
        }
        ev.stderr_tail.push_back(clean);
    }
    return ev;
}

// Order is the contract. Each test is reached only when every test above it
// has been ruled out, so later branches may assume the earlier conditions are
// false — that is what keeps the classes disjoint.
inline Class classify_target(const TargetFacts& f,
                             const ChildEvidence& ev,
                             int per_job_timeout_s) {
    // 1. Healthy: exactly the conjunction benchmark_runtime_exit_code requires.
    if (f.docking_completed && f.docking_exit_code == 0 &&
        f.num_poses > 0 && !f.stuck) {
        return Class::None;
    }
    // 2. No child ran. DatasetRunner launches via `sh -c "... 2>DIR/stderr.log"`,
    //    so the shell creates stderr.log before the engine gets control: its
    //    absence means the fork/exec never happened for this target.
    if (!ev.out_dir_exists || !ev.stderr_log_exists) {
        return Class::EngineNotLaunched;
    }
    // 3. The engine terminated through its own top-level handler. Whether that
    //    counts as "startup" is decided by output, not by timing: zero poses
    //    means it died before producing any, which is the incident signature.
    if (ev.engine_terminal_handler) {
        return (f.num_poses > 0) ? Class::EngineCrashMidrun
                                 : Class::EngineStartupAbort;
    }
    // 4. SubprocessGuard::wait_with_timeout (DatasetRunner.cpp:451) returns
    //    WIFEXITED(status) ? WEXITSTATUS(status) : -1, and -1 again on timeout
    //    and on an unwaitable child. A negative code therefore means "did not
    //    exit normally" and NOTHING finer. Wall time against the configured
    //    budget is the only evidence that separates the timeout case.
    if (f.docking_exit_code < 0) {
        const bool hit_budget =
            per_job_timeout_s > 0 &&
            f.wall_time_s >= 0.95 * static_cast<double>(per_job_timeout_s);
        return hit_budget ? Class::Timeout : Class::EngineKilledSignal;
    }
    // 5. Clean exit, nonzero status, no handler line: the engine chose to fail
    //    and said so somewhere other than its terminal handler.
    if (f.docking_exit_code > 0) {
        return Class::EngineNonzeroExit;
    }
    // 6-8. Exit 0 from here down.
    if (f.stuck)           return Class::GaStuckClashes;
    if (f.num_poses <= 0)  return Class::ZeroPosesWritten;
    // Poses, clean exit, not stuck, yet not completed: the only remaining
    // producer is the skip-cache branch replaying a cached row whose
    // docking_completed cell was false (DatasetRunner.cpp:8245).
    return Class::CachedRowIncomplete;
}

struct TargetDiagnosis {
    TargetFacts   facts;
    Class         cls{Class::None};
    std::string   child_stderr_log;
    std::vector<std::string> stderr_tail;
};

struct RunDiagnosis {
    std::string dataset_name;
    Class       dominant{Class::None};
    int         process_exit_code{0};
    int         n_targets{0};
    std::vector<TargetDiagnosis> targets;
    std::vector<std::string>     startup_abort_targets;
    // Populated from the first startup abort seen, which is the one worth
    // quoting: in the 2026-09-20 pattern every later casualty carries the same
    // environmental message.
    std::string first_abort_pdb_id;
    std::string first_abort_log_path;
    std::vector<std::string> first_abort_stderr_tail;
};

inline std::string safe_dataset_name(std::string name) {
    // Mirrors DatasetRunner::write_report so the sidecar lands beside
    // <safe_name>_results.csv rather than in a name of its own invention.
    std::replace(name.begin(), name.end(), ' ', '_');
    std::replace(name.begin(), name.end(), '-', '_');
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name;
}

inline std::string csv_escape(const std::string& v) {
    bool needs = v.find_first_of(",\"\n") != std::string::npos;
    if (!needs) return v;
    std::string out = "\"";
    for (char c : v) { if (c == '"') out += "\"\""; else out += c; }
    out += "\"";
    return out;
}

inline std::string json_escape(const std::string& v) {
    std::string out;
    out.reserve(v.size() + 8);
    for (unsigned char c : v) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); out += b; }
                else out += static_cast<char>(c);
        }
    }
    return out;
}

// ── Sidecar schema ───────────────────────────────────────────────────────
// APPEND-ONLY. New columns go at the END, header and row changed in the same
// edit, and schema_version increments. Never reorder, never rename, never
// insert: a column added in the middle silently shifts every later named field
// in every downstream reader, which this project has already been bitten by.
// tests/test_startup_diagnosis_contract.py pins this prefix.
inline constexpr int kDiagnosisSchemaVersion = 1;

inline const char* diagnosis_csv_header() {
    return "pdb_id,failure_reason,failure_exit_code,docking_exit_code,num_poses,"
           "docking_completed,stuck,wall_time_s,child_stderr_log,schema_version";
}

inline std::string diagnosis_csv_row(const TargetDiagnosis& t) {
    std::ostringstream os;
    os << csv_escape(t.facts.pdb_id) << ','
       << to_string(t.cls) << ','
       << exit_code_for(t.cls) << ','
       << t.facts.docking_exit_code << ','
       << t.facts.num_poses << ','
       << (t.facts.docking_completed ? 1 : 0) << ','
       << (t.facts.stuck ? 1 : 0) << ','
       << t.facts.wall_time_s << ','
       << csv_escape(t.child_stderr_log) << ','
       << kDiagnosisSchemaVersion;
    return os.str();
}

// Atomic publish: write beside the destination, then rename. A reader that
// opens the path mid-write sees either the previous file or the new one, never
// a half-written header. Returns false rather than throwing.
inline bool write_text_atomic(const std::string& path, const std::string& body) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream ofs(tmp, std::ios::binary | std::ios::trunc);
        if (!ofs) return false;
        ofs << body;
        ofs.flush();
        if (!ofs) return false;
    }
    std::error_code ec;
    fs::rename(fs::path(tmp), fs::path(path), ec);
    if (ec) { std::error_code rm; fs::remove(fs::path(tmp), rm); return false; }
    return true;
}

inline std::string render_diagnosis_csv(const RunDiagnosis& run) {
    std::ostringstream os;
    os << diagnosis_csv_header() << '\n';
    for (const auto& t : run.targets) os << diagnosis_csv_row(t) << '\n';
    return os.str();
}

inline std::string render_diagnosis_json(const RunDiagnosis& run) {
    // Counts over the closed vocabulary, emitted only for classes present.
    std::vector<std::pair<std::string, int>> counts;
    for (const auto& t : run.targets) {
        const std::string key = to_string(t.cls);
        bool found = false;
        for (auto& kv : counts) if (kv.first == key) { kv.second++; found = true; break; }
        if (!found) counts.emplace_back(key, 1);
    }

    std::ostringstream os;
    os << "{\n";
    os << "  \"schema_version\": " << kDiagnosisSchemaVersion << ",\n";
    os << "  \"generated_by\": \"LIB/benchmark_datasets.cpp\",\n";
    os << "  \"dataset\": \"" << json_escape(run.dataset_name) << "\",\n";
    os << "  \"run_failure_reason\": \"" << to_string(run.dominant) << "\",\n";
    os << "  \"process_exit_code\": " << run.process_exit_code << ",\n";
    os << "  \"n_targets\": " << run.n_targets << ",\n";
    os << "  \"counts\": {";
    for (std::size_t i = 0; i < counts.size(); ++i) {
        os << (i ? ", " : " ") << "\"" << json_escape(counts[i].first) << "\": " << counts[i].second;
    }
    os << (counts.empty() ? "}" : " }") << ",\n";
    os << "  \"startup_abort_targets\": [";
    for (std::size_t i = 0; i < run.startup_abort_targets.size(); ++i) {
        os << (i ? ", " : "") << "\"" << json_escape(run.startup_abort_targets[i]) << "\"";
    }
    os << "],\n";
    os << "  \"first_startup_abort\": ";
    if (run.first_abort_pdb_id.empty()) {
        os << "null\n";
    } else {
        os << "{\n";
        os << "    \"pdb_id\": \"" << json_escape(run.first_abort_pdb_id) << "\",\n";
        os << "    \"child_stderr_log\": \"" << json_escape(run.first_abort_log_path) << "\",\n";
        os << "    \"stderr_tail\": [";
        for (std::size_t i = 0; i < run.first_abort_stderr_tail.size(); ++i) {
            os << (i ? ", " : "") << "\"" << json_escape(run.first_abort_stderr_tail[i]) << "\"";
        }
        os << "]\n  }\n";
    }
    os << "}\n";
    return os.str();
}

// The human message. It keeps the historical "Incomplete docking" wording so a
// log grep for it still hits, but the class name and the child log path are on
// the same line, and a startup abort quotes the child's own last words. The
// nine-cell incident would have been over in one minute with this output.
inline void print_human_summary(const RunDiagnosis& run, std::ostream& os) {
    if (run.dominant == Class::None) return;

    os << "ERROR: Incomplete docking [failure_reason=" << to_string(run.dominant)
       << " exit_code=" << run.process_exit_code << "]\n";

    if (run.dominant == Class::EngineStartupAbort ||
        run.dominant == Class::EngineCrashMidrun) {
        os << "  The ENGINE terminated through its own top-level handler. This is "
              "not a docking\n"
              "  result — the targets below were never given a chance to dock. "
              "Check the\n"
              "  environment (TMPDIR, data dir, input paths) before re-queuing "
              "anything.\n";
    }
    if (!run.startup_abort_targets.empty()) {
        os << "  startup-aborted targets (" << run.startup_abort_targets.size()
           << "):";
        std::size_t shown = 0;
        for (const auto& id : run.startup_abort_targets) {
            if (shown++ == 12) { os << " …"; break; }
            os << ' ' << id;
        }
        os << "\n";
    }
    if (!run.first_abort_stderr_tail.empty()) {
        os << "  child stderr — " << run.first_abort_log_path << "\n";
        for (const auto& l : run.first_abort_stderr_tail) os << "    | " << l << "\n";
    }
    os << "  per-target classes: see <output_dir>/"
       << safe_dataset_name(run.dataset_name)
       << "_failure_diagnosis.csv (column failure_reason)\n";
}

}  // namespace startup_diag

#ifndef FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY

static void print_usage(const char* progname) {
    printf("FlexAIDdS Benchmark Dataset Runner\n\n");
    printf("Usage:\n");
    printf("  %s --benchmark <dataset> [options]\n\n", progname);
    printf("Datasets:\n");
    printf("  astex            Astex Diverse 85 (Hartshorn et al. 2007)\n");
    printf("  astex_nonnative  Astex Non-Native 1112 (Verdonk et al. 2008)\n");
    printf("  hap2             HAP2 59 targets (Gaudreault & Najmanovich 2015)\n");
    printf("  casf2016         CASF-2016 285 complexes (PDBbind core set)\n");
    printf("  posebusters      PoseBusters 308 (Buttenschoen et al. 2024)\n");
    printf("  dude             DUD-E 102 targets (Mysinger et al. 2012)\n");
    printf("  bindingdb_itc    BindingDB ITC thermodynamic data\n");
    printf("  sampl6           SAMPL6 host-guest 27 (ITC: dG, dH, TdS)\n");
    printf("  sampl7           SAMPL7 host-guest ~30 (ITC: dG, dH, TdS)\n");
    printf("  pdbbind          PDBbind Refined 5316 (v2020)\n");
    printf("  all              Run all standard benchmarks\n");
    printf("  doi:<DOI>              Parse PDB codes from a DOI\n");
    printf("  pdb_list:<file>        Load PDB codes from a text file\n");
    printf("  crossdock_json:<file>  Cross-docking pairs from a JSON spec file\n\n");
    printf("Options:\n");
    printf("  --output <dir>        Output directory (default: benchmark_results/)\n");
    printf("  --threads <N>         Concurrent FlexAIDdS workers (default: 1)\n");
    printf("  --omp-threads <N>     OMP threads per worker (default: auto = hw_cores/workers)\n");
    printf("                        Rule: workers × omp-threads ≤ physical P-cores\n");
    printf("                        M3 Pro optimal: --threads 1 --omp-threads 6\n");
    printf("  --gpu <backend>       Enable GPU (cuda or metal)\n");
    printf("  --cache <dir>         Cache directory (default: ~/.flexaidds/benchmarks/)\n");
    printf("  --force               Re-run even if results already exist\n");
    printf("  --prepare-only        Download and prepare only (no docking)\n");
    printf("  --list-codes          List PDB codes for a dataset and exit\n");
    printf("  --only-codes <list>   Restrict dataset run to comma/space-separated PDB codes, or a file\n");
    printf("  --ga-generations <N>  GA generations (default: 500)\n");
    printf("  --ga-population <N>   GA population size (default: 1000)\n");
    printf("  --grid-spacing <F>    Grid spacing in Å (default: 0.375; use 0.5 for coarse pass)\n");
    printf("  --job-timeout-seconds <N>  Per-complex timeout in s (default: 3600)\n");
    printf("  --fleet               Enable immutable Fleet chunk-result output\n");
    printf("  --campaign-id <ID>    Fleet campaign identifier (required with --fleet)\n");
    printf("  --chunk-id <ID>       Fleet chunk identifier (required with --fleet)\n");
    printf("  --attempt-id <ID>     Fleet attempt identifier (required with --fleet)\n");
    printf("  --worker-id <ID>      Fleet worker identifier (required with --fleet)\n");
    printf("  --manifest-sha256 <H> SHA-256 of immutable campaign manifest\n");
    printf("  --runner-sha256 <H>   SHA-256 of this benchmark runner\n");
    printf("  --engine-sha256 <H>   SHA-256 of the pinned FlexAIDdS engine\n");
    printf("  --output-json <path>  Immutable Fleet chunk JSON destination\n");
    printf("  --mode <mode>         Benchmark protocol (Layer 1):\n");
    printf("                        oracle-ceiling  seed_elitism=ON,  blinding=OFF (ceiling)\n");
    printf("                        defined-cleft-redock  seed_elitism=OFF, blinding=ON, known cleft/site\n");
    printf("                        autonomous      seed_elitism=OFF, blinding=ON  (thesis number)\n");
    printf("                        (default: unset — reads FLEXAIDDS_SEED_ELITISM env var)\n");
    printf("  -h, --help            Show this help\n\n");
    printf("Thread sizing guide (M3 Pro, 6 P-cores):\n");
    printf("  --threads 1 --omp-threads 6   → 6 min/complex, optimal throughput\n");
    printf("  --threads 2 --omp-threads 3   → 10 min/complex, 2× dataset parallelism\n");
    printf("  --threads 2                   → AUTO: 6/2=3 OMP threads/worker (same as above)\n");
    printf("  --threads 2 (no --omp-threads, OMP_NUM_THREADS=6) → BUG: 12 threads on 6 cores\n\n");
    printf("Examples:\n");
    printf("  %s --benchmark astex --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark astex --prepare-only\n", progname);
    printf("  %s --benchmark casf2016 --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark all --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark doi:10.1021/acs.jcim.3c00817\n", progname);
    printf("  %s --benchmark astex --list-codes\n", progname);
}

static void print_publication_table(const dataset::BenchmarkReport& report) {
    // Print a publication-ready summary table matching manuscript format
    printf("\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("  FlexAIDdS Benchmark: %s\n", report.dataset_name.c_str());
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  ┌─────────────────────────────┬────────────────────┐\n");
    printf("  │ Metric                      │ Value              │\n");
    printf("  ├─────────────────────────────┼────────────────────┤\n");
    printf("  │ Total systems               │ %18d │\n", report.total_systems);
    printf("  │ Successful (RMSD <= 2.0 A)  │ %18d │\n", report.successful);
    printf("  │ Success rate                │ %17.1f%% │\n", report.success_rate * 100.0);
    printf("  │ Valid RMSDs                 │ %18d │\n", report.valid_rmsd_count);
    if (report.valid_rmsd_count > 0 && std::isfinite(report.mean_rmsd) && std::isfinite(report.median_rmsd)) {
        printf("  │ Mean RMSD (Å)               │ %18.2f │\n", report.mean_rmsd);
        printf("  │ Median RMSD (Å)             │ %18.2f │\n", report.median_rmsd);
    } else {
        printf("  │ Mean RMSD (Å)               │ %18s │\n", "NA");
        printf("  │ Median RMSD (Å)             │ %18s │\n", "NA");
    }
    printf("  │ Affinity pairs              │ %18d │\n", report.affinity_pairs);
    if (report.affinity_pairs >= 3 &&
        std::isfinite(report.pearson_r) &&
        std::isfinite(report.spearman_rho) &&
        std::isfinite(report.kendall_tau)) {
        printf("  │ Pearson r                   │ %18.3f │\n", report.pearson_r);
        printf("  │ Spearman ρ                  │ %18.3f │\n", report.spearman_rho);
        printf("  │ Kendall τ                   │ %18.3f │\n", report.kendall_tau);
    } else {
        printf("  │ Pearson r                   │ %18s │\n", "NA");
        printf("  │ Spearman ρ                  │ %18s │\n", "NA");
        printf("  │ Kendall τ                   │ %18s │\n", "NA");
    }
    printf("  └─────────────────────────────┴────────────────────┘\n");
    printf("\n");
}

static void print_itc_table(const dataset::BenchmarkReport& report,
                             const std::vector<dataset::DatasetEntry>& entries) {
    // Print ITC-specific thermodynamic comparison table
    bool has_itc = false;
    for (const auto& e : entries) {
        if (e.has_enthalpy()) { has_itc = true; break; }
    }
    if (!has_itc) return;

    // Claim firewall. The `predicted_*` side of every pair below comes from the
    // uncalibrated contact-function optimizer (EnergyDomain =
    // ContactFunctionArbitraryUnits, EnsembleMeasure = OptimizerSamples, no
    // sha256 receipt ⇒ ProxyOnly per LIB/statmech.h). Correlating a proxy
    // against experiment is a legitimate diagnostic — the correlation
    // coefficients are dimensionless and unaffected by the missing calibration
    // — but the proxy columns themselves must not be presented as physical
    // thermodynamics. Only the labels/units change here; the statistics are
    // computed from exactly the same inputs as before.
    printf("\n");
    printf("  ITC correlation diagnostic — CF-proxy predictions vs experiment\n");
    printf("  claim_validity: proxy_only. Predicted values are uncalibrated\n");
    printf("  contact-function diagnostics in arbitrary units, not kcal/mol; only\n");
    printf("  the dimensionless correlations below are interpretable. No ΔG, ΔH,\n");
    printf("  TΔS, Kd, Ki or affinity value is claimed for the predicted side.\n");
    printf("  ─────────────────────────────────────────────────────────────\n");

    // Collect ITC pairs
    std::vector<double> exp_dG, pred_dG;
    std::vector<double> exp_dH, pred_dH;
    std::vector<double> exp_TdS, pred_TdS;

    for (size_t i = 0; i < entries.size() && i < report.results.size(); ++i) {
        const auto& entry = entries[i];
        const auto& result = report.results[i];

        if (entry.has_affinity() && result.predicted_dG != 0.0f) {
            exp_dG.push_back(-entry.experimental_affinity * 1.3636);
            pred_dG.push_back(result.predicted_dG);
        }
        if (entry.has_enthalpy() && result.predicted_dH != 0.0f) {
            exp_dH.push_back(entry.experimental_dH);
            pred_dH.push_back(result.predicted_dH);
        }
        if (entry.has_entropy() && result.predicted_TdS != 0.0f) {
            exp_TdS.push_back(entry.experimental_TdS);
            pred_TdS.push_back(result.predicted_TdS);
        }
    }

    // Row labels name the PROXY channel that was correlated against the
    // experimental ITC column, e.g. "ΔG-like proxy" = predicted_dG (CF units)
    // vs experimental ΔG (kcal/mol). Cell width is 16 columns; keep the label
    // ≤ 14 display columns so the box stays aligned.
    printf("  ┌────────────────┬──────────┬──────────┬──────────┐\n");
    printf("  │ Proxy channel  │ Pearson  │ Spearman │ Kendall  │\n");
    printf("  ├────────────────┼──────────┼──────────┼──────────┤\n");

    if (exp_dG.size() >= 3) {
        printf("  │ ΔG-like proxy  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dG, exp_dG),
               dataset::compute_spearman_rho(pred_dG, exp_dG),
               dataset::compute_kendall_tau(pred_dG, exp_dG));
    }
    if (exp_dH.size() >= 3) {
        printf("  │ ΔH-like proxy  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dH, exp_dH),
               dataset::compute_spearman_rho(pred_dH, exp_dH),
               dataset::compute_kendall_tau(pred_dH, exp_dH));
    }
    if (exp_TdS.size() >= 3) {
        printf("  │ TΔS-like proxy │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_TdS, exp_TdS),
               dataset::compute_spearman_rho(pred_TdS, exp_TdS),
               dataset::compute_kendall_tau(pred_TdS, exp_TdS));
    }

    printf("  └────────────────┴──────────┴──────────┴──────────┘\n");
    printf("  Correlation only — proxy columns are in arbitrary CF units.\n");
    printf("\n");
}

static void list_pdb_codes(dataset::BenchmarkSet set) {
    std::vector<std::string> codes;
    switch (set) {
        case dataset::BenchmarkSet::ASTEX_DIVERSE:
            codes = dataset::DatasetRunner::astex_diverse_codes();
            break;
        case dataset::BenchmarkSet::CASF_2016:
            codes = dataset::DatasetRunner::casf2016_codes();
            break;
        case dataset::BenchmarkSet::DUD_E:
            codes = dataset::DatasetRunner::dude_targets();
            break;
        case dataset::BenchmarkSet::HAP2:
            codes = dataset::DatasetRunner::hap2_codes();
            break;
        default:
            printf("No hardcoded PDB list for this dataset. Use --prepare-only to fetch.\n");
            return;
    }

    printf("%s — %zu entries:\n", dataset::benchmark_set_name(set).c_str(), codes.size());
    int col = 0;
    for (const auto& code : codes) {
        printf("%-6s", code.c_str());
        if (++col % 12 == 0) printf("\n");
    }
    if (col % 12 != 0) printf("\n");
}

static std::string uppercase_code(std::string code) {
    std::transform(code.begin(), code.end(), code.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return code;
}

static std::vector<std::string> parse_only_codes(const std::string& spec) {
    std::vector<std::string> codes;
    if (spec.empty()) return codes;

    std::stringstream input;
    if (fs::exists(spec)) {
        std::ifstream ifs(spec);
        if (!ifs) {
            std::cerr << "ERROR: cannot open --only-codes file: " << spec << "\n";
            return codes;
        }
        input << ifs.rdbuf();
    } else {
        input << spec;
    }

    std::string text = input.str();
    for (char& ch : text) {
        if (ch == ',' || ch == ';' || ch == '\n' || ch == '\r' || ch == '\t') ch = ' ';
    }

    std::istringstream iss(text);
    std::string code;
    while (iss >> code) {
        auto hash = code.find('#');
        if (hash == 0) {
            std::string ignored;
            std::getline(iss, ignored);
            continue;
        }
        if (hash != std::string::npos) code.erase(hash);
        code = uppercase_code(code);
        if (code.size() == 4) codes.push_back(code);
    }
    return codes;
}

static void filter_entries_by_code(std::vector<dataset::DatasetEntry>& entries,
                                   const std::vector<std::string>& only_codes) {
    if (only_codes.empty()) return;

    std::set<std::string> wanted;
    for (const auto& code : only_codes) wanted.insert(uppercase_code(code));

    std::vector<dataset::DatasetEntry> filtered;
    filtered.reserve(entries.size());
    for (auto& entry : entries) {
        if (wanted.count(uppercase_code(entry.pdb_id))) {
            filtered.push_back(std::move(entry));
        }
    }

    std::cout << "  --only-codes selected " << filtered.size()
              << " / " << wanted.size() << " requested entries\n";
    entries = std::move(filtered);
}

static dataset::BenchmarkReport run_single_benchmark(const std::string& name,
                                  dataset::DatasetRunner& runner,
                                  const dataset::DockingConfig& config,
                                  bool prepare_only,
                                  bool list_codes_only,
                                  const std::vector<std::string>& only_codes) {
    using BS = dataset::BenchmarkSet;

    // Check for special prefixes: doi: and pdb_list:
    if (name.substr(0, 4) == "doi:") {
        std::string doi = name.substr(4);
        auto entries = runner.prepare_from_doi(doi);
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return {};
    }
    if (name.substr(0, 9) == "pdb_list:") {
        std::string file_path = name.substr(9);
        auto entries = runner.prepare_from_pdb_list(file_path);
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return {};
    }
    if (name.substr(0, 15) == "crossdock_json:") {
        // Cross-docking from a JSON pairs file.
        // JSON format (benchmark_crossdock_85.json):
        //   { "pairs": [ { "receptor_id": "1G9V", "ligand_id": "1GM8",
        //                   "receptor_pdb": "/abs/path/1G9V_apo.pdb",
        //                   "ligand_sdf":   "/abs/path/1GM8_ligand.sdf",
        //                   "oracle_site_pdb": "/abs/path/1G9V_binding_site.pdb",
        //                   "cleft_sphere_file": "/abs/path/1G9V_sph_1.pdb" }, ... ] }
        std::string json_file = name.substr(15);
        // Expand leading ~
        if (!json_file.empty() && json_file[0] == '~') {
            if (const char* home = std::getenv("HOME"))
                json_file = std::string(home) + json_file.substr(1);
        }
        std::ifstream ifs(json_file);
        if (!ifs) {
            fprintf(stderr, "ERROR: cannot open crossdock_json file: %s\n", json_file.c_str());
            return {};
        }
        std::string content((std::istreambuf_iterator<char>(ifs)),
                             std::istreambuf_iterator<char>());
        const fs::path json_base =
            fs::path(flexaids::fs_safe::absolute_or(json_file)).parent_path();

        // Minimal JSON string-field extractor (no external JSON lib needed).
        auto extract_str = [&](const std::string& obj, const std::string& key) -> std::string {
            std::string needle = "\"" + key + "\": \"";
            auto pos = obj.find(needle);
            if (pos == std::string::npos) return "";
            pos += needle.size();
            auto end = obj.find('"', pos);
            return (end != std::string::npos) ? obj.substr(pos, end - pos) : "";
        };
        auto resolve_input_path = [&](const std::string& value) -> std::string {
            if (value.empty()) return value;
            fs::path path(value);
            if (path.is_relative()) path = json_base / path;
            std::error_code ec;
            fs::path normalized = fs::weakly_canonical(path, ec);
            return ec ? path.lexically_normal().string() : normalized.string();
        };

        // Walk JSON locating each pair object by scanning for "receptor_id" keys.
        std::vector<dataset::DatasetEntry> entries;
        std::size_t scan = 0;
        while (true) {
            auto kpos = content.find("\"receptor_id\"", scan);
            if (kpos == std::string::npos) break;
            auto obj_start = content.rfind('{', kpos);
            auto obj_end   = content.find('}', kpos);
            if (obj_start == std::string::npos || obj_end == std::string::npos) break;
            std::string obj = content.substr(obj_start, obj_end - obj_start + 1);

            dataset::DatasetEntry entry;
            entry.pdb_id            = extract_str(obj, "receptor_id");
            entry.receptor_path     = resolve_input_path(extract_str(obj, "receptor_pdb"));
            entry.ligand_path       = resolve_input_path(extract_str(obj, "ligand_sdf"));
            entry.rmsd_reference_path = resolve_input_path(extract_str(obj, "rmsd_ref_sdf"));
            entry.binding_site_path = resolve_input_path(extract_str(obj, "oracle_site_pdb"));
            entry.cleft_sphere_path = resolve_input_path(extract_str(obj, "cleft_sphere_file"));
            entry.source            = "astex_crossdock_85";

            if (!entry.pdb_id.empty() &&
                !entry.receptor_path.empty() &&
                !entry.ligand_path.empty()) {
                entries.push_back(std::move(entry));
            }
            scan = obj_end + 1;
        }

        printf("  crossdock_json: %s — loaded %zu pairs\n",
               json_file.c_str(), entries.size());
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
            return report;
        }
        return {};
    }

    auto bs = dataset::parse_benchmark_set(name);
    if (!bs.has_value()) {
        fprintf(stderr, "ERROR: Unknown benchmark: '%s'\n", name.c_str());
        fprintf(stderr, "Use --help for available datasets.\n");
        return {};
    }

    if (list_codes_only) {
        list_pdb_codes(*bs);
        return {};
    }

    auto entries = runner.prepare(*bs);
    filter_entries_by_code(entries, only_codes);
    printf("  → %zu entries prepared\n", entries.size());
    // P1 diagnostic (additive, non-behavior): explicit T + progress for early diagnosis of best BindingMode runs
    printf("  [P1] Docking phase starting for best BindingMode search at temperature from config (exact 298/310 K fidelity required). Live progress + run_status sidecar will be emitted during run(entries).\n");

    if (prepare_only) {
        printf("  [prepare-only mode] Skipping docking.\n");
        return {};
    }

    if (!entries.empty()) {
        auto report = runner.run(entries, config);
        print_publication_table(report);
        print_itc_table(report, entries);
        runner.write_report(report, config.output_dir);
        return report;
    }
    return {};
}

// Bridge from BenchmarkReport to the classifier, and publisher of the sidecar.
// Deliberately NOT inside namespace startup_diag: keeping the classifier free
// of DatasetRunner.h is what lets the unit test compile it against libc++
// alone, and that property is easy to lose by accident.
//
// It publishes two files next to the report write_report() just wrote:
//   <output_dir>/<safe_dataset>_failure_diagnosis.csv   one row per target
//   <output_dir>/<safe_dataset>_failure_diagnosis.json  run-level + stderr tail
// Both are NEW files. No existing CSV is opened, re-read or rewritten: adding
// failure_reason to <dataset>_results.csv would mean a read-modify-write of a
// file a live campaign is appending to, and the column belongs in the writer
// (DatasetRunner::write_report) rather than in a post-pass here.
static startup_diag::RunDiagnosis diagnose_and_publish(
        const dataset::BenchmarkReport& report,
        const dataset::DockingConfig& config,
        const std::string& fallback_dataset_name,
        bool no_docking) {
    using namespace startup_diag;

    RunDiagnosis run;
    run.dataset_name = report.dataset_name.empty() ? fallback_dataset_name
                                                   : report.dataset_name;
    run.n_targets = static_cast<int>(report.results.size());
    if (no_docking) return run;   // --prepare-only / --list-codes: nothing ran

    // Report-level conditions first, in the order benchmark_runtime_exit_code
    // tests them (LIB/DatasetRunner.h:427). Both short-circuit: with no results
    // at all, or a count that disagrees with total_systems, no per-target class
    // would be trustworthy.
    if (report.total_systems <= 0) {
        run.dominant = Class::NoTargetResults;
        run.process_exit_code = exit_code_for(run.dominant);
        return run;
    }
    if (report.results.size() != static_cast<std::size_t>(report.total_systems)) {
        run.dominant = Class::ResultCountMismatch;
        run.process_exit_code = exit_code_for(run.dominant);
        return run;
    }

    for (const auto& r : report.results) {
        TargetFacts f;
        f.pdb_id            = r.pdb_id;
        f.docking_exit_code = r.docking_exit_code;
        f.num_poses         = r.num_poses;
        f.docking_completed = r.docking_completed;
        f.stuck             = r.stuck;
        f.wall_time_s       = r.wall_time_s;

        // Per-target directory layout is fixed at DatasetRunner.cpp:6330.
        const ChildEvidence ev =
            read_child_evidence(config.output_dir + "/" + r.pdb_id);

        TargetDiagnosis td;
        td.facts            = f;
        td.cls              = classify_target(f, ev, config.per_job_timeout_s);
        td.child_stderr_log = ev.stderr_log_path;
        td.stderr_tail      = ev.stderr_tail;

        if (td.cls == Class::EngineStartupAbort) {
            run.startup_abort_targets.push_back(f.pdb_id);
            if (run.first_abort_pdb_id.empty()) {
                run.first_abort_pdb_id      = f.pdb_id;
                run.first_abort_log_path    = ev.stderr_log_path;
                run.first_abort_stderr_tail = ev.stderr_tail;
            }
        }
        if (severity_rank(td.cls) > severity_rank(run.dominant)) run.dominant = td.cls;
        run.targets.push_back(std::move(td));
    }
    run.process_exit_code = exit_code_for(run.dominant);

    std::error_code ec;
    fs::create_directories(fs::path(config.output_dir), ec);
    const std::string base = config.output_dir + "/" +
                             safe_dataset_name(run.dataset_name) + "_failure_diagnosis";
    // A sidecar that cannot be written must not change the run's outcome.
    if (!write_text_atomic(base + ".csv", render_diagnosis_csv(run)))
        std::cerr << "  [WARN] could not publish " << base << ".csv\n";
    if (!write_text_atomic(base + ".json", render_diagnosis_json(run)))
        std::cerr << "  [WARN] could not publish " << base << ".json\n";
    return run;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Keep the benchmark runner and docking engine from the same build tree.
    // An explicit FLEXAIDDS_BINARY still wins, but absent that override the
    // sibling FlexAIDdS executable is authoritative. This prevents a
    // build_lto/benchmark_datasets invocation from silently using build/FlexAIDdS.
    if (std::getenv("FLEXAIDDS_BINARY") == nullptr) {
        std::error_code exe_ec;
        fs::path runner_path = fs::weakly_canonical(flexaids::fs_safe::absolute_or(argv[0]), exe_ec);
        if (!exe_ec) {
            const fs::path sibling = runner_path.parent_path() / "FlexAIDdS";
            if (fs::is_regular_file(sibling, exe_ec) && !exe_ec) {
#ifdef _WIN32
                _putenv_s("FLEXAIDDS_BINARY", sibling.string().c_str());
#else
                setenv("FLEXAIDDS_BINARY", sibling.string().c_str(), 0);
#endif
            }
        }
    }

    // Parse arguments
    std::string benchmark_name;
    std::string output_dir = "benchmark_results";
    std::string cache_dir;
    int threads = 1;
    int omp_threads = 0;
    int job_timeout_s = 3600;
    bool use_gpu = false;
    std::string gpu_backend = "cuda";
    bool prepare_only = false;
    bool list_codes_only = false;
    bool force_rerun = false;
    int ga_generations = 0;
    int ga_population = 0;
    double temperature = 0.0;
    double grid_spacing = 0.0;
    std::string clustering;
    std::vector<std::string> only_codes;
    // Fleet mode options
    bool fleet_mode = false;
    std::string campaign_id;
    std::string chunk_id;
    std::string attempt_id;
    std::string worker_id;
    std::string manifest_sha256;
    std::string runner_sha256;
    std::string engine_sha256;
    std::string output_json;
    // Layer 1: benchmark protocol mode
    std::string mode_str;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--benchmark" && i + 1 < argc) {
            benchmark_name = argv[++i];
            continue;
        }
        if (arg == "--output" && i + 1 < argc) {
            output_dir = argv[++i];
            continue;
        }
        if (arg == "--cache" && i + 1 < argc) {
            cache_dir = argv[++i];
            continue;
        }
        if (arg == "--threads" && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--omp-threads" && i + 1 < argc) {
            omp_threads = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--gpu" && i + 1 < argc) {
            use_gpu = true;
            gpu_backend = argv[++i];
            continue;
        }
        if (arg == "--prepare-only") {
            prepare_only = true;
            continue;
        }
        if (arg == "--list-codes") {
            list_codes_only = true;
            continue;
        }
        if (arg == "--only-codes" && i + 1 < argc) {
            only_codes = parse_only_codes(argv[++i]);
            continue;
        }
        if (arg == "--force") {
            force_rerun = true;
            continue;
        }
        if (arg == "--ga-generations" && i + 1 < argc) {
            ga_generations = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--ga-population" && i + 1 < argc) {
            ga_population = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--grid-spacing" && i + 1 < argc) {
            grid_spacing = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--temperature" && i + 1 < argc) {
            temperature = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--clustering" && i + 1 < argc) {
            clustering = argv[++i];
            continue;
        }
        if (arg == "--job-timeout-seconds" && i + 1 < argc) {
            job_timeout_s = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--fleet") {
            fleet_mode = true;
            continue;
        }
        if (arg == "--campaign-id" && i + 1 < argc) {
            campaign_id = argv[++i];
            continue;
        }
        if (arg == "--chunk-id" && i + 1 < argc) {
            chunk_id = argv[++i];
            continue;
        }
        if (arg == "--attempt-id" && i + 1 < argc) {
            attempt_id = argv[++i];
            continue;
        }
        if (arg == "--worker-id" && i + 1 < argc) {
            worker_id = argv[++i];
            continue;
        }
        if (arg == "--manifest-sha256" && i + 1 < argc) {
            manifest_sha256 = argv[++i];
            continue;
        }
        if (arg == "--runner-sha256" && i + 1 < argc) {
            runner_sha256 = argv[++i];
            continue;
        }
        if (arg == "--engine-sha256" && i + 1 < argc) {
            engine_sha256 = argv[++i];
            continue;
        }
        if (arg == "--output-json" && i + 1 < argc) {
            output_json = argv[++i];
            continue;
        }
        if (arg == "--mode" && i + 1 < argc) {
            mode_str = argv[++i];
            continue;
        }

        // Fallback: if first positional arg, treat as benchmark name
        if (benchmark_name.empty()) {
            benchmark_name = arg;
        }
    }

    if (benchmark_name.empty()) {
        fprintf(stderr, "ERROR: No benchmark specified. Use --benchmark <name>\n");
        print_usage(argv[0]);
        return 1;
    }

    if (fleet_mode) {
        const bool missing_metadata = campaign_id.empty() || chunk_id.empty() ||
            attempt_id.empty() || worker_id.empty() || manifest_sha256.empty() ||
            runner_sha256.empty() || engine_sha256.empty() || output_json.empty();
        if (missing_metadata) {
            fprintf(stderr, "ERROR: --fleet requires campaign/chunk/attempt/worker IDs, "
                    "manifest/runner/engine SHA-256 values, and --output-json\n");
            return 1;
        }
        if (benchmark_name == "all" || prepare_only || list_codes_only ||
            only_codes.empty() || mode_str.empty() || output_json == "-") {
            fprintf(stderr, "ERROR: Fleet mode requires one explicit benchmark, --only-codes, "
                    "an explicit --mode, and a file output; prepare/list/all modes are unsupported\n");
            return 1;
        }
    }

    // Create runner and config
    dataset::DatasetRunner runner(cache_dir);

    dataset::DockingConfig config;
    config.num_threads            = threads;
    config.omp_threads_per_worker = omp_threads;   // 0 → auto-detect in DatasetRunner
    config.use_gpu                = use_gpu;
    config.gpu_backend            = gpu_backend;
    config.output_dir             = output_dir;
    config.skip_completed         = !force_rerun;
    if (ga_generations > 0)       config.ga_generations    = ga_generations;
    if (ga_population  > 0)       config.ga_population     = ga_population;
    if (grid_spacing   > 0.0)     config.grid_spacing      = static_cast<float>(grid_spacing);
    if (temperature    > 0.0)     config.temperature       = static_cast<float>(temperature);
    if (job_timeout_s  > 0)       config.per_job_timeout_s = job_timeout_s;
    if (!clustering.empty())      config.clustering_algorithm = clustering;

    // Layer 1: explicit benchmark protocol mode
    if (mode_str == "oracle-ceiling") {
        config.mode = dataset::BenchmarkMode::ORACLE_CEILING;
    } else if (mode_str == "defined-cleft-redock" || mode_str == "cognate-redock") {
        config.mode = dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK;
    } else if (mode_str == "autonomous") {
        config.mode = dataset::BenchmarkMode::AUTONOMOUS;
    } else if (!mode_str.empty()) {
        fprintf(stderr, "ERROR: Unknown --mode '%s'. Use 'oracle-ceiling', 'defined-cleft-redock', or 'autonomous'\n",
                mode_str.c_str());
        return 1;
    }

    // Ablation hook: FLEXAIDDS_FORCE_RIGID re-pins legacy rigid-body docking
    // (DatasetRunner writes flexibility.intramolecular=false → engine builds a
    // 4-gene chromosome: translation + rotation only, no ligand torsional DoF).
    // Isolates whether a flexible-docking regression is eval-budget dilution
    // (more genes, same generations) vs the oracle-confinement changes.
    if (const char* fr = std::getenv("FLEXAIDDS_FORCE_RIGID")) {
        if (fr[0] && fr[0] != '0') {
            config.force_rigid = true;
            std::cout << "  FORCE_RIGID:  ON (intramolecular=false, num_genes=4)\n";
        }
    }

    // Override receptor rotamer prep gate (default true since v44).
    // Set FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0 to disable for cross-docking
    // benchmarks where the native ligand is absent and sidechain prep against
    // a ghost occupancy makes no physical sense.
    if (const char* rrp = std::getenv("FLEXAIDDS_RECEPTOR_ROTAMER_PREP")) {
        bool enable = (rrp[0] && rrp[0] != '0' &&
                       std::string(rrp) != "false" && std::string(rrp) != "off");
        config.receptor_rotamer_prep = enable;
        std::cout << "  RECEPTOR_ROTAMER_PREP: " << (enable ? "ON" : "OFF")
                  << " (override via env)\n";
    }

    // Compute effective OMP threads for display (mirrors DatasetRunner logic)
    int effective_omp = config.omp_threads_per_worker;
    if (effective_omp <= 0) {
        const char* env_omp = std::getenv("OMP_NUM_THREADS");
        int base = (env_omp && std::atoi(env_omp) > 0)
            ? std::atoi(env_omp)
            : static_cast<int>(std::thread::hardware_concurrency());
        effective_omp = std::max(1, base / std::max(1, config.num_threads));
    }

    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  FlexAIDdS Benchmark Dataset Runner\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    std::cout << "  Cache:        " << runner.cache_dir() << "\n";
    std::cout << "  Output:       " << output_dir << "\n";
    std::cout << "  Workers:      " << threads << " concurrent FlexAIDdS process(es)\n";
    std::cout << "  OMP/worker:   " << effective_omp << " thread(s)"
              << (config.omp_threads_per_worker > 0 ? " (explicit)" : " (auto)") << "\n";
    std::cout << "  Total threads:" << (threads * effective_omp) << " across "
              << std::thread::hardware_concurrency() << " logical cores\n";
    std::cout << "  Skip done:    " << (config.skip_completed ? "yes (--force to override)" : "no") << "\n";
    if (use_gpu) {
        std::cout << "  GPU:          " << gpu_backend << "\n";
    }
    std::cout << "  GA:           pop=" << config.ga_population
              << "  gen=" << config.ga_generations
              << "  (" << (config.ga_population * config.ga_generations / 1000) << "k evals/complex)\n";
    std::cout << "  Temp:         " << config.temperature << " K\n";
    std::cout << "  Cluster:      " << config.clustering_algorithm << "\n";
    std::cout << "  Timeout/job:  " << config.per_job_timeout_s << " s\n";
    // Layer 1: mode
    {
        const char* mode_label =
            (config.mode == dataset::BenchmarkMode::ORACLE_CEILING)       ? "oracle-ceiling" :
            (config.mode == dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK) ? "defined-cleft-redock" :
            (config.mode == dataset::BenchmarkMode::AUTONOMOUS)           ? "autonomous" :
                                                                            "unset (env-var)";
        std::cout << "  Mode:         " << mode_label << "\n";
    }
    if (fleet_mode) {
        std::cout << "  Fleet:        enabled\n";
        std::cout << "  Campaign:     " << campaign_id << "\n";
        std::cout << "  Chunk/attempt:" << chunk_id << "/" << attempt_id << "\n";
        std::cout << "  JSON:         " << output_json << "\n";
    }
    std::cout << "\n";

    int runtime_exit_code = 0;
    // benchmark_runtime_exit_code() stays the sole authority on WHETHER the run
    // failed — this change renames failures, it does not reclassify any run as
    // passing. `worst` only supplies the name and the granular code.
    startup_diag::RunDiagnosis worst;
    // Handle "all" benchmark
    if (benchmark_name == "all") {
        std::vector<std::string> all_benchmarks = {
            "astex", "astex_nonnative", "hap2", "casf2016",
            "posebusters", "dude", "bindingdb_itc",
            "sampl6", "sampl7"
        };

        for (const auto& name : all_benchmarks) {
            std::cout << "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
            std::cout << "  Running: " << name << "\n";
            std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n";

            auto report = run_single_benchmark(name, runner, config, prepare_only, list_codes_only, only_codes);
            const bool no_docking = prepare_only || list_codes_only;
            const int legacy_rc =
                dataset::benchmark_runtime_exit_code(report, no_docking);
            const auto diag = diagnose_and_publish(report, config, name, no_docking);
            if (legacy_rc != 0) {
                // Worst ACROSS benchmarks is chosen by diagnostic urgency, not
                // by numeric code: 30 (result_count_mismatch) is a larger
                // integer than 20 (engine_startup_abort) and a far less
                // actionable finding, so std::max over the codes would bury the
                // one class that means "stop and fix the environment".
                if (startup_diag::severity_rank(diag.dominant) >
                    startup_diag::severity_rank(worst.dominant)) {
                    worst = diag;
                }
                runtime_exit_code = (worst.dominant != startup_diag::Class::None)
                    ? worst.process_exit_code
                    : std::max(runtime_exit_code, legacy_rc);
            }
        }

        // Print combined summary
        std::cout << "\n\n═══════════════════════════════════════════════════════════════\n";
        std::cout << "  All benchmark attempts finished. Results in: " << output_dir << "\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
    } else {
        const auto fleet_started = std::chrono::steady_clock::now();
        auto report = run_single_benchmark(benchmark_name, runner, config, prepare_only, list_codes_only, only_codes);

        const bool no_docking = prepare_only || list_codes_only;
        runtime_exit_code = dataset::benchmark_runtime_exit_code(report, no_docking);
        worst = diagnose_and_publish(report, config, benchmark_name, no_docking);
        // FAIL-CLOSED COMPOSITION. The legacy predicate decides pass/fail; the
        // classifier only supplies a name and a finer code. If the predicate
        // says "failed" and the classifier found no class, the generic 2 is
        // kept — an unnamed failure must never be promoted to success, and a
        // disagreement between the two is itself worth seeing in the message.
        if (runtime_exit_code != 0 && worst.dominant != startup_diag::Class::None)
            runtime_exit_code = worst.process_exit_code;
        if (fleet_mode) {
            const double duration_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - fleet_started).count();
            std::ostringstream command;
            for (int i = 0; i < argc; ++i) {
                if (i > 0) command << ' ';
                command << argv[i];
            }
            std::error_code path_error;
            fs::path runner_path = fs::weakly_canonical(flexaids::fs_safe::absolute_or(argv[0]), path_error);
            if (path_error) runner_path = flexaids::fs_safe::absolute_or(argv[0]);
            const char* engine_env = std::getenv("FLEXAIDDS_BINARY");
            fleet::ChunkMetadata metadata{
                campaign_id,
                chunk_id,
                attempt_id,
                worker_id,
                benchmark_name,
                command.str(),
                runner_path.string(),
                runner_sha256,
                engine_env ? engine_env : "",
                engine_sha256,
                manifest_sha256,
            };
            const std::string payload = fleet::FleetRunner::serialize_chunk_result(
                metadata, report, config, duration_s);
            std::string write_error;
            if (!fleet::FleetRunner::write_chunk_result_atomic(
                    output_json, payload, &write_error)) {
                std::cerr << "ERROR: Fleet result publication failed: " << write_error << "\n";
                return 3;
            }
            std::cout << "  [Fleet] immutable chunk result: " << output_json << "\n";
            if (report.total_systems == 0 || report.results.empty()) {
                std::cerr << "ERROR: Fleet chunk produced no target results\n";
                return 2;
            }
        }
    }

    // The literal "ERROR: Incomplete docking" is retained in both arms so an
    // existing log grep still matches; what changed is that the class name, the
    // granular exit code and the child log path now travel with it.
    if (runtime_exit_code != 0) {
        if (worst.dominant != startup_diag::Class::None) {
            startup_diag::print_human_summary(worst, std::cerr);
        } else {
            // benchmark_runtime_exit_code() failed the run while the classifier
            // named nothing. That is a gap in the vocabulary, not a clean run,
            // so it keeps the legacy code and says which state it was in.
            std::cerr << "ERROR: Incomplete docking [failure_reason=unclassified"
                         " exit_code=" << runtime_exit_code << "]; the runtime"
                         " predicate failed but no failure class matched."
                         " Inspect per-target runtime fields and child logs\n";
        }
    }
    return runtime_exit_code;
}

#endif  // FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY
