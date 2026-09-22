// =============================================================================
// test_startup_diagnosis.cpp — a startup abort must not look like a docking
// failure.
//
// WHAT THIS PROTECTS
// ------------------
// On 2026-09-20 nine cells of an 85-target Astex campaign died because $TMPDIR
// named a directory an agent harness had swept. The engine threw out of input
// parsing, LIB/top.cpp's catch(std::exception) printed
//
//     Fatal error: filesystem error: in temp_directory_path: path "/tmp/…" is
//     not a directory: Not a directory
//
// to the child's stderr.log and returned 1. The parent printed "ERROR:
// Incomplete docking" — the same sentence, and the same exit code, it prints
// when the engine runs perfectly and the GA simply finds no pose. Those nine
// targets looked like nine hard targets for six hours.
//
// The central assertion in this file is therefore a DIFFERENCE, not a value:
// a simulated startup abort and a genuine zero-pose docking result must not
// classify the same and must not exit the same. A test that only checked
// "startup abort is classified as engine_startup_abort" would still pass if
// every other class were also engine_startup_abort.
//
// The stderr fixture below is the verbatim first line of
// /Users/lp.more/flexaidds_results/tempdir_e2e/old_ghost/1G9V/stderr.log, the
// project's own reproduction of the incident, so the string being matched is
// one the engine has actually emitted rather than one invented here.
//
// WHAT IT DOES NOT COVER
// ----------------------
// Stated per test, at the end of each test body.
// =============================================================================

#define FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY 1
#include "../LIB/benchmark_datasets.cpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace {

namespace sd = startup_diag;

// Verbatim from the project's own ghost-TMPDIR reproduction.
constexpr const char* kIncidentStderr =
    "Fatal error: filesystem error: in temp_directory_path: "
    "path \"/tmp/ghost_old_13910\" is not a directory: Not a directory\n";

// A realistic tail from a run that reached the GA (new_ghost/1G9V/stderr.log).
constexpr const char* kHealthyStderr =
    "[DIRECT-IC] GPA topology atoms=90005,90004,90006 local=5,4,6\n"
    "[SMFREE] gen=0  beta_sel=0.003333  T=300.0  F=-45.039  <E>=-44.986\n"
    "[ELITE] GA-internal elitism active: protecting 1 lowest-CF individual(s)\n";

class StartupDiagnosis : public ::testing::Test {
protected:
    std::filesystem::path root_;

    void SetUp() override {
        // error_code overload deliberately: this suite is about not dying on a
        // temp path, so it must not use the throwing form itself.
        std::error_code ec;
        std::filesystem::path base = std::filesystem::temp_directory_path(ec);
        if (ec || base.empty()) base = "/tmp";
        root_ = base / ("flexaids_startup_diag_" +
                        std::to_string(::testing::UnitTest::GetInstance()
                                           ->random_seed()) +
                        "_" + std::to_string(reinterpret_cast<uintptr_t>(this)));
        std::filesystem::create_directories(root_, ec);
    }

    void TearDown() override {
        // Fixtures only ever live under this process's own scratch root.
        std::error_code ec;
        std::filesystem::remove_all(root_, ec);
    }

    // Create <root>/<pdb_id>/stderr.log with the given body.
    std::string make_target(const std::string& pdb_id, const char* stderr_body) {
        std::error_code ec;
        const std::filesystem::path dir = root_ / pdb_id;
        std::filesystem::create_directories(dir, ec);
        if (stderr_body) {
            std::ofstream ofs(dir / "stderr.log", std::ios::binary);
            ofs << stderr_body;
        }
        return dir.string();
    }
};

sd::TargetFacts facts(const std::string& id, int exit_code, int poses,
                      bool completed, bool stuck, double wall_s) {
    sd::TargetFacts f;
    f.pdb_id            = id;
    f.docking_exit_code = exit_code;
    f.num_poses         = poses;
    f.docking_completed = completed;
    f.stuck             = stuck;
    f.wall_time_s       = wall_s;
    return f;
}

// ── THE CENTRAL TEST ─────────────────────────────────────────────────────
TEST_F(StartupDiagnosis, StartupAbortIsNotAZeroPoseDockingResult) {
    // Case A: the 2026-09-20 signature. Engine's terminal handler fired,
    // nothing was written. Exit 1 is what LIB/top.cpp returns from
    // catch(std::exception).
    const std::string abort_dir = make_target("1G9V", kIncidentStderr);
    const auto abort_ev = sd::read_child_evidence(abort_dir);
    const auto abort_cls = sd::classify_target(
        facts("1G9V", /*exit=*/1, /*poses=*/0, /*completed=*/false,
              /*stuck=*/false, /*wall_s=*/0.4),
        abort_ev, /*per_job_timeout_s=*/3600);

    // Case B: a genuine docking failure. The engine ran to completion, exited
    // cleanly, and produced no pose. This is a RESULT, not a fault.
    const std::string zero_dir = make_target("1HWI", kHealthyStderr);
    const auto zero_ev = sd::read_child_evidence(zero_dir);
    const auto zero_cls = sd::classify_target(
        facts("1HWI", /*exit=*/0, /*poses=*/0, /*completed=*/false,
              /*stuck=*/false, /*wall_s=*/812.0),
        zero_ev, /*per_job_timeout_s=*/3600);

    EXPECT_EQ(abort_cls, sd::Class::EngineStartupAbort);
    EXPECT_EQ(zero_cls, sd::Class::ZeroPosesWritten);

    // The assertions that actually matter: distinguishable by name AND by the
    // integer a shell driver branches on.
    EXPECT_NE(abort_cls, zero_cls);
    EXPECT_STRNE(sd::to_string(abort_cls), sd::to_string(zero_cls));
    EXPECT_NE(sd::exit_code_for(abort_cls), sd::exit_code_for(zero_cls));

    // Both must still be failures. Naming a class must never turn a failed
    // target into a passing one.
    EXPECT_NE(0, sd::exit_code_for(abort_cls));
    EXPECT_NE(0, sd::exit_code_for(zero_cls));

    // NOT COVERED: that the ENGINE exits 1 on a ghost TMPDIR — that is
    // LIB/top.cpp's behaviour, asserted by tests/test_temp_dir.cpp, and this
    // test would pass unchanged if the engine's exit code changed.
}

TEST_F(StartupDiagnosis, HealthyTargetIsNotClassifiedAtAll) {
    const std::string dir = make_target("1G9V", kHealthyStderr);
    const auto cls = sd::classify_target(
        facts("1G9V", 0, 34, /*completed=*/true, false, 812.0),
        sd::read_child_evidence(dir), 3600);
    EXPECT_EQ(cls, sd::Class::None);
    EXPECT_EQ(0, sd::exit_code_for(cls));

    // NOT COVERED: pose quality. A completed target with a 9 Å RMSD is
    // Class::None here by design — RMSD lives on DockingResult::rmsd_fail_reason.
}

TEST_F(StartupDiagnosis, TerminalHandlerWithPosesIsMidrunNotStartup) {
    const std::string dir = make_target("1OF1", kIncidentStderr);
    const auto cls = sd::classify_target(
        facts("1OF1", 1, /*poses=*/12, false, false, 400.0),
        sd::read_child_evidence(dir), 3600);
    EXPECT_EQ(cls, sd::Class::EngineCrashMidrun);
    EXPECT_NE(sd::exit_code_for(cls),
              sd::exit_code_for(sd::Class::EngineStartupAbort));

    // NOT COVERED: whether the crash was actually mid-GA. "Poses exist" is the
    // discriminator, so a crash after the first pose is written but before the
    // GA starts would also land here.
}

TEST_F(StartupDiagnosis, MissingOutputDirMeansNoChildWasLaunched) {
    // No make_target call: the directory never existed.
    const auto ev = sd::read_child_evidence((root_ / "9XYZ").string());
    EXPECT_FALSE(ev.out_dir_exists);
    EXPECT_FALSE(ev.stderr_log_exists);
    const auto cls = sd::classify_target(facts("9XYZ", -1, 0, false, false, 0.0),
                                         ev, 3600);
    EXPECT_EQ(cls, sd::Class::EngineNotLaunched);

    // A directory with no stderr.log is the same finding: DatasetRunner's
    // `sh -c "... 2>DIR/stderr.log"` creates the log before the engine runs.
    std::error_code ec;
    std::filesystem::create_directories(root_ / "8XYZ", ec);
    const auto ev2 = sd::read_child_evidence((root_ / "8XYZ").string());
    EXPECT_TRUE(ev2.out_dir_exists);
    EXPECT_FALSE(ev2.stderr_log_exists);
    EXPECT_EQ(sd::Class::EngineNotLaunched,
              sd::classify_target(facts("8XYZ", -1, 0, false, false, 0.0), ev2, 3600));

    // NOT COVERED: WHY no child ran. Missing receptor, missing ligand and a
    // fork failure are indistinguishable here, which is exactly why the class
    // is not called input_missing.
}

TEST_F(StartupDiagnosis, NegativeExitSplitsOnWallTimeAgainstTheBudget) {
    // wait_with_timeout returns -1 for BOTH a signal death and a timeout
    // (DatasetRunner.cpp:451). Wall time against per_job_timeout_s is the only
    // available discriminator.
    const std::string dir = make_target("1P2Y", kHealthyStderr);
    const auto ev = sd::read_child_evidence(dir);

    EXPECT_EQ(sd::Class::Timeout,
              sd::classify_target(facts("1P2Y", -1, 0, false, false, 3600.0), ev, 3600));
    EXPECT_EQ(sd::Class::EngineKilledSignal,
              sd::classify_target(facts("1P2Y", -1, 0, false, false, 12.0), ev, 3600));
    EXPECT_NE(sd::exit_code_for(sd::Class::Timeout),
              sd::exit_code_for(sd::Class::EngineKilledSignal));

    // With no budget configured the split cannot be made and the honest answer
    // is the weaker claim.
    EXPECT_EQ(sd::Class::EngineKilledSignal,
              sd::classify_target(facts("1P2Y", -1, 0, false, false, 9e9), ev, 0));

    // NOT COVERED: a child SIGKILLed by the OOM killer one second before its
    // deadline classifies as timeout. Separating those needs the signal number,
    // which wait_with_timeout discards.
}

TEST_F(StartupDiagnosis, StuckAndCachedRowsAreTheirOwnClasses) {
    const std::string dir = make_target("1R1H", kHealthyStderr);
    const auto ev = sd::read_child_evidence(dir);

    EXPECT_EQ(sd::Class::GaStuckClashes,
              sd::classify_target(facts("1R1H", 0, 0, false, /*stuck=*/true, 900.0), ev, 3600));
    // Poses, clean exit, not stuck, yet not completed: only the skip-cache
    // branch produces this (DatasetRunner.cpp:8245).
    EXPECT_EQ(sd::Class::CachedRowIncomplete,
              sd::classify_target(facts("1R1H", 0, 20, /*completed=*/false, false, 0.0), ev, 3600));
    EXPECT_EQ(sd::Class::EngineNonzeroExit,
              sd::classify_target(facts("1R1H", 3, 0, false, false, 30.0), ev, 3600));

    // NOT COVERED: whether `stuck` is a correct GA diagnosis. This trusts
    // DockingResult::stuck; the clash-rate threshold is DatasetRunner's.
}

TEST_F(StartupDiagnosis, EveryClassHasADistinctExitCode) {
    const sd::Class all[] = {
        sd::Class::None, sd::Class::NoTargetResults, sd::Class::ResultCountMismatch,
        sd::Class::EngineNotLaunched, sd::Class::EngineStartupAbort,
        sd::Class::EngineCrashMidrun, sd::Class::Timeout,
        sd::Class::EngineKilledSignal, sd::Class::EngineNonzeroExit,
        sd::Class::GaStuckClashes, sd::Class::ZeroPosesWritten,
        sd::Class::CachedRowIncomplete,
    };
    for (const auto a : all) {
        for (const auto b : all) {
            if (a == b) continue;
            EXPECT_NE(sd::exit_code_for(a), sd::exit_code_for(b))
                << sd::to_string(a) << " and " << sd::to_string(b)
                << " share an exit code; a driver cannot branch";
            EXPECT_STRNE(sd::to_string(a), sd::to_string(b));
        }
    }
    // The four pre-existing codes must not be reassigned to a new class.
    for (const auto c : all) {
        if (c == sd::Class::None) continue;
        EXPECT_GT(sd::exit_code_for(c), 3)
            << sd::to_string(c) << " collides with a pre-existing exit code (1/2/3)";
    }
}

TEST_F(StartupDiagnosis, StartupAbortQuotesTheChildStderrInTheParentMessage) {
    // This is the property that would have ended the incident in one minute.
    const std::string dir = make_target("1G9V", kIncidentStderr);
    const auto ev = sd::read_child_evidence(dir);
    ASSERT_TRUE(ev.engine_terminal_handler);
    ASSERT_FALSE(ev.stderr_tail.empty());

    sd::RunDiagnosis run;
    run.dataset_name = "Astex Diverse";
    run.dominant = sd::Class::EngineStartupAbort;
    run.process_exit_code = sd::exit_code_for(run.dominant);
    run.n_targets = 1;
    run.startup_abort_targets.push_back("1G9V");
    run.first_abort_pdb_id = "1G9V";
    run.first_abort_log_path = ev.stderr_log_path;
    run.first_abort_stderr_tail = ev.stderr_tail;

    std::ostringstream os;
    sd::print_human_summary(run, os);
    const std::string msg = os.str();

    EXPECT_NE(std::string::npos, msg.find("engine_startup_abort"));
    EXPECT_NE(std::string::npos, msg.find("temp_directory_path"));
    EXPECT_NE(std::string::npos, msg.find("1G9V"));
    EXPECT_NE(std::string::npos, msg.find(ev.stderr_log_path));
    // Historical wording retained so an existing log grep still hits.
    EXPECT_NE(std::string::npos, msg.find("ERROR: Incomplete docking"));

    // A healthy run must print nothing at all.
    sd::RunDiagnosis clean;
    std::ostringstream quiet;
    sd::print_human_summary(clean, quiet);
    EXPECT_TRUE(quiet.str().empty());

    // NOT COVERED: that the parent's stderr actually reaches the operator's
    // log. The 2026-09-20 report was that stderr WAS captured and unread.
}

TEST_F(StartupDiagnosis, SidecarHeaderAndRowAgreeAndNewColumnsGoLast) {
    const std::string header = sd::diagnosis_csv_header();

    // The v1 prefix is frozen. A future column is APPENDED, so the header must
    // still START with this exact string — an inserted column shifts every
    // later named field in every downstream reader.
    EXPECT_EQ(0u, header.find("pdb_id,failure_reason,failure_exit_code,"
                              "docking_exit_code,num_poses,docking_completed,"
                              "stuck,wall_time_s,child_stderr_log"));
    EXPECT_EQ(0u, header.find("pdb_id,"))
        << "pdb_id must stay the first column; it is the join key";

    const auto count_commas = [](const std::string& s) {
        std::size_t n = 0, depth = 0;
        for (char c : s) {
            if (c == '"') depth ^= 1u;
            if (c == ',' && !depth) ++n;
        }
        return n;
    };

    sd::TargetDiagnosis td;
    td.facts = facts("1G9V", 1, 0, false, false, 0.4);
    td.cls = sd::Class::EngineStartupAbort;
    td.child_stderr_log = "/out/1G9V/stderr.log";
    const std::string row = sd::diagnosis_csv_row(td);

    EXPECT_EQ(count_commas(header), count_commas(row))
        << "header/row arity drift:\n  " << header << "\n  " << row;
    EXPECT_NE(std::string::npos, row.find("engine_startup_abort"));

    // A path containing a comma must not change the arity.
    td.child_stderr_log = "/out/odd,name/1G9V/stderr.log";
    EXPECT_EQ(count_commas(header), count_commas(sd::diagnosis_csv_row(td)));

    // NOT COVERED: the columns of <dataset>_results.csv and result.csv. Those
    // are written by DatasetRunner::write_report and are untouched by this
    // change; a guard for them belongs beside that writer.
}

TEST_F(StartupDiagnosis, JsonCarriesTheClassAndTheTail) {
    sd::RunDiagnosis run;
    run.dataset_name = "Astex Diverse";
    run.dominant = sd::Class::EngineStartupAbort;
    run.process_exit_code = sd::exit_code_for(run.dominant);
    run.n_targets = 2;
    run.startup_abort_targets = {"1G9V"};
    run.first_abort_pdb_id = "1G9V";
    run.first_abort_log_path = "/out/1G9V/stderr.log";
    run.first_abort_stderr_tail = {"Fatal error: filesystem error: in "
                                   "temp_directory_path: path \"/tmp/x\""};

    sd::TargetDiagnosis a;
    a.facts = facts("1G9V", 1, 0, false, false, 0.4);
    a.cls = sd::Class::EngineStartupAbort;
    sd::TargetDiagnosis b;
    b.facts = facts("1HWI", 0, 34, true, false, 800.0);
    b.cls = sd::Class::None;
    run.targets = {a, b};

    const std::string js = sd::render_diagnosis_json(run);
    EXPECT_NE(std::string::npos, js.find("\"run_failure_reason\": \"engine_startup_abort\""));
    EXPECT_NE(std::string::npos, js.find("\"process_exit_code\": 20"));
    EXPECT_NE(std::string::npos, js.find("\"engine_startup_abort\": 1"));
    EXPECT_NE(std::string::npos, js.find("\"none\": 1"));
    // The embedded quotes in the stderr tail must be escaped, not emitted raw.
    EXPECT_NE(std::string::npos, js.find("\\\"/tmp/x\\\""));
    EXPECT_EQ(std::string::npos, js.find("path \"/tmp/x\""));

    // NOT COVERED: full JSON well-formedness. That is asserted by
    // tests/test_startup_diagnosis_contract.py, which parses it with json.loads.
}

TEST_F(StartupDiagnosis, TailReadIsBoundedAndSurvivesAHugeLog) {
    // A healthy stderr.log is megabytes; reading it whole in the parent would
    // be a new failure mode in the code meant to diagnose failure modes.
    std::error_code ec;
    const std::filesystem::path dir = root_ / "BIGL";
    std::filesystem::create_directories(dir, ec);
    {
        std::ofstream ofs(dir / "stderr.log", std::ios::binary);
        for (int i = 0; i < 200000; ++i) ofs << "[SMFREE] gen=" << i << " filler line\n";
        ofs << kIncidentStderr;   // terminal handler, as it always is: last
    }
    const auto ev = sd::read_child_evidence(dir.string());
    EXPECT_TRUE(ev.engine_terminal_handler);
    EXPECT_LE(ev.stderr_tail.size(), 6u);
    EXPECT_NE(std::string::npos, ev.stderr_tail.back().find("temp_directory_path"));

    // An unreadable log degrades to "no evidence" instead of throwing.
    const auto missing = sd::read_child_evidence((root_ / "NOPE").string());
    EXPECT_FALSE(missing.engine_terminal_handler);
    EXPECT_TRUE(missing.stderr_tail.empty());

    // NOT COVERED: a "Fatal error: " emitted more than 8 KB before the end of
    // the log. The handler runs immediately before exit, so a marker that far
    // back is not the terminal one — but a future engine that logged after its
    // handler would defeat this window.
}

}  // namespace
