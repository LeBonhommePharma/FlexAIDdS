#include <gtest/gtest.h>

#include "GaPopulationReceipt.h"
#include "json_value.h"

#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Population {
    std::array<std::array<gene, 2>, 2> genes{};
    std::array<chromosome, 2> chrom{};

    Population() {
        for (std::size_t i = 0; i < chrom.size(); ++i) {
            chrom[i].genes = genes[i].data();
            chrom[i].status = 'n';
        }
        genes[0][0] = {INT32_MIN, -0.0};
        genes[0][1] = {INT32_MAX, std::bit_cast<double>(UINT64_C(1))};
        genes[1][0] = {-7, 1.0};
        genes[1][1] = {42, -2.0};
    }
};

class TempDirectory {
public:
    TempDirectory() {
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int attempt = 0; attempt < 100; ++attempt) {
            path = std::filesystem::temp_directory_path() /
                   ("ga-population-receipt-test-" + std::to_string(tick) + "-" +
                    std::to_string(attempt));
            if (std::filesystem::create_directory(path)) return;
        }
        throw std::runtime_error("Cannot create test directory");
    }
    ~TempDirectory() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
    std::filesystem::path path;
};

class ReceiptEnvironment {
public:
    ReceiptEnvironment() {
        const char* previous = std::getenv("FLEXAIDDS_GEN0_RECEIPT");
        present_ = previous != nullptr;
        if (previous) previous_ = previous;
    }
    ~ReceiptEnvironment() { set(present_ ? previous_.c_str() : nullptr); }
    static void set(const char* value) {
#ifdef _WIN32
        _putenv_s("FLEXAIDDS_GEN0_RECEIPT", value ? value : "");
#else
        if (value) setenv("FLEXAIDDS_GEN0_RECEIPT", value, 1);
        else unsetenv("FLEXAIDDS_GEN0_RECEIPT");
#endif
    }
private:
    bool present_;
    std::string previous_;
};

std::string read_file(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

std::uint64_t encoded_bits(const json::Value& value) {
    const auto text = value.as_string();
    if (text.size() != 18 || text.substr(0, 2) != "0x")
        throw std::runtime_error("Not an exact IEEE binary64 encoding");
    std::size_t end = 0;
    const auto result = std::stoull(text, &end, 16);
    if (end != text.size()) throw std::runtime_error("Trailing bit encoding content");
    return result;
}

}  // namespace

TEST(GaPopulationReceipt, ExactGeneIdentityAndEveryStoredScoreComponent) {
    Population p;
    const std::array<std::pair<const char*, double cfstr::*>, 14> fields{{
        {"com_bits", &cfstr::com}, {"con_bits", &cfstr::con},
        {"wal_bits", &cfstr::wal}, {"sas_bits", &cfstr::sas},
        {"elec_bits", &cfstr::elec}, {"gist_bits", &cfstr::gist},
        {"hbond_bits", &cfstr::hbond}, {"gist_desolv_bits", &cfstr::gist_desolv},
        {"metal_coord_bits", &cfstr::metal_coord}, {"h_rep_bits", &cfstr::h_rep},
        {"entropy_bits", &cfstr::entropy}, {"pb_clash_bits", &cfstr::pb_clash},
        {"totsas_bits", &cfstr::totsas}, {"nor_bits", &cfstr::nor}
    }};
    for (std::size_t f = 0; f < fields.size(); ++f)
        p.chrom[0].cf.*fields[f].second =
            std::bit_cast<double>(UINT64_C(0x3ff0000000000000) + f);
    p.chrom[0].cf.rclash = 1;
    p.chrom[0].evalue = -1.0;
    p.chrom[0].app_evalue = 2.0;
    p.chrom[0].fitnes = 4.0;
    p.chrom[0].boltzmann_weight = 0.5;
    p.chrom[0].free_energy = -2.0;
    p.chrom[0].ring_phases[0] = -0.0f;
    p.chrom[0].ring_phases[MAX_RING_FLEX - 1] = 1.0f;
    p.chrom[0].ring_six[MAX_RING_FLEX - 1] = 255;
    p.chrom[0].ring_five[0] = 17;

    const auto receipt = json::parse(flexaids::serialize_ga_population_receipt(
        p.chrom, 2, UINT64_MAX));
    EXPECT_EQ(receipt["schema"].as_string(), "flexaidds.ga_population_receipt.v1");
    EXPECT_EQ(receipt["boundary"].as_string(),
              "initial_population_complete_before_reproduction");
    EXPECT_EQ(receipt["generation"].as_int(-1), 0);
    EXPECT_EQ(receipt["population_count"].as_int(), 2);
    EXPECT_EQ(receipt["n_genes"].as_int(), 2);
    EXPECT_EQ(receipt["seed"].as_string(), "18446744073709551615");
    EXPECT_TRUE(receipt["complete"].as_bool());
    const auto& records = receipt["records"].as_array();
    ASSERT_EQ(records.size(), 2U);
    const auto& first = records[0];
    EXPECT_EQ(first["index"].as_int(-1), 0);
    EXPECT_EQ(records[1]["index"].as_int(-1), 1);
    EXPECT_EQ(first["status"].as_int(), static_cast<int>('n'));
    const auto& genes = first["genes"].as_array();
    ASSERT_EQ(genes.size(), 2U);
    EXPECT_EQ(genes[0]["to_int32"].as_int(), INT32_MIN);
    EXPECT_EQ(genes[1]["to_int32"].as_int(), INT32_MAX);
    EXPECT_EQ(encoded_bits(genes[0]["to_ic_bits"]), UINT64_C(0x8000000000000000));
    EXPECT_EQ(encoded_bits(genes[1]["to_ic_bits"]), UINT64_C(1));
    ASSERT_EQ(first["cf"].size(), fields.size() + 1);
    for (std::size_t f = 0; f < fields.size(); ++f)
        EXPECT_EQ(encoded_bits(first["cf"][fields[f].first]),
                  UINT64_C(0x3ff0000000000000) + f) << fields[f].first;
    EXPECT_EQ(first["cf"]["rclash"].as_int(), 1);
    EXPECT_EQ(encoded_bits(first["evalue_bits"]), UINT64_C(0xbff0000000000000));
    EXPECT_EQ(encoded_bits(first["app_evalue_bits"]), UINT64_C(0x4000000000000000));
    EXPECT_EQ(encoded_bits(first["fitnes_bits"]), UINT64_C(0x4010000000000000));
    EXPECT_EQ(encoded_bits(first["boltzmann_weight_bits"]), UINT64_C(0x3fe0000000000000));
    EXPECT_EQ(encoded_bits(first["free_energy_bits"]), UINT64_C(0xc000000000000000));
    ASSERT_EQ(first["ring_phases_bits"].size(), static_cast<std::size_t>(MAX_RING_FLEX));
    EXPECT_EQ(first["ring_phases_bits"].as_array().front().as_string(), "0x80000000");
    EXPECT_EQ(first["ring_phases_bits"].as_array().back().as_string(), "0x3f800000");
    ASSERT_EQ(first["ring_six"].size(), static_cast<std::size_t>(MAX_RING_FLEX));
    ASSERT_EQ(first["ring_five"].size(), static_cast<std::size_t>(MAX_RING_FLEX));
    EXPECT_EQ(first["ring_six"].as_array().back().as_int(), 255);
    EXPECT_EQ(first["ring_five"].as_array().front().as_int(), 17);
}

TEST(GaPopulationReceipt, PreservesNonfinitePayloadsAndDoesNotMutatePopulation) {
    Population p;
    const auto signaling_nan = UINT64_C(0x7ff0000000001234);
    std::memcpy(&p.chrom[0].cf.com, &signaling_nan, sizeof(signaling_nan));
    p.chrom[0].cf.wal = std::numeric_limits<double>::infinity();
    std::array<unsigned char, sizeof(p.chrom)> before_chrom;
    std::array<unsigned char, sizeof(p.genes)> before_genes;
    std::memcpy(before_chrom.data(), p.chrom.data(), sizeof(p.chrom));
    std::memcpy(before_genes.data(), p.genes.data(), sizeof(p.genes));
    const auto payload = flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345);
    const auto receipt = json::parse(payload);
    EXPECT_EQ(encoded_bits(receipt["records"].as_array()[0]["cf"]["com_bits"]), signaling_nan);
    EXPECT_EQ(encoded_bits(receipt["records"].as_array()[0]["cf"]["wal_bits"]),
              UINT64_C(0x7ff0000000000000));
    EXPECT_EQ(std::memcmp(before_chrom.data(), p.chrom.data(), sizeof(p.chrom)), 0);
    EXPECT_EQ(std::memcmp(before_genes.data(), p.genes.data(), sizeof(p.genes)), 0);
    EXPECT_EQ(flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345), payload);
}

TEST(GaPopulationReceipt, DisabledObserverDoesNotInspectInvalidInput) {
    ReceiptEnvironment environment;
    environment.set(nullptr);
    EXPECT_NO_THROW(flexaids::write_ga_population_receipt_if_requested({}, -1, 0));
    environment.set("");
    EXPECT_NO_THROW(flexaids::write_ga_population_receipt_if_requested({}, -1, 0));
}

TEST(GaPopulationReceipt, RequestedObserverWritesCompleteReceiptAndRefusesOverwrite) {
    Population p;
    TempDirectory directory;
    ReceiptEnvironment environment;
    const auto path = (directory.path / "initial.json").string();
    environment.set(path.c_str());
    flexaids::write_ga_population_receipt_if_requested(p.chrom, 2, 12345);
    const auto original = read_file(path);
    EXPECT_EQ(original, flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345));
    p.genes[0][0].to_int32 = 19;
    EXPECT_THROW(flexaids::write_ga_population_receipt_if_requested(p.chrom, 2, 12345),
                 FlexAIDException);
    EXPECT_EQ(read_file(path), original);
}

TEST(GaPopulationReceipt, ConcurrentInvocationsCannotOverwriteEachOther) {
    Population p;
    TempDirectory directory;
    const auto path = (directory.path / "initial.json").string();
    std::atomic<int> successes{0}, failures{0};
    auto write = [&] {
        try {
            flexaids::write_ga_population_receipt(path.c_str(), p.chrom, 2, 12345);
            ++successes;
        } catch (const FlexAIDException&) { ++failures; }
    };
    std::thread first(write), second(write);
    first.join();
    second.join();
    EXPECT_EQ(successes.load(), 1);
    EXPECT_EQ(failures.load(), 1);
    EXPECT_EQ(read_file(path), flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345));
}

TEST(GaPopulationReceipt, InvalidInputAndUnwritableDestinationsFailClosed) {
    Population p;
    TempDirectory directory;
    const auto path = (directory.path / "initial.json").string();
    EXPECT_THROW(flexaids::write_ga_population_receipt(path.c_str(), {}, 2, 0), FlexAIDException);
    EXPECT_THROW(flexaids::write_ga_population_receipt(path.c_str(), p.chrom, 0, 0), FlexAIDException);
    p.chrom[1].genes = nullptr;
    EXPECT_THROW(flexaids::write_ga_population_receipt(path.c_str(), p.chrom, 2, 0), FlexAIDException);
    EXPECT_FALSE(std::filesystem::exists(path));
    p.chrom[1].genes = p.genes[1].data();
    EXPECT_THROW(flexaids::write_ga_population_receipt(nullptr, p.chrom, 2, 0), FlexAIDException);
    EXPECT_THROW(flexaids::write_ga_population_receipt("", p.chrom, 2, 0), FlexAIDException);
    const auto absent_parent = (directory.path / "absent" / "initial.json").string();
    EXPECT_THROW(flexaids::write_ga_population_receipt(absent_parent.c_str(), p.chrom, 2, 0),
                 FlexAIDException);
    const auto directory_path = directory.path.string();
    EXPECT_THROW(flexaids::write_ga_population_receipt(directory_path.c_str(), p.chrom, 2, 0),
                 FlexAIDException);
}

TEST(GaPopulationReceipt, ActualStreamWriteFailureIsReportedWithoutReplacingEvidence) {
    Population p;
    TempDirectory directory;
    const auto path = (directory.path / "existing.json").string();
    const std::string evidence = "preserved evidence\n";
    { std::ofstream output(path, std::ios::binary); output << evidence; }
    FILE* read_only = std::fopen(path.c_str(), "r");
    ASSERT_NE(read_only, nullptr);
    const auto payload = flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345);
    // The production finalizer owns/closes this stream on both success and error.
    EXPECT_THROW(flexaids::population_receipt_detail::finish_output(
                     read_only, payload, path.c_str()), FlexAIDException);
    EXPECT_EQ(read_file(path), evidence);
}

TEST(GaPopulationReceipt, CountsActualEvaluationWorkersAndTeamSizes) {
    ReceiptEnvironment environment;
    environment.set("observation-only-no-file");
    Population p;
    constexpr int jobs = 64;
    std::vector<chromosome> population(jobs, p.chrom[0]);
#ifdef _OPENMP
    const int was_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    const std::array<int, 2> teams{1, 4};
#else
    const std::array<int, 1> teams{1};
#endif
    for (const int requested : teams) {
        flexaids::GaPopulationObservation observation;
        ASSERT_TRUE(flexaids::ga_population_observation_active());
        std::vector<flexaids::GaPopulationWorkerReceipt> workers(requested);
#ifdef _OPENMP
#pragma omp parallel for num_threads(requested) schedule(static)
#endif
        for (int i = 0; i < jobs; ++i) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int team = omp_get_num_threads();
#else
            const int tid = 0, team = 1;
#endif
            flexaids::record_ga_population_worker(workers[tid], team);
        }
        flexaids::observe_ga_population_workers("populate_chromosomes", jobs, 0, workers);
        const auto receipt = json::parse(flexaids::serialize_ga_population_receipt(population, 2, 12345));
        const auto& execution = receipt["execution"];
#ifdef _OPENMP
        EXPECT_TRUE(execution["openmp_compiled"].as_bool());
#else
        EXPECT_FALSE(execution["openmp_compiled"].as_bool(true));
#endif
        const auto& batches = execution["evaluation_batches"].as_array();
        ASSERT_EQ(batches.size(), 1U);
        EXPECT_EQ(batches[0]["region"].as_string(), "populate_chromosomes");
        EXPECT_EQ(batches[0]["population_count"].as_int(), jobs);
        EXPECT_EQ(batches[0]["popoffset"].as_int(-1), 0);
        EXPECT_EQ(batches[0]["workspace_slots"].as_int(), requested);
        const auto& observed_workers = batches[0]["workers"].as_array();
        ASSERT_EQ(observed_workers.size(), static_cast<std::size_t>(requested));
        int count = 0;
        for (int tid = 0; tid < requested; ++tid) {
            EXPECT_EQ(observed_workers[tid]["worker_id"].as_int(-1), tid);
            EXPECT_EQ(observed_workers[tid]["team_size"].as_int(), requested);
            EXPECT_GT(observed_workers[tid]["evaluated_chromosomes"].as_int(), 0);
            count += observed_workers[tid]["evaluated_chromosomes"].as_int();
        }
        EXPECT_EQ(count, jobs);
    }
#ifdef _OPENMP
    omp_set_dynamic(was_dynamic);
#endif
    EXPECT_FALSE(flexaids::ga_population_observation_active());
}

TEST(GaPopulationReceipt, RecursiveInitialBatchesAccumulateButLaterRepopulationIsNotObserved) {
    ReceiptEnvironment environment;
    environment.set("observation-only-no-file");
    Population p;
    const std::array<flexaids::GaPopulationWorkerReceipt, 1> one{{{1, 1}}};
    {
        flexaids::GaPopulationObservation observation;
        // Model IPFILE: remaining random chromosome first, loaded chromosome
        // evaluated by initial calculate_fitness second. Retain both witnesses.
        flexaids::observe_ga_population_workers("populate_chromosomes", 2, 1, one);
        flexaids::observe_ga_population_workers("calculate_fitness", 2, 0, one);
        const std::array<flexaids::GaPopulationWorkerReceipt, 1> no_evaluations{};
        flexaids::observe_ga_population_workers("calculate_fitness", 2, 0, no_evaluations);
        const auto receipt = json::parse(flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345));
        const auto& batches = receipt["execution"]["evaluation_batches"].as_array();
        ASSERT_EQ(batches.size(), 2U);
        EXPECT_EQ(batches[0]["popoffset"].as_int(), 1);
        EXPECT_EQ(batches[1]["region"].as_string(), "calculate_fitness");
        EXPECT_THROW((flexaids::GaPopulationObservation{}), FlexAIDException);
    }
    EXPECT_FALSE(flexaids::ga_population_observation_active());
    flexaids::observe_ga_population_workers("populate_chromosomes", 2, 0, one);
    const auto after = json::parse(flexaids::serialize_ga_population_receipt(p.chrom, 2, 12345));
    EXPECT_EQ(after["execution"]["evaluation_batches"].size(), 0U);
}

TEST(GaPopulationReceipt, CallingThreadCollectorsAreIsolatedAndRejectInvalidWorkerMetadata) {
    ReceiptEnvironment environment;
    environment.set("observation-only-no-file");
    flexaids::GaPopulationObservation parent;
    std::atomic<int> observed{0};
    auto child = [&] {
        if (flexaids::ga_population_observation_active()) return;
        flexaids::GaPopulationObservation own;
        const std::array<flexaids::GaPopulationWorkerReceipt, 1> worker{{{1, 1}}};
        flexaids::observe_ga_population_workers("populate_chromosomes", 1, 0, worker);
        if (own.batches.size() == 1) ++observed;
    };
    std::thread first(child), second(child);
    first.join();
    second.join();
    EXPECT_EQ(observed.load(), 2);
    EXPECT_TRUE(parent.batches.empty());
    const std::array<flexaids::GaPopulationWorkerReceipt, 1> oversized_team{{{4, 1}}};
    const std::array<flexaids::GaPopulationWorkerReceipt, 1> too_many{{{1, 3}}};
    const std::array<flexaids::GaPopulationWorkerReceipt, 2> too_many_combined{{{2, 2}, {2, 2}}};
    EXPECT_THROW(flexaids::observe_ga_population_workers("populate_chromosomes", 2, 0, oversized_team), FlexAIDException);
    EXPECT_THROW(flexaids::observe_ga_population_workers("populate_chromosomes", 2, 0, too_many), FlexAIDException);
    EXPECT_THROW(flexaids::observe_ga_population_workers("populate_chromosomes", 2, 0, too_many_combined), FlexAIDException);
    EXPECT_THROW(flexaids::observe_ga_population_workers("unknown", 2, 0, oversized_team), FlexAIDException);
    EXPECT_TRUE(parent.batches.empty());
}
