// test_fleet_runner.cpp - Fleet chunk contract tests.
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "FleetRunner.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

namespace {

fleet::ChunkMetadata metadata() {
    return {
        "campaign-1",
        "chunk-0001",
        "attempt-1",
        "worker\"one",
        "astex",
        "benchmark_datasets --fleet",
        "/build/benchmark_datasets",
        std::string(64, 'a'),
        "/build/FlexAIDdS",
        std::string(64, 'b'),
        std::string(64, 'c'),
    };
}

dataset::DockingResult successful_result() {
    dataset::DockingResult result;
    result.pdb_id = "1GPK";
    result.num_poses = 20;
    result.docking_completed = true;
    result.docking_exit_code = 0;
    result.elected_pose_path = "1GPK/elected_pose.pdb";
    result.elected_pose_source = "r2/1GPK_0.pdb";
    result.rmsd_hungarian = 1.25f;
    result.success_rmsd = true;
    result.pb_ran = true;
    result.pb_pass = true;
    result.success_pb = true;
    result.protocol_claim_eligible = true;
    result.claim_ready = true;
    result.tencom_status = "ok";
    result.eigen_status = "ok";
    result.eigen_n_modes = 12;
    result.pose_sha256 = std::string(64, 'd');
    result.rmsd_pose_sha256 = result.pose_sha256;
    result.posebusters_pose_sha256 = result.pose_sha256;
    result.tencom_pose_sha256 = result.pose_sha256;
    return result;
}

} // namespace

TEST(FleetRunner, SeparatesExecutionFromScientificSuccess) {
    dataset::BenchmarkReport report;
    report.dataset_name = "Astex Diverse";
    report.total_systems = 2;
    report.results.push_back(successful_result());

    dataset::DockingResult miss;
    miss.pdb_id = "1HNN";
    miss.num_poses = 20;
    miss.docking_completed = true;
    miss.docking_exit_code = 0;
    miss.elected_pose_path = "1HNN/elected_pose.pdb";
    miss.rmsd_hungarian = 8.0f;
    miss.pb_ran = true;
    miss.pb_pass = true;
    miss.tencom_status = "ok";
    miss.eigen_status = "ok";
    miss.pose_sha256 = std::string(64, 'e');
    miss.rmsd_pose_sha256 = miss.pose_sha256;
    miss.posebusters_pose_sha256 = miss.pose_sha256;
    miss.tencom_pose_sha256 = miss.pose_sha256;
    report.results.push_back(miss);

    dataset::DockingConfig config;
    config.mode = dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK;
    const std::string json = fleet::FleetRunner::serialize_chunk_result(
        metadata(), report, config, 12.5);

    EXPECT_NE(json.find("\"execution_completed\": 2"), std::string::npos);
    EXPECT_NE(json.find("\"execution_failed\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"success_rmsd\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"success_pb\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"claim_ready\": 1"), std::string::npos);
    EXPECT_EQ(json.find("\"dG_mean\""), std::string::npos);
}

TEST(FleetRunner, EscapesMetadataAndEmitsBinaryHashes) {
    dataset::BenchmarkReport report;
    report.total_systems = 1;
    report.results.push_back(successful_result());
    dataset::DockingConfig config;

    const std::string json = fleet::FleetRunner::serialize_chunk_result(
        metadata(), report, config, 0.5);

    EXPECT_NE(json.find("worker\\\"one"), std::string::npos);
    EXPECT_NE(json.find(std::string(64, 'a')), std::string::npos);
    EXPECT_NE(json.find(std::string(64, 'b')), std::string::npos);
    EXPECT_NE(json.find(std::string(64, 'c')), std::string::npos);
}

TEST(FleetRunner, RefusesToOverwritePublishedResult) {
    const fs::path root = fs::temp_directory_path() /
        ("fleet_runner_test_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    const fs::path result = root / "result.json";
    std::string error;

    ASSERT_TRUE(fleet::FleetRunner::write_chunk_result_atomic(
        result.string(), "{\"attempt\":1}\n", &error)) << error;
    EXPECT_FALSE(fleet::FleetRunner::write_chunk_result_atomic(
        result.string(), "{\"attempt\":2}\n", &error));

    std::ifstream input(result);
    std::string contents((std::istreambuf_iterator<char>(input)),
                         std::istreambuf_iterator<char>());
    EXPECT_EQ(contents, "{\"attempt\":1}\n");
    fs::remove_all(root);
}

TEST(FleetRunner, RetainedPoseDoesNotHideChildFailure) {
    dataset::BenchmarkReport report;
    report.total_systems = 1;
    auto result = successful_result();
    result.docking_exit_code = 6;
    result.docking_completed = false;
    result.claim_ready = false;
    result.matrix_md5 = std::string(32, 'a');
    report.results = {result};
    dataset::DockingConfig config;
    const auto json = fleet::FleetRunner::serialize_chunk_result(metadata(), report, config, 0.5);
    EXPECT_NE(json.find("\"execution_completed\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"execution_failed\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"docking_exit_code\": 6"), std::string::npos);
    EXPECT_NE(json.find("\"matrix_md5\": \"" + std::string(32, 'a') + "\""), std::string::npos);
}

TEST(FleetRunner, CompletedDockWithoutReferenceIsStillRuntimeCompleted) {
    dataset::BenchmarkReport report;
    report.total_systems = 1;
    dataset::DockingResult result;
    result.pdb_id = "1GPK";
    result.num_poses = 20;
    result.docking_completed = true;
    result.docking_exit_code = 0;
    // DatasetRunner cannot elect/measure against an unavailable RMSD reference,
    // but that does not undo the witnessed child exit and its produced poses.
    result.rmsd_fail_reason = "input_missing";
    EXPECT_TRUE(result.elected_pose_path.empty());
    EXPECT_FALSE(result.claim_ready);
    report.results = {result};
    EXPECT_EQ(dataset::benchmark_runtime_exit_code(report), 0);
    dataset::DockingConfig config;
    const auto json = fleet::FleetRunner::serialize_chunk_result(metadata(), report, config, 0.5);
    EXPECT_NE(json.find("\"execution_completed\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"execution_failed\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"execution_completed\": true"), std::string::npos);
    EXPECT_NE(json.find("\"validators_complete\": false"), std::string::npos);
    EXPECT_NE(json.find("\"claim_ready\": false"), std::string::npos);
}
