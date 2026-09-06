#ifndef FLEXAIDS_GA_POPULATION_RECEIPT_H
#define FLEXAIDS_GA_POPULATION_RECEIPT_H

#include "gaboom.h"
#include "flexaid_exception.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <locale>
#include <span>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace flexaids {

struct GaPopulationWorkerReceipt {
    int team_size = 0;
    std::uint64_t evaluated_chromosomes = 0;
};

struct GaPopulationEvaluationBatch {
    std::string region;
    int population_count;
    int popoffset;
    std::vector<GaPopulationWorkerReceipt> workers;
};

inline bool ga_population_receipt_requested() {
    const char* path = std::getenv("FLEXAIDDS_GEN0_RECEIPT");
    return path && *path;
}

// One collector per calling GA thread, scoped around initial population creation
// only. Recursive IPFILE completion uses the same collector; later adaptive
// repopulation/offspring evaluation runs after the scope has ended.
class GaPopulationObservation {
public:
    GaPopulationObservation() : enabled_(ga_population_receipt_requested()) {
        if (enabled_) {
            if (current)
                throw FlexAIDException("Generation-zero receipt: nested GA observation");
            current = this;
        }
    }
    ~GaPopulationObservation() { if (enabled_) current = nullptr; }
    GaPopulationObservation(const GaPopulationObservation&) = delete;
    GaPopulationObservation& operator=(const GaPopulationObservation&) = delete;

    std::vector<GaPopulationEvaluationBatch> batches;
    inline static thread_local GaPopulationObservation* current = nullptr;
private:
    bool enabled_;
};

inline bool ga_population_observation_active() {
    return GaPopulationObservation::current != nullptr;
}

// Called only by the OpenMP worker owning this slot, after a completed evaluation.
// The caller must wait for the parallel region before reading the slots.
inline void record_ga_population_worker(GaPopulationWorkerReceipt& slot, int team_size) {
    slot.team_size = team_size;
    ++slot.evaluated_chromosomes;
}

inline void observe_ga_population_workers(
        const char* region, int population_count, int popoffset,
        std::span<const GaPopulationWorkerReceipt> workers) {
    auto* observation = GaPopulationObservation::current;
    if (!observation) return;
    if (!region || (std::strcmp(region, "populate_chromosomes") != 0 &&
                    std::strcmp(region, "calculate_fitness") != 0))
        throw FlexAIDException("Generation-zero receipt: unknown evaluation region");
    if (population_count <= 0 || popoffset < 0 || popoffset > population_count)
        throw FlexAIDException("Generation-zero receipt: invalid evaluation range");
    std::uint64_t total = 0;
    for (std::size_t tid = 0; tid < workers.size(); ++tid) {
        const auto& worker = workers[tid];
        if (!worker.evaluated_chromosomes) continue;
        if (worker.team_size <= 0 || tid >= static_cast<std::size_t>(worker.team_size) ||
            static_cast<std::size_t>(worker.team_size) > workers.size() ||
            worker.evaluated_chromosomes > static_cast<std::uint64_t>(population_count - popoffset))
            throw FlexAIDException("Generation-zero receipt: invalid observed worker participation");
        if (worker.evaluated_chromosomes >
                static_cast<std::uint64_t>(population_count - popoffset) - total)
            throw FlexAIDException("Generation-zero receipt: evaluation count exceeds observed range");
        total += worker.evaluated_chromosomes;
    }
    if (total > static_cast<std::uint64_t>(population_count - popoffset))
        throw FlexAIDException("Generation-zero receipt: evaluation count exceeds observed range");
    if (total)
        observation->batches.push_back({region, population_count, popoffset,
                                       {workers.begin(), workers.end()}});
}

namespace population_receipt_detail {

// Read object representations, never convert or perform arithmetic on a score.
// This preserves signed zero, subnormals, and even NaN payloads for diagnosis.
inline std::string hex_integer(std::uint64_t bits, int digits) {
    constexpr char hex[] = "0123456789abcdef";
    std::string result(static_cast<std::size_t>(digits) + 2, '0');
    result[1] = 'x';
    for (int i = digits + 1; i >= 2; --i) {
        result[static_cast<std::size_t>(i)] = hex[bits & 15];
        bits >>= 4;
    }
    return result;
}

inline std::string bits(const double& value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559);
    std::uint64_t representation;
    std::memcpy(&representation, &value, sizeof(representation));
    return hex_integer(representation, 16);
}

inline std::string bits(const float& value) {
    static_assert(sizeof(float) == sizeof(std::uint32_t) &&
                  std::numeric_limits<float>::is_iec559);
    std::uint32_t representation;
    std::memcpy(&representation, &value, sizeof(representation));
    return hex_integer(representation, 8);
}

inline void field(std::ostream& out, const char* name, const double& value) {
    out << ",\"" << name << "_bits\":\"" << bits(value) << '\"';
}

inline FlexAIDException io_error(const char* operation, const char* path, int code) {
    return FlexAIDException(std::string("Generation-zero receipt: ") + operation +
                           " '" + path + "': " +
                           std::error_code(code, std::generic_category()).message());
}

// Always close the owned stream, even after a short write or flush error. A
// failed write may leave a partial receipt; successful process completion is
// mandatory in addition to the JSON completion marker when admitting a run.
inline void finish_output(FILE* output, const std::string& payload, const char* path) {
    int error = 0;
    if (std::fwrite(payload.data(), 1, payload.size(), output) != payload.size())
        error = errno ? errno : EIO;
    if (std::fflush(output) != 0 && !error) error = errno ? errno : EIO;
    if (std::fclose(output) != 0 && !error) error = errno ? errno : EIO;
    if (error) throw io_error("cannot finish", path, error);
}

}  // namespace population_receipt_detail

// Called only after initial populate_chromosomes() has returned, including its
// existing fitness calculation and sort, and before any generation/reproduction.
// Array order is retained; consumers may compare multisets while keeping duplicate
// multiplicities. Sums belong in the external comparator, not the docking engine.
inline std::string serialize_ga_population_receipt(
        std::span<const chromosome> population, int n_genes, std::uint64_t seed) {
    using namespace population_receipt_detail;
    if (population.empty() || n_genes <= 0)
        throw FlexAIDException("Generation-zero receipt: empty population or invalid gene count");
    for (const auto& chrom : population) {
        if (!chrom.genes)
            throw FlexAIDException("Generation-zero receipt: null chromosome genes");
    }

    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "{\"schema\":\"flexaidds.ga_population_receipt.v1\","
           "\"boundary\":\"initial_population_complete_before_reproduction\","
           "\"generation\":0,\"population_count\":" << population.size()
        << ",\"n_genes\":" << n_genes << ",\"seed\":\"" << seed
        << "\",\"execution\":{\"openmp_compiled\":"
#ifdef _OPENMP
        << "true"
#else
        << "false"
#endif
        << ",\"deterministic_compile\":"
#ifdef FLEXAID_DETERMINISTIC
        << "true"
#else
        << "false"
#endif
        << ",\"evaluation_batches\":[";
    if (const auto* observation = GaPopulationObservation::current) {
        for (std::size_t b = 0; b < observation->batches.size(); ++b) {
            const auto& batch = observation->batches[b];
            if (b) out << ',';
            out << "{\"region\":\"" << batch.region << "\",\"population_count\":"
                << batch.population_count << ",\"popoffset\":" << batch.popoffset
                << ",\"workspace_slots\":" << batch.workers.size() << ",\"workers\":[";
            bool first = true;
            for (std::size_t tid = 0; tid < batch.workers.size(); ++tid) {
                const auto& worker = batch.workers[tid];
                if (!worker.evaluated_chromosomes) continue;
                if (!first) out << ',';
                first = false;
                out << "{\"worker_id\":" << tid << ",\"team_size\":" << worker.team_size
                    << ",\"evaluated_chromosomes\":" << worker.evaluated_chromosomes << '}';
            }
            out << "]}";
        }
    }
    out << "]},\"records\":[\n";
    for (std::size_t i = 0; i < population.size(); ++i) {
        const auto& chrom = population[i];
        if (i) out << ",\n";
        out << "{\"index\":" << i << ",\"status\":"
            << static_cast<unsigned int>(static_cast<unsigned char>(chrom.status))
            << ",\"genes\":[";
        for (int g = 0; g < n_genes; ++g) {
            if (g) out << ',';
            out << "{\"to_int32\":" << chrom.genes[g].to_int32
                << ",\"to_ic_bits\":\"" << bits(chrom.genes[g].to_ic) << "\"}";
        }
        out << "],\"cf\":{\"rclash\":" << chrom.cf.rclash;
        field(out, "com", chrom.cf.com);
        field(out, "con", chrom.cf.con);
        field(out, "wal", chrom.cf.wal);
        field(out, "sas", chrom.cf.sas);
        field(out, "elec", chrom.cf.elec);
        field(out, "gist", chrom.cf.gist);
        field(out, "hbond", chrom.cf.hbond);
        field(out, "gist_desolv", chrom.cf.gist_desolv);
        field(out, "metal_coord", chrom.cf.metal_coord);
        field(out, "h_rep", chrom.cf.h_rep);
        field(out, "entropy", chrom.cf.entropy);
        field(out, "pb_clash", chrom.cf.pb_clash);
        field(out, "totsas", chrom.cf.totsas);
        field(out, "nor", chrom.cf.nor);
        out << '}';
        field(out, "evalue", chrom.evalue);
        field(out, "app_evalue", chrom.app_evalue);
        field(out, "fitnes", chrom.fitnes);
        field(out, "boltzmann_weight", chrom.boltzmann_weight);
        field(out, "free_energy", chrom.free_energy);
        out << ",\"ring_phases_bits\":[";
        for (int r = 0; r < MAX_RING_FLEX; ++r) {
            if (r) out << ',';
            out << '\"' << bits(chrom.ring_phases[r]) << '\"';
        }
        out << "],\"ring_six\":[";
        for (int r = 0; r < MAX_RING_FLEX; ++r) {
            if (r) out << ',';
            out << static_cast<unsigned int>(chrom.ring_six[r]);
        }
        out << "],\"ring_five\":[";
        for (int r = 0; r < MAX_RING_FLEX; ++r) {
            if (r) out << ',';
            out << static_cast<unsigned int>(chrom.ring_five[r]);
        }
        out << "]}";
    }
    out << "\n],\"complete\":true}\n";
    if (!out) throw FlexAIDException("Generation-zero receipt: serialization failed");
    return out.str();
}

// Exclusive creation deliberately rejects an existing receipt, including another
// GA invocation/restart using this path. Adaptive repopulations do not call this
// observer. Each independently launched GA needs a distinct receipt path.
inline void write_ga_population_receipt(
        const char* path, std::span<const chromosome> population,
        int n_genes, std::uint64_t seed) {
    using namespace population_receipt_detail;
    if (!path || !*path)
        throw FlexAIDException("Generation-zero receipt: empty output path");
    const std::string payload = serialize_ga_population_receipt(population, n_genes, seed);
    FILE* output = std::fopen(path, "wbx");
    if (!output) throw io_error("cannot create", path, errno);
    finish_output(output, payload, path);
}

// OFF unless a nonempty path is explicitly supplied. No serialization, scoring,
// RNG calls, or input inspection occurs when disabled.
inline void write_ga_population_receipt_if_requested(
        std::span<const chromosome> population, int n_genes, std::uint64_t seed) {
    const char* path = std::getenv("FLEXAIDDS_GEN0_RECEIPT");
    if (path && *path) write_ga_population_receipt(path, population, n_genes, seed);
}

}  // namespace flexaids

#endif  // FLEXAIDS_GA_POPULATION_RECEIPT_H
