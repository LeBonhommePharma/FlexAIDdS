// GENERATED -- DO NOT EDIT, run scripts/gen_dataset_codes.py
//
// Every PDB code list the engine executes, generated from ONE declared canonical
// source per dataset. Hand-editing this file is what the generator exists to stop:
// scripts/gen_dataset_codes.py --check regenerates and diffs, and fails with both
// sha256 values if the committed file disagrees with its source.
//
// Sources (path, sha256 of the source AT GENERATION TIME, provenance):
//   astex_diverse  benchmarks/astex_diverse/astex_diverse_set.csv
//                  sha256=a6384bdda9ba5a73d58284149287f584452805b549b0ea665be944f6e6497708
//                  n=85 provenance=ROSTER_CSV primary_source=10.1021/jm061277y
//   hap2           benchmarks/datasets/hap2.yaml
//                  sha256=8dcf34e6a025093587d066dd8ecc110bcf02d7d48df7ef05b329d8b4f52cbda1
//                  n=59 provenance=DEFINITION_YAML primary_source=10.1021/acs.jcim.5b00078
//   casf2016       benchmarks/datasets/casf2016_core_codes.csv
//                  sha256=c90f591e09fb3c2edeb808a19ba2d9b4f96b6f5b9bd726009de91d17c88a1b1a
//                  n=283 provenance=EXTRACTED_FROM_CODE primary_source=UNVERIFIED
//
// provenance=EXTRACTED_FROM_CODE means the source was lifted out of the C++ literal
// because no definition file held the list. That is CIRCULAR: it proves only what
// the code already said. scripts/check_dataset_identity.py reports it UNVERIFIED.
// =============================================================================

#ifndef FLEXAIDDS_GENERATED_DATASET_CODES_INC
#define FLEXAIDDS_GENERATED_DATASET_CODES_INC

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace flexaids {
namespace generated {

// astex_diverse -- benchmarks/astex_diverse/astex_diverse_set.csv (ROSTER_CSV)
inline const std::vector<std::string>& astex_diverse_codes() {
    static const std::vector<std::string> v = {
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
        "2D3U", "2GBP", "2HB1", "2HR7", "2J62",
    };
    return v;
}

// hap2 -- benchmarks/datasets/hap2.yaml (DEFINITION_YAML)
inline const std::vector<std::string>& hap2_codes() {
    static const std::vector<std::string> v = {
        "1A28", "1A4Q", "1A9M", "1ADB", "1AI5", "1B6M", "1B9V", "1BMA",
        "1C1B", "1C5C", "1C83", "1CBX", "1CIL", "1D3H", "1D4P", "1DBB",
        "1DWD", "1EBY", "1EED", "1ETA", "1ETR", "1F0R", "1F0S", "1FCX",
        "1FEN", "1FKI", "1FL3", "1FPC", "1GKC", "1HPV", "1HTF", "1HWI",
        "1IDA", "1IGJ", "1IMB", "1IVC", "1K1J", "1KZK", "1LAM", "1LPM",
        "1MEH", "1MLD", "1MMV", "1MRK", "1MTS", "1N2V", "1OKL", "1OPK",
        "1OWE", "1PHD", "1POC", "1QPJ", "1RBP", "1STP", "1TLP", "1TMN",
        "1TNI", "1ULB", "1UNL",
    };
    return v;
}

// casf2016 -- benchmarks/datasets/casf2016_core_codes.csv (EXTRACTED_FROM_CODE)
inline const std::vector<std::string>& casf2016_codes() {
    static const std::vector<std::string> v = {
        "1A30", "1B6J", "1B6K", "1BMA", "1C5Z", "1E66", "1EBY", "1F8B",
        "1F8D", "1FEN", "1FKI", "1G2K", "1GKC", "1GNI", "1GNM", "1GPK",
        "1HFS", "1HNN", "1HP0", "1HQ2", "1IA1", "1J3J", "1J4R", "1JD0",
        "1JJE", "1K1J", "1K3U", "1KZK", "1L2S", "1L7F", "1LPZ", "1M2Z",
        "1MQ6", "1N1M", "1N2J", "1N2V", "1N46", "1NAV", "1OF1", "1OF6",
        "1OPK", "1OQ5", "1OWE", "1OYT", "1P2Y", "1P62", "1PMN", "1PSO",
        "1Q1G", "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O", "1S19",
        "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T49", "1T9B",
        "1TT1", "1TW6", "1TZ8", "1U1C", "1U4D", "1UML", "1UNL", "1UOU",
        "1V0P", "1V48", "1V4S", "1VCJ", "1W1P", "1W2G", "1X8X", "1XM6",
        "1XOZ", "1Y6B", "1Y6R", "1YGC", "1YQY", "1YV3", "1YVF", "1YWR",
        "1Z95", "2AL5", "2BM2", "2BR1", "2BSM", "2BYS", "2C3I", "2CET",
        "2CGR", "2D3U", "2FVD", "2G70", "2GBP", "2GQG", "2HB1", "2HR7",
        "2IW1", "2J62", "2J78", "2JDM", "2JDY", "2OBF", "2P4Y", "2PQ9",
        "2QBP", "2QBQ", "2QBR", "2QBS", "2R9W", "2V00", "2VO5", "2VVN",
        "2VW5", "2W66", "2W97", "2WBG", "2WCA", "2WER", "2WHB", "2WN9",
        "2WT2", "2WTV", "2WYG", "2X00", "2X0Y", "2XB8", "2XBV", "2XDL",
        "2XHM", "2XJ7", "2XJJ", "2XNB", "2XYS", "2Y5H", "2YFE", "2YGE",
        "2YLB", "2YMD", "2YPL", "2ZB1", "2ZXD", "3AO4", "3AGN", "3BL1",
        "3BV9", "3CJ4", "3CJ2", "3CKZ", "3CYU", "3D4Z", "3DD0", "3DDQ",
        "3DXG", "3EBP", "3EIG", "3EL1", "3F3A", "3F3C", "3F3D", "3F3E",
        "3FV1", "3FV2", "3GBB", "3GEN", "3GI5", "3GP0", "3GQL", "3GV9",
        "3GVU", "3HUC", "3IAR", "3JVR", "3JVS", "3JY0", "3K5V", "3KGP",
        "3KMZ", "3KR8", "3KWA", "3L3N", "3L4U", "3L4W", "3L7B", "3LKA",
        "3MFV", "3MNA", "3MUZ", "3MY5", "3N7A", "3N86", "3NW9", "3NZK",
        "3OAF", "3OOF", "3OUP", "3OZS", "3OZT", "3P3G", "3P5O", "3PCG",
        "3PE2", "3PFQ", "3PRS", "3PWW", "3QAA", "3QBH", "3QGW", "3QGY",
        "3QQK", "3QTI", "3R88", "3RLQ", "3RT4", "3RUX", "3RYJ", "3S8O",
        "3SXR", "3SYR", "3U5J", "3U5L", "3UAH", "3UAJ", "3UIB", "3UP2",
        "3UPV", "3UTU", "3UWK", "3VD4", "3VF5", "3VHE", "3VRI", "3WMC",
        "3ZSO", "3ZYX", "4AGM", "4AGN", "4AGQ", "4BKT", "4CIG", "4CRA",
        "4CRC", "4DE1", "4DE2", "4DJP", "4DLI", "4E5W", "4EA2", "4EOR",
        "4F09", "4F2W", "4F3C", "4GAM", "4GFM", "4GID", "4GIH", "4GKM",
        "4GR0", "4HGE", "4IQJ", "4IVB", "4IVC", "4IVD", "4J21", "4J28",
        "4JFS", "4JIA", "4JSZ", "4JXS", "4K18", "4K77", "4KAW", "4KEL",
        "4KNE", "4KZ6", "4KZQ",
    };
    return v;
}

/// What a generated list claims to implement: the canonical source it was read
/// from, that source's sha256 at generation time, and how trustworthy the
/// provenance is. Recorded into RUN_RECEIPT.json so a run's executed roster can
/// be checked against its declared source instead of taken on faith.
struct DatasetCodeSource {
    const char* slug;
    const char* cpp_function;
    const char* source_path;
    const char* source_sha256;
    const char* provenance;
    const char* primary_source;
    std::size_t count;
};

inline const std::vector<DatasetCodeSource>& dataset_code_sources() {
    static const std::vector<DatasetCodeSource> v = {
        {"astex_diverse", "DatasetRunner::astex_diverse_codes",
         "benchmarks/astex_diverse/astex_diverse_set.csv",
         "a6384bdda9ba5a73d58284149287f584452805b549b0ea665be944f6e6497708",
         "ROSTER_CSV", "10.1021/jm061277y", 85},
        {"hap2", "DatasetRunner::hap2_codes",
         "benchmarks/datasets/hap2.yaml",
         "8dcf34e6a025093587d066dd8ecc110bcf02d7d48df7ef05b329d8b4f52cbda1",
         "DEFINITION_YAML", "10.1021/acs.jcim.5b00078", 59},
        {"casf2016", "DatasetRunner::casf2016_codes",
         "benchmarks/datasets/casf2016_core_codes.csv",
         "c90f591e09fb3c2edeb808a19ba2d9b4f96b6f5b9bd726009de91d17c88a1b1a",
         "EXTRACTED_FROM_CODE", "UNVERIFIED", 283},
    };
    return v;
}

/// Codes for a dataset slug. THROWS on an unknown slug: an engine that cannot
/// name the declaration it is executing must not execute a list at all.
inline const std::vector<std::string>& dataset_codes(const std::string& slug) {
    if (slug == "astex_diverse") return astex_diverse_codes();
    if (slug == "hap2") return hap2_codes();
    if (slug == "casf2016") return casf2016_codes();
    throw std::runtime_error("dataset_codes: no generated code list declared for slug '"
                             + slug + "' (see scripts/gen_dataset_codes.py REGISTRY)");
}

/// Source record for a dataset slug. THROWS on an unknown slug, same rationale.
inline const DatasetCodeSource& dataset_code_source(const std::string& slug) {
    for (const auto& s : dataset_code_sources()) {
        if (slug == s.slug) return s;
    }
    throw std::runtime_error("dataset_code_source: no generated code list declared for "
                             "slug '" + slug + "'");
}

}  // namespace generated
}  // namespace flexaids

#endif  // FLEXAIDDS_GENERATED_DATASET_CODES_INC
