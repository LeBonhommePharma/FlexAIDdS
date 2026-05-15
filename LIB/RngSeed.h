#pragma once

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <random>

namespace flexaids_rng {

inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline bool env_seed(std::uint64_t& seed)
{
    const char* raw = std::getenv("FLEXAID_SEED");
    if (!raw || !*raw) return false;

    char* end = nullptr;
    unsigned long long parsed = std::strtoull(raw, &end, 10);
    if (end == raw) return false;

    seed = static_cast<std::uint64_t>(parsed);
    return true;
}

inline std::uint32_t seed_from_env_or_random(std::uint64_t stream = 0)
{
    std::uint64_t base = 0;
    if (env_seed(base)) {
        return static_cast<std::uint32_t>(splitmix64(base ^ stream));
    }

    std::random_device rd;
    std::uint64_t random_base =
        (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());
    return static_cast<std::uint32_t>(splitmix64(random_base ^ stream));
}

inline std::uint64_t next_thread_stream()
{
    static std::atomic<std::uint64_t> counter{0};
    return counter.fetch_add(1, std::memory_order_relaxed);
}

inline std::mt19937 make_rng(std::uint64_t stream = 0)
{
    return std::mt19937(seed_from_env_or_random(stream));
}

inline std::mt19937 make_thread_rng(std::uint64_t stream = 0)
{
    return std::mt19937(seed_from_env_or_random(stream ^ next_thread_stream()));
}

} // namespace flexaids_rng
