#pragma once

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace flexaids_rng {

inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline std::atomic<std::uint64_t> g_master_seed{0};
inline std::atomic<bool> g_has_master_seed{false};
inline std::atomic<std::uint64_t> g_seed_epoch{0};

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

inline void set_master_seed(std::uint64_t seed)
{
    g_master_seed.store(seed, std::memory_order_release);
    g_has_master_seed.store(true, std::memory_order_release);
    g_seed_epoch.fetch_add(1, std::memory_order_acq_rel);
}

inline bool has_master_seed()
{
    return g_has_master_seed.load(std::memory_order_acquire);
}

inline std::uint64_t master_seed()
{
    return g_master_seed.load(std::memory_order_acquire);
}

inline void init_from_env()
{
    std::uint64_t base = 0;
    if (env_seed(base)) set_master_seed(base);
}

inline std::uint64_t next_thread_stream()
{
    static std::atomic<std::uint64_t> counter{0};
    return counter.fetch_add(1, std::memory_order_relaxed);
}

inline std::uint64_t deterministic_thread_salt()
{
#ifdef _OPENMP
    return static_cast<std::uint64_t>(omp_get_thread_num());
#else
    thread_local std::uint64_t salt = next_thread_stream();
    return salt;
#endif
}

inline bool seeding_is_deterministic()
{
    if (g_has_master_seed.load(std::memory_order_acquire)) return true;
    std::uint64_t ignored = 0;
    return env_seed(ignored);
}

inline std::uint32_t seed_from_env_or_random(std::uint64_t stream = 0)
{
    if (g_has_master_seed.load(std::memory_order_acquire)) {
        return static_cast<std::uint32_t>(
            splitmix64(g_master_seed.load(std::memory_order_relaxed) ^ stream));
    }

    std::uint64_t base = 0;
    if (env_seed(base)) {
        return static_cast<std::uint32_t>(splitmix64(base ^ stream));
    }

    // Non-deterministic fallback only when no explicit seed was provided.
    std::random_device rd;
    std::uint64_t random_base =
        (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());
    return static_cast<std::uint32_t>(splitmix64(random_base ^ stream));
}

inline std::mt19937 make_rng(std::uint64_t stream = 0)
{
    return std::mt19937(seed_from_env_or_random(stream));
}

inline std::mt19937 make_thread_rng(std::uint64_t stream = 0)
{
    const std::uint64_t salt = seeding_is_deterministic()
        ? (deterministic_thread_salt() << 32)
        : next_thread_stream();
    return std::mt19937(seed_from_env_or_random(stream ^ salt));
}

// One generator per logical stream on this thread. A single cached generator
// that re-seeds on stream switch collapses interleaved streams to a repeating
// first draw (F1). Epoch bump (set_master_seed) still rebuilds the map.
inline std::mt19937& lazy_thread_rng(std::uint64_t stream)
{
    thread_local std::uint64_t cached_epoch = ~0ULL;
    // std::map: insert of a new stream must not invalidate references held
    // by callers (auto& rng = lazy_thread_rng(id)). unordered_map rehash can.
    // Do not hold a reference across set_master_seed(): that bumps g_seed_epoch
    // and this function clears the map.
    thread_local std::map<std::uint64_t, std::mt19937> rngs;

    const std::uint64_t epoch = g_seed_epoch.load(std::memory_order_acquire);
    if (cached_epoch != epoch) {
        rngs.clear();
        cached_epoch = epoch;
    }
    auto it = rngs.find(stream);
    if (it == rngs.end()) {
        it = rngs.emplace(stream, make_thread_rng(stream)).first;
    }
    return it->second;
}

} // namespace flexaids_rng