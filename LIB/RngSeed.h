#pragma once

#include "EnvFlags.h"  // flexaids::env_bool — one parser for FLEXAIDDS_* switches

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
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

// FLEXAIDDS_RNG_STREAM_FIX — DEFAULT OFF.
//
// OFF reproduces the historical single-generator behaviour bit-for-bit. That
// behaviour is defective: three streams interleave on one thread (GA 0x9A800D,
// Vcontacts 0x0C0A11 inside chromosome evaluation, FOPTICS 0xF0701C5) and the
// generator re-seeds on every stream switch, so each stream collapses to its
// first draw forever. It is nonetheless the behaviour every frozen reference
// number in this repo was produced under, so it remains the default until a
// baseline re-run retires it. Enabling this flag changes the draw sequence for
// a given FLEXAID_SEED and therefore invalidates comparison against those
// numbers.
//
// Parsed by flexaids::env_bool (EnvFlags.h), the repo's one parser for these
// switches. A hand-rolled "not 0/n/f" test would read FLEXAIDDS_RNG_STREAM_FIX=off
// as ON — a knob documented OFF silently being ON is the exact failure this
// gate exists to prevent.
//
// The flag is re-read on seed-epoch change rather than on every call: the hot
// path stays an atomic load plus a compare, and a test can flip the variable
// and call set_master_seed() to pick it up.
/// Voronoi hull-failure jitter source. DEFAULT OFF, independent of
/// FLEXAIDDS_RNG_STREAM_FIX so an A/B of either gate stays single-variable.
inline bool voronoi_keyed_jitter_enabled() noexcept
{
    return flexaids::env_bool("FLEXAIDDS_VORONOI_KEYED_JITTER", false);
}

inline bool rng_stream_fix_enabled()
{
    thread_local std::uint64_t flag_epoch = ~0ULL;
    thread_local bool flag = false;

    const std::uint64_t epoch = g_seed_epoch.load(std::memory_order_acquire);
    if (flag_epoch != epoch) {
        const char* raw = std::getenv("FLEXAIDDS_RNG_STREAM_FIX"); // DEFAULT OFF
        flag = (raw != nullptr && raw[0] != '\0')
            ? flexaids::env_bool("FLEXAIDDS_RNG_STREAM_FIX", false)
            : false;
        flag_epoch = epoch;
    }
    return flag;
}

inline std::mt19937& lazy_thread_rng(std::uint64_t stream)
{
    const std::uint64_t epoch = g_seed_epoch.load(std::memory_order_acquire);

    if (!rng_stream_fix_enabled()) {
        // ---- LEGACY PATH (default) — byte-identical to the pre-fix code. ----
        // One generator per thread, re-seeded whenever the requested stream
        // differs from the cached one. Do not "improve" this branch: its exact
        // draw sequence is the project baseline.
        thread_local std::uint64_t cached_stream = ~0ULL;
        thread_local std::uint64_t cached_epoch  = ~0ULL;
        thread_local std::mt19937 rng = make_thread_rng(stream);

        if (cached_stream != stream || cached_epoch != epoch) {
            rng = make_thread_rng(stream);
            cached_stream = stream;
            cached_epoch = epoch;
        }
        return rng;
    }

    // ---- STREAM-FIX PATH (opt-in) — one generator per logical stream. ----
    // std::map: insert of a new stream must not invalidate references held
    // by callers (auto& rng = lazy_thread_rng(id)). unordered_map rehash can.
    // Do not hold a reference across set_master_seed(): that bumps g_seed_epoch
    // and this function clears the map.
    thread_local std::uint64_t fix_cached_epoch = ~0ULL;
    thread_local std::map<std::uint64_t, std::mt19937> rngs;

    if (fix_cached_epoch != epoch) {
        rngs.clear();
        fix_cached_epoch = epoch;
    }
    auto it = rngs.find(stream);
    if (it == rngs.end()) {
        it = rngs.emplace(stream, make_thread_rng(stream)).first;
    }
    return it->second;
}

// Pose/atom identity for Voronoi hull-failure jitter (F2).
// Quantized from the atom's integer PDB number and the three coordinates
// *before* perturbation so two threads that fail the hull on the same pose
// apply the same displacement. Independent of omp_get_thread_num().
inline std::uint64_t pose_atom_identity(int atom_number,
                                        float x, float y, float z)
{
    std::uint32_t bx = 0, by = 0, bz = 0;
    std::memcpy(&bx, &x, sizeof(bx));
    std::memcpy(&by, &y, sizeof(by));
    std::memcpy(&bz, &z, sizeof(bz));
    const std::uint64_t id =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(atom_number)) << 32)
        ^ static_cast<std::uint64_t>(bx)
        ^ (static_cast<std::uint64_t>(by) << 11)
        ^ (static_cast<std::uint64_t>(bz) << 22);
    return splitmix64(id);
}

// Deterministic jitter in [-amp, amp] for axis 0/1/2. Uses the master /
// FLEXAID_SEED when set; does not consume lazy_thread_rng (so it cannot
// race with GA/FOPTICS streams or depend on draw order).
inline float keyed_jitter(std::uint64_t identity, int axis, float amp = 0.005f)
{
    std::uint64_t seed = 0;
    if (g_has_master_seed.load(std::memory_order_acquire))
        seed = g_master_seed.load(std::memory_order_relaxed);
    else
        (void)env_seed(seed);
    const std::uint64_t h = splitmix64(
        seed ^ identity ^ (0xA11E0000ULL + static_cast<std::uint64_t>(axis)));
    const float u = static_cast<float>(
        static_cast<std::uint32_t>(h >> 40) * (1.0 / 16777216.0));
    return (2.f * u - 1.f) * amp;
}

} // namespace flexaids_rng