/// @file SharedPosePool.cpp
/// @brief Thread-safe best-pose sharing pool implementation.
///
/// Copyright 2026 Le Bonhomme Pharma
/// SPDX-License-Identifier: Apache-2.0

#include "SharedPosePool.h"
#include <cstring>
#include <cmath>
#include <algorithm>

SharedPosePool::SharedPosePool(int pool_size)
    : capacity_(pool_size), used_(0)
{
    if (pool_size <= 0) {
        throw std::invalid_argument(
            "SharedPosePool: pool_size must be > 0 (got " + std::to_string(pool_size) + ")");
    }
    pool_.resize(static_cast<size_t>(capacity_));
}

SharedPosePool::SharedPosePool(SharedPosePool&& other) noexcept
    : capacity_(0), used_(0)
{
    std::lock_guard<std::mutex> lock(other.mtx_);
    capacity_ = other.capacity_;
    pool_ = std::move(other.pool_);
    used_ = other.used_;
    other.capacity_ = 0;
    other.used_ = 0;
}

SharedPosePool& SharedPosePool::operator=(SharedPosePool&& other) noexcept
{
    if (this == &other) return *this;
    // Lock both mutexes in address order to prevent deadlock.
    std::mutex* first = &mtx_;
    std::mutex* second = &other.mtx_;
    if (first > second) std::swap(first, second);
    std::lock_guard<std::mutex> lock1(*first);
    std::lock_guard<std::mutex> lock2(*second);

    capacity_ = other.capacity_;
    pool_ = std::move(other.pool_);
    used_ = other.used_;
    other.capacity_ = 0;
    other.used_ = 0;
    return *this;
}

// ── Core operations ────────────────────────────────────────────────────

bool SharedPosePool::publish(const SharedPose& pose)
{
    // Reject NaN energies — they would corrupt sort order.
    if (!std::isfinite(pose.energy)) return false;

    std::lock_guard<std::mutex> lock(mtx_);

    if (used_ < capacity_) {
        insert_sorted_unlocked(pose);
        return true;
    }

    // Pool full: only accept if strictly better than worst.
    assert(used_ > 0);
    if (pose.energy < pool_[static_cast<size_t>(used_ - 1)].energy) {
        // Overwrite worst entry, then bubble it into sorted position.
        pool_[static_cast<size_t>(used_ - 1)] = pose;
        for (int i = used_ - 1; i > 0; --i) {
            const auto ci = static_cast<size_t>(i);
            if (pool_[ci].energy < pool_[ci - 1].energy) {
                std::swap(pool_[ci], pool_[ci - 1]);
            } else {
                break;
            }
        }
        return true;
    }

    return false;  // pose not good enough
}

std::vector<SharedPose> SharedPosePool::get_top(int k) const
{
    std::lock_guard<std::mutex> lock(mtx_);
    const int n = std::clamp(k, 0, used_);
    return {pool_.begin(), pool_.begin() + n};
}

int SharedPosePool::count() const noexcept
{
    std::lock_guard<std::mutex> lock(mtx_);
    return used_;
}

bool SharedPosePool::is_full() const noexcept
{
    std::lock_guard<std::mutex> lock(mtx_);
    return used_ >= capacity_;
}

double SharedPosePool::best_energy() const noexcept
{
    std::lock_guard<std::mutex> lock(mtx_);
    if (used_ == 0) return std::numeric_limits<double>::infinity();
    return pool_[0].energy;
}

double SharedPosePool::worst_energy() const noexcept
{
    std::lock_guard<std::mutex> lock(mtx_);
    if (used_ == 0) return std::numeric_limits<double>::infinity();
    return pool_[static_cast<size_t>(used_ - 1)].energy;
}

void SharedPosePool::clear()
{
    std::lock_guard<std::mutex> lock(mtx_);
    used_ = 0;
    // Reset entries to default state for cleanliness.
    for (auto& p : pool_) p = SharedPose();
}

// ── Sorted insertion (binary search) ───────────────────────────────────

void SharedPosePool::insert_sorted_unlocked(const SharedPose& pose)
{
    assert(used_ < capacity_);

    // Binary search for the insertion point (ascending energy order).
    int lo = 0, hi = used_;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (pool_[static_cast<size_t>(mid)].energy <= pose.energy) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    // lo is the insertion index.

    // Shift elements right to make room.
    for (int i = used_; i > lo; --i) {
        pool_[static_cast<size_t>(i)] = pool_[static_cast<size_t>(i - 1)];
    }
    pool_[static_cast<size_t>(lo)] = pose;
    ++used_;
}

// ── Serialization for MPI ──────────────────────────────────────────────

std::vector<char> SharedPosePool::serialize() const
{
    std::lock_guard<std::mutex> lock(mtx_);

    // Format: [int32 used_count][SharedPose × used_count]
    const size_t payload = static_cast<size_t>(used_) * sizeof(SharedPose);
    const size_t total = sizeof(int) + payload;
    std::vector<char> buf(total);
    char* ptr = buf.data();

    std::memcpy(ptr, &used_, sizeof(int));
    ptr += sizeof(int);

    if (used_ > 0) {
        std::memcpy(ptr, pool_.data(), payload);
    }

    return buf;
}

void SharedPosePool::deserialize_merge(const char* data, size_t len)
{
    if (!data || len < sizeof(int)) return;

    int remote_count = 0;
    std::memcpy(&remote_count, data, sizeof(int));
    data += sizeof(int);

    // Validate: reject negative counts and check for integer overflow.
    if (remote_count <= 0) return;

    const size_t expected_payload = static_cast<size_t>(remote_count) * sizeof(SharedPose);
    if (expected_payload / sizeof(SharedPose) != static_cast<size_t>(remote_count)) {
        return;  // integer overflow in size computation
    }
    if (len < sizeof(int) + expected_payload) {
        return;  // truncated buffer
    }

    // Merge each remote pose into this pool via publish().
    for (int i = 0; i < remote_count; ++i) {
        SharedPose pose;
        std::memcpy(&pose, data + static_cast<size_t>(i) * sizeof(SharedPose),
                     sizeof(SharedPose));
        publish(pose);
    }
}
