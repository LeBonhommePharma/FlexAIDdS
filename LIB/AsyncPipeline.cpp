// =============================================================================
// AsyncPipeline.cpp — Async I/O pipeline implementation for FlexAIDdS
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "AsyncPipeline.h"

#include <iostream>

namespace dataset {

// =============================================================================
// Construction / destruction
// =============================================================================

AsyncPipeline::AsyncPipeline(size_t max_queue_depth, size_t num_workers)
    : max_queue_depth_(max_queue_depth > 0 ? max_queue_depth : 1)
    , num_workers_(num_workers > 0 ? num_workers : 1)
{
}

AsyncPipeline::~AsyncPipeline()
{
    stop();
}

// =============================================================================
// Lifecycle
// =============================================================================

void AsyncPipeline::start()
{
    bool expected = false;
    if (!running_.compare_exchange_strong(expected, true)) {
        return;  // already running
    }

    shutdown_requested_.store(false, std::memory_order_release);
    error_count_.store(0, std::memory_order_release);

    workers_.reserve(num_workers_);
    for (size_t i = 0; i < num_workers_; ++i) {
        workers_.emplace_back(&AsyncPipeline::worker_loop, this);
    }
}

void AsyncPipeline::stop()
{
    // Set shutdown flag so workers drain the queue then exit.
    shutdown_requested_.store(true, std::memory_order_release);
    cv_not_empty_.notify_all();
    cv_not_full_.notify_all();

    // Join all worker threads.
    for (auto& w : workers_) {
        if (w.joinable()) {
            w.join();
        }
    }
    workers_.clear();

    running_.store(false, std::memory_order_release);

    if (error_count_.load(std::memory_order_acquire) > 0) {
        std::cerr << "[AsyncPipeline] Completed with "
                  << error_count_.load() << " I/O error(s).\n";
    }
}

// =============================================================================
// Producer interface
// =============================================================================

void AsyncPipeline::enqueue(std::function<void()> task)
{
    if (!task) return;
    if (shutdown_requested_.load(std::memory_order_acquire)) return;

    {
        std::unique_lock<std::mutex> lock(mutex_);
        // Backpressure: block while the queue is full and we are still running.
        cv_not_full_.wait(lock, [&] {
            return queue_.size() < max_queue_depth_ ||
                   shutdown_requested_.load(std::memory_order_acquire);
        });

        if (shutdown_requested_.load(std::memory_order_acquire)) {
            return;  // dropped — pipeline shutting down
        }

        queue_.push(std::move(task));
    }

    cv_not_empty_.notify_one();
}

size_t AsyncPipeline::pending() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return queue_.size();
}

// =============================================================================
// Worker implementation
// =============================================================================

void AsyncPipeline::worker_loop()
{
    for (;;) {
        std::function<void()> task;

        {
            std::unique_lock<std::mutex> lock(mutex_);
            cv_not_empty_.wait(lock, [&] {
                return !queue_.empty() ||
                       shutdown_requested_.load(std::memory_order_acquire);
            });

            // Drain: keep processing even after shutdown is requested.
            if (queue_.empty() && shutdown_requested_.load(std::memory_order_acquire)) {
                return;  // no more work — worker exits
            }

            if (!queue_.empty()) {
                task = std::move(queue_.front());
                queue_.pop();
            }
        }

        if (task) {
            cv_not_full_.notify_one();  // unblock a blocked producer

            try {
                task();
            } catch (const std::exception& e) {
                ++error_count_;
                std::cerr << "[AsyncPipeline] I/O task exception: "
                          << e.what() << "\n";
            } catch (...) {
                ++error_count_;
                std::cerr << "[AsyncPipeline] Unknown I/O task exception.\n";
            }
        }
    }
}

} // namespace dataset
