// =============================================================================
// AsyncPipeline.h — Async I/O pipeline for FlexAIDdS benchmark runner
//
// Producer-consumer pipeline: docking (producer) enqueues I/O tasks that are
// executed by a background worker pool.  This overlaps result writing for
// complex N with docking of complex N+1.
//
// Thread-safe, RAII, backpressure via bounded queue.
// C++20, no GPL dependencies.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace dataset {

/// Async I/O pipeline with bounded queue and graceful shutdown.
///
/// Usage:
///   AsyncPipeline pipe(/*max_queue_depth=*/4, /*num_workers=*/2);
///   pipe.start();
///   for (auto& complex : complexes) {
///       dock(complex);                           // synchronous compute
///       pipe.enqueue([&]() { write_results(...); }); // async I/O
///   }
///   pipe.stop();  // blocks until all queued work completes
///
/// Thread-safe: enqueue() may be called from any thread.
class AsyncPipeline {
public:
    /// Construct a pipeline with given queue capacity and worker count.
    /// @param max_queue_depth  Maximum pending tasks before enqueue() blocks.
    /// @param num_workers      Number of background I/O threads (>= 1).
    explicit AsyncPipeline(size_t max_queue_depth = 4,
                           size_t num_workers = 2);

    /// Destructor: waits for all pending work to complete (RAII).
    ~AsyncPipeline();

    // Non-copyable, non-movable (owns threads + mutexes)
    AsyncPipeline(const AsyncPipeline&) = delete;
    AsyncPipeline& operator=(const AsyncPipeline&) = delete;
    AsyncPipeline(AsyncPipeline&&) = delete;
    AsyncPipeline& operator=(AsyncPipeline&&) = delete;

    /// Start the worker threads.  Must be called before enqueue().
    /// Calling start() on an already-running pipeline is a no-op.
    void start();

    /// Signal shutdown and wait for all pending tasks to finish.
    /// After stop(), enqueue() becomes a no-op and pending() returns 0.
    /// Safe to call multiple times.  Also called by the destructor.
    void stop();

    /// Enqueue an I/O task for background execution.
    /// Blocks if the queue is full (backpressure).
    /// No-op if the pipeline is stopped.
    /// @param task  A callable that performs I/O.  Captured data must be
    ///              owned by value or shared_ptr (no dangling references).
    void enqueue(std::function<void()> task);

    /// Return the number of tasks currently queued (approximate).
    size_t pending() const;

    /// Return true if the pipeline is running (workers active).
    bool running() const { return running_.load(std::memory_order_acquire); }

private:
    void worker_loop();

    size_t max_queue_depth_;
    size_t num_workers_;

    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cv_not_empty_;  // signaled when a task is enqueued
    std::condition_variable cv_not_full_;   // signaled when a task is consumed
    std::atomic<bool> running_{false};
    std::atomic<bool> shutdown_requested_{false};

    // Exception handling: if an I/O task throws, we catch and log.
    // The pipeline does not propagate exceptions across threads; instead
    // it increments an error counter for diagnostics.
    std::atomic<size_t> error_count_{0};
};

} // namespace dataset
