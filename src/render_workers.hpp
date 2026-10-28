#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

// A serialized render stage reuses its helpers across frames. Work owns its
// row counter and cancellation policy; run() never outlives its captures.
namespace renderworkers {
class Pool {
    std::mutex serial_, mutex_;
    std::condition_variable ready_, done_;
    std::vector<std::thread> threads_;
    std::function<void()> task_;
    std::exception_ptr error_;
    unsigned active_ = 1, remaining_ = 0;
    std::atomic<unsigned> lastActive_{1};
    unsigned long long epoch_ = 0;
    bool stopping_ = false;

    void invoke(const std::function<void()>& task) noexcept {
        try { task(); }
        catch (...) {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!error_) error_ = std::current_exception();
        }
    }
    void worker(unsigned index, unsigned long long seen) {
        std::unique_lock<std::mutex> lock(mutex_);
        for (;;) {
            ready_.wait(lock, [&] { return stopping_ || epoch_ != seen; });
            if (stopping_) return;
            seen = epoch_;
            if (index >= active_) continue;
            // task_ is cleared only after every participating helper returns.
            lock.unlock();
            invoke(task_);
            lock.lock();
            if (--remaining_ == 0) done_.notify_one();
        }
    }

public:
    Pool() = default;
    Pool(const Pool&) = delete;
    Pool& operator=(const Pool&) = delete;
    unsigned lastWorkerCount() const noexcept { return lastActive_.load(std::memory_order_relaxed); }
    ~Pool() {
        std::lock_guard<std::mutex> serial(serial_);
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (auto& thread : threads_) thread.join();
    }

    // Calls from different owners serialize. Do not recursively call the same
    // pool from a work function. Failed helper creation reduces parallelism.
    template<class Work> void run(unsigned totalWorkers, Work work) {
        std::lock_guard<std::mutex> serial(serial_);
        const unsigned requested = std::clamp(totalWorkers, 1u, 64u);
        std::function<void()> task(std::move(work));
        {
            std::lock_guard<std::mutex> lock(mutex_);
            while (threads_.size() + 1 < requested) {
                const unsigned index = unsigned(threads_.size()) + 1;
                const auto seen = epoch_;
                try { threads_.emplace_back([this, index, seen] { worker(index, seen); }); }
                catch (const std::system_error&) { break; }
            }
            active_ = std::min(requested, unsigned(threads_.size()) + 1);
            lastActive_.store(active_, std::memory_order_relaxed);
            remaining_ = active_ - 1;
            error_ = nullptr;
            task_ = std::move(task);
            ++epoch_;
        }
        ready_.notify_all();
        invoke(task_);
        std::exception_ptr failure;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            done_.wait(lock, [&] { return remaining_ == 0; });
            task_ = {};
            failure = error_;
        }
        if (failure) std::rethrow_exception(failure);
    }
};
} // namespace renderworkers
