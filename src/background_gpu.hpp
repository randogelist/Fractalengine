#pragma once

#include "amp_gpu.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// Shader compilation and device-only resource creation may run alongside an
// independent CPU renderer. The background thread never submits GPU commands.
// The successful release/acquire publication transfers the prepared renderer
// to its one render owner; that owner alone may use its immediate context.
namespace rhbackground {
using Clock = std::chrono::steady_clock;

class BackgroundGpu {
public:
    using Completion = std::function<void(bool)>;

private:
    enum class Status : unsigned char { Idle, Loading, Ready, Failed, Disabled };
    std::atomic<Status> status_{Status::Idle};
    std::unique_ptr<amp::Renderer> engine_;
    std::thread worker_;
    Completion completion_;
    Clock::time_point startedAt_{}, finishedAt_{};
    double elapsedMs_ = 0.;
    std::string error_;

    void rememberError(const char* message) noexcept {
        try { error_ = message ? message : "Unknown GPU preparation failure"; }
        catch (...) { /* Readiness must still be published if allocation fails. */ }
    }

    void publish(Status result) noexcept {
        finishedAt_ = Clock::now();
        elapsedMs_ = std::chrono::duration<double, std::milli>(finishedAt_ - startedAt_).count();
        status_.store(result, std::memory_order_release);
        // This callback runs on the preparation thread. It must only signal
        // its owner (for example PostMessage), never use or destroy this object.
        try { if (completion_) completion_(result == Status::Ready); }
        catch (...) { /* A notification cannot invalidate a prepared backend. */ }
    }

    void prepare(unsigned threads) noexcept {
        try {
            engine_ = std::make_unique<amp::Renderer>();
            if (!engine_->init(false, threads)) {
                rememberError(engine_->diagnostic.c_str());
                publish(Status::Failed);
                return;
            }
            // prepareRhGpu creates/compiles on ID3D11Device only. It must not
            // render, flush, map, query or otherwise use the immediate context.
            engine_->prepareRhGpu();
            publish(Status::Ready);
        } catch (const std::exception& error) {
            rememberError(error.what());
            publish(Status::Failed);
        } catch (...) {
            rememberError("Unknown GPU preparation failure");
            publish(Status::Failed);
        }
    }

public:
    BackgroundGpu() = default;
    BackgroundGpu(const BackgroundGpu&) = delete;
    BackgroundGpu& operator=(const BackgroundGpu&) = delete;
    BackgroundGpu(BackgroundGpu&&) = delete;
    BackgroundGpu& operator=(BackgroundGpu&&) = delete;
    ~BackgroundGpu() {
        // BackgroundGpu outlives its completion callback. In particular, the
        // callback must not delete it from worker_ itself.
        if (worker_.joinable()) worker_.join();
    }

    // Single-owner, one-shot start. False means CPU-only, launch failure or an
    // already-started object. CPU-only never creates a preparation thread and
    // does not invoke the asynchronous completion callback.
    bool start(bool forceCpu, unsigned threads, Completion completion = {}) noexcept {
        if (status_.load(std::memory_order_acquire) != Status::Idle) return false;
        startedAt_ = Clock::now();
        if (forceCpu) {
            rememberError("CPU requested");
            finishedAt_ = startedAt_;
            status_.store(Status::Disabled, std::memory_order_release);
            return false;
        }
        try {
            completion_ = std::move(completion);
            status_.store(Status::Loading, std::memory_order_release);
            worker_ = std::thread([this, threads] { prepare(threads); });
            return true;
        } catch (const std::exception& error) {
            rememberError(error.what());
        } catch (...) {
            rememberError("Cannot start GPU preparation thread");
        }
        finishedAt_ = Clock::now();
        elapsedMs_ = std::chrono::duration<double, std::milli>(finishedAt_ - startedAt_).count();
        status_.store(Status::Failed, std::memory_order_release);
        return false;
    }

    bool ready() const noexcept {
        const auto value = status_.load(std::memory_order_acquire);
        return value == Status::Ready || value == Status::Failed || value == Status::Disabled;
    }
    bool succeeded() const noexcept {
        return status_.load(std::memory_order_acquire) == Status::Ready;
    }
    bool loading() const noexcept {
        return status_.load(std::memory_order_acquire) == Status::Loading;
    }
    bool skipped() const noexcept {
        return status_.load(std::memory_order_acquire) == Status::Disabled;
    }
    amp::Renderer* renderer() noexcept { return succeeded() ? engine_.get() : nullptr; }
    const amp::Renderer* renderer() const noexcept { return succeeded() ? engine_.get() : nullptr; }

    // Result fields are immutable after terminal publication. Returning empty
    // values before readiness also makes accidental early status polling safe.
    double elapsedMs() const noexcept { return ready() ? elapsedMs_ : 0.; }
    std::string error() const { return ready() ? error_ : std::string{}; }
    Clock::time_point startedAt() const noexcept { return ready() ? startedAt_ : Clock::time_point{}; }
    Clock::time_point finishedAt() const noexcept { return ready() ? finishedAt_ : Clock::time_point{}; }
    bool wait() {
        if (worker_.joinable()) worker_.join();
        return succeeded();
    }
};

struct CompileOverlapResult {
    bool gpuReady = false, cpuCancelled = false;
    unsigned cpuThreads = 1, cpuFrames = 0, cpuFramesBeforeReady = 0;
    double compileMs = 0., cpuWorkMs = 0., overlapMs = 0.;
    std::uint64_t digest = 0;
    std::string error;
};

// Standalone diagnostic. A larger probe should instantiate BackgroundGpu
// directly, then reuse *renderer() for its later GPU runs to avoid recompiling.
// The limit bounds additional CPU frames, not D3DCompile: compiler termination
// is cooperative only on return, and the owned preparation thread is joined.
inline CompileOverlapResult runCompileOverlapProbe(std::ostream& out, const amp::State& source,
                                                   unsigned availableThreads,
                                                   double maxSeconds = 30.) {
    CompileOverlapResult result;
    const unsigned available = availableThreads ? availableThreads :
        std::max(1u, std::thread::hardware_concurrency());
    result.cpuThreads = std::max(1u, std::min(8u, available > 2 ? available - 2 : 1u));
    BackgroundGpu background;
    background.start(false, availableThreads);
    amp::Renderer cpu;
    cpu.init(true, result.cpuThreads);
    auto state = source;
    state.iterations = std::clamp(state.iterations, 1, 32);
    const auto scene = amp::scene(state);
    constexpr int width = 80, height = 48;
    std::vector<std::uint32_t> pixels;
    std::vector<std::pair<Clock::time_point, Clock::time_point>> work;
    std::atomic<std::uint64_t> revision{0};
    const auto deadline = Clock::now() + std::chrono::duration_cast<Clock::duration>(
        std::chrono::duration<double>(std::max(0., maxSeconds)));
    while (!background.ready() && Clock::now() < deadline) {
        const auto begin = Clock::now();
        if (!cpu.render(scene, width, height, double(height) / width, 1, pixels, &revision, 0)) {
            result.cpuCancelled = true;
            break;
        }
        const auto end = Clock::now();
        work.emplace_back(begin, end);
        result.cpuWorkMs += std::chrono::duration<double, std::milli>(end - begin).count();
        ++result.cpuFrames;
        for (auto pixel : pixels) result.digest = result.digest * 1099511628211ull ^ pixel;
    }
    revision.fetch_add(1, std::memory_order_relaxed);
    result.gpuReady = background.wait();
    result.compileMs = background.elapsedMs();
    result.error = background.error();
    const auto compileBegin = background.startedAt(), compileEnd = background.finishedAt();
    for (const auto& interval : work) {
        if (interval.second <= compileEnd) ++result.cpuFramesBeforeReady;
        const auto begin = std::max(interval.first, compileBegin);
        const auto end = std::min(interval.second, compileEnd);
        if (end > begin)
            result.overlapMs += std::chrono::duration<double, std::milli>(end - begin).count();
    }
    out << "BACKGROUND_GPU_READY=" << (result.gpuReady ? "YES" : "NO")
        << " COMPILE_MS=" << result.compileMs << " CPU_WORKERS=" << result.cpuThreads
        << " CPU_FRAMES=" << result.cpuFrames
        << " CPU_FRAMES_BEFORE_READY=" << result.cpuFramesBeforeReady
        << " CPU_WORK_MS=" << result.cpuWorkMs << " HOST_INTERVAL_OVERLAP_MS=" << result.overlapMs
        << "\nBACKGROUND_OVERLAP_MEASURE=HOST_INTERVALS_NOT_GPU_DEVICE_TIME\n";
    if (!result.error.empty()) out << "BACKGROUND_GPU_ERROR=" << result.error << '\n';
    return result;
}
} // namespace rhbackground
