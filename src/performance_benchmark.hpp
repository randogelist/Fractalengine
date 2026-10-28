#pragma once
#include "amp_gpu.hpp"
#include "focus_presets.hpp"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace performance_benchmark {
inline int run(std::ostream& out, bool forceCpu, unsigned threads,
               int width, int height, const std::filesystem::path& destination,
               const rhgpu::EngineTuning& tuning=rhgpu::EngineTuning{}) {
    if (width < 1 || height < 1 || std::uint64_t(width) * height > 16000000ULL)
        throw std::runtime_error("Benchmark dimensions must be positive and at most 16 million pixels");
    amp::Renderer renderer;
    renderer.init(forceCpu, threads);
    renderer.setGpuTuning(tuning);
    std::ofstream csv(destination);
    if (!csv) throw std::runtime_error("Cannot write performance.csv");
    csv << "preset,pass,width,height,iterations,total_ms,prepare_ms,dispatch_ms,readback_ms,repair_ms,dispatches,cache_hit,cpu_repairs,gpu,surface_allocations\n";
    out << "BENCHMARK=" << width << 'x' << height << " ITERATIONS=64\n"
        << "PRESET PASS TOTAL_MS DISPATCHES CACHE_HIT CPU_REPAIRS\n";
    for (int bank : {1, 2}) {
        auto state = amp::defaults(20, 1);
        focus_presets::apply(state, bank, 0);
        state.iterations = 64;
        const char* name = bank == 1 ? "Cayley_currents" : "Tricorn_weave";
        for (int pass = 0; pass < 4; ++pass) {
            const char* label = pass == 0 ? "orbit" : pass == 1 ? "repeat" : pass == 2 ? "recolour" : "dynamics";
            if (pass == 2) {
                state.hue += .17;
                state.rhExperiment.exposure = .8;
                state.rhExperiment.frequency += .4;
            }
            if (pass == 3) state.rhExperiment.gain = .85;
            const auto start = std::chrono::steady_clock::now();
            const auto scene = amp::scene(state);
            std::vector<std::uint32_t> pixels;
            if (!renderer.render(scene, width, height, double(height) / width, 1, pixels))
                throw std::runtime_error("Benchmark render was cancelled");
            const double milliseconds = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
            const auto metrics = renderer.rhGpuMetrics();
            csv << name << ',' << label << ',' << width << ',' << height << ',' << state.iterations << ','
                << std::fixed << std::setprecision(3) << milliseconds << ',' << metrics.prepareMs << ','
                << metrics.dispatchMs << ',' << metrics.readbackMs << ',' << metrics.repairMs << ','
                << metrics.dispatches << ',' << metrics.cacheHit << ',' << renderer.rhCpuRepairs() << ','
                << renderer.rhUsedGpu() << ',' << metrics.surfaceAllocations << '\n';
            out << name << ' ' << label << ' ' << std::fixed << std::setprecision(2) << milliseconds << ' '
                << metrics.dispatches << ' ' << (metrics.cacheHit ? "YES" : "NO") << ' '
                << renderer.rhCpuRepairs() << '\n';
        }
    }
    csv.flush();
    if (!csv) throw std::runtime_error("Cannot finish performance.csv");
    out << "BENCHMARK_SAVED=performance.csv\n";
    return 0;
}
} // namespace performance_benchmark
