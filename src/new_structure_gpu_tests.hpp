#pragma once

#include "amp_gpu.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <ostream>
#include <string>
#include <vector>

// Native hardware checks for the five spectral/projective mathematical drivers.
// Correct precision repairs are permitted at individual float-conditioned
// points; every driver must also demonstrate real unrepaired GPU coverage.
namespace new_structure_gpu_tests {

inline int run(std::ostream& out) {
    try {
        amp::Renderer warm, cold, cpu;
        if (!warm.init(false, 1)) {
            const bool unavailable = warm.diagnostic.find("Direct3D 11 device") != std::string::npos;
            out << "NEW_STRUCTURE_GPU_TEST=" << (unavailable ? "UNAVAILABLE" : "FAIL")
                << " STAGE=INITIALIZE REASON=" << warm.diagnostic << '\n';
            return unavailable ? 2 : 1;
        }
        if (!cold.init(false, 1)) {
            out << "NEW_STRUCTURE_GPU_TEST=FAIL STAGE=REFERENCE_INITIALIZE REASON="
                << cold.diagnostic << '\n';
            return 1;
        }
        cpu.init(true, 1);
        constexpr int width = 8, height = 4, tolerance = 8;
        constexpr std::size_t pixels = std::size_t(width) * height;
        int cases = 0, failures = 0, worstError = 0, worstCacheError = 0;
        std::array<std::size_t, 5> rawByCandidate{{0, 0, 0, 0, 0}};
        std::size_t repairsTotal = 0;
        auto base = [](int generation, int order = 2) {
            auto s = amp::defaults(20, generation);
            s.plane = 1;
            s.iterations = 1;
            s.spectralOrder = order;
            s.c = amp::C(-.3, .1);
            s.alpha = amp::C(.15, .03);
            s.palette = 0;
            s.hue = .2;
            s.bands = 1.;
            s.rhExperiment = rhmaster::Experiment{};
            s.rhExperiment.map = 1;
            s.rhExperiment.coloring = 1;
            for (auto& camera : s.cameras) camera = amp::Camera{.31, .27, .02};
            return s;
        };
        auto error = [](const std::vector<std::uint32_t>& a,
                        const std::vector<std::uint32_t>& b) {
            if (a.size() != b.size()) return 256;
            int maximum = 0;
            for (std::size_t i = 0; i < a.size(); ++i)
                for (int shift : {0, 8, 16}) {
                    const int x = int((a[i] >> shift) & 255u), y = int((b[i] >> shift) & 255u);
                    maximum = std::max(maximum, x > y ? x - y : y - x);
                }
            return maximum;
        };
        auto compare = [&](amp::Renderer& renderer, const amp::State& state,
                           const std::string& label, int expectedHit,
                           std::vector<std::uint32_t>& image,
                           const std::vector<std::uint32_t>* cached = nullptr) {
            ++cases;
            const auto scene = amp::scene(state);
            const bool gpuOk = renderer.render(scene, width, height, .5, 1, image);
            std::vector<std::uint32_t> cpuImage;
            const bool cpuOk = cpu.render(scene, width, height, .5, 1, cpuImage);
            const auto metrics = renderer.rhGpuMetrics();
            const std::size_t repairs = renderer.rhCpuRepairs();
            const bool shape = image.size() == pixels && cpuImage.size() == pixels;
            const bool sentinel = std::find(image.begin(), image.end(), 0xffffffffu) != image.end()
                || std::find(cpuImage.begin(), cpuImage.end(), 0xffffffffu) != cpuImage.end();
            const int cpuError = shape ? error(image, cpuImage) : 256;
            const int cacheError = cached ? error(image, *cached) : 0;
            worstError = std::max(worstError, cpuError);
            worstCacheError = std::max(worstCacheError, cacheError);
            const bool good = gpuOk && cpuOk && shape && !sentinel && repairs <= pixels
                && renderer.rhUsedGpu() == (repairs < pixels)
                && metrics.dispatches > 0 && metrics.metricsCached
                && (expectedHit < 0 || metrics.cacheHit == bool(expectedHit))
                && cpuError <= tolerance && (!cached || cacheError <= 1);
            if (good) {
                rawByCandidate[std::size_t(state.generation - 14)] += pixels - repairs;
                repairsTotal += repairs;
            } else {
                ++failures;
                out << "NEW_STRUCTURE_GPU_TEST_FAIL=" << label
                    << " CANDIDATE=" << state.generation
                    << " DEGREE=" << (8 << (std::clamp(state.spectralOrder, 1, 3) - 1))
                    << " MAP=" << state.rhExperiment.map << " CHART=" << state.rhExperiment.chart
                    << " USED_GPU=" << renderer.rhUsedGpu() << " REPAIRS=" << repairs
                    << " GPU_RENDER=" << gpuOk << " CPU_RENDER=" << cpuOk
                    << " GPU_PIXELS=" << image.size() << " CPU_PIXELS=" << cpuImage.size()
                    << " CACHE=" << metrics.cacheHit << " EXPECTED_CACHE=" << expectedHit
                    << " METRICS_CACHED=" << metrics.metricsCached
                    << " MAX_RGB_ERROR=" << cpuError << " TOLERANCE=" << tolerance
                    << " CACHE_RGB_ERROR=" << cacheError
                    << " DIAGNOSTIC=" << renderer.rhGpuDiagnostic() << '\n';
            }
        };
        std::vector<std::uint32_t> image, reference;
        for (int generation = 14; generation <= 18; ++generation)
            for (int order = 1; order <= 3; ++order)
                compare(warm, base(generation, order), "DRIVER_DEGREE", 0, image);
        for (int generation = 14; generation <= 18; ++generation) {
            for (int map : {2, 5}) {
                auto s = base(generation);
                s.rhExperiment.map = map;
                s.rhExperiment.relaxation = .7;
                compare(warm, s, map == 2 ? "TRICORN" : "CAYLEY", 0, image);
            }
            auto s = base(generation);
            s.rhExperiment.chart = 1;
            for (auto& camera : s.cameras) camera = amp::Camera{.7, .4, .02};
            compare(warm, s, "RECIPROCAL_CHART", 0, image);
        }
        // Appearance remains outside the orbit key. A freshly initialized
        // second renderer must reproduce the warm cache within quantisation.
        for (int generation = 14; generation <= 18; ++generation) {
            auto s = base(generation);
            s.rhExperiment.map = 5;
            s.rhExperiment.gain = .8;
            compare(warm, s, "CACHE_WARMUP", 0, image);
            s.palette = 2;
            s.hue = .41;
            s.bands = 1.3;
            s.rhExperiment.exposure = 1.1;
            s.rhExperiment.frequency = 1.2;
            compare(warm, s, "APPEARANCE_CACHE_HIT", 1, image);
            compare(cold, s, "COLD_MATCHES_RECOLOUR", 0, reference, &image);
        }
        for (int generation = 14; generation <= 18; ++generation) {
            ++cases;
            const std::size_t raw = rawByCandidate[std::size_t(generation - 14)];
            out << "NEW_STRUCTURE_GPU_CANDIDATE=" << generation << " RAW_GPU_PIXELS=" << raw << '\n';
            if (raw == 0) {
                ++failures;
                out << "NEW_STRUCTURE_GPU_TEST_FAIL=NO_RAW_GPU_COVERAGE CANDIDATE=" << generation << '\n';
            }
        }
        out << "NEW_STRUCTURE_GPU_TEST_CASES=" << cases << " FAILURES=" << failures
            << " MAX_RGB_ERROR=" << worstError << " TOLERANCE=" << tolerance
            << " MAX_CACHE_RGB_ERROR=" << worstCacheError << " CPU_REPAIRS=" << repairsTotal << '\n';
        out << "NEW_STRUCTURE_GPU_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& exception) {
        out << "NEW_STRUCTURE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=" << exception.what() << '\n';
        return 1;
    } catch (...) {
        out << "NEW_STRUCTURE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return 1;
    }
}
} // namespace new_structure_gpu_tests
