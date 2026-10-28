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

// Native, user-side gate for the optional curvature polynomial evaluator.
// These are short smooth orbits, not a pixel-equality claim about long
// chaotic trajectories. The default logarithmic-jet evaluator remains off
// this experimental path unless EngineTuning::stableCurvature is enabled.
namespace stable_curvature_gpu_tests {

inline int run(std::ostream& out) {
    try {
        amp::Renderer warm, cold, cpu;
        if (!warm.init(false, 1)) {
            const bool unavailable = warm.diagnostic.find("Direct3D 11 device") != std::string::npos;
            out << "STABLE_CURVATURE_GPU_TEST=" << (unavailable ? "UNAVAILABLE" : "FAIL")
                << " STAGE=INITIALIZE REASON=" << warm.diagnostic << '\n';
            return unavailable ? 2 : 1;
        }
        if (!cold.init(false, 1)) {
            out << "STABLE_CURVATURE_GPU_TEST=FAIL STAGE=REFERENCE_INITIALIZE REASON="
                << cold.diagnostic << '\n';
            return 1;
        }
        cpu.init(true, 1);
        auto stableTuning = warm.gpuTuning();
        stableTuning.stableCurvature = true;
        warm.setGpuTuning(stableTuning);
        cold.setGpuTuning(stableTuning);
        constexpr int width = 8, height = 4, tolerance = 8;
        constexpr std::size_t pixelCount = std::size_t(width) * height;
        int cases = 0, failures = 0, worstError = 0, worstPeerError = 0;
        std::size_t rawGpuPixels = 0, repairedPixels = 0;

        auto base = [] {
            auto s = amp::defaults(20, 1);
            s.plane = 1;
            s.iterations = 1;
            s.spectralOrder = 1;
            // Set y=z explicitly so cameras below exercise |y|/4, rather
            // than silently moving the evaluator boundary via model scale.
            s.seed.imag(1.);
            s.spectralD = 0.;
            s.c = amp::C(-.3, .1);
            s.alpha = amp::C(.15, .03);
            s.palette = 0;
            s.hue = .2;
            s.bands = 1.;
            s.rhExperiment = rhmaster::Experiment{};
            s.rhExperiment.coloring = 1;
            for (auto& camera : s.cameras) camera = amp::Camera{.3, .2, .015625};
            return s;
        };

        std::array<bool, 3> ready{{false, false, false}};
        for (int order = 1; order <= 3; ++order) {
            const int degree = 8 << (order - 1);
            const auto payload = rhgpu::makeStableCurvatureData(degree);
            ready[std::size_t(order - 1)] = payload.valid;
            out << "STABLE_CURVATURE_GPU_PAYLOAD_DEGREE=" << degree
                << " READY=" << (payload.valid ? "YES" : "NO");
            if (!payload.valid) out << " REASON=" << payload.reason;
            out << '\n';
        }

        auto rgbError = [](const std::vector<std::uint32_t>& a,
                           const std::vector<std::uint32_t>& b) {
            if (a.size() != b.size()) return 256;
            int error = 0;
            for (std::size_t i = 0; i < a.size(); ++i)
                for (int shift : {0, 8, 16}) {
                    const int x = int((a[i] >> shift) & 255u);
                    const int y = int((b[i] >> shift) & 255u);
                    error = std::max(error, x > y ? x - y : y - x);
                }
            return error;
        };
        auto compare = [&](amp::Renderer& renderer, const amp::State& state,
                           const std::string& label, int samples, int expectedHit,
                           std::vector<std::uint32_t>& image,
                           const std::vector<std::uint32_t>* peer = nullptr,
                           int peerTolerance = 1, int expectedRepairs = -1) {
            ++cases;
            const auto scene = amp::scene(state);
            const bool gpuOk = renderer.render(scene, width, height, .5, samples, image);
            std::vector<std::uint32_t> cpuImage;
            const bool cpuOk = cpu.render(scene, width, height, .5, samples, cpuImage);
            const auto metrics = renderer.rhGpuMetrics();
            const std::size_t repairs = renderer.rhCpuRepairs();
            const int order = std::clamp(state.spectralOrder, 1, 3);
            const bool stable = renderer.gpuTuning().stableCurvature;
            const auto data = rhgpu::makeData(state);
            const std::size_t metadata = std::size_t(rhgpu::CurvatureMetaRow)
                * rhgpu::CoefficientWidth + rhgpu::CurvatureMetaColumn;
            const bool markerReady = data.valid && data.texels.size() > metadata
                && data.texels[metadata][0] == 1.f;
            const bool payloadReady = ready[std::size_t(order - 1)];
            const bool payloadConsistent = data.valid && markerReady == payloadReady;
            // Unsupported optional payloads must use actual CPU repairs,
            // rather than reporting an untested all-CPU image as GPU work.
            const bool rawRequired = !stable || payloadReady || state.alpha == amp::C{};
            // Only a fixture explicitly requiring every precision marker
            // may accept a fully repaired frame with a ready payload. All
            // ordinary fixtures still require actual unrepaired GPU pixels.
            const bool allRepairsExpected = expectedRepairs == int(pixelCount);
            const bool rawValid = allRepairsExpected
                ? stable && payloadReady && !renderer.rhUsedGpu() && repairs == pixelCount
                : (rawRequired ? renderer.rhUsedGpu() && repairs < pixelCount
                    : !renderer.rhUsedGpu() && repairs == pixelCount);
            const bool shape = image.size() == pixelCount && cpuImage.size() == pixelCount;
            const bool sentinel = std::find(image.begin(), image.end(), 0xffffffffu) != image.end()
                || std::find(cpuImage.begin(), cpuImage.end(), 0xffffffffu) != cpuImage.end();
            const int error = shape ? rgbError(image, cpuImage) : 256;
            const int peerError = peer ? rgbError(image, *peer) : 0;
            worstError = std::max(worstError, error);
            worstPeerError = std::max(worstPeerError, peerError);
            const bool cacheEligible = samples == 1;
            const bool good = gpuOk && cpuOk && shape && !sentinel && payloadConsistent
                && rawValid && metrics.dispatches > 0 && error <= tolerance
                && (!peer || peerError <= peerTolerance)
                && (expectedHit < 0 || metrics.cacheHit == bool(expectedHit))
                && metrics.metricsCached == cacheEligible
                && (expectedRepairs < 0 || repairs == std::size_t(expectedRepairs));
            if (good) {
                rawGpuPixels += pixelCount - repairs;
                repairedPixels += repairs;
            } else {
                ++failures;
                out << "STABLE_CURVATURE_GPU_TEST_FAIL=" << label
                    << " DEGREE=" << (8 << (order - 1))
                    << " STABLE=" << stable << " PAYLOAD_READY=" << payloadReady
                    << " PAYLOAD_CONSISTENT=" << payloadConsistent
                    << " MAP=" << state.rhExperiment.map
                    << " CHART=" << state.rhExperiment.chart
                    << " SAMPLES_PER_AXIS=" << samples
                    << " GPU_RENDER=" << gpuOk << " CPU_RENDER=" << cpuOk
                    << " USED_GPU=" << renderer.rhUsedGpu()
                    << " REPAIRS=" << repairs << " EXPECTED_REPAIRS=" << expectedRepairs
                    << " EXPECTED_ALL_REPAIRS=" << allRepairsExpected
                    << " CACHE=" << metrics.cacheHit << " EXPECTED_CACHE=" << expectedHit
                    << " METRICS_CACHED=" << metrics.metricsCached
                    << " GPU_PIXELS=" << image.size() << " CPU_PIXELS=" << cpuImage.size()
                    << " MAX_RGB_ERROR=" << error << " PEER_RGB_ERROR=" << peerError
                    << " PEER_TOLERANCE=" << peerTolerance
                    << " DIAGNOSTIC=" << renderer.rhGpuDiagnostic() << '\n';
            }
        };

        std::vector<std::uint32_t> image, reference;
        // Below/above points are separated from |y|=4 by far more than
        // float coordinate error; the narrow cameras do not straddle it.
        const std::array<amp::Camera, 4> cameras{{
            {.3, .2, .015625}, {3.96875, .015625, .015625},
            {4.03125, .015625, .015625}, {6., .25, .015625}
        }};
        const std::array<const char*, 4> branches{{
            "NORMAL", "BELOW_SCALED_BOUNDARY", "ABOVE_SCALED_BOUNDARY", "REVERSED"
        }};
        for (int order = 1; order <= 3; ++order)
            for (std::size_t branch = 0; branch < cameras.size(); ++branch) {
                auto s = base();
                s.spectralOrder = order;
                // Newton-style carrier avoids an escape-radius transition
                // while checking the evaluator's own reversal transition.
                s.rhExperiment.map = 1;
                for (auto& camera : s.cameras) camera = cameras[branch];
                compare(warm, s, branches[branch], 1, 0, image);
            }

        // Both favourite maps retain nonzero complex feedback. Two steps
        // also exercise conjugation and composition after the first update.
        for (int map : {2, 5}) {
            auto s = base();
            s.spectralOrder = 2;
            s.iterations = 2;
            s.rhExperiment.map = map;
            s.rhExperiment.relaxation = .7;
            compare(warm, s, map == 2 ? "TRICORN" : "CAYLEY", 1, 0, image);
        }
        for (int chart = 1; chart <= 3; ++chart) {
            auto s = base();
            s.spectralOrder = 2;
            s.rhExperiment.chart = chart;
            for (auto& camera : s.cameras) camera = amp::Camera{.7, .4, .015625};
            compare(warm, s, "CHART_" + std::to_string(chart), 1, 0, image);
        }
        compare(warm, base(), "SUPERSAMPLES_NO_METRIC_CACHE", 2, 0, image);

        auto tuning = stableTuning;
        auto cacheState = base();
        for (auto& camera : cacheState.cameras) camera = amp::Camera{.31, .27, .015625};
        tuning.stableCurvature = false;
        warm.setGpuTuning(tuning);
        std::vector<std::uint32_t> original, stableFirst, stableCached, stableCold;
        compare(warm, cacheState, "FLAG_FALSE_WARMUP", 1, 0, original);
        compare(warm, cacheState, "FLAG_FALSE_CACHE_HIT", 1, 1, image, &original, 1);
        warm.setGpuTuning(stableTuning);
        compare(warm, cacheState, "FLAG_TRUE_INVALIDATES", 1, 0, stableFirst, &original, tolerance);
        cacheState.hue = .41;
        cacheState.bands = 1.3;
        cacheState.rhExperiment.exposure = 1.1;
        compare(warm, cacheState, "STABLE_RECOLOUR_HIT", 1, 1, stableCached);
        compare(cold, cacheState, "STABLE_COLD_MATCHES_CACHE", 1, 0, stableCold, &stableCached, 1);
        warm.setGpuTuning(tuning);
        compare(warm, cacheState, "FLAG_FALSE_INVALIDATES_AGAIN", 1, 0, image, &stableCold, tolerance);

        // The real log-polar precision guard supplies exactly 16 repair
        // pixels without requiring a pixel to round to a mathematical pole.
        // x=-95,-93,-91,-89 require repair; x=-87,-85,-83,-81 are normal.
        auto repairState = base();
        repairState.rhExperiment.chart = 2;
        repairState.alpha = amp::C{};
        repairState.c = amp::C(.05, .03);
        repairState.iterations = 8;
        for (auto& camera : repairState.cameras) camera = amp::Camera{-88., 0., 16.};
        warm.setGpuTuning(stableTuning);
        compare(warm, repairState, "PRECISION_GUARD_WARMUP", 1, 0, image, nullptr, 1, 16);
        repairState.hue = .43;
        compare(warm, repairState, "PRECISION_MARKERS_REPAIRED_ON_CACHE_HIT", 1, 1, image, nullptr, 1, 16);
        compare(cold, repairState, "PRECISION_COLD_MATCHES_CACHE", 1, 0, reference, &image, 1, 16);

        // The degree-32 expanded J polynomial is severely cancelled near
        // y=4.2i. Its Horner majorant must trigger the real shader guard,
        // rather than accepting a finite but unreliable Newton correction.
        // Nonzero feedback ensures this evaluates the candidate; a narrow
        // camera keeps all 32 samples inside the cancelled neighbourhood.
        auto cancellationState = base();
        cancellationState.spectralOrder = 3;
        cancellationState.rhExperiment.map = 1;
        for (auto& camera : cancellationState.cameras)
            camera = amp::Camera{0., 4.2, .015625};
        compare(warm, cancellationState, "CANCELLED_J_GUARD_WARMUP", 1, 0,
            image, nullptr, 1, int(pixelCount));
        cancellationState.hue = .39;
        compare(warm, cancellationState, "CANCELLED_J_MARKERS_REPAIRED_ON_CACHE_HIT", 1, 1,
            image, nullptr, 1, int(pixelCount));
        compare(cold, cancellationState, "CANCELLED_J_COLD_MATCHES_CACHE", 1, 0,
            reference, &image, 1, int(pixelCount));

        out << "STABLE_CURVATURE_GPU_TEST_CASES=" << cases << " FAILURES=" << failures
            << " MAX_RGB_ERROR=" << worstError << " TOLERANCE=" << tolerance
            << " MAX_PEER_RGB_ERROR=" << worstPeerError
            << " RAW_GPU_PIXELS=" << rawGpuPixels << " CPU_REPAIRS=" << repairedPixels << '\n';
        out << "STABLE_CURVATURE_GPU_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        out << "STABLE_CURVATURE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=" << error.what() << '\n';
        return 1;
    } catch (...) {
        out << "STABLE_CURVATURE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return 1;
    }
}

} // namespace stable_curvature_gpu_tests
