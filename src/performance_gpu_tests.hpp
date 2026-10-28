#pragma once

#include "amp_gpu.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// User-side D3D11 checks for the performance paths. These deliberately use
// short, smooth orbits: resource/cache correctness must not depend on whether
// different precisions follow the same long chaotic boundary trajectory.
namespace performance_gpu_tests {

inline int run(std::ostream& out) {
    try {
        amp::Renderer warm, cold, cpu;
        if (!warm.init(false, 1)) {
            const bool unavailable = warm.diagnostic.find("Direct3D 11 device") != std::string::npos;
            out << "PERFORMANCE_GPU_TEST=" << (unavailable ? "UNAVAILABLE" : "FAIL")
                << " STAGE=INITIALIZE REASON=" << warm.diagnostic << '\n';
            return unavailable ? 2 : 1;
        }
        if (!cold.init(false, 1)) {
            out << "PERFORMANCE_GPU_TEST=FAIL STAGE=REFERENCE_INITIALIZE REASON="
                << cold.diagnostic << '\n';
            return 1;
        }
        cpu.init(true, 1);
        int cases = 0, failures = 0, worstCpuError = 0;
        constexpr std::size_t poolLimit = 128u * 1024u * 1024u;
        constexpr int width = 8, height = 4;
        constexpr int cpuTolerance = 8;

        auto check = [&](bool ok, const std::string& name) {
            ++cases;
            if (!ok) {
                ++failures;
                out << "PERFORMANCE_GPU_TEST_FAIL=" << name << '\n';
            }
            return ok;
        };
        auto base = [] {
            auto s = amp::defaults(20, 1);
            s.plane = 1;
            s.iterations = 8;
            s.spectralOrder = 1;
            s.c = amp::C(.05, .03);
            s.alpha = amp::C{};
            s.palette = 0;
            s.hue = .2;
            s.bands = 1.;
            s.rhExperiment = rhmaster::Experiment{};
            s.rhExperiment.coloring = 1;
            for (auto& camera : s.cameras) camera = amp::Camera{.1, .1, .02};
            return s;
        };
        auto gpuRender = [&](amp::Renderer& renderer, const amp::State& s,
                             int w, int h, double aspect, int samples,
                             std::vector<std::uint32_t>& image, int expectedHit,
                             const std::string& label) {
            const bool ok = renderer.render(amp::scene(s), w, h, aspect, samples, image);
            const auto metrics = renderer.rhGpuMetrics();
            bool valid = ok && renderer.rhUsedGpu() && image.size() == std::size_t(w) * h
                && std::find(image.begin(), image.end(), 0xffffffffu) == image.end()
                && metrics.dispatches > 0 && metrics.pooledBytes <= poolLimit;
            for (double value : {metrics.prepareMs, metrics.dispatchMs,
                                 metrics.readbackMs, metrics.repairMs})
                valid = valid && std::isfinite(value) && value >= 0.;
            if (expectedHit >= 0) valid = valid && metrics.cacheHit == bool(expectedHit);
            const bool eligible = samples == 1 && s.rhExperiment.coloring >= 1
                && s.rhExperiment.coloring <= 7;
            valid = valid && metrics.metricsCached == eligible;
            if (!valid) {
                out << "PERFORMANCE_GPU_RENDER_DETAIL=" << label
                    << " GPU=" << renderer.rhUsedGpu() << " CACHE=" << metrics.cacheHit
                    << " METRICS=" << metrics.metricsCached
                    << " DISPATCHES=" << metrics.dispatches
                    << " POOLED_BYTES=" << metrics.pooledBytes
                    << " REPAIRS=" << renderer.rhCpuRepairs()
                    << " DIAGNOSTIC=" << renderer.rhGpuDiagnostic() << '\n';
            }
            return check(valid, label);
        };
        auto exact = [&](const std::vector<std::uint32_t>& a,
                         const std::vector<std::uint32_t>& b, const std::string& label) {
            check(a == b, label);
        };
        auto nearCold = [&](const std::vector<std::uint32_t>& a,
                            const std::vector<std::uint32_t>& b, const std::string& label) {
            int error = 0;
            if (a.size() != b.size()) error = 256;
            else for (std::size_t i = 0; i < a.size(); ++i)
                for (int shift : {0, 8, 16})
                    error = std::max(error, std::abs(int((a[i] >> shift) & 255u)
                        - int((b[i] >> shift) & 255u)));
            // The two entrypoints have different compiler optimisation
            // contexts. Allow one final quantisation unit, while keeping a
            // much tighter check than the independent CPU/float32 tolerance.
            if (error > 1) out << "PERFORMANCE_GPU_CACHE_RGB_ERROR=" << error << '\n';
            check(error <= 1, label);
        };
        auto cpuCompare = [&](const amp::State& s, const std::vector<std::uint32_t>& gpuImage,
                              const std::string& label) {
            std::vector<std::uint32_t> cpuImage;
            const bool ok = cpu.render(amp::scene(s), width, height,
                                       double(height) / width, 1, cpuImage);
            int error = 0;
            if (ok && cpuImage.size() == gpuImage.size()) {
                for (std::size_t i = 0; i < cpuImage.size(); ++i)
                    for (int shift : {0, 8, 16})
                        error = std::max(error, std::abs(int((gpuImage[i] >> shift) & 255u)
                            - int((cpuImage[i] >> shift) & 255u)));
            } else error = 256;
            worstCpuError = std::max(worstCpuError, error);
            check(error <= cpuTolerance, label);
        };

        // Every cache-eligible colouring must preserve the raw-GPU
        // result through appearance edits. The independent renderer sees a
        // new colouring on each iteration and therefore takes the cold path.
        for (int coloring = 1; coloring <= 7; ++coloring) {
            auto s = base();
            s.rhExperiment.coloring = coloring;
            if (coloring == 6) s.iterations = 64;
            if (coloring == 3 || coloring == 5) {
                // Exercise the favourite Tricorn/Cayley paths with actual
                // nonzero rational feedback rather than an empty statistic.
                s.spectralOrder = 2; s.iterations = 2;
                s.c = amp::C(-.3, .1); s.alpha = amp::C(.2, .11);
                s.rhExperiment.map = coloring == 3 ? 2 : 5;
                s.rhExperiment.power = 3; s.rhExperiment.gain = .8;
                s.rhExperiment.relaxation = .7;
                for (auto& camera : s.cameras) camera = amp::Camera{.12, -.09, .02};
            }
            std::vector<std::uint32_t> first, cached, fresh;
            const std::string tag = "COLOUR_" + std::to_string(coloring);
            gpuRender(warm, s, width, height, .5, 1, first, 0, tag + "_WARMUP");
            s.palette = 2; s.hue = .37; s.bands = 1.7;
            s.rhExperiment.exposure = 1.2; s.rhExperiment.frequency = 1.4;
            gpuRender(warm, s, width, height, .5, 1, cached, 1, tag + "_RECOLOUR");
            gpuRender(cold, s, width, height, .5, 1, fresh, 0, tag + "_COLD");
            nearCold(cached, fresh, tag + "_MATCH_COLD");
            cpuCompare(s, cached, tag + "_CPU_REFERENCE");
        }

        // Individual appearance fields are checked separately so an
        // accidentally overbroad key cannot hide behind a combined edit.
        auto appearance = base();
        std::vector<std::uint32_t> image, reference;
        gpuRender(warm, appearance, width, height, .5, 1, image, 0, "PAINT_WARMUP");
        const std::vector<std::pair<const char*, std::function<void(amp::State&)>>> paintEdits{
            {"PALETTE", [](amp::State& s) { s.palette = 1; }},
            {"HUE", [](amp::State& s) { s.hue = .43; }},
            {"BANDS", [](amp::State& s) { s.bands = 2.1; }},
            {"EXPOSURE", [](amp::State& s) { s.rhExperiment.exposure = 1.1; }},
            {"FREQUENCY", [](amp::State& s) { s.rhExperiment.frequency = 1.3; }}
        };
        for (const auto& edit : paintEdits) {
            edit.second(appearance);
            const std::string tag = std::string("PAINT_") + edit.first;
            gpuRender(warm, appearance, width, height, .5, 1, image, 1, tag);
            // Invalidate the reference renderer explicitly; its final render
            // must really evaluate orbits rather than use another cache hit.
            auto different = appearance;
            ++different.iterations;
            gpuRender(cold, different, width, height, .5, 1, reference, 0, tag + "_RESET");
            gpuRender(cold, appearance, width, height, .5, 1, reference, 0, tag + "_COLD");
            nearCold(image, reference, tag + "_MATCH_COLD");
        }

        // Check each key field from a known warm baseline. Although some
        // fields are inactive in this smooth fixture, conservative misses
        // are required because other maps/colourings do depend on them.
        const std::vector<std::pair<const char*, std::function<void(amp::State&)>>> orbitEdits{
            {"ITERATIONS", [](amp::State& s) { ++s.iterations; }},
            {"CANDIDATE", [](amp::State& s) { s.generation = 2; }},
            {"DEGREE", [](amp::State& s) { s.spectralOrder = 2; }},
            {"RADIUS", [](amp::State& s) { s.seed.real(s.seed.real() + .125); }},
            {"SCALE", [](amp::State& s) { s.seed.imag(s.seed.imag() + .125); }},
            {"SHIFT", [](amp::State& s) { s.spectralD += .01; }},
            {"LEVEL", [](amp::State& s) { s.rhLevel += .125; }},
            {"CAMERA_X", [](amp::State& s) { s.cameras[1].x += .01; }},
            {"CAMERA_Y", [](amp::State& s) { s.cameras[1].y += .01; }},
            {"CAMERA_SPAN", [](amp::State& s) { s.cameras[1].span *= 1.1; }},
            {"PARAMETER_REAL", [](amp::State& s) { s.c.real(.06); }},
            {"PARAMETER_IMAG", [](amp::State& s) { s.c.imag(.04); }},
            {"STRENGTH_REAL", [](amp::State& s) { s.alpha.real(.01); }},
            {"STRENGTH_IMAG", [](amp::State& s) { s.alpha.imag(.01); }},
            {"PLANE", [](amp::State& s) { s.plane = 0; }},
            {"MAP", [](amp::State& s) { s.rhExperiment.map = 2; }},
            {"CHART", [](amp::State& s) { s.rhExperiment.chart = 2; }},
            {"COLOURING", [](amp::State& s) { s.rhExperiment.coloring = 2; }},
            {"POWER", [](amp::State& s) { s.rhExperiment.power = 3; }},
            {"GAIN", [](amp::State& s) { s.rhExperiment.gain = .8; }},
            {"RELAXATION", [](amp::State& s) { s.rhExperiment.relaxation = .8; }},
            {"MEMORY", [](amp::State& s) { s.rhExperiment.memory = .1; }},
            {"TRAP_RADIUS", [](amp::State& s) { s.rhExperiment.trapRadius += .1; }},
            {"TRAP_ANGLE", [](amp::State& s) { s.rhExperiment.trapAngle += .1; }},
            {"TWIST", [](amp::State& s) { s.rhExperiment.twist = .2; }},
            {"BLEND", [](amp::State& s) { s.rhExperiment.blend = .6; }}
        };
        for (const auto& edit : orbitEdits) {
            auto s = base();
            s.iterations = 1;
            gpuRender(warm, s, width, height, .5, 1, image, -1,
                      std::string("KEY_") + edit.first + "_WARMUP");
            edit.second(s);
            gpuRender(warm, s, width, height, .5, 1, image, 0,
                      std::string("KEY_") + edit.first + "_MISS");
        }

        // Supersamples are shaded separately before averaging. Recolouring
        // an averaged OrbitResult would be mathematically incorrect, so that
        // path and the legacy colour mode intentionally never cache metrics.
        auto s = base();
        gpuRender(warm, s, width, height, .5, 2, image, 0, "SUPERSAMPLES_NO_CACHE");
        gpuRender(warm, s, width, height, .5, 2, image, 0, "SUPERSAMPLES_REPEAT_NO_CACHE");
        s.rhExperiment.coloring = 0;
        gpuRender(warm, s, width, height, .5, 1, image, 0, "CLASSIC_NO_CACHE");
        gpuRender(warm, s, width, height, .5, 1, image, 0, "CLASSIC_REPEAT_NO_CACHE");

        // A non-multiple of the shader group shape tests pooled texture
        // boundaries. A -> B -> A must retain A's actual allocation, not
        // merely keep a constant aggregate memory footprint while replacing it.
        s = base();
        gpuRender(warm, s, width, height, .5, 1, image, -1, "POOL_A");
        const auto allocationA = warm.rhGpuMetrics().surfaceAllocations;
        gpuRender(warm, s, 19, 7, 7. / 19., 1, reference, 0, "POOL_B");
        const auto allocationB = warm.rhGpuMetrics().surfaceAllocations;
        const auto bytesB = warm.rhGpuMetrics().pooledBytes;
        gpuRender(warm, s, width, height, .5, 1, reference, 1, "POOL_A_RETURN");
        check(allocationB == allocationA + 1
              && warm.rhGpuMetrics().surfaceAllocations == allocationB
              && warm.rhGpuMetrics().pooledBytes == bytesB, "POOL_A_B_A_REUSES_SURFACE");
        nearCold(image, reference, "POOL_A_B_A_MATCH");
        gpuRender(warm, s, width, height, .6, 1, image, 0, "ASPECT_INVALIDATES");
        gpuRender(warm, s, width, height, .6, 1, reference, 1, "ASPECT_REPEAT_HITS");
        nearCold(image, reference, "ASPECT_REPEAT_MATCH");
        gpuRender(warm, s, width, height, .6, 1, image, 1, "RECOLOUR_REPEAT_HITS");
        exact(image, reference, "SAME_RECOLOUR_KERNEL_EXACT_REPEAT");

        // Exercise the real log-polar precision guard without requiring an
        // exact zero after shader division by odd canvas dimensions. This
        // dyadic grid has real coordinates -95,-93,-91,-89,-87,-85,-83,-81.
        // log(FLT_MIN) is between -89 and -87: the first four columns require
        // CPU precision, and the remaining four stay on the GPU. Their tiny
        // initial states enter the smooth zero-feedback fixture immediately.
        auto repairState = base();
        repairState.rhExperiment.chart = 2;
        for (auto& camera : repairState.cameras) camera = amp::Camera{-88., 0., 16.};
        constexpr std::size_t expectedRepairs = 16;
        gpuRender(warm, repairState, 8, 4, .5, 1, image, 0, "REPAIR_WARMUP");
        out << "REPAIR_FIXTURE_STAGE=WARMUP EXPECTED=" << expectedRepairs
            << " ACTUAL=" << warm.rhCpuRepairs() << '\n';
        check(warm.rhCpuRepairs() == expectedRepairs,
              "REPAIR_WARMUP_MIXED_GPU_CPU");
        repairState.hue = .41;
        repairState.rhExperiment.exposure = 1.1;
        gpuRender(warm, repairState, 8, 4, .5, 1, image, 1, "REPAIR_RECOLOUR_HIT");
        out << "REPAIR_FIXTURE_STAGE=RECOLOUR EXPECTED=" << expectedRepairs
            << " ACTUAL=" << warm.rhCpuRepairs() << '\n';
        check(warm.rhCpuRepairs() == expectedRepairs,
              "REPAIR_MARKER_REPAIRED_AGAIN");
        gpuRender(cold, repairState, 8, 4, .5, 1, reference, 0, "REPAIR_COLD");
        out << "REPAIR_FIXTURE_STAGE=COLD EXPECTED=" << expectedRepairs
            << " ACTUAL=" << cold.rhCpuRepairs() << '\n';
        check(cold.rhCpuRepairs() == expectedRepairs, "REPAIR_COLD_MIXED_GPU_CPU");
        nearCold(image, reference, "REPAIR_RECOLOUR_MATCH_COLD");
        cpuCompare(repairState, image, "REPAIR_RECOLOUR_CPU_REFERENCE");

        // Pre-cancellation is deterministic and checks the contract without
        // making driver scheduling or a wall-clock timeout a pass criterion.
        std::atomic<std::uint64_t> revision{2};
        const std::vector<std::uint32_t> unpublished{0x00123456u, 0x00654321u};
        image = unpublished;
        const bool cancelled = warm.render(amp::scene(s), width, height, .5, 1,
                                           image, &revision, 1);
        check(!cancelled && image == unpublished, "PRE_CANCEL_NO_OUTPUT_MUTATION");
        revision.store(3);
        check(warm.render(amp::scene(s), width, height, .5, 1, image, &revision, 3),
              "CANCEL_RECOVERY_COMPLETES");
        auto resetCold = s;
        ++resetCold.iterations;
        gpuRender(cold, resetCold, width, height, .5, 1, reference, 0, "CANCEL_RECOVERY_RESET_COLD");
        gpuRender(cold, s, width, height, .5, 1, reference, 0, "CANCEL_RECOVERY_COLD");
        nearCold(image, reference, "CANCEL_RECOVERY_MATCH");

        out << "PERFORMANCE_GPU_TEST_CASES=" << cases << " FAILURES=" << failures
            << " MAX_CPU_RGB_ERROR=" << worstCpuError << " CPU_TOLERANCE=" << cpuTolerance
            << " SURFACE_ALLOCATIONS=" << warm.rhGpuMetrics().surfaceAllocations
            << " POOLED_BYTES=" << warm.rhGpuMetrics().pooledBytes << '\n';
        out << "PERFORMANCE_GPU_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        out << "PERFORMANCE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=" << error.what() << '\n';
        return 1;
    } catch (...) {
        out << "PERFORMANCE_GPU_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return 1;
    }
}

} // namespace performance_gpu_tests
