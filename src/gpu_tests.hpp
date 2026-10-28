#pragma once

#include "amp_gpu.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <ostream>
#include <string>
#include <vector>

// A short, smooth float32 preview comparison on the user's D3D11 hardware.
// This is a numerical integration check, not a proof or a requirement that
// long chaotic orbits agree with the CPU at individual fractal boundaries.
namespace gpu_tests {

inline int run(std::ostream& out) {
    try {
        amp::Renderer gpu, cpu;
        if (!gpu.init(false, 1)) {
            if (gpu.diagnostic.find("Direct3D 11 device") != std::string::npos) {
                out << "GPU_SELF_TEST=UNAVAILABLE REASON=" << gpu.diagnostic << '\n';
                return 2;
            }
            out << "GPU_SELF_TEST=FAIL STAGE=INITIALIZE REASON=" << gpu.diagnostic << '\n';
            return 1;
        }
        cpu.init(true, 1);
        constexpr int width = 8, height = 4;
        constexpr int tolerance = 8;
        constexpr std::size_t pixels = std::size_t(width) * height;
        int cases = 0, failures = 0, worstError = 0;
        std::size_t rawGpuPixels = 0, repairedPixels = 0;

        auto base = [](int candidate = 1) {
            auto s = amp::defaults(20, candidate);
            s.plane = 1;
            s.iterations = 1;
            s.palette = 0;
            s.hue = .2;
            s.bands = 1.;
            s.c = amp::C(-.3, .1);
            s.alpha = amp::C(.15, .03);
            s.rhExperiment = rhmaster::Experiment{};
            s.rhExperiment.coloring = 1;
            for (auto& camera : s.cameras) camera = amp::Camera{.31, .27, .02};
            return s;
        };

        auto compare = [&](const std::string& name, const amp::State& state, int samplesPerAxis = 1) {
            ++cases;
            const int degree = 8 << (std::clamp(state.spectralOrder, 1, 3) - 1);
            std::vector<std::uint32_t> gpuImage, cpuImage;
            bool gpuOk = false, cpuOk = false;
            std::string exception;
            try {
                const auto scene = amp::scene(state);
                gpuOk = gpu.render(scene, width, height, double(height) / width, samplesPerAxis, gpuImage);
                cpuOk = cpu.render(scene, width, height, double(height) / width, samplesPerAxis, cpuImage);
            } catch (const std::exception& e) {
                exception = e.what();
            } catch (...) {
                exception = "Unknown render exception";
            }
            const std::size_t repairs = gpu.rhCpuRepairs();
            const bool shape = gpuImage.size() == pixels && cpuImage.size() == pixels;
            const bool usedGpu = gpu.rhUsedGpu();
            const bool rawPixels = repairs < pixels;
            int maxError = 0;
            bool sentinel = false;
            if (shape) {
                for (std::size_t i = 0; i < pixels; ++i) {
                    sentinel = sentinel || gpuImage[i] == 0xffffffffu || cpuImage[i] == 0xffffffffu;
                    for (int shift : {0, 8, 16}) {
                        const int a = int((gpuImage[i] >> shift) & 255u);
                        const int b = int((cpuImage[i] >> shift) & 255u);
                        maxError = std::max(maxError, a > b ? a - b : b - a);
                    }
                }
            }
            const bool good = exception.empty() && gpuOk && cpuOk && shape && usedGpu
                && rawPixels && !sentinel && maxError <= tolerance;
            worstError = std::max(worstError, maxError);
            if (good) {
                rawGpuPixels += pixels - repairs;
                repairedPixels += repairs;
            } else {
                ++failures;
                out << "GPU_TEST_FAIL=" << name
                    << " CANDIDATE=" << state.generation
                    << " DEGREE=" << degree
                    << " MAP=" << state.rhExperiment.map
                    << " CHART=" << state.rhExperiment.chart
                    << " COLOURING=" << state.rhExperiment.coloring
                    << " PLANE=" << state.plane
                    << " ITERATIONS=" << state.iterations
                    << " SAMPLES_PER_AXIS=" << samplesPerAxis
                    << " MAX_RGB_ERROR=" << maxError
                    << " TOLERANCE=" << tolerance
                    << " REPAIRS=" << repairs
                    << " USED_GPU=" << (usedGpu ? "YES" : "NO")
                    << " GPU_PIXELS=" << gpuImage.size()
                    << " CPU_PIXELS=" << cpuImage.size()
                    << " GPU_RENDER=" << (gpuOk ? "YES" : "NO")
                    << " CPU_RENDER=" << (cpuOk ? "YES" : "NO")
                    << " SENTINEL=" << (sentinel ? "YES" : "NO")
                    << " DIAGNOSTIC=" << gpu.rhGpuDiagnostic();
                if (!exception.empty()) out << " EXCEPTION=" << exception;
                out << '\n';
            }
        };

        // Nonzero complex feedback exercises every rational target and the
        // coefficient/derivative payload at each supported Taylor degree.
        for (int candidate = 1; candidate <= rhfresh::OriginalCandidateCount; ++candidate) {
            for (int order = 1; order <= 3; ++order) {
                auto s = base(candidate);
                s.spectralOrder = order;
                compare("CANDIDATE_DEGREE", s);
            }
        }

        // Two steps expose alternating sign and nonzero inertial memory.
        // One-step cases also check the common relaxation/power path.
        for (int map = 0; map < rhmaster::MapCount; ++map) {
            for (int iterations = 1; iterations <= 2; ++iterations) {
                auto s = base();
                s.alpha = amp::C(.2, .11);
                s.iterations = iterations;
                s.rhExperiment.map = map;
                s.rhExperiment.power = 3;
                s.rhExperiment.gain = .8;
                s.rhExperiment.relaxation = .7;
                s.rhExperiment.memory = .15;
                if(map>=6) {s.rhExperiment.twist=.4;s.rhExperiment.blend=.7;}
                for (auto& camera : s.cameras) camera = amp::Camera{.12, -.09, .02};
                compare("MAP_ONE_TWO_STEPS", s);
            }
        }

        // Charts affect the displayed Julia seed. The zero-gain smooth map
        // isolates coordinate composition from sensitive rational feedback.
        for (int chart = 0; chart < rhmaster::ChartCount; ++chart) {
            auto s = base();
            s.rhExperiment.chart = chart;
            s.rhExperiment.gain = 0.;
            for (auto& camera : s.cameras) camera = amp::Camera{.7, .4, .02};
            compare("CHART_PULLBACK", s);
        }

        // A bounded polynomial fixed point gives stable traps, winding and
        // radius statistics. A small nonzero c avoids the float32 underflow
        // of repeated squaring toward zero, which would erase endpoint phase
        // well before the corresponding CPU orbit. The period colour gets
        // enough tail samples to detect the fixed point numerically.
        for (int colouring = 0; colouring < rhmaster::ColouringCount; ++colouring) {
            auto s = base();
            s.c = amp::C(.05, .03);
            s.alpha = amp::C{};
            s.iterations = colouring == 6 ? 64 : 8;
            s.rhExperiment.coloring = colouring;
            for (auto& camera : s.cameras) camera = amp::Camera{.1, .1, .02};
            compare("ORBIT_COLOURING", s);
        }
        {
            auto s = base();
            s.iterations = 2;
            s.rhExperiment.coloring = 5;
            compare("NONZERO_FEEDBACK_COLOURING", s);
        }

        for (int plane = 0; plane < 3; ++plane) {
            auto s = base();
            s.plane = plane;
            compare("PARAMETER_JULIA_LAMBDA_PLANES", s);
        }
        for (int palette = 0; palette < 3; ++palette) {
            auto s = base();
            s.palette = palette;
            s.alpha = amp::C{};
            s.rhExperiment.coloring = 0;
            // Escaped colours use each palette's hue/smoothing branch.
            // An unresolved one-step seed would return the same dark colour
            // for all three palettes and fail to cover those branches.
            for (auto& camera : s.cameras) camera = amp::Camera{18., .2, .02};
            compare("ORIGINAL_PALETTE", s);
        }
        compare("SUPERSAMPLING", base(), 2);
        // The original lambda view uses finite-period diagnostic colours.
        {
            auto s = base();
            s.plane = 2;
            s.c = amp::C{};
            s.alpha = amp::C{};
            s.iterations = 64;
            s.rhExperiment.gain = 0.;
            s.rhExperiment.coloring = 0;
            for (auto& camera : s.cameras) camera = amp::Camera{.1, .1, .02};
            compare("LAMBDA_PERIOD_DIAGNOSTIC", s);
        }
        out << "GPU_SELF_TEST_CASES=" << cases
            << " FAILURES=" << failures
            << " MAX_RGB_ERROR=" << worstError
            << " TOLERANCE=" << tolerance
            << " RAW_GPU_PIXELS=" << rawGpuPixels
            << " CPU_REPAIRS=" << repairedPixels << '\n';
        out << "GPU_SELF_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& e) {
        out << "GPU_SELF_TEST=FAIL STAGE=EXCEPTION REASON=" << e.what() << '\n';
        return 1;
    } catch (...) {
        out << "GPU_SELF_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return 1;
    }
}

} // namespace gpu_tests
