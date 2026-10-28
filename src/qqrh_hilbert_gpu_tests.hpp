#pragma once

#include "amp_gpu.hpp"
#include "frk_extra_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <ostream>
#include <string>
#include <vector>

// Native D3D11 execution checks. The raw test entry calls the exact evaluator
// used by the production shader, avoiding the uninformative comparison of two
// dark unresolved images. Only short, well-conditioned orbits are compared.
// Separately exercise the production family-25 field dispatch and colour path.
namespace qqrh_hilbert_gpu_tests {

inline int run(std::ostream& out) {
    try {
        amp::Com<ID3D11Device> device;
        amp::Com<ID3D11DeviceContext> context;
        D3D_FEATURE_LEVEL wanted = D3D_FEATURE_LEVEL_11_0, obtained{};
        const HRESULT created = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE,
            nullptr, 0, &wanted, 1, D3D11_SDK_VERSION, device.put(), &obtained, context.put());
        if (FAILED(created)) {
            out << "QQRH_HILBERT_GPU_TEST=UNAVAILABLE STAGE=DEVICE HRESULT=0x"
                << std::hex << uint32_t(created) << std::dec << '\n';
            return 2;
        }

        constexpr unsigned rawCount = 72;
        const std::string source = std::string(frkExtraShader) + R"QQRH_TEST(
RWStructuredBuffer<float4> qqrhTestOutput : register(u1);
[numthreads(16,1,1)]
void qqrhTestMain(uint3 id : SV_DispatchThreadID) {
    uint i = id.x;
    if(i >= 72u) return;
    uint group = i / 9u, radiusIndex = (i % 9u) / 3u, seedIndex = i % 3u;
    int driver = group == 0u ? 26 : group == 2u ? 28 : group == 4u ? 29 : group == 5u ? 30 : group >= 6u ? 31 : 27;
    float radius = radiusIndex == 0u ? .2 : radiusIndex == 1u ? .6 : 1.1;
    float2 seed = float2(.2 + .25 * seedIndex, .12 - .1 * seedIndex);
    FrkxResult r = frkx_evaluate(driver, seed, radius, group != 3u && group != 7u,
        float2(-.16,.19), (driver == 28 || driver == 30) ? 1 : 3);
    qqrhTestOutput[i] = float4(r.z, r.value, (float)r.iterations);
}
)QQRH_TEST";
        amp::Com<ID3DBlob> code, errors;
        HRESULT compiled = D3DCompile(source.data(), source.size(), "qqrh_hilbert_dispatch_test",
            nullptr, nullptr, "qqrhTestMain", "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3,
            0, code.put(), errors.put());
        if (FAILED(compiled)) {
            out << "QQRH_HILBERT_GPU_TEST=FAIL STAGE=TEST_SHADER "
                << (errors.p ? std::string(static_cast<char*>(errors->GetBufferPointer()),
                    errors->GetBufferSize()) : "No compiler diagnostic") << '\n';
            return 1;
        }
        amp::Com<ID3D11ComputeShader> shader;
        amp::dxcheck(device->CreateComputeShader(code->GetBufferPointer(), code->GetBufferSize(),
            nullptr, shader.put()), "QQRH Hilbert test compute shader");
        amp::Com<ID3D11Buffer> output, staging;
        D3D11_BUFFER_DESC descriptor{};
        descriptor.ByteWidth = rawCount * 4u * sizeof(float);
        descriptor.Usage = D3D11_USAGE_DEFAULT;
        descriptor.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
        descriptor.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
        descriptor.StructureByteStride = 4u * sizeof(float);
        amp::dxcheck(device->CreateBuffer(&descriptor, nullptr, output.put()), "QQRH Hilbert output");
        D3D11_UNORDERED_ACCESS_VIEW_DESC view{};
        view.Format = DXGI_FORMAT_UNKNOWN;
        view.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        view.Buffer.NumElements = rawCount;
        amp::Com<ID3D11UnorderedAccessView> uav;
        amp::dxcheck(device->CreateUnorderedAccessView(output.p, &view, uav.put()),
            "QQRH Hilbert output UAV");
        descriptor.Usage = D3D11_USAGE_STAGING;
        descriptor.BindFlags = 0;
        descriptor.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        descriptor.MiscFlags = 0;
        descriptor.StructureByteStride = 0;
        amp::dxcheck(device->CreateBuffer(&descriptor, nullptr, staging.put()),
            "QQRH Hilbert readback");
        context->CSSetShader(shader.p, nullptr, 0);
        context->CSSetUnorderedAccessViews(1, 1, &uav.p, nullptr);
        context->Dispatch((rawCount + 15u) / 16u, 1, 1);
        ID3D11UnorderedAccessView* empty = nullptr;
        context->CSSetUnorderedAccessViews(1, 1, &empty, nullptr);
        context->CopyResource(staging.p, output.p);
        D3D11_MAPPED_SUBRESOURCE mapped{};
        amp::dxcheck(context->Map(staging.p, 0, D3D11_MAP_READ, 0, &mapped), "QQRH Hilbert readback map");
        std::array<std::array<float,4>,rawCount> values{};
        std::memcpy(values.data(), mapped.pData, sizeof(values));
        context->Unmap(staging.p, 0);

        int failures = 0;
        double worstEndpointError = 0, worstFieldError = 0;
        for (unsigned i = 0; i < rawCount; ++i) {
            const unsigned group = i / 9u, radiusIndex = (i % 9u) / 3u, seedIndex = i % 3u;
            const int driver = group == 0u ? 26 : group == 2u ? 28 : group == 4u ? 29 : group == 5u ? 30 : group >= 6u ? 31 : 27;
            const double radius = radiusIndex == 0u ? .2 : radiusIndex == 1u ? .6 : 1.1;
            const frk_extra::C seed(.2 + .25 * seedIndex, .12 - .1 * seedIndex);
            frk_extra::Params parameters;
            parameters.lambda = radius;
            parameters.julia = group != 3u && group != 7u;
            parameters.julia_c = {-.16,.19};
            parameters.iterations = (driver == 28 || driver == 30) ? 1 : 3;
            const auto expected = frk_extra::evaluate(driver, seed, parameters);
            const auto& actual = values[i];
            const double endpointError = std::abs(frk_extra::C(actual[0],actual[1]) - expected.final_z)
                / (1. + std::abs(expected.final_z));
            const double fieldError = std::abs(double(actual[2]) - expected.value);
            worstEndpointError = std::max(worstEndpointError, endpointError);
            worstFieldError = std::max(worstFieldError, fieldError);
            const bool finite = std::all_of(actual.begin(), actual.end(), [](float x) { return std::isfinite(x); });
            if (!finite || endpointError > 2e-4 || fieldError > 2e-5 ||
                actual[3] != float(expected.iterations)) {
                ++failures;
                out << "QQRH_HILBERT_GPU_FAIL=RAW_EVALUATOR CASE=" << i << " DRIVER=" << driver
                    << " RADIUS=" << radius << " ENDPOINT_REL_ERROR=" << endpointError
                    << " FIELD_ABS_ERROR=" << fieldError << " GPU_STEPS=" << actual[3]
                    << " CPU_STEPS=" << expected.iterations << '\n';
            }
        }

        // Validate the registered renderer and production compute entry as well
        // as the evaluator. No fallback can satisfy this check: gpu must remain
        // active, and the camera span exceeds the renderer's float-zoom cutoff.
        amp::Renderer gpu, cpu;
        if (!gpu.init(false, 1)) {
            out << "QQRH_HILBERT_GPU_TEST=FAIL STAGE=PRODUCTION_INITIALIZE REASON="
                << gpu.diagnostic << '\n';
            return 1;
        }
        cpu.init(true, 1);
        constexpr int width = 24, height = 12, rgbTolerance = 4;
        int imageCases = 0, worstRgbError = 0;
        for (int driver : {26,29}) for (double radius : {.2,.6,1.1}) {
            for (int palette = 0; palette < 3; ++palette) {
                ++imageCases;
                auto state = amp::defaults(25, driver-17);
                state.plane = 1;
                state.alpha = {radius,0};
                state.palette = palette;
                for (auto& camera : state.cameras) camera = amp::Camera{.5,0.,3.};
                const auto scene = amp::scene(state);
                std::vector<uint32_t> actual, expected;
                const bool gpuOk = gpu.render(scene, width, height, double(height)/width, 1, actual);
                const bool cpuOk = cpu.render(scene, width, height, double(height)/width, 1, expected);
                int error = 0;
                const bool shape = actual.size() == size_t(width*height) && actual.size() == expected.size();
                if (shape) for (size_t j = 0; j < actual.size(); ++j)
                    for (int shift : {0,8,16})
                        error = std::max(error, std::abs(int((actual[j] >> shift)&255u)
                            - int((expected[j] >> shift)&255u)));
                const bool nonconstant = shape && std::any_of(actual.begin(), actual.end(),
                    [&](uint32_t pixel) { return (pixel & 0xffffffu) != (actual.front() & 0xffffffu); });
                const bool actualGpu = gpu.gpu && gpu.cpuWorkerCount() == 0 &&
                    gpu.currentBackend(scene,width).find("Direct3D 11") != std::string::npos;
                worstRgbError = std::max(worstRgbError,error);
                if (!gpuOk || !cpuOk || !shape || !actualGpu || !nonconstant || error > rgbTolerance) {
                    ++failures;
                    out << "QQRH_HILBERT_GPU_FAIL=PRODUCTION_FIELD DRIVER=" << driver << " RADIUS=" << radius
                        << " PALETTE=" << palette << " GPU_ACTIVE=" << actualGpu
                        << " NONCONSTANT=" << nonconstant << " MAX_RGB_ERROR=" << error
                        << " TOLERANCE=" << rgbTolerance << " DIAGNOSTIC=" << gpu.diagnostic << '\n';
                }
            }
        }
        out << "QQRH_HILBERT_GPU_RAW_CASES=" << rawCount << " PRODUCTION_FIELD_CASES=" << imageCases
            << " MAX_ENDPOINT_REL_ERROR=" << worstEndpointError << " MAX_FIELD_ABS_ERROR=" << worstFieldError
            << " MAX_RGB_ERROR=" << worstRgbError << " FAILURES=" << failures << '\n';
        out << "QQRH_HILBERT_GPU_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        out << "QQRH_HILBERT_GPU_TEST=FAIL STAGE=EXCEPTION REASON=" << error.what() << '\n';
        return 1;
    } catch (...) {
        out << "QQRH_HILBERT_GPU_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return 1;
    }
}
} // namespace qqrh_hilbert_gpu_tests
