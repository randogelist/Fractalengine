#pragma once
#define NOMINMAX
#include "amp_core.hpp"
#include "amp_shader.hpp"
#include "garden_shader.hpp"
#include "spectral_shader.hpp"
#include "quantum_render.hpp"
#include "spherical_render.hpp"
#include "spin_render.hpp"
#include "rh_lab_gpu.hpp"
#include "frk_extra_shader.hpp"
#include "render_workers.hpp"
#include <cstring>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <iostream>
#include <thread>
#include <windows.h>

namespace amp {
template <class T> struct Com {
    T *p = nullptr;
    Com() = default;
    Com(const Com &) = delete;
    Com &operator=(const Com &) = delete;
    ~Com() { reset(); }
    void reset() {
        if (p) {
            p->Release();
            p = nullptr;
        }
    }
    T *operator->() const { return p; }
    T **put() {
        reset();
        return &p;
    }
};
inline std::string utf8(const std::wstring &s) {
    int n = WideCharToMultiByte(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0, nullptr, nullptr);
    std::string r(n, '\0');
    WideCharToMultiByte(CP_UTF8, 0, s.data(), int(s.size()), r.data(), n, nullptr, nullptr);
    return r;
}
inline std::wstring wide(const std::string &s) {
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring r(n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), r.data(), n);
    return r;
}
inline void dxcheck(HRESULT hr, const char *what) {
    if (FAILED(hr)) {
        std::ostringstream o;
        o << what << " (0x" << std::hex << uint32_t(hr) << ")";
        throw std::runtime_error(o.str());
    }
}
struct alignas(16) Uniforms {
    uint32_t sizeRow[4], mode[4];
    float view[4], param[4], option[4], refA[4], refB[4];
    uint32_t flags[4], garden[4], record[4];
};
static_assert(sizeof(Uniforms) == 160, "Constant buffer layout");
class Renderer {
    Com<ID3D11Device> device;
    Com<ID3D11DeviceContext> context;
    Com<ID3D11ComputeShader> shader, candidateShader, spectralComputeShader, extraComputeShader;
    Com<ID3D11Texture2D> coefficientBuffer;
    Com<ID3D11ShaderResourceView> coefficientView;
    int coefficientKey = -1;
    C coefficientQ{};
    Com<ID3D11Buffer> constants;
    Com<ID3D11Texture2D> buffer, staging;
    Com<ID3D11UnorderedAccessView> uav;
    int allocatedW = 0, allocatedH = 0;
    unsigned cpuThreads = 0;
    unsigned lastCpuWorkers = 0;
    quantum::CloudRenderer quantumRenderer;
    spherical::CloudRenderer sphericalRenderer;
    spin::CloudRenderer spinRenderer;
    rhgpu::Engine rhLabGpu;
    renderworkers::Pool cpuWorkers;
    struct BkqrCachedScan {
        bkqrlambda::Settings settings;
        C center{};
        double span=0;
        int gridCap=0,seedCap=0,stepCap=0;
        std::shared_ptr<const bkqrlambda::Grid> grid;
    };
    std::array<BkqrCachedScan,3> bkqrCache;

    static bool sameBkqrScan(const bkqrlambda::Settings& a,const bkqrlambda::Settings& b) {
        // All metrics are computed by one scan; changing the displayed metric
        // must not rerun the ensemble or alter its percentile reference grid.
        return a.c==b.c && a.host==b.host && a.branch==b.branch && a.coupling==b.coupling &&
            a.power==b.power && a.quality==b.quality && a.precision==b.precision && a.iterations==b.iterations &&
            a.rhChannel==b.rhChannel && a.rhN==b.rhN && a.rhX==b.rhX && a.rhB==b.rhB &&
            a.rhH==b.rhH && a.eta==b.eta &&
            a.qqrhChannel==b.qqrhChannel && a.qqrhSlope==b.qqrhSlope;
    }

  public:
    bool gpu = false;
    std::string backend = "CPU / double precision", diagnostic;
    void setThreadLimit(unsigned count) { cpuThreads=count?std::clamp(count,1u,64u):0; }
    unsigned threadLimit() const { return cpuThreads; }
    unsigned cpuWorkerCount() const { return lastCpuWorkers; }
    void setGpuTuning(const rhgpu::EngineTuning& tuning) { rhLabGpu.setTuning(tuning); }
    void setRhTuning(const rhgpu::EngineTuning& tuning) { setGpuTuning(tuning); }
    const rhgpu::EngineTuning& gpuTuning() const { return rhLabGpu.tuning(); }
    // Background initialization owns this renderer until publication. Only
    // device methods and CPU shader compilation run here; no context calls.
    void prepareRhGpu() {
        if(!gpu||!device.p)throw std::runtime_error("No GPU device for RH shader preparation");
        rhLabGpu.prepare(device.p);
    }
    double drainGpu() { return gpu&&context.p?rhLabGpu.drain(context.p):0.; }
    double drainRhGpu() { return drainGpu(); }
    bool init(bool forceCPU = false, unsigned requestedThreads = 0) {
        cpuThreads = requestedThreads;
        if (forceCPU) {
            diagnostic = "CPU requested";
            return false;
        }
        try {
            D3D_FEATURE_LEVEL wanted = D3D_FEATURE_LEVEL_11_0, got;
            dxcheck(D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 0, &wanted, 1,
                                      D3D11_SDK_VERSION, device.put(), &got, context.put()),
                    "Direct3D 11 device");
            Com<IDXGIDevice> dxgi;
            Com<IDXGIAdapter> adapter;
            dxcheck(device->QueryInterface(__uuidof(IDXGIDevice), reinterpret_cast<void **>(dxgi.put())),
                    "DXGI device");
            dxcheck(dxgi->GetAdapter(adapter.put()), "DXGI adapter");
            DXGI_ADAPTER_DESC desc{};
            adapter->GetDesc(&desc);
            Com<ID3DBlob> blob, errors;
            HRESULT hr = D3DCompile(ampShader, std::strlen(ampShader), "amp_shader", nullptr, nullptr, "main",
                                    "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, blob.put(), errors.put());
            if (FAILED(hr)) {
                std::string msg = errors.p ? std::string(static_cast<char *>(errors->GetBufferPointer()),
                                                         errors->GetBufferSize())
                                           : "Unknown shader error";
                throw std::runtime_error("Shader compilation: " + msg);
            }
            dxcheck(device->CreateComputeShader(blob->GetBufferPointer(), blob->GetBufferSize(), nullptr,
                                                shader.put()),
                    "Compute shader");
            D3D11_BUFFER_DESC b{};
            b.ByteWidth = sizeof(Uniforms);
            b.Usage = D3D11_USAGE_DEFAULT;
            b.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
            dxcheck(device->CreateBuffer(&b, nullptr, constants.put()), "Constant buffer");
            backend = "Direct3D 11 / " + utf8(desc.Description);
            gpu = true;
            diagnostic = "Shader Model 5 / float32";
            return true;
        } catch (const std::exception &e) {
            diagnostic = e.what();
            gpu = false;
            return false;
        }
    }
    bool render(const Scene &scene, int w, int h, double aspect, int samples, std::vector<uint32_t> &out,
                const std::atomic<uint64_t> *revision = nullptr, uint64_t wanted = 0, int historyPreviewPower = 0) {
        auto cancelled = [&] { return revision && revision->load() != wanted; };
        if (cancelled())
            return false;
        lastCpuWorkers=0;
        if (scene.state.family == 20 && gpu) {
            if (rhLabGpu.render(device.p, context.p, scene, w, h, aspect, samples, out,
                                revision, wanted, cpuThreads))
                return true;
            if (cancelled())return false;
        }
        if (scene.state.family == 21 || scene.state.family == 22 || scene.state.family == 23) {
            if(w<1 || h<1) return false;
            const auto& q=scene.bkqr;
            const auto& cam=scene.state.cameras[2];
            const C center(cam.x,cam.y);
            const int stage=historyPreviewPower==0 ? 2 : historyPreviewPower<=17 ? 0 : 1;
            const int gridCap=stage==0 ? std::min(33,bkqrlambda::gridSize(q))
                : stage==1 ? std::min(65,bkqrlambda::gridSize(q)) : 0;
            const int seedCap=stage==0 ? 4 : stage==1 ? 8 : 0;
            const int stepCap=stage==0 ? 32 : stage==1 ? 64 : 0;
            auto& cached=bkqrCache[stage];
            if(!cached.grid || !sameBkqrScan(q,cached.settings) || cached.center!=center ||
               cached.span!=cam.span || cached.gridCap!=gridCap || cached.seedCap!=seedCap || cached.stepCap!=stepCap) {
                auto grid=std::make_shared<bkqrlambda::Grid>(
                    bkqrlambda::scanGrid(q,center,cam.span,gridCap,cancelled,seedCap,stepCap,cpuThreads));
                if(cancelled() || grid->cancelled) return false;
                if(!grid->valid || grid->side<1) throw std::runtime_error("BKQR lambda scan: "+grid->diagnostic);
                cached.settings=q;cached.center=center;cached.span=cam.span;
                cached.gridCap=gridCap;cached.seedCap=seedCap;cached.stepCap=stepCap;
                cached.grid=std::move(grid);
            }
            const auto& grid=*cached.grid;
            out.assign(size_t(w)*h,0xff101119u);
            const int side=std::min(w,h),left=(w-side)/2,top=(h-side)/2;
            for(int y=0;y<side;++y) {
                if(cancelled()) return false;
                // The legacy grid starts at lower Im; screen rows start at
                // upper Im.  Preserve its nearest-sample projection exactly.
                const int iy=std::clamp(int(std::llround((1.0-(y+.5)/side)*(grid.side-1))),0,grid.side-1);
                for(int x=0;x<side;++x) {
                    const int ix=std::clamp(int(std::llround((x+.5)/side*(grid.side-1))),0,grid.side-1);
                    const auto rgb=bkqrlambda::metricCell(grid,ix,iy,q.metric);
                    const auto byte=[](double v) {return uint32_t(std::clamp(v,0.0,1.0)*255.0+.5);};
                    out[size_t(top+y)*w+left+x]=0xff000000u|(byte(rgb[0])<<16)|(byte(rgb[1])<<8)|byte(rgb[2]);
                }
            }
            return !cancelled();
        }
        if (scene.state.family == 14) {
            const auto& s = scene.state;
            return spinRenderer.render(s.generation, s.spin, s.quantum, s.plane, s.palette, s.hue, s.method != 0,
                                       w, h, samples, cpuThreads, out, revision, wanted, historyPreviewPower);
        }
        if (scene.state.family == 13) {
            const auto& s = scene.state;
            return sphericalRenderer.render(s.generation, s.sphere, s.quantum, s.plane, s.palette, s.hue, s.method != 0,
                                            w, h, samples, cpuThreads, out, revision, wanted, historyPreviewPower);
        }
        if (scene.state.family == 12) {
            const auto& s = scene.state;
            return quantumRenderer.render(s.quantum, s.plane, s.palette, s.hue, s.method != 0,
                                          w, h, samples, cpuThreads, out, revision, wanted, historyPreviewPower);
        }
        out.resize(size_t(w) * h);
        // Float coordinates must not be allowed to collapse silently at deep zoom.
        const auto &cam = scene.state.cameras[scene.state.plane];
        double gpuLimit = std::max(.0002, 8.0 * std::numeric_limits<float>::epsilon() *
                                              std::max({1., std::abs(cam.x), std::abs(cam.y)}) * w);
        bool useGpu = gpu && cam.span >= gpuLimit && scene.state.family != 6 && scene.state.family != 7 && scene.state.family != 8 && scene.state.family != 9 && scene.state.family != 10 && scene.state.family != 15 && scene.state.family != 16 && scene.state.family != 17 && scene.state.family != 19 && scene.state.family != 18 && scene.state.family != 20 && scene.state.family != 24;
        if (useGpu) {
            try {

                if (w != allocatedW || h != allocatedH) {
                    uav.reset();
                    buffer.reset();
                    staging.reset();
                    D3D11_TEXTURE2D_DESC b{};
                    b.Width = w;
                    b.Height = h;
                    b.MipLevels = 1;
                    b.ArraySize = 1;
                    b.Format = DXGI_FORMAT_R32_UINT;
                    b.SampleDesc.Count = 1;
                    b.Usage = D3D11_USAGE_DEFAULT;
                    b.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
                    dxcheck(device->CreateTexture2D(&b, nullptr, buffer.put()), "Pixel texture");
                    D3D11_UNORDERED_ACCESS_VIEW_DESC u{};
                    u.Format = DXGI_FORMAT_R32_UINT;
                    u.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
                    dxcheck(device->CreateUnorderedAccessView(buffer.p, &u, uav.put()), "Pixel UAV");
                    b.Usage = D3D11_USAGE_STAGING;
                    b.BindFlags = 0;
                    b.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
                    dxcheck(device->CreateTexture2D(&b, nullptr, staging.put()), "Readback texture");
                    allocatedW = w;
                    allocatedH = h;
                }
                const State &s = scene.state;
                if (s.family == 25 && !extraComputeShader.p) {
                    Com<ID3DBlob> blob,errors;
                    HRESULT hr=D3DCompile(frkExtraShader,std::strlen(frkExtraShader),"frk_extra_shader",
                                          nullptr,nullptr,"main","cs_5_0",D3DCOMPILE_OPTIMIZATION_LEVEL3,
                                          0,blob.put(),errors.put());
                    if(FAILED(hr))
                        throw std::runtime_error(errors.p ? std::string(
                            static_cast<char*>(errors->GetBufferPointer()),errors->GetBufferSize())
                            : "Additional fractal compute shader compilation failed");
                    dxcheck(device->CreateComputeShader(blob->GetBufferPointer(),blob->GetBufferSize(),
                                                        nullptr,extraComputeShader.put()),
                            "Additional fractal compute shader");
                }
                if (s.family == 11 && !spectralComputeShader.p) {
                    // Compiled lazily on the user's computer when this family is selected.
                    Com<ID3DBlob> blob, errors;
                    HRESULT hr = D3DCompile(spectralShader, std::strlen(spectralShader), "spectral_shader",
                                            nullptr, nullptr, "main", "cs_5_0",
                                            D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, blob.put(), errors.put());
                    if (FAILED(hr))
                        throw std::runtime_error(
                            errors.p ? std::string(static_cast<char *>(errors->GetBufferPointer()),
                                                   errors->GetBufferSize())
                                     : "Spectral shader compilation failed");
                    dxcheck(device->CreateComputeShader(blob->GetBufferPointer(), blob->GetBufferSize(),
                                                        nullptr, spectralComputeShader.put()),
                            "Spectral compute shader");
                }
                if (s.family && s.family != 11 && s.family != 25 && !candidateShader.p) {
                    Com<ID3DBlob> blob, errors;
                    HRESULT hr = D3DCompile(gardenShader, std::strlen(gardenShader), "garden_shader", nullptr,
                                            nullptr, "main", "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3, 0,
                                            blob.put(), errors.put());
                    if (FAILED(hr))
                        throw std::runtime_error(
                            errors.p ? std::string(static_cast<char *>(errors->GetBufferPointer()),
                                                   errors->GetBufferSize())
                                     : "Garden shader compilation failed");
                    dxcheck(device->CreateComputeShader(blob->GetBufferPointer(), blob->GetBufferSize(),
                                                        nullptr, candidateShader.put()),
                            "Garden compute shader");
                }
                int key = s.family * 16 + s.generation + (s.family == 3 && s.plane == 1 ? 256 : 0);
                if (s.family && s.family != 11 && s.family != 25 &&
                    (coefficientKey != key || (s.family == 3 && coefficientQ != s.c))) {
                    std::vector<std::array<float, 4>> records;
                    for (size_t i = 0; i < scene.gardenData->records.size(); i++) {
                        auto row = scene.gardenData->records[i];
                        if (s.family == 3) {
                            C root = scene.crownRoots.empty() ? C{} : scene.crownRoots[i];
                            records.push_back(
                                {float(root.real()), float(root.imag()), float(row[0]), float(row[1])});
                        } else
                            records.push_back({float(row[0]), float(row[1]), float(row[2]), float(row[3])});
                    }
                    D3D11_TEXTURE2D_DESC b{};
                    b.Width = UINT(records.size());
                    b.Height = 1;
                    b.MipLevels = 1;
                    b.ArraySize = 1;
                    b.Format = DXGI_FORMAT_R32G32B32A32_FLOAT;
                    b.SampleDesc.Count = 1;
                    b.Usage = D3D11_USAGE_IMMUTABLE;
                    b.BindFlags = D3D11_BIND_SHADER_RESOURCE;
                    D3D11_SUBRESOURCE_DATA initial{};
                    initial.pSysMem = records.data();
                    initial.SysMemPitch = UINT(records.size() * 16);
                    coefficientView.reset();
                    dxcheck(device->CreateTexture2D(&b, &initial, coefficientBuffer.put()),
                            "Garden coefficient buffer");
                    D3D11_SHADER_RESOURCE_VIEW_DESC v{};
                    v.Format = DXGI_FORMAT_R32G32B32A32_FLOAT;
                    v.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
                    v.Texture2D.MipLevels = 1;
                    dxcheck(device->CreateShaderResourceView(coefficientBuffer.p, &v, coefficientView.put()),
                            "Garden coefficient view");
                    coefficientKey = key;
                    coefficientQ = s.c;
                }
                Uniforms u{};
                u.sizeRow[0] = w;
                u.sizeRow[1] = h;
                u.mode[0] = s.plane;
                u.mode[1] = s.orbitMode;
                u.mode[2] = s.palette;
                u.mode[3] = s.iterations;
                u.view[0] = float(cam.x);
                u.view[1] = float(cam.y);
                u.view[2] = float(cam.span);
                u.view[3] = float(aspect);
                u.param[0] = float(s.c.real());
                u.param[1] = float(s.c.imag());
                u.param[2] = float(s.alpha.real());
                u.param[3] = float(s.alpha.imag());
                u.option[0] = float(s.seed.real());
                u.option[1] = float(s.seed.imag());
                u.option[2] = float(s.hue);
                u.option[3] = float(s.bands);
                auto copyQ = [](float *a, Q q) {
                    a[0] = float(q.x.real());
                    a[1] = float(q.x.imag());
                    a[2] = float(q.y.real());
                    a[3] = float(q.y.imag());
                };
                copyQ(u.refA, scene.a.last);
                copyQ(u.refB, scene.b.last);
                if (s.family == 11) {
                    u.refA[0] = float(s.spectralD);
                    u.refA[1] = float(s.spectralWeight);
                    u.record[1] = uint32_t(s.spectralOrder);
                }
                u.flags[0] = scene.a.period ? 1 : 0;
                u.flags[1] = scene.b.period && !scene.same ? 1 : 0;
                u.flags[2] = samples;
                u.garden[0] = s.family;
                u.garden[1] = s.generation;
                u.garden[2] = s.method;
                u.garden[3] = scene.gardenData ? scene.gardenData->degree : 0;
                u.record[0] = scene.gardenData ? UINT(scene.gardenData->records.size()) : 0;
                context->CSSetShader(s.family == 25 ? extraComputeShader.p :
                                     s.family == 11 ? spectralComputeShader.p :
                                     s.family ? candidateShader.p : shader.p, nullptr, 0);
                if (s.family && s.family != 11 && s.family != 25)
                    context->CSSetShaderResources(0, 1, &coefficientView.p);
                else {
                    ID3D11ShaderResourceView *emptyView = nullptr;
                    context->CSSetShaderResources(0, 1, &emptyView);
                }
                context->CSSetConstantBuffers(0, 1, &constants.p);
                context->CSSetUnorderedAccessViews(0, 1, &uav.p, nullptr);
                // Small dispatch batches keep cancellation responsive and avoid long
                // single dispatches.
                for (int row = 0; row < h; row += 16) {
                    if (cancelled())
                        return false;
                    u.sizeRow[2] = row;
                    u.sizeRow[3] = std::min(16, h - row);
                    context->UpdateSubresource(constants.p, 0, nullptr, &u, 0, 0);
                    context->Dispatch((w + 15) / 16, (u.sizeRow[3] + 7) / 8, 1);
                    context->Flush();
                }
                ID3D11UnorderedAccessView *empty = nullptr;
                context->CSSetUnorderedAccessViews(0, 1, &empty, nullptr);
                context->CopyResource(staging.p, buffer.p);
                D3D11_MAPPED_SUBRESOURCE mapped{};
                dxcheck(context->Map(staging.p, 0, D3D11_MAP_READ, 0, &mapped), "GPU readback");
                for (int y = 0; y < h; y++)
                    std::memcpy(out.data() + size_t(y) * w,
                                static_cast<const char *>(mapped.pData) + size_t(y) * mapped.RowPitch,
                                size_t(w) * 4);
                context->Unmap(staging.p, 0);
                return !cancelled();
            } catch (const std::exception &e) {
                diagnostic = e.what();
                std::cerr << "GPU_FALLBACK=" << diagnostic << std::endl;
                gpu = false;
                backend = "CPU / double precision (GPU fallback)";
            }
        }
        // Family 9 is a synchronized three-panel diagnostic.  Evaluate each
        // (log2(a),logit(r)) sample once, then fan the same result out to all
        // three panels.  This avoids tripling the completed-L work.
        if (scene.state.family == 9) {
            int panelW = std::max(1, w / 3);
            int commonW = panelW;
            double localAspect = double(h) / panelW;
            std::atomic<int> prow{0};
            std::atomic<bool> pfailed{false};
            auto pwork = [&] {
                try {
                    for (int y; (y = prow.fetch_add(1)) < h;) {
                        if (cancelled()) return;
                        for (int x = 0; x < commonW; x++) {
                            std::array<Color,3> colors{};
                            for (int sy=0; sy<samples; sy++) for (int sx=0; sx<samples; sx++) {
                                C p{cam.x + ((x + (sx + .5) / samples) / panelW - .5) * cam.span,
                                    cam.y + (.5 - (y + (sy + .5) / samples) / h) * cam.span * localAspect};
                                auto cc=grhObservableColors(p,scene.state);
                                for(int panel=0;panel<3;panel++) for(int k=0;k<3;k++)
                                    colors[panel][k]+=cc[panel][k]/double(samples*samples);
                            }
                            for (int panel=0; panel<3; panel++) {
                                int xx=x+panel*panelW;
                                if (xx<w) out[size_t(y)*w+xx]=pack(colors[panel]);
                            }
                        }
                        // At most two remainder columns when w is not divisible by 3.
                        for (int x=3*panelW; x<w; x++) {
                            int local=x-2*panelW;
                            Color col{};
                            for (int sy=0;sy<samples;sy++) for(int sx=0;sx<samples;sx++) {
                                C p{cam.x + ((local + (sx + .5) / samples) / (w-2*panelW) - .5) * cam.span,
                                    cam.y + (.5 - (y + (sy + .5) / samples) / h) * cam.span * (double(h)/(w-2*panelW))};
                                auto cc=grhObservableColors(p,scene.state)[2];
                                for(int k=0;k<3;k++) col[k]+=cc[k]/double(samples*samples);
                            }
                            out[size_t(y)*w+x]=pack(col);
                        }
                    }
                } catch (...) { pfailed=true; }
            };
            unsigned available=std::max(1u,std::thread::hardware_concurrency());
            unsigned requested=cpuThreads?cpuThreads:std::min(16u,available);
            unsigned count=std::min({64u,std::max(1u,requested),unsigned(h)});
            cpuWorkers.run(count, pwork);
            lastCpuWorkers=cpuWorkers.lastWorkerCount();
            if(pfailed) throw std::runtime_error("CPU GRH observable render failed");
            return !cancelled();
        }
        std::atomic<int> row{0};
        std::atomic<bool> failed{false};
        auto work = [&] {
            try {
                for (int y; (y = row.fetch_add(1)) < h;) {
                    if (cancelled())
                        return;
                    for (int x = 0; x < w; x++) {
                        if((x&15)==0 && cancelled())return;
                        Color col{};
                        for (int sy = 0; sy < samples; sy++)
                            for (int sx = 0; sx < samples; sx++) {
                                C p{cam.x + ((x + (sx + .5) / samples) / w - .5) * cam.span,
                                    cam.y + (.5 - (y + (sy + .5) / samples) / h) * cam.span * aspect};
                                auto cc = renderColor(p, scene);
                                for (int k = 0; k < 3; k++)
                                    col[k] += cc[k] / double(samples * samples);
                            }
                        out[size_t(y) * w + x] = pack(col);
                    }
                }
            } catch (...) {
                failed = true;
            }
        };
        unsigned available = std::max(1u, std::thread::hardware_concurrency());
        unsigned requested = cpuThreads ? cpuThreads : std::min(16u, available);
        unsigned count = std::min({64u, std::max(1u, requested), unsigned(h)});
        cpuWorkers.run(count, work);
        lastCpuWorkers=cpuWorkers.lastWorkerCount();
        if (failed)
            throw std::runtime_error("CPU render failed");
        return !cancelled();
    }
    std::string currentBackend(const Scene &s, int w) const {
        if (s.state.family == 24) return "CPU / parallel double-precision direct QQRH dynamics";
        if (s.state.family == 23) return "CPU / BKQR lambda ensemble with QQRH carrier";
        if (s.state.family == 22) return "CPU / BKQR lambda ensemble with connected RH host";
        if (s.state.family == 21) return "CPU / original BKQR ensemble grid and percentile metrics";
        if (s.state.family == 20) {
            if (gpu && rhLabGpu.lastUsedGpu)
                return "Direct3D 11 / RH Lab float32" +
                    (rhLabGpu.repairedPixels ? " + CPU repairs="+std::to_string(rhLabGpu.repairedPixels) : "");
            return "CPU / RH Lab double + long-double jets" +
                (gpu && !rhLabGpu.diagnostic.empty() ? " / "+rhLabGpu.diagnostic : "");
        }
        if (s.state.family == 19) return "CPU / parallel double-precision moment dynamics";
        if (s.state.family == 18) return "CPU / parallel double-precision rational Julia dynamics";
        if (s.state.family == 17) return "CPU / parallel double-precision cubic Julia dynamics";
        if (s.state.family == 16) return "CPU / parallel double-precision Riemann witness fractals";
        if (s.state.family == 15) return "CPU / parallel double-precision project fractals";
        if (s.state.family == 14) return s.state.generation==3 && s.state.plane==0
            ? "CPU / parallel double-precision spin parameter atlas"
            : "CPU / parallel double-precision spin trajectories and analytic tangents";
        if (s.state.family == 13) return s.state.generation==3 ? "CPU / parallel double-precision Julia inverse branches" : "CPU / parallel double-precision circle reflections";
        if (s.state.family == 12) return "CPU / parallel double-precision quantum trajectories";
        if (s.state.family == 6) return "CPU / double precision (analytic L jets)";
        if (s.state.family == 7) return "CPU / double precision (primitive GRH Cayley jets)";
        if (s.state.family == 8) return "CPU / double precision (primitive character-pair jets)";
        if (s.state.family == 9) return "CPU / double precision (literal GRH Cayley observable)";
        if (s.state.family == 10) return "CPU / double precision (GRH phase-feedback jets)";
        auto v = s.state.cameras[s.state.plane];
        double min = std::max(.0002, 8.0 * std::numeric_limits<float>::epsilon() *
                                         std::max({1., std::abs(v.x), std::abs(v.y)}) * w);
        return gpu && v.span < min ? "CPU / double precision (deep zoom)" : backend;
    }
    quantum::Stats quantumStats() const { return quantumRenderer.stats(); }
    bool rhUsedGpu() const { return gpu && rhLabGpu.lastUsedGpu; }
    std::size_t rhCpuRepairs() const { return rhLabGpu.repairedPixels; }
    const std::vector<std::uint8_t>& rhRepairMask() const { return rhLabGpu.repairMask(); }
    std::string rhGpuDiagnostic() const { return rhLabGpu.diagnostic; }
    const rhgpu::RenderMetrics& rhGpuMetrics() const { return rhLabGpu.lastMetrics; }
    spherical::Stats sphericalStats() const { return sphericalRenderer.stats(); }
    spin::Stats spinStats() const { return spinRenderer.stats(); }
};
} // namespace amp
