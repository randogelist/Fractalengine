#pragma once
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include "amp_core.hpp"
#include "rh_gpu_data.hpp"
#include "rh_lab_shader.hpp"
#include "render_workers.hpp"
#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace rhgpu {

// These resources belong only to the RH compute path. The device and immediate
// context are borrowed from the existing renderer, which serializes rendering.
template<class T> class GpuCom {
public:
    T* p = nullptr;
    GpuCom() = default;
    GpuCom(const GpuCom&) = delete;
    GpuCom& operator=(const GpuCom&) = delete;
    ~GpuCom() { reset(); }
    void reset() { if(p) { p->Release(); p=nullptr; } }
    T* operator->() const { return p; }
    T** put() { reset(); return &p; }
};

class GpuFailure : public std::runtime_error {
public:
    explicit GpuFailure(const std::string& message) : std::runtime_error(message) {}
};
inline void gpuCheck(HRESULT hr, const char* operation) {
    if(FAILED(hr)) {
        std::ostringstream message;
        message << operation << " (0x" << std::hex << std::uint32_t(hr) << ')';
        throw GpuFailure(message.str());
    }
}

// No compute binding may leak into the legacy renderer, including when the
// user changes a parameter while a tiled dispatch or readback is in progress.
class ComputeBindings {
    ID3D11DeviceContext* context_;
public:
    explicit ComputeBindings(ID3D11DeviceContext* context) : context_(context) {}
    ComputeBindings(const ComputeBindings&) = delete;
    ComputeBindings& operator=(const ComputeBindings&) = delete;
    ~ComputeBindings() { clear(); }
    void clear() const {
        if(!context_) return;
        std::array<ID3D11ShaderResourceView*,D3D11_COMMONSHADER_INPUT_RESOURCE_SLOT_COUNT> views{};
        std::array<ID3D11UnorderedAccessView*,D3D11_PS_CS_UAV_REGISTER_COUNT> outputs{};
        std::array<ID3D11Buffer*,D3D11_COMMONSHADER_CONSTANT_BUFFER_API_SLOT_COUNT> constants{};
        context_->CSSetShaderResources(0,UINT(views.size()),views.data());
        context_->CSSetUnorderedAccessViews(0,UINT(outputs.size()),outputs.data(),nullptr);
        context_->CSSetConstantBuffers(0,UINT(constants.size()),constants.data());
        context_->CSSetShader(nullptr,nullptr,0);
    }
};

// Scheduling and waiting are independent: the diagnostic can reproduce the
// old Sleep(1) wait with either schedule while normal use selects short waits.
enum class Schedule { Adaptive, Legacy };
enum class WaitPolicy { LegacySleep, Responsive };
struct EngineTuning {
    Schedule schedule=Schedule::Adaptive;
    double targetBatchMs=4.;
    std::size_t minimumBatchPixels=32;
    bool metricCache=true;
    bool gpuTimestamps=false;
    WaitPolicy waitPolicy=WaitPolicy::Responsive;
    double waitSliceMs=.35;
    // Zero retains the work-predicted target. Explicit targets are still
    // clipped to the workload cap, independently of minimumBatchPixels.
    std::size_t batchPixels=0;
    bool captureRepairMask=false;
    bool stableCurvature=false;
};
inline bool validTuning(const EngineTuning& tuning) {
    return (tuning.schedule==Schedule::Adaptive||tuning.schedule==Schedule::Legacy)&&
        (tuning.waitPolicy==WaitPolicy::LegacySleep||tuning.waitPolicy==WaitPolicy::Responsive)&&
        std::isfinite(tuning.targetBatchMs)&&tuning.targetBatchMs>0.&&tuning.targetBatchMs<=100.&&
        std::isfinite(tuning.waitSliceMs)&&tuning.waitSliceMs>=.25&&tuning.waitSliceMs<=.5&&
        tuning.minimumBatchPixels>=32&&tuning.minimumBatchPixels<=65536&&
        (tuning.batchPixels==0||(tuning.batchPixels>=32&&tuning.batchPixels<=65536));
}

// Host wall-clock stages: dispatchMs includes bounded completion waits. A cache
// hit still submits a cheap recolour kernel. The optional GPU timestamp span
// covers the entire dispatch timeline (including gaps between host batches),
// not a sum of shader execution durations. It must be read with its valid bit.
struct RenderMetrics {
    double prepareMs=0.,dispatchMs=0.,readbackMs=0.,repairMs=0.;
    double waitMs=0.,maxBatchMs=0.,gpuComputeMs=0.;
    double maxIndividualWaitMs=0.;
    double meanBatchPixels=0.;
    std::size_t minBatchPixels=0,maxBatchPixels=0;
    std::size_t requestedBatchPixels=0,targetBatchPixels=0;
    std::uint64_t batchPixels=0;
    unsigned dispatches=0,completionBatches=0,queryPolls=0,sleepPolls=0,flushes=0;
    unsigned completedBatches=0;
    unsigned timerWaits=0,timerFailures=0,timerFallbacks=0,spinPolls=0,yieldPolls=0;
    unsigned timerTimeouts=0,timerRecreations=0;
    std::uint32_t timerLastError=0,timerLastStatus=0;
    unsigned repairWorkers=0;
    bool cacheHit=false,metricsCached=false;
    bool gpuTimingValid=false,gpuTimingDisjoint=false;
    std::size_t pooledBytes=0,surfaceAllocations=0;
};

class Engine {
    using Clock=std::chrono::steady_clock;
    static double milliseconds(Clock::time_point start) {
        return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    }
    // Lazy creation keeps the timer on the serialized render owner. No global
    // timer-resolution setting is changed. Older Windows versions which lack
    // high-resolution waitable timers use bounded yielding instead of a coarse
    // timer: otherwise the original ~15ms wake floor would return unnoticed.
    class ShortWaitTimer {
        HANDLE timer_=nullptr;
        DWORD ownerThread_=0;
        Clock::time_point retryAt_{};
        unsigned failures_=0;
        bool createdBefore_=false;
        DWORD lastError_=0,lastStatus_=0;
        void reset() {
            if(timer_){CancelWaitableTimer(timer_);CloseHandle(timer_);timer_=nullptr;}
        }
        void fail(DWORD error,DWORD status,RenderMetrics* metrics) {
            lastError_=error;lastStatus_=status;
            if(metrics) {
                ++metrics->timerFailures;
                metrics->timerLastError=error;metrics->timerLastStatus=status;
            }
            reset();
            failures_=std::min(7u,failures_+1);
            // A real API error cannot disable timing for the remaining run.
            // Retry after 4..250ms, avoiding a hot CreateHandle failure loop.
            const unsigned delay=std::min(250u,4u<<(failures_-1));
            retryAt_=Clock::now()+std::chrono::milliseconds(delay);
        }
        bool ensure(RenderMetrics* metrics) {
            const DWORD thread=GetCurrentThreadId();
            if(ownerThread_!=thread) {
                reset();retryAt_={};failures_=0;ownerThread_=thread;
            }
            if(timer_)return true;
            if(Clock::now()<retryAt_) {
                if(metrics){metrics->timerLastError=lastError_;metrics->timerLastStatus=lastStatus_;}
                return false;
            }
            using CreateTimer=HANDLE(WINAPI*)(LPSECURITY_ATTRIBUTES,LPCWSTR,DWORD,DWORD);
            const HMODULE kernel=GetModuleHandleW(L"kernel32.dll");
            const auto create=kernel?reinterpret_cast<CreateTimer>(
                GetProcAddress(kernel,"CreateWaitableTimerExW")):nullptr;
            constexpr DWORD highResolution=0x00000002u;
            if(create)timer_=create(nullptr,nullptr,highResolution,TIMER_MODIFY_STATE|SYNCHRONIZE);
            if(!timer_) {
                const DWORD error=create?GetLastError():ERROR_PROC_NOT_FOUND;
                fail(error,WAIT_FAILED,metrics);
                return false;
            }
            if(createdBefore_&&metrics)++metrics->timerRecreations;
            createdBefore_=true;failures_=0;retryAt_={};
            return true;
        }
    public:
        ShortWaitTimer()=default;
        ShortWaitTimer(const ShortWaitTimer&)=delete;
        ShortWaitTimer& operator=(const ShortWaitTimer&)=delete;
        ~ShortWaitTimer(){reset();}
        template<class Cancelled> bool wait(double sliceMs,Cancelled cancelled,RenderMetrics* metrics) {
            if(cancelled())return false;
            const auto start=Clock::now();
            struct RecordIndividualWait {
                RenderMetrics* metrics;Clock::time_point start;
                ~RecordIndividualWait(){if(metrics)metrics->maxIndividualWaitMs=
                    std::max(metrics->maxIndividualWaitMs,Engine::milliseconds(start));}
            } record{metrics,start};
            if(ensure(metrics)) {
                LARGE_INTEGER due{};
                due.QuadPart=-std::max<LONGLONG>(1,LONGLONG(std::llround(sliceMs*10000.)));
                if(SetWaitableTimer(timer_,&due,0,nullptr,nullptr,FALSE)) {
                    if(metrics)++metrics->timerWaits;
                    // The timer is armed for 0.25-0.5ms. Cancellation is tested
                    // immediately before and after each slice; queued GPU work
                    // remains drainable by the same context owner.
                    // A separate short watchdog prevents an indefinitely
                    // blocked wait if a timer/backend malfunctions. It is not
                    // the requested timer duration or a latency guarantee.
                    const DWORD status=WaitForSingleObject(timer_,2);
                    lastStatus_=status;
                    if(metrics)metrics->timerLastStatus=status;
                    if(status==WAIT_OBJECT_0)return !cancelled();
                    if(status==WAIT_TIMEOUT) {
                        // Scheduling contention can delay a valid short timer.
                        // Preserve the handle; the next slice rearms it. A
                        // watchdog timeout is not an API/handle failure.
                        if(metrics)++metrics->timerTimeouts;
                        CancelWaitableTimer(timer_);
                        return !cancelled();
                    }
                    const DWORD error=status==WAIT_FAILED?GetLastError():ERROR_INVALID_DATA;
                    fail(error,status,metrics);
                } else {
                    const DWORD error=GetLastError();
                    fail(error,WAIT_FAILED,metrics);
                }
            }
            if(metrics)++metrics->timerFallbacks;
            const auto deadline=start+std::chrono::duration<double,std::milli>(sliceMs);
            do {
                if(cancelled())return false;
                if(metrics)++metrics->yieldPolls;
                SwitchToThread();
            } while(Clock::now()<deadline);
            return !cancelled();
        }
    };
    // Each invocation owns new queries, so cancellation cannot cause a later
    // request to reuse a timestamp query still queued on the immediate context.
    // Destruction balances Begin/End even when a cancelled render returns early.
    struct TimestampSpan {
        ID3D11DeviceContext* context=nullptr;
        GpuCom<ID3D11Query> disjoint,start,finish;
        bool begun=false,ended=false;
        bool begin(ID3D11Device* device,ID3D11DeviceContext* owner,bool enabled) {
            if(!enabled)return false;
            D3D11_QUERY_DESC description{};
            description.Query=D3D11_QUERY_TIMESTAMP_DISJOINT;
            if(FAILED(device->CreateQuery(&description,disjoint.put())))return false;
            description.Query=D3D11_QUERY_TIMESTAMP;
            if(FAILED(device->CreateQuery(&description,start.put()))||
               FAILED(device->CreateQuery(&description,finish.put())))return false;
            context=owner;context->Begin(disjoint.p);context->End(start.p);
            begun=true;return true;
        }
        void end() {
            if(!begun||ended)return;
            context->End(finish.p);context->End(disjoint.p);ended=true;
        }
        void read(RenderMetrics& metrics) {
            if(!ended)return;
            D3D11_QUERY_DATA_TIMESTAMP_DISJOINT frequency{};
            UINT64 first=0,last=0;
            // Readback has completed all earlier commands. Never add a second
            // polling loop or flush just to obtain an optional diagnostic.
            if(context->GetData(disjoint.p,&frequency,sizeof(frequency),D3D11_ASYNC_GETDATA_DONOTFLUSH)!=S_OK)
                return;
            metrics.gpuTimingDisjoint=frequency.Disjoint!=FALSE;
            if(frequency.Disjoint||!frequency.Frequency)return;
            if(context->GetData(start.p,&first,sizeof(first),D3D11_ASYNC_GETDATA_DONOTFLUSH)!=S_OK||
               context->GetData(finish.p,&last,sizeof(last),D3D11_ASYNC_GETDATA_DONOTFLUSH)!=S_OK||last<first)
                return;
            metrics.gpuComputeMs=double(last-first)*1000./double(frequency.Frequency);
            metrics.gpuTimingValid=true;
        }
        ~TimestampSpan(){end();}
    };
    // Cache excludes ONLY palette/hue/bands/exposure/frequency. The exact
    // double inputs supplement GPU constants so every changed geometry,
    // dynamic, trap or finite model invalidates, even before float rounding.
    struct MetricKey {
        Constants constants{};
        std::array<double,20> exact{};
        bool operator==(const MetricKey& other) const {
            return std::memcmp(&constants,&other.constants,sizeof(constants))==0&&exact==other.exact;
        }
    };
    static MetricKey metricKey(const amp::Scene& scene,const Constants& uniforms,double aspect) {
        MetricKey key;
        key.constants=uniforms;
        key.constants.style[0]=0;key.constants.style[2]=0;
        key.constants.paint[1]=0.;key.constants.paint[2]=0.;key.constants.paint[3]=0.;
        key.constants.extra[0]=0.;
        const auto& s=scene.state;
        const auto& camera=s.cameras[s.plane];const auto& e=s.rhExperiment;
        key.exact={{camera.x,camera.y,camera.span,aspect,
            s.c.real(),s.c.imag(),s.alpha.real(),s.alpha.imag(),
            s.seed.real(),s.seed.imag(),s.spectralD,scene.masterModel->level,
            e.gain,e.relaxation,e.memory,e.trapRadius,e.trapAngle,e.twist,e.blend,0.}};
        return key;
    }
    struct Surface {
        int width=0,height=0;
        std::uint64_t stamp=0;
        GpuCom<ID3D11Texture2D> output,staging,metrics;
        GpuCom<ID3D11UnorderedAccessView> outputView,metricOutput;
        GpuCom<ID3D11ShaderResourceView> metricInput;
        bool cacheValid=false;
        MetricKey key{};
        std::size_t bytes() const {
            return output.p?std::size_t(width)*height*(8u+(metrics.p?16u:0u)):0u;
        }
        void reset() {
            outputView.reset();metricOutput.reset();metricInput.reset();
            output.reset();staging.reset();metrics.reset();
            width=height=0;stamp=0;cacheValid=false;
        }
    };
    // Three preview dimensions remain resident. Output+staging are8 bytes
    // per pixel, optional selected diagnostic float4 adds16. Large exports
    // may omit metrics, and requests exceeding the entire budget use CPU.
    static constexpr std::size_t PoolBudget=128u*1024u*1024u;
    std::array<Surface,3> surfaces_{};
    std::uint64_t surfaceClock_=0;
    std::size_t surfaceAllocations_=0;
    renderworkers::Pool repairWorkers_;
    GpuCom<ID3D11ComputeShader> shader_,recolourShader_;
    GpuCom<ID3D11Buffer> constants_;
    GpuCom<ID3D11Texture2D> coefficientTexture_;
    GpuCom<ID3D11ShaderResourceView> coefficientView_;
    GpuCom<ID3D11Query> completion_;
    ShortWaitTimer waitTimer_;
    std::vector<std::uint8_t> repairMask_;
    EngineTuning tuning_{};
    ID3D11Device* resourceDevice_ = nullptr;
    bool haveDataKey_ = false, dataReady_ = false;
    int dataGeneration_ = 0, dataOrder_ = 0;
    double dataRadius_ = 0.;
    std::string dataReason_;

    std::size_t pooledBytes() const {
        std::size_t bytes=0;for(const auto& surface:surfaces_)bytes+=surface.bytes();return bytes;
    }
    void invalidateMetricCaches() {
        for(auto& surface:surfaces_)surface.cacheValid=false;
    }
    bool sameDataKey(const amp::State& s) const {
        return haveDataKey_ && dataGeneration_==s.generation &&
            dataOrder_==s.spectralOrder && dataRadius_==s.seed.real();
    }
    void resetResources(ID3D11Device* device) {
        shader_.reset();recolourShader_.reset();constants_.reset();
        coefficientView_.reset();coefficientTexture_.reset();completion_.reset();
        for(auto& surface:surfaces_)surface.reset();
        surfaceClock_=surfaceAllocations_=0;haveDataKey_=dataReady_=false;
        dataReason_.clear();resourceDevice_=device;disabled=false;
    }
    std::string unavailableReason(const amp::Scene& scene, int width) const {
        if(disabled) return diagnostic.empty() ? "RH compute path disabled after a GPU error" : diagnostic;
        if(scene.state.family!=20) return "RH compute path is for the registered RH candidates";
        if(width<1 || scene.state.plane<0 || scene.state.plane>2)
            return "Invalid RH render geometry";
        if(!scene.masterModel || !scene.masterModel->validModel)
            return "Invalid RH model; using CPU";
        const auto& camera=scene.state.cameras[scene.state.plane];
        const double limit=std::max(.0002,8.*std::numeric_limits<float>::epsilon()*
            std::max({1.,std::abs(camera.x),std::abs(camera.y)})*width);
        if(!std::isfinite(camera.x) || !std::isfinite(camera.y) ||
           !std::isfinite(camera.span) || camera.span<limit)
            return "Deep zoom requires double-precision CPU coordinates";
        if(sameDataKey(scene.state) && !dataReady_) return dataReason_;
        return {};
    }
    void initialize(ID3D11Device* device) {
        if(shader_.p&&recolourShader_.p&&constants_.p&&completion_.p) return;
        auto compile=[&](const char* entry,GpuCom<ID3D11ComputeShader>& target) {
            if(target.p)return;
            GpuCom<ID3DBlob> blob,errors;
            const HRESULT compiled=D3DCompile(rhLabShader,std::strlen(rhLabShader),"rh_lab_shader",
                nullptr,nullptr,entry,"cs_5_0",
                D3DCOMPILE_OPTIMIZATION_LEVEL3|D3DCOMPILE_ENABLE_STRICTNESS,0,blob.put(),errors.put());
            if(FAILED(compiled)) {
                const std::string details=errors.p ?
                    std::string(static_cast<const char*>(errors->GetBufferPointer()),errors->GetBufferSize()) :
                    "No compiler diagnostic was returned.";
                throw GpuFailure(std::string("RH compute shader ")+entry+": "+details);
            }
            gpuCheck(device->CreateComputeShader(blob->GetBufferPointer(),blob->GetBufferSize(),
                nullptr,target.put()),"RH compute shader creation");
        };
        compile("main",shader_);compile("recolour",recolourShader_);
        D3D11_BUFFER_DESC buffer{};
        buffer.ByteWidth=UINT(sizeof(Constants));
        buffer.Usage=D3D11_USAGE_DEFAULT;
        buffer.BindFlags=D3D11_BIND_CONSTANT_BUFFER;
        gpuCheck(device->CreateBuffer(&buffer,nullptr,constants_.put()),"RH compute constants");
        D3D11_QUERY_DESC query{};
        query.Query=D3D11_QUERY_EVENT;
        gpuCheck(device->CreateQuery(&query,completion_.put()),"RH completion query");
    }
    bool prepareData(ID3D11Device* device, const amp::State& state) {
        if(sameDataKey(state)) return dataReady_;
        const Data prepared=makeData(state);
        dataGeneration_=state.generation; dataOrder_=state.spectralOrder;
        dataRadius_=state.seed.real(); haveDataKey_=true; dataReady_=false;
        dataReason_=prepared.reason;
        coefficientView_.reset(); coefficientTexture_.reset();
        if(!prepared.valid) return false;
        if(prepared.texels.size()!=std::size_t(CoefficientWidth)*CoefficientHeight) {
            dataReason_="Unexpected RH coefficient texture dimensions; using CPU";
            return false;
        }
        D3D11_TEXTURE2D_DESC texture{};
        texture.Width=CoefficientWidth; texture.Height=CoefficientHeight;
        texture.MipLevels=1; texture.ArraySize=1;
        texture.Format=DXGI_FORMAT_R32G32B32A32_FLOAT;
        texture.SampleDesc.Count=1;
        texture.Usage=D3D11_USAGE_IMMUTABLE;
        texture.BindFlags=D3D11_BIND_SHADER_RESOURCE;
        D3D11_SUBRESOURCE_DATA initial{};
        initial.pSysMem=prepared.texels.data();
        initial.SysMemPitch=UINT(CoefficientWidth*sizeof(std::array<float,4>));
        gpuCheck(device->CreateTexture2D(&texture,&initial,coefficientTexture_.put()),
            "RH coefficient texture");
        D3D11_SHADER_RESOURCE_VIEW_DESC view{};
        view.Format=texture.Format;
        view.ViewDimension=D3D11_SRV_DIMENSION_TEXTURE2D;
        view.Texture2D.MipLevels=1;
        gpuCheck(device->CreateShaderResourceView(coefficientTexture_.p,&view,coefficientView_.put()),
            "RH coefficient texture view");
        dataReady_=true;
        return true;
    }
    bool reserveBytes(std::size_t added,Surface* keep) {
        if(added>PoolBudget||(keep&&keep->bytes()>PoolBudget-added))return false;
        while(pooledBytes()+added>PoolBudget) {
            Surface* victim=nullptr;
            for(auto& candidate:surfaces_)
                if(&candidate!=keep&&candidate.output.p&&(!victim||candidate.stamp<victim->stamp))
                    victim=&candidate;
            if(!victim)return false;
            victim->reset();
        }
        return true;
    }
    Surface* allocate(ID3D11Device* device,int width,int height,bool wantMetrics) {
        Surface* surface=nullptr;
        for(auto& candidate:surfaces_)
            if(candidate.output.p&&candidate.width==width&&candidate.height==height){surface=&candidate;break;}
        if(!surface) {
            for(auto& candidate:surfaces_)
                if(!candidate.output.p){surface=&candidate;break;}
            if(!surface) {
                surface=&surfaces_[0];
                for(auto& candidate:surfaces_)if(candidate.stamp<surface->stamp)surface=&candidate;
            }
            surface->reset();
            const std::size_t bytes=std::size_t(width)*height*8u;
            if(!reserveBytes(bytes,surface))return nullptr;
            D3D11_TEXTURE2D_DESC texture{};
            texture.Width=UINT(width);texture.Height=UINT(height);
            texture.MipLevels=1;texture.ArraySize=1;
            texture.Format=DXGI_FORMAT_R32_UINT;texture.SampleDesc.Count=1;
            texture.Usage=D3D11_USAGE_DEFAULT;texture.BindFlags=D3D11_BIND_UNORDERED_ACCESS;
            gpuCheck(device->CreateTexture2D(&texture,nullptr,surface->output.put()),"RH pooled output texture");
            D3D11_UNORDERED_ACCESS_VIEW_DESC view{};
            view.Format=texture.Format;view.ViewDimension=D3D11_UAV_DIMENSION_TEXTURE2D;
            gpuCheck(device->CreateUnorderedAccessView(surface->output.p,&view,surface->outputView.put()),"RH pooled output view");
            texture.Usage=D3D11_USAGE_STAGING;texture.BindFlags=0;texture.CPUAccessFlags=D3D11_CPU_ACCESS_READ;
            gpuCheck(device->CreateTexture2D(&texture,nullptr,surface->staging.put()),"RH pooled staging texture");
            surface->width=width;surface->height=height;++surfaceAllocations_;
        }
        surface->stamp=++surfaceClock_;
        if(wantMetrics&&!surface->metrics.p&&reserveBytes(std::size_t(width)*height*16u,surface)) {
            D3D11_TEXTURE2D_DESC texture{};
            texture.Width=UINT(width);texture.Height=UINT(height);texture.MipLevels=1;texture.ArraySize=1;
            texture.Format=DXGI_FORMAT_R32G32B32A32_FLOAT;texture.SampleDesc.Count=1;
            texture.Usage=D3D11_USAGE_DEFAULT;
            texture.BindFlags=D3D11_BIND_UNORDERED_ACCESS|D3D11_BIND_SHADER_RESOURCE;
            gpuCheck(device->CreateTexture2D(&texture,nullptr,surface->metrics.put()),"RH compact orbit metrics");
            D3D11_UNORDERED_ACCESS_VIEW_DESC output{};
            output.Format=texture.Format;output.ViewDimension=D3D11_UAV_DIMENSION_TEXTURE2D;
            gpuCheck(device->CreateUnorderedAccessView(surface->metrics.p,&output,surface->metricOutput.put()),"RH orbit metrics output");
            D3D11_SHADER_RESOURCE_VIEW_DESC input{};
            input.Format=texture.Format;input.ViewDimension=D3D11_SRV_DIMENSION_TEXTURE2D;input.Texture2D.MipLevels=1;
            gpuCheck(device->CreateShaderResourceView(surface->metrics.p,&input,surface->metricInput.put()),"RH orbit metrics input");
        }
        return surface;
    }
    template<class Cancelled> bool finishBatch(ID3D11DeviceContext* context, Cancelled cancelled,
                                               bool record=true,bool forceResponsive=false) {
        const auto waitingStart=Clock::now();
        struct RecordWait {
            RenderMetrics* metrics;Clock::time_point start;
            ~RecordWait(){if(metrics)metrics->waitMs+=Engine::milliseconds(start);}
        } recordWait{record?&lastMetrics:nullptr,waitingStart};
        if(record)++lastMetrics.completionBatches;
        context->End(completion_.p);
        context->Flush();
        if(record)++lastMetrics.flushes;
        const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(15);
        const auto spinDeadline=Clock::now()+std::chrono::microseconds(75);
        const WaitPolicy policy=forceResponsive?WaitPolicy::Responsive:tuning_.waitPolicy;
        unsigned polls=0;
        for(;;) {
            if(cancelled()) return false;
            BOOL finished=FALSE;
            const HRESULT hr=context->GetData(completion_.p,&finished,sizeof(finished),
                D3D11_ASYNC_GETDATA_DONOTFLUSH);
            if(record)++lastMetrics.queryPolls;
            gpuCheck(hr,"RH dispatch completion");
            if(hr==S_OK && finished) return true;
            if(std::chrono::steady_clock::now()>=deadline)
                throw GpuFailure("RH GPU dispatch did not finish within 15 seconds; using CPU");
            ++polls;
            if(policy==WaitPolicy::LegacySleep) {
                if(polls<8) {
                    if(record)++lastMetrics.yieldPolls;
                    std::this_thread::yield();
                } else {
                    if(record)++lastMetrics.sleepPolls;
                    const auto start=Clock::now();
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    if(record)lastMetrics.maxIndividualWaitMs=std::max(
                        lastMetrics.maxIndividualWaitMs,milliseconds(start));
                }
            } else if(polls<=32&&Clock::now()<spinDeadline) {
                if((polls&3)==0) {
                    if(record)++lastMetrics.yieldPolls;
                    SwitchToThread();
                } else {
                    if(record)++lastMetrics.spinPolls;
                    YieldProcessor();
                }
            } else if(!waitTimer_.wait(tuning_.waitSliceMs,cancelled,record?&lastMetrics:nullptr)) {
                return false;
            }
        }
    }
public:
    RenderMetrics lastMetrics;
    std::string diagnostic;
    bool lastUsedGpu=false;
    std::size_t repairedPixels=0;
    bool disabled=false;

    // Run once on a background owner before publishing the renderer. Shader
    // compilation and resource creation use the device only; no immediate
    // context is touched and no placeholder frame/model is evaluated. The
    // caller must serialize this with rendering, tuning and resource teardown.
    void prepare(ID3D11Device* device) {
        if(!device)throw GpuFailure("No Direct3D 11 device for RH shader preparation");
        if(resourceDevice_!=device)resetResources(device);
        initialize(device);
    }

    const EngineTuning& tuning() const {return tuning_;}
    const std::vector<std::uint8_t>& repairMask() const {return repairMask_;}
    void setTuning(const EngineTuning& tuning) {
        if(!validTuning(tuning))throw std::invalid_argument("Invalid RH engine tuning");
        // Invoke only on the owner of the serialized immediate context.
        if(!tuning.metricCache||tuning.metricCache!=tuning_.metricCache||
           tuning.stableCurvature!=tuning_.stableCurvature)
            invalidateMetricCaches();
        if(!tuning.captureRepairMask)repairMask_.clear();
        tuning_=tuning;
    }
    // A cancelled request may leave already submitted work on the GPU. Engine
    // probes drain it outside their timed region; regular rendering keeps its
    // existing queue behavior. The caller owns the same immediate-context lock.
    double drain(ID3D11DeviceContext* context) {
        if(!context||!completion_.p)return 0.;
        const auto start=Clock::now();
        finishBatch(context,[]{return false;},false,true);
        return milliseconds(start);
    }

    bool availableFor(const amp::Scene& scene, int width) const {
        return unavailableReason(scene,width).empty();
    }
    std::string status(const amp::Scene& scene, int width) const {
        const std::string reason=unavailableReason(scene,width);
        return reason.empty() ? "Direct3D 11 / RH float32 preview with CPU precision repair" :
            "CPU / "+reason;
    }
    bool render(ID3D11Device* device, ID3D11DeviceContext* context,
                const amp::Scene& scene, int width, int height, double aspect, int samples,
                std::vector<std::uint32_t>& out,
                const std::atomic<std::uint64_t>* revision, std::uint64_t wanted,
                unsigned cpuThreads) {
        lastUsedGpu=false;repairedPixels=0;lastMetrics={};repairMask_.clear();
        const auto preparationStart=Clock::now();
        // Cancellation is read concurrently by repair workers. Only the
        // rendering thread mutates cache validity, after workers have joined.
        struct CacheCommit {
            Engine* engine;bool committed=false;
            ~CacheCommit(){if(!committed)engine->invalidateMetricCaches();}
        } cacheCommit{this};
        auto cancelled=[&] {
            return revision&&revision->load(std::memory_order_relaxed)!=wanted;
        };
        if(!device || !context) { diagnostic="No Direct3D 11 device; using CPU"; return false; }
        ComputeBindings bindings(context);
        if(resourceDevice_!=device) resetResources(device);
        if(cancelled()) return false;
        const std::string reason=unavailableReason(scene,width);
        if(!reason.empty()) { diagnostic=reason; return false; }
        if(height<1 || width>D3D11_REQ_TEXTURE2D_U_OR_V_DIMENSION ||
           height>D3D11_REQ_TEXTURE2D_U_OR_V_DIMENSION || samples<1 || samples>16 ||
           !std::isfinite(aspect) || aspect<=0.) {
            diagnostic="RH render dimensions or sample count require the CPU path";
            return false;
        }
        try {
            const auto& state=scene.state;
            const auto& camera=state.cameras[state.plane];
            const auto& experiment=state.rhExperiment;
            const double originals[]={camera.x,camera.y,camera.span,aspect,
                state.c.real(),state.c.imag(),state.alpha.real(),state.alpha.imag(),
                state.seed.imag(),state.spectralD,scene.masterModel->level,
                experiment.gain,experiment.relaxation,experiment.memory,experiment.trapRadius,
                experiment.trapAngle,experiment.exposure,experiment.frequency,
                experiment.twist,experiment.blend,state.hue,state.bands};
            for(double value:originals)
                if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<float>::max() ||
                   (value!=0. && std::abs(value)<std::numeric_limits<float>::min())) {
                    diagnostic="RH parameters exceed normal float32 arithmetic; using CPU";
                    return false;
                }
            if(!prepareData(device,scene.state)) { diagnostic=dataReason_; return false; }
            if(cancelled()) return false;
            initialize(device);
            if(std::size_t(width)*height*8u>PoolBudget) {
                diagnostic="RH GPU texture pool budget exceeded; using CPU";return false;
            }
            const bool cacheEligible=tuning_.metricCache&&samples==1&&
                experiment.coloring>=1&&experiment.coloring<=7;
            Surface* surface=allocate(device,width,height,cacheEligible);
            if(!surface){diagnostic="RH GPU texture pool budget exceeded; using CPU";return false;}
            out.resize(std::size_t(width)*height);
            Constants uniforms=makeConstants(scene,width,height,aspect,samples);
            uniforms.style[3]=tuning_.stableCurvature?1u:0u;
            const MetricKey key=metricKey(scene,uniforms,aspect);
            const bool cacheAvailable=cacheEligible&&surface->metricOutput.p&&surface->metricInput.p;
            const bool cacheHit=cacheAvailable&&surface->cacheValid&&surface->key==key;
            surface->cacheValid=false; // Commit only after readback and CPU repair succeed.
            uniforms.style[2]=cacheAvailable&&!cacheHit?1u:0u;
            lastMetrics.cacheHit=cacheHit;
            lastMetrics.pooledBytes=pooledBytes();
            lastMetrics.surfaceAllocations=surfaceAllocations_;
            const float* floatRows[]={uniforms.view,uniforms.parameters,uniforms.model,
                uniforms.dynamics,uniforms.paint,uniforms.extra};
            for(const float* values:floatRows)
                for(int i=0;i<4;++i)
                    if(!std::isfinite(values[i]) ||
                       (values[i]!=0.f && std::abs(values[i])<std::numeric_limits<float>::min())) {
                        diagnostic="RH shader constants exceed normal float32 arithmetic; using CPU";
                        return false;
                    }
            lastMetrics.prepareMs=milliseconds(preparationStart);
            const auto dispatchStart=Clock::now();
            TimestampSpan timestamps;
            timestamps.begin(device,context,tuning_.gpuTimestamps);
            context->CSSetShader(cacheHit?recolourShader_.p:shader_.p,nullptr,0);
            ID3D11ShaderResourceView* inputs[]={coefficientView_.p,cacheHit?surface->metricInput.p:nullptr};
            ID3D11UnorderedAccessView* outputs[]={surface->outputView.p,
                uniforms.style[2]?surface->metricOutput.p:nullptr};
            context->CSSetShaderResources(0,2,inputs);
            context->CSSetConstantBuffers(0,1,&constants_.p);
            context->CSSetUnorderedAccessViews(0,2,outputs,nullptr);
            // Predict a bounded batch from the workload. The target is a size
            // hint, never a promised GPU deadline. Do not feed host wake time
            // into this estimate: Sleep(1) can take ~15ms and formerly shrank
            // this scheduler to 32 pixels, creating hundreds of extra waits.
            const double work=cacheHit?1.:std::max(1.,double(uniforms.selection[2])/16.)*
                std::max(1.,double(uniforms.selection[3])/160.)*samples*samples*
                (uniforms.selection[1]==9?2.:1.);
            const std::size_t maximumBatch=cacheHit?65536u:
                std::max(std::size_t(32),std::size_t(65536./work));
            // The diagnostic floor cannot override the work-scaled queue cap.
            // High-degree or supersampled views retain their conservative cap.
            const std::size_t minimumBatch=std::clamp(tuning_.minimumBatchPixels,
                std::size_t(32),maximumBatch);
            const std::size_t targetPixels=tuning_.batchPixels?
                std::clamp(tuning_.batchPixels,std::size_t(32),maximumBatch):
                cacheHit?65536u:std::clamp(std::size_t(2048.*(tuning_.targetBatchMs/4.)/work),
                                         minimumBatch,maximumBatch);
            lastMetrics.requestedBatchPixels=tuning_.batchPixels;
            // Legacy uses fixed 16-dispatch completion groups, not this target.
            lastMetrics.targetBatchPixels=tuning_.schedule==Schedule::Legacy?0:targetPixels;
            std::size_t pendingPixels=0;
            auto batchStart=Clock::now();
            auto complete=[&]() {
                if(!pendingPixels)return true;
                if(!finishBatch(context,cancelled))return false;
                ++lastMetrics.completedBatches;
                if(lastMetrics.completedBatches==1)lastMetrics.minBatchPixels=pendingPixels;
                else lastMetrics.minBatchPixels=std::min(lastMetrics.minBatchPixels,pendingPixels);
                lastMetrics.maxBatchPixels=std::max(lastMetrics.maxBatchPixels,pendingPixels);
                lastMetrics.batchPixels+=pendingPixels;
                lastMetrics.meanBatchPixels=double(lastMetrics.batchPixels)/lastMetrics.completedBatches;
                const double elapsed=std::max(.05,milliseconds(batchStart));
                lastMetrics.maxBatchMs=std::max(lastMetrics.maxBatchMs,elapsed);
                pendingPixels=0;batchStart=Clock::now();
                return true;
            };
            const bool legacy=tuning_.schedule==Schedule::Legacy;
            unsigned legacyPending=0;
            for(int row=0;row<height;) {
                const int tileWidth=legacy?std::min(width,128):
                    std::min(width,int(std::min<std::size_t>(256,targetPixels)));
                const int rows=legacy?std::min(height-row,4):
                    std::min(height-row,std::max(1,std::min(32,int(targetPixels/tileWidth))));
                uniforms.dispatch[2]=UINT(row);uniforms.dispatch[3]=UINT(rows);
                for(int column=0;column<width;column+=tileWidth) {
                    if(cancelled())return false;
                    const int count=std::min(tileWidth,width-column);
                    uniforms.extra[1]=float(column);uniforms.extra[2]=float(count);
                    context->UpdateSubresource(constants_.p,0,nullptr,&uniforms,0,0);
                    context->Dispatch(UINT((count+7)/8),UINT((rows+3)/4),1);
                    ++lastMetrics.dispatches;pendingPixels+=std::size_t(count)*rows;
                    if(legacy) {
                        context->Flush();++lastMetrics.flushes;
                        if(++legacyPending==16) {
                            if(!complete())return false;
                            legacyPending=0;
                        }
                    } else if(pendingPixels>=targetPixels&&!complete())return false;
                }
                row+=rows;
            }
            timestamps.end();
            if(!complete())return false;
            lastMetrics.dispatchMs=milliseconds(dispatchStart);
            if(cancelled())return false;
            bindings.clear();
            const auto readbackStart=Clock::now();
            context->CopyResource(surface->staging.p,surface->output.p);
            D3D11_MAPPED_SUBRESOURCE mapped{};
            gpuCheck(context->Map(surface->staging.p,0,D3D11_MAP_READ,0,&mapped),"RH pixel readback");
            // Nothing in the copy can throw; unmap before cancellation or CPU work.
            for(int y=0;y<height;++y)
                std::memcpy(out.data()+std::size_t(y)*width,
                    static_cast<const char*>(mapped.pData)+std::size_t(y)*mapped.RowPitch,
                    std::size_t(width)*sizeof(std::uint32_t));
            context->Unmap(surface->staging.p,0);
            lastMetrics.readbackMs=milliseconds(readbackStart);
            timestamps.read(lastMetrics);
            const auto repairStart=Clock::now();
            if(cancelled()) return false;
            constexpr std::uint32_t repairMarker=0xffffffffu;
            const std::size_t repairCount=std::size_t(std::count(out.begin(),out.end(),repairMarker));
            if(tuning_.captureRepairMask) {
                repairMask_.resize(out.size());
                std::transform(out.begin(),out.end(),repairMask_.begin(),[](std::uint32_t value) {
                    return std::uint8_t(value==repairMarker);
                });
            }
            if(repairCount) {
                std::atomic<int> nextRow{0};
                std::atomic<bool> failed{false};
                std::atomic<std::size_t> finished{0};
                auto repair=[&] {
                    try {
                        std::size_t local=0;
                        for(int y;(y=nextRow.fetch_add(1,std::memory_order_relaxed))<height;) {
                            if(cancelled()) break;
                            for(int x=0;x<width;++x) {
                                if((x&31)==0 && cancelled()) break;
                                auto& packed=out[std::size_t(y)*width+x];
                                if(packed==repairMarker) {
                                    if(aspect==double(height)/width)
                                        packed=amp::pixel(x,y,width,height,scene,samples);
                                    else {
                                        amp::Color color{};
                                        for(int sy=0;sy<samples;++sy)
                                            for(int sx=0;sx<samples;++sx) {
                                                const amp::C point{
                                                    camera.x+((x+(sx+.5)/samples)/width-.5)*camera.span,
                                                    camera.y+(.5-(y+(sy+.5)/samples)/height)*camera.span*aspect};
                                                const auto sample=amp::renderColor(point,scene);
                                                for(int k=0;k<3;++k)
                                                    color[k]+=sample[k]/double(samples*samples);
                                            }
                                        packed=amp::pack(color);
                                    }
                                    ++local;
                                }
                            }
                        }
                        finished.fetch_add(local,std::memory_order_relaxed);
                    } catch(...) { failed=true; }
                };
                const unsigned available=std::max(1u,std::thread::hardware_concurrency());
                const unsigned requested=cpuThreads ? cpuThreads : std::min(16u,available);
                const unsigned count=std::min({64u,std::max(1u,requested),unsigned(height)});
                repairWorkers_.run(count,repair);
                lastMetrics.repairWorkers=repairWorkers_.lastWorkerCount();
                repairedPixels=finished.load();
                if(cancelled()) return false;
                if(failed || repairedPixels!=repairCount) {
                    diagnostic="RH CPU precision repair failed; using full CPU rendering";
                    return false;
                }
            }
            if(cancelled()) return false;
            lastMetrics.repairMs=milliseconds(repairStart);
            if(cacheAvailable){surface->key=key;surface->cacheValid=true;}
            lastMetrics.metricsCached=surface->cacheValid;
            lastUsedGpu=repairedPixels<out.size();
            diagnostic=lastUsedGpu ? "Direct3D 11 / RH float32 preview; CPU pixel repairs=" :
                "CPU / all RH GPU pixels required precision repair; repaired=";
            diagnostic+=std::to_string(repairedPixels)+"/"+std::to_string(out.size());
            if(cacheHit)diagnostic+="; cached recolour";
            cacheCommit.committed=true;
            return true;
        } catch(const GpuFailure& error) {
            invalidateMetricCaches();diagnostic=error.what();disabled=true;
        } catch(const std::exception& error) {
            invalidateMetricCaches();diagnostic=std::string("RH GPU preparation: ")+error.what();
        } catch(...) {
            invalidateMetricCaches();diagnostic="Unknown RH GPU preparation failure; using CPU";
        }
        return false;
    }
};

} // namespace rhgpu
