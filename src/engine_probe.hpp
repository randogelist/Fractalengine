#pragma once
#include "amp_gpu.hpp"
#include "background_gpu.hpp"
#include "focus_presets.hpp"
#include <dxgi.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <optional>
#include <sstream>
#include <streambuf>
#include <thread>
#include <vector>
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
#include <cpuid.h>
#endif

// Native, headless measurements. Never changes a user's rendering settings and
// never scores a cancelled or partial frame as a completed performance sample.
namespace engine_probe {
using Clock=std::chrono::steady_clock;
inline double ms(Clock::time_point start,Clock::time_point end=Clock::now()) {
    return std::chrono::duration<double,std::milli>(end-start).count();
}
struct Options {
    int width=320,height=180,iterations=64,repeats=3;
    double maxSeconds=120.,rowSeconds=5.;
    bool full=false;
    std::optional<amp::State> customState;
};
class TeeBuffer final:public std::streambuf {
    std::streambuf *first_,*second_;
protected:
    int_type overflow(int_type c) override {
        if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);
        if(traits_type::eq_int_type(first_->sputc(traits_type::to_char_type(c)),traits_type::eof()) ||
           traits_type::eq_int_type(second_->sputc(traits_type::to_char_type(c)),traits_type::eof()))
            return traits_type::eof();
        return c;
    }
    int sync() override {return first_->pubsync()==0&&second_->pubsync()==0?0:-1;}
public:
    TeeBuffer(std::streambuf* first,std::streambuf* second):first_(first),second_(second){}
};
inline unsigned bits(ULONG_PTR mask) {
    unsigned n=0;for(;mask;mask&=mask-1)++n;return n;
}
struct Host {
    unsigned hardware=std::max(1u,std::thread::hardware_concurrency());
    unsigned logical=hardware,physical=0,affinity=hardware;
    DWORD_PTR processMask=0,systemMask=0;
    std::string brand;
};
inline Host host() {
    Host h;
    if(GetProcessAffinityMask(GetCurrentProcess(),&h.processMask,&h.systemMask))
        h.affinity=std::max(1u,bits(h.processMask));
    DWORD bytes=0;
    GetLogicalProcessorInformationEx(RelationProcessorCore,nullptr,&bytes);
    if(bytes) {
        std::vector<unsigned char> storage(bytes);
        auto* info=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(storage.data());
        if(GetLogicalProcessorInformationEx(RelationProcessorCore,info,&bytes)) {
            unsigned logical=0,physical=0;
            for(DWORD offset=0;offset+sizeof(DWORD)*2<=bytes;) {
                const auto* item=reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(storage.data()+offset);
                if(!item->Size || item->Size>bytes-offset)break;
                if(item->Relationship==RelationProcessorCore) {
                    ++physical;
                    for(WORD group=0;group<item->Processor.GroupCount;++group)
                        logical+=bits(item->Processor.GroupMask[group].Mask);
                }
                offset+=item->Size;
            }
            if(logical)h.logical=logical;
            h.physical=physical;
        }
    }
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
    unsigned a=0,b=0,c=0,d=0;
    if(__get_cpuid_max(0x80000000u,nullptr)>=0x80000004u) {
        char brand[49]{};
        for(unsigned leaf=0;leaf<3;++leaf) {
            __cpuid(0x80000002u+leaf,a,b,c,d);
            const unsigned values[]={a,b,c,d};
            std::memcpy(brand+leaf*16,values,16);
        }
        h.brand=brand;
    }
#endif
    if(h.brand.empty()) {
        const char* value=std::getenv("PROCESSOR_IDENTIFIER");
        h.brand=value?value:"unavailable";
    }
    return h;
}
inline double processCpuMs() {
    FILETIME creation{},exit{},kernel{},user{};
    if(!GetProcessTimes(GetCurrentProcess(),&creation,&exit,&kernel,&user))return -1.;
    ULARGE_INTEGER k{},u{};k.LowPart=kernel.dwLowDateTime;k.HighPart=kernel.dwHighDateTime;
    u.LowPart=user.dwLowDateTime;u.HighPart=user.dwHighDateTime;
    return double(k.QuadPart+u.QuadPart)/10000.;
}
inline void describeAdapters(std::ostream& out) {
    amp::Com<IDXGIFactory1> factory;
    if(FAILED(CreateDXGIFactory1(__uuidof(IDXGIFactory1),reinterpret_cast<void**>(factory.put())))) {
        out<<"DXGI_ADAPTER_INVENTORY=UNAVAILABLE\n";return;
    }
    for(UINT index=0;;++index) {
        amp::Com<IDXGIAdapter1> adapter;
        const HRESULT result=factory->EnumAdapters1(index,adapter.put());
        if(result==DXGI_ERROR_NOT_FOUND)break;
        if(FAILED(result))break;
        DXGI_ADAPTER_DESC1 description{};
        if(FAILED(adapter->GetDesc1(&description)))continue;
        out<<"DXGI_ADAPTER="<<index<<" NAME="<<amp::utf8(description.Description)
           <<" VRAM_MB="<<description.DedicatedVideoMemory/(1024u*1024u)
           <<" SOFTWARE="<<((description.Flags&DXGI_ADAPTER_FLAG_SOFTWARE)!=0)<<'\n';
    }
}
inline std::uint64_t hashPixels(const std::vector<std::uint32_t>& pixels) {
    std::uint64_t hash=14695981039346656037ull;
    for(auto pixel:pixels)for(int byte=0;byte<4;++byte) {
        hash^=(pixel>>(byte*8))&255u;hash*=1099511628211ull;
    }
    return hash;
}
inline std::string hexHash(std::uint64_t hash) {
    std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;return out.str();
}
inline double median(std::vector<double> values) {
    if(values.empty())return 0.;
    std::sort(values.begin(),values.end());
    const auto middle=values.size()/2;
    return values.size()%2?values[middle]:(values[middle-1]+values[middle])*.5;
}
class Budget {
    std::mutex mutex_;std::condition_variable changed_;bool finished_=false;
    std::thread thread_;
    Clock::time_point deadline_=Clock::time_point::max();
public:
    std::atomic<std::uint64_t> revision{0};
    explicit Budget(double seconds) {
        if(seconds<=0.)return;
        deadline_=Clock::now()+std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(seconds));
        thread_=std::thread([this,seconds] {
            std::unique_lock<std::mutex> lock(mutex_);
            if(!changed_.wait_for(lock,std::chrono::duration<double>(seconds),[&]{return finished_;}))
                revision.store(1,std::memory_order_relaxed);
        });
    }
    ~Budget() {
        {std::lock_guard<std::mutex> lock(mutex_);finished_=true;}
        changed_.notify_one();if(thread_.joinable())thread_.join();
    }
    bool available() const {return revision.load(std::memory_order_relaxed)==0&&Clock::now()<deadline_;}
    double remaining() const {return deadline_==Clock::time_point::max()?1.e9:std::max(.000001,std::chrono::duration<double>(deadline_-Clock::now()).count());}
};
struct Sample {
    std::string suite,scene,backend,policy,path,status="OK";
    unsigned threads=1,effectiveWorkers=0;
    int width=0,height=0,iterations=0,repeat=0;
    unsigned frames=1;
    double total=0.,processCpu=0.,cpuWall=0.,gpuWall=0.,overlap=0.;
    double windowMs=0.,cpuWindowMs=-1.;
    double cancelReturn=-1.,drain=0.,recovery=0.;
    rhgpu::RenderMetrics gpu{};
    bool usedGpu=false,stableAlgebra=false;
    std::size_t repairs=0,rawPixels=0;
    std::uint64_t hash=0,cpuHash=0;
};
inline void csvHeader(std::ostream& out) {
    out<<"suite,scene,backend,policy,path,stable_algebra,threads,effective_workers,width,height,iterations,repeat,frames,ms_per_frame,process_cpu_ms_per_frame,window_ms,process_cpu_window_ms,cpu_counter_status,cpu_cores_equiv,cpu_wall_ms,gpu_wall_ms,prepare_ms,dispatch_ms,wait_ms,readback_ms,repair_ms,dispatches,batches,query_polls,sleep_polls,timer_waits,timer_failures,timer_timeouts,timer_recreations,timer_last_error,timer_last_status,timer_fallbacks,spin_polls,yield_polls,max_individual_wait_ms,flushes,max_batch_ms,requested_batch_pixels,target_batch_pixels,min_batch_pixels,mean_batch_pixels,max_batch_pixels,repair_workers,cpu_repairs,raw_gpu_pixels,gpu_used,cache_hit,gpu_timeline_ms,gpu_timing_valid,gpu_disjoint,hash,cpu_hash,overlap_ms,cancel_return_ms,drain_ms,recovery_ms,status\n";
}
inline void csvSample(std::ostream& out,const Sample& r) {
    const auto& g=r.gpu;
    out<<r.suite<<','<<r.scene<<','<<r.backend<<','<<r.policy<<','<<r.path<<','<<r.stableAlgebra<<','<<r.threads<<','<<r.effectiveWorkers<<','
       <<r.width<<','<<r.height<<','<<r.iterations<<','<<r.repeat<<','<<r.frames<<','
       <<std::fixed<<std::setprecision(4)<<r.total<<','<<r.processCpu<<','
       <<r.windowMs<<','<<r.cpuWindowMs<<','
       <<(r.windowMs>=100.&&r.cpuWindowMs>=15.625?"AGGREGATED":"COARSE")<<','
       <<(r.windowMs>=100.&&r.cpuWindowMs>=15.625?r.cpuWindowMs/r.windowMs:-1.)<<','<<r.cpuWall<<','<<r.gpuWall<<','
       <<g.prepareMs<<','<<g.dispatchMs<<','<<g.waitMs<<','<<g.readbackMs<<','<<g.repairMs<<','
       <<g.dispatches<<','<<g.completionBatches<<','<<g.queryPolls<<','<<g.sleepPolls<<','
       <<g.timerWaits<<','<<g.timerFailures<<','<<g.timerTimeouts<<','<<g.timerRecreations<<','<<g.timerLastError<<','<<g.timerLastStatus<<','<<g.timerFallbacks<<','<<g.spinPolls<<','<<g.yieldPolls<<','<<g.maxIndividualWaitMs<<','
       <<g.flushes<<','<<g.maxBatchMs<<','<<g.requestedBatchPixels<<','<<g.targetBatchPixels<<','<<g.minBatchPixels<<','<<g.meanBatchPixels<<','<<g.maxBatchPixels<<','<<g.repairWorkers<<','<<r.repairs<<','
       <<r.rawPixels<<','<<r.usedGpu<<','<<g.cacheHit<<','
       <<g.gpuComputeMs<<','<<g.gpuTimingValid<<','<<g.gpuTimingDisjoint<<','<<hexHash(r.hash)<<','
       <<hexHash(r.cpuHash)<<','<<r.overlap<<','<<r.cancelReturn<<','<<r.drain<<','<<r.recovery<<','<<r.status<<'\n';
    out.flush();
}
struct NamedState {std::string name; amp::State state;};
inline rhgpu::EngineTuning tuning(int profile,bool cache=false) {
    rhgpu::EngineTuning value;
    value.schedule=profile==2?rhgpu::Schedule::Adaptive:rhgpu::Schedule::Legacy;
    value.waitPolicy=profile==0?rhgpu::WaitPolicy::LegacySleep:rhgpu::WaitPolicy::Responsive;
    value.waitSliceMs=.35;
    value.targetBatchMs=4.;value.minimumBatchPixels=2048u;
    value.metricCache=cache;value.gpuTimestamps=true;return value;
}
inline const char* profileName(int profile) {
    return profile==0?"sleep_legacy":profile==1?"responsive_legacy":"responsive_bounded2048";
}
inline rhgpu::EngineTuning batchTuning(std::size_t pixels,bool stable=false,bool cache=false) {
    auto value=tuning(2,cache);value.batchPixels=pixels;value.stableCurvature=stable;
    value.captureRepairMask=true;return value;
}
inline std::string batchName(std::size_t pixels,bool stable=false) {
    return std::string(stable?"stable_batch":"reference_batch")+std::to_string(pixels);
}
inline Sample measured(amp::Renderer& engine,const NamedState& named,int width,int height,
                       unsigned threads,const std::string& suite,const std::string& backend,
                       const std::string& policy,int repeat,Budget& budget,double rowSeconds,double renderAspect=0.,
                       double minimumWindowMs=0.,std::vector<std::uint32_t>* captured=nullptr,
                       std::vector<std::uint8_t>* capturedMask=nullptr) {
    Sample r;r.suite=suite;r.scene=named.name;r.backend=backend;r.policy=policy;r.repeat=repeat;
    r.threads=threads;r.width=width;r.height=height;r.iterations=named.state.iterations;
    r.stableAlgebra=backend!="CPU"&&engine.gpuTuning().stableCurvature;
    engine.setThreadLimit(threads);
    const auto scene=amp::scene(named.state);
    std::vector<std::uint32_t> image;
    Budget rowBudget(std::min(rowSeconds,budget.remaining()));
    const double beforeCpu=processCpuMs();const auto start=Clock::now();
    bool ok=false;r.frames=0;
    do {
        ok=engine.render(scene,width,height,renderAspect>0.?renderAspect:double(height)/width,1,image,&rowBudget.revision,0);
        if(ok)++r.frames;
    } while(ok&&minimumWindowMs>0.&&ms(start)<minimumWindowMs&&rowBudget.available());
    r.windowMs=ms(start);const double afterCpu=processCpuMs();
    r.cpuWindowMs=beforeCpu>=0.&&afterCpu>=0.?std::max(0.,afterCpu-beforeCpu):-1.;
    r.total=r.windowMs/std::max(1u,r.frames);
    r.processCpu=r.cpuWindowMs<0.?-1.:r.cpuWindowMs/std::max(1u,r.frames);
    r.gpu=engine.rhGpuMetrics();r.usedGpu=engine.rhUsedGpu();r.repairs=engine.rhCpuRepairs();
    r.rawPixels=r.usedGpu&&image.size()>=r.repairs?image.size()-r.repairs:0;
    r.effectiveWorkers=backend=="CPU"?engine.cpuWorkerCount():r.gpu.repairWorkers;
    r.path=backend=="CPU"?"CPU":r.usedGpu?"GPU":
        (r.repairs==image.size()&&r.gpu.dispatches?"GPU_ALL_CPU_REPAIRS":"CPU_FALLBACK");
    if(r.path=="CPU_FALLBACK")r.effectiveWorkers=engine.cpuWorkerCount();
    if(!ok || image.size()!=std::size_t(width)*height) {
        r.status=!budget.available()?"TIME_BUDGET":!rowBudget.available()?"ROW_TIMEOUT":"INCOMPLETE";
        r.drain=engine.drainGpu();
    }
    else {
        r.hash=hashPixels(image);
        if(capturedMask)*capturedMask=backend=="GPU"?engine.rhRepairMask():std::vector<std::uint8_t>{};
        if(captured)*captured=std::move(image);
    }
    if(backend=="CPU")r.cpuWall=r.total;else r.gpuWall=r.total;
    return r;
}
inline bool warmup(amp::Renderer& engine,const NamedState& named,int width,int height,
                   unsigned threads,Budget& budget,double rowSeconds) {
    engine.setThreadLimit(threads);
    std::vector<std::uint32_t> image;
    Budget rowBudget(std::min(rowSeconds,budget.remaining()));
    const bool ok=engine.render(amp::scene(named.state),width,height,double(height)/width,1,image,&rowBudget.revision,0);
    if(!ok)engine.drainGpu();
    return ok;
}
inline void describeState(std::ostream& out,const NamedState& named) {
    const auto& s=named.state;const auto& c=s.cameras[s.plane];const auto& e=s.rhExperiment;
    out<<"VIEW="<<named.name<<" GENERATION="<<s.generation<<" DEGREE_SELECTOR="<<s.spectralOrder
       <<" ITERATIONS="<<s.iterations<<" MAP="<<e.map<<" COLOUR="<<e.coloring
       <<" CHART="<<e.chart<<" CAMERA="<<std::setprecision(17)<<c.x<<','<<c.y<<','<<c.span
       <<" C="<<s.c.real()<<','<<s.c.imag()<<" ALPHA="<<s.alpha.real()<<','<<s.alpha.imag()
       <<" GAIN="<<e.gain<<" RELAX="<<e.relaxation<<" MEMORY="<<e.memory
       <<" TWIST="<<e.twist<<" BLEND="<<e.blend<<" SEED="<<s.seed.real()<<','<<s.seed.imag()
       <<" SPECTRAL_D="<<s.spectralD<<" RH_LEVEL="<<s.rhLevel<<" POWER="<<e.power
       <<" TRAP="<<e.trapRadius<<','<<e.trapAngle<<" EXPOSURE="<<e.exposure
       <<" FREQUENCY="<<e.frequency<<" PALETTE="<<s.palette<<" HUE="<<s.hue<<" BANDS="<<s.bands<<'\n';
}
inline void summary(std::ostream& out,const std::vector<Sample>& rows,double reference=0.,const char* referenceName="reference") {
    if(rows.empty())return;
    std::vector<double> totals,cpu,wait,readback,repair,polls,dispatches,gpuTimes,overlap,cpuWall,gpuWall;
    bool hashesSame=true,gpuValid=true,cpuValid=true;const auto firstHash=rows.front().hash;
    double cpuWindow=0.,wallWindow=0.;unsigned frames=0;
    std::uint64_t timerWaits=0,timerFailures=0,timerTimeouts=0,timerRecreations=0,timerFallbacks=0,spins=0,yields=0;
    std::uint32_t lastTimerError=0,lastTimerStatus=0;double maximumWait=0.;
    for(const auto& r:rows) {
        totals.push_back(r.total);cpu.push_back(r.processCpu);wait.push_back(r.gpu.waitMs);
        readback.push_back(r.gpu.readbackMs);repair.push_back(r.gpu.repairMs);
        polls.push_back(r.gpu.queryPolls);dispatches.push_back(r.gpu.dispatches);
        gpuTimes.push_back(r.gpu.gpuComputeMs);overlap.push_back(r.overlap);
        cpuWall.push_back(r.cpuWall);gpuWall.push_back(r.gpuWall);
        cpuWindow+=r.cpuWindowMs;wallWindow+=r.windowMs;frames+=r.frames;
        cpuValid=cpuValid&&r.cpuWindowMs>=0.;
        timerWaits+=r.gpu.timerWaits;timerFailures+=r.gpu.timerFailures;timerTimeouts+=r.gpu.timerTimeouts;
        timerRecreations+=r.gpu.timerRecreations;timerFallbacks+=r.gpu.timerFallbacks;
        spins+=r.gpu.spinPolls;yields+=r.gpu.yieldPolls;maximumWait=std::max(maximumWait,r.gpu.maxIndividualWaitMs);
        if(r.gpu.timerLastError)lastTimerError=r.gpu.timerLastError;
        if(r.gpu.timerLastStatus)lastTimerStatus=r.gpu.timerLastStatus;
        hashesSame=hashesSame&&r.hash==firstHash;gpuValid=gpuValid&&r.gpu.gpuTimingValid;
    }
    const auto& r=rows.front();const double total=median(totals),cpuTime=median(cpu);
    out<<std::fixed<<std::setprecision(3)<<"ROW="<<r.suite<<'/'<<r.scene<<'/'<<r.policy
       <<" PATH="<<r.path<<" STABLE="<<r.stableAlgebra<<" T="<<r.threads<<" WORKERS="<<r.effectiveWorkers<<" SIZE="<<r.width<<'x'<<r.height
       <<" MS_PER_FRAME="<<total<<" FRAMES="<<frames<<" WINDOW_MS="<<wallWindow<<" CPU_MS_PER_FRAME="<<cpuTime;
    if(cpuValid&&wallWindow>=100.&&cpuWindow>=15.625)
        out<<" CPU_COUNTER=AGGREGATED CPU_WINDOW_MS="<<cpuWindow<<" CORES="<<cpuWindow/wallWindow;
    else out<<" CPU_COUNTER=COARSE CORES=NA";
    if(reference>0.)out<<" VS_"<<referenceName<<"_SPEEDUP="<<reference/std::max(.000001,total);
    if(r.backend!="CPU") {
        const auto fallbacks=std::count_if(rows.begin(),rows.end(),[](const Sample& row){return !row.usedGpu;});
        out<<" GPU="<<r.usedGpu<<" DISPATCH="<<median(dispatches)<<" WAIT_MS="<<median(wait)
           <<" POLLS="<<median(polls)<<" BATCHES="<<r.gpu.completionBatches<<" SLEEPS="<<r.gpu.sleepPolls
           <<" TIMER_WAITS_SUM="<<timerWaits<<" TIMER_FAILURES_SUM="<<timerFailures
           <<" TIMER_TIMEOUTS_SUM="<<timerTimeouts<<" TIMER_RECREATIONS_SUM="<<timerRecreations
           <<" TIMER_LAST_ERROR=0x"<<std::hex<<lastTimerError<<" TIMER_LAST_STATUS=0x"<<lastTimerStatus<<std::dec
           <<" TIMER_FALLBACKS_SUM="<<timerFallbacks<<" SPINS_SUM="<<spins
           <<" YIELDS_SUM="<<yields<<" MAX_WAIT_MS="<<maximumWait
           <<" REQUESTED_BATCH="<<r.gpu.requestedBatchPixels<<" TARGET_BATCH="<<r.gpu.targetBatchPixels
           <<" BATCH_PIXELS="<<r.gpu.minBatchPixels<<','<<r.gpu.meanBatchPixels<<','<<r.gpu.maxBatchPixels
           <<" FLUSHES="<<r.gpu.flushes<<" READ_MS="<<median(readback)<<" REPAIR_MS="<<median(repair)
           <<" REPAIRS="<<r.repairs<<" RAW_GPU_PIXELS="<<r.rawPixels<<" REPAIR_T="<<r.gpu.repairWorkers
           <<" GPU_TIMELINE_MS="<<(gpuValid?median(gpuTimes):-1.)
           <<" TIMESTAMP_VALID="<<gpuValid<<" CACHE="<<r.gpu.cacheHit
           <<" FALLBACK_ROWS="<<fallbacks<<" GPU_THROUGHPUT_VALID="<<(fallbacks==0);
    }
    if(r.suite=="split")out<<" CPU_WALL="<<median(cpuWall)<<" GPU_WALL="<<median(gpuWall)<<" OVERLAP_MS="<<median(overlap);
    out<<" HASH="<<hexHash(firstHash)<<" REPEAT_HASH="<<(hashesSame?"SAME":"VARIES")
       <<" N="<<rows.size()<<'\n'<<std::flush;
}
inline Sample splitMeasured(amp::Renderer& cpu,amp::Renderer& gpu,const NamedState& named,
                            int width,int height,int percent,unsigned cpuThreads,
                            int repeat,Budget& budget,double rowSeconds,unsigned gpuRepairThreads=16) {
    Sample r;r.suite="split";r.scene=named.name;r.backend="CPU_GPU";
    r.stableAlgebra=gpu.gpuTuning().stableCurvature;
    r.policy="cpu_"+std::to_string(percent)+"pct";r.width=width;r.height=height;
    r.iterations=named.state.iterations;r.repeat=repeat;r.threads=cpuThreads;
    const int left=std::clamp(width*percent/100,1,width-1),right=width-left;
    auto a=named.state,b=named.state;
    const auto source=named.state.cameras[named.state.plane];
    auto& ca=a.cameras[a.plane];auto& cb=b.cameras[b.plane];
    ca.span=source.span*double(left)/width;cb.span=source.span*double(right)/width;
    ca.x=source.x-source.span*.5+ca.span*.5;
    cb.x=source.x+source.span*.5-cb.span*.5;
    const auto cpuScene=amp::scene(a),gpuScene=amp::scene(b);
    cpu.setThreadLimit(cpuThreads);gpu.setThreadLimit(gpuRepairThreads);
    Budget rowBudget(std::min(rowSeconds,budget.remaining()));
    std::vector<std::uint32_t> cpuImage,gpuImage;
    bool cpuOk=false,gpuOk=false;std::exception_ptr cpuFailure,gpuFailure;
    Clock::time_point cpuStart{},cpuEnd{},gpuStart{},gpuEnd{};
    const double beforeCpu=processCpuMs();const auto start=Clock::now();
    std::thread worker([&] {
        cpuStart=Clock::now();
        try {cpuOk=cpu.render(cpuScene,left,height,double(height)/left,1,cpuImage,&rowBudget.revision,0);}
        catch(...) {cpuFailure=std::current_exception();}
        cpuEnd=Clock::now();
    });
    gpuStart=Clock::now();
    try {gpuOk=gpu.render(gpuScene,right,height,double(height)/right,1,gpuImage,&rowBudget.revision,0);}
    catch(...) {gpuFailure=std::current_exception();}
    gpuEnd=Clock::now();worker.join();r.total=ms(start);
    const double afterCpu=processCpuMs();
    r.processCpu=beforeCpu>=0.&&afterCpu>=0.?std::max(0.,afterCpu-beforeCpu):-1.;
    r.windowMs=r.total;r.cpuWindowMs=r.processCpu;r.frames=cpuOk&&gpuOk?1u:0u;
    if(cpuFailure)std::rethrow_exception(cpuFailure);
    if(gpuFailure)std::rethrow_exception(gpuFailure);
    r.cpuWall=ms(cpuStart,cpuEnd);r.gpuWall=ms(gpuStart,gpuEnd);
    r.overlap=std::max(0.,ms(std::max(cpuStart,gpuStart),std::min(cpuEnd,gpuEnd)));
    r.gpu=gpu.rhGpuMetrics();r.usedGpu=gpu.rhUsedGpu();r.repairs=gpu.rhCpuRepairs();r.effectiveWorkers=cpu.cpuWorkerCount();
    r.rawPixels=r.usedGpu&&gpuImage.size()>=r.repairs?gpuImage.size()-r.repairs:0;
    r.path=r.usedGpu?"CPU_GPU":"CPU_GPU_FALLBACK";
    if(!cpuOk||!gpuOk) {r.status=!budget.available()?"TIME_BUDGET":!rowBudget.available()?"ROW_TIMEOUT":"INCOMPLETE";r.drain=gpu.drainGpu();}
    else {r.hash=hashPixels(gpuImage);r.cpuHash=hashPixels(cpuImage);}
    return r;
}
struct Difference {
    std::size_t pixels=0,different=0,aboveTolerance=0;
    std::uint64_t channelError=0;int maximum=0;
    void add(int maximumError,unsigned errorSum) {
        ++pixels;different+=maximumError!=0;aboveTolerance+=maximumError>8;
        maximum=std::max(maximum,maximumError);channelError+=errorSum;
    }
    double mean() const {return pixels?double(channelError)/(3.*pixels):0.;}
};
struct Comparison {Difference full,repaired,raw;bool valid=false,maskValid=false;};
inline void printDifference(std::ostream& out,const char* prefix,const Difference& difference) {
    out<<' '<<prefix<<"MAX_RGB="<<difference.maximum<<' '<<prefix<<"MEAN_CHANNEL="<<difference.mean()
       <<' '<<prefix<<"DIFFERENT="<<difference.different<<'/'<<difference.pixels
       <<' '<<prefix<<"ABOVE_TOL8="<<difference.aboveTolerance;
}
inline Comparison compareImages(std::ostream& out,const NamedState& named,int width,int height,
                                const std::string& policy,const std::vector<std::uint32_t>& reference,
                                const std::vector<std::uint32_t>& image,
                                const std::vector<std::uint8_t>& mask={}) {
    Comparison comparison;
    if(reference.size()!=std::size_t(width)*height||image.size()!=reference.size()) {
        out<<"IMAGE_COMPARE="<<named.name<<'/'<<policy<<" SIZE="<<width<<'x'<<height
           <<" STATUS=UNAVAILABLE\n";return comparison;
    }
    comparison.valid=true;comparison.maskValid=mask.size()==image.size();
    for(std::size_t i=0;i<image.size();++i) {
        int maximum=0;unsigned sum=0;
        for(int shift:{0,8,16}) {
            const int error=std::abs(int((reference[i]>>shift)&255u)-int((image[i]>>shift)&255u));
            sum+=unsigned(error);maximum=std::max(maximum,error);
        }
        comparison.full.add(maximum,sum);
        if(comparison.maskValid)(mask[i]?comparison.repaired:comparison.raw).add(maximum,sum);
    }
    out<<std::fixed<<std::setprecision(5)<<"IMAGE_COMPARE="<<named.name<<'/'<<policy
       <<" SIZE="<<width<<'x'<<height<<" MASK="<<(comparison.maskValid?"AVAILABLE":"UNAVAILABLE");
    printDifference(out,"",comparison.full);
    if(comparison.maskValid) {
        printDifference(out,"REPAIRED_",comparison.repaired);printDifference(out,"RAW_",comparison.raw);
    }
    out<<" REFERENCE_HASH="<<hexHash(hashPixels(reference))<<" IMAGE_HASH="<<hexHash(hashPixels(image))
       <<" PURPOSE=FULL_VIEW_DIAGNOSTIC\n"<<std::flush;return comparison;
}
inline void cancellationProbe(amp::Renderer& engine,const NamedState& named,int width,int height,
                              unsigned threads,const std::string& backend,const std::string& policy,
                              std::ostream& out,std::ostream& csv,double rowSeconds,Budget* outerBudget=nullptr) {
    auto heavy=named.state;
    heavy.iterations=std::max(128,heavy.iterations);
    // A new camera deliberately invalidates cached orbit data.
    heavy.cameras[heavy.plane].x+=.00001;
    const auto scene=amp::scene(heavy);
    engine.setThreadLimit(threads);
    std::atomic<std::uint64_t> revision{0};
    std::mutex mutex;std::condition_variable changed;bool finished=false;
    Clock::time_point requested{};
    bool sent=false;
    std::thread canceller([&] {
        std::unique_lock<std::mutex> lock(mutex);
        if(!changed.wait_for(lock,std::chrono::milliseconds(2),[&]{return finished;})) {
            requested=Clock::now();sent=true;revision.store(1,std::memory_order_relaxed);
        }
    });
    Sample r;r.suite="cancel";r.scene=named.name;r.backend=backend;r.policy=policy;
    r.stableAlgebra=backend=="GPU"&&engine.gpuTuning().stableCurvature;
    r.threads=threads;r.width=width;r.height=height;r.iterations=heavy.iterations;
    std::vector<std::uint32_t> image;
    const double beforeCpu=processCpuMs();const auto start=Clock::now();
    bool ok=false;std::exception_ptr failure;
    try {ok=engine.render(scene,width,height,double(height)/width,1,image,&revision,0);}
    catch(...) {failure=std::current_exception();}
    const auto returned=Clock::now();r.total=ms(start,returned);
    const double afterCpu=processCpuMs();
    r.processCpu=beforeCpu>=0.&&afterCpu>=0.?std::max(0.,afterCpu-beforeCpu):-1.;
    r.windowMs=r.total;r.cpuWindowMs=r.processCpu;r.frames=ok?1u:0u;
    {std::lock_guard<std::mutex> lock(mutex);finished=true;}
    changed.notify_one();canceller.join();
    if(failure)std::rethrow_exception(failure);
    r.gpu=engine.rhGpuMetrics();r.usedGpu=engine.rhUsedGpu();r.repairs=engine.rhCpuRepairs();
    r.effectiveWorkers=backend=="CPU"?engine.cpuWorkerCount():r.gpu.repairWorkers;
    r.path=backend=="CPU"?"CPU":"GPU_CANCEL_PROBE";
    if(sent&&requested>=start&&requested<=returned&&!ok) {
        r.status="CANCELLED";r.cancelReturn=ms(requested,returned);
    } else r.status=ok?"COMPLETED_BEFORE_CANCEL":"NO_INFLIGHT_CANCEL";
    // Only the calling thread ever operates this immediate D3D context.
    if(backend=="GPU")r.drain=engine.drainGpu();
    auto recovery=named.state;recovery.iterations=32;
    recovery.cameras[recovery.plane]=amp::Camera{0.,0.,4.};
    const auto recoveryScene=amp::scene(recovery);
    if(outerBudget&&!outerBudget->available())r.status+="_RECOVERY_SKIPPED_BUDGET";
    else {
        Budget recoveryBudget(outerBudget?std::min(rowSeconds,outerBudget->remaining()):rowSeconds);
        const auto recoveryStart=Clock::now();
        const bool recovered=engine.render(recoveryScene,80,48,48./80.,1,image,&recoveryBudget.revision,0);
        r.recovery=ms(recoveryStart);
        if(!recovered) {r.status+="_RECOVERY_INCOMPLETE";r.drain+=engine.drainGpu();}
    }
    csvSample(csv,r);
    out<<std::fixed<<std::setprecision(3)<<"CANCEL="<<backend<<'/'<<policy<<" T="<<threads
       <<" STATUS="<<r.status<<" REQUEST_TO_RETURN_MS="<<r.cancelReturn
       <<" GPU_DRAIN_MS="<<r.drain<<" RECOVERY_AFTER_DRAIN_MS="<<r.recovery
       <<" SUBMITTED_DISPATCHES="<<r.gpu.dispatches<<'\n'<<std::flush;
}
inline int run(std::ostream& output,const std::filesystem::path& directory,const Options& options=Options{}) {
    if(options.width<8||options.height<4||options.width>4096||options.height>4096||
       std::uint64_t(options.width)*options.height>4000000ull||options.iterations<1||
       options.iterations>1024||options.repeats<1||options.repeats>9||
       !std::isfinite(options.maxSeconds)||options.maxSeconds<0.||
       !std::isfinite(options.rowSeconds)||options.rowSeconds<=0.)
        throw std::runtime_error("Engine probe: invalid size, iterations, repeats or time budget");
    std::filesystem::create_directories(directory);
    const auto reportPath=directory/"engine-report.txt",csvPath=directory/"engine-samples.csv";
    std::ofstream reportFile(reportPath),csv(csvPath);
    if(!reportFile||!csv)throw std::runtime_error("Cannot write engine diagnostic output files");
    TeeBuffer tee(output.rdbuf(),reportFile.rdbuf());std::ostream out(&tee);csvHeader(csv);
    try {
        const Host machine=host();const unsigned available=std::min({64u,machine.hardware,machine.affinity});
        const unsigned cpuThreads=std::min(32u,available),hybridGpu=std::min(16u,std::max(1u,available/2));
        const unsigned hybridCpu=available>hybridGpu+2?available-hybridGpu-2:1u;
        out<<"ENGINE_PROBE=1.8.0 MODE="<<(options.full?"FULL":"QUICK")
           <<" BASE_SIZE="<<options.width<<'x'<<options.height<<" REPEATS="<<options.repeats
           <<" TIME_BUDGET_S="<<options.maxSeconds<<" ROW_LIMIT_S="<<options.rowSeconds<<'\n'
           <<"CPU="<<machine.brand<<'\n'
           <<"CPU_LOGICAL_SYSTEM="<<machine.logical<<" CPU_PHYSICAL_SYSTEM="<<machine.physical
           <<" CPU_HARDWARE_CONCURRENCY="<<machine.hardware<<" PROCESS_AFFINITY_LOGICAL="<<machine.affinity
           <<" PROCESS_AFFINITY_MASK=0x"<<std::hex<<machine.processMask<<std::dec<<'\n'
           <<"CPU_THREADS="<<cpuThreads<<" GPU_REPAIR_BUDGET="<<cpuThreads
           <<" HYBRID_CPU_BUDGET="<<hybridCpu<<" HYBRID_GPU_REPAIR_BUDGET="<<hybridGpu<<'\n'
           <<"POLICY=Targeted third sweep: CPU32, explicit GPU queue sizes, reference versus experimental stable curvature. No GUI autotuning.\n"
           <<"CPU_ACCOUNTING=Complete CPU frames repeat for >=100ms. MS_PER_FRAME is batch mean; coarse CPU counters suppress CORES.\n"
           <<"GPU_TIMELINE=Dispatch timestamp span includes gaps between host batches; valid flag required.\n"
           <<"BATCH_POLICY=Requested queue pixels are explicit; workload safety cap and edge batches remain visible as target/min/mean/max.\n"
           <<"STABLE_CURVATURE=Opt-in experimental algebra only; reference path remains available.\n"
           <<"MASK_CAPTURE=Enabled for GPU probe frames; its byte-mask capture cost is included equally in reference/stable timing.\n"
           <<"QUALITY=RGB comparisons occur after timing. Repair mask1 means CPU repaired; mask0 means raw GPU. Hybrid requires full-frame MAX_RGB<=8 against CPU and a measured speed advantage.\n"
           <<"CANCELLED_SAMPLES_SCORED=0\n"<<std::flush;
        describeAdapters(out);out<<"GPU_GROUP_SHAPE=8x4 GPU_THREADS_PER_GROUP=32\n";
        amp::Renderer cpu,unavailableGpu;
        std::vector<NamedState> scenes;
        for(int bank:{1,2}) {
            auto state=amp::defaults(20,1);focus_presets::apply(state,bank,0);state.iterations=options.iterations;
            scenes.push_back({bank==1?"Cayley_currents":"Tricorn_weave",state});
        }
        if(options.customState) {
            if(options.customState->family!=20||options.customState->iterations<1||options.customState->iterations>1024)
                throw std::runtime_error("Custom probe view must use RH Lab and 1..1024 iterations");
            scenes.push_back({"User_view",*options.customState});
        }
        for(const auto& named:scenes)describeState(out,named);
        const unsigned startupCpuThreads=std::min(8u,available>2?available-2:1u);
        cpu.init(true,startupCpuThreads);
        auto startupState=scenes.front().state;startupState.iterations=32;
        const auto startupScene=amp::scene(startupState);
        rhbackground::BackgroundGpu gpuPreparation;
        out<<"STAGE=BACKGROUND_GPU_PREPARATION CPU_LOGIC_CONTINUES=1 CPU_WORKER_LIMIT="<<startupCpuThreads
           <<" CPU_FRAME_SIZE=160x90 CPU_ITERATIONS="<<startupState.iterations<<" CPU_WINDOW_LIMIT_S=30\n"<<std::flush;
        const auto startupBegin=Clock::now();const double startupCpuBegin=processCpuMs();
        gpuPreparation.start(false,cpuThreads);
        std::vector<std::pair<Clock::time_point,Clock::time_point>> startupWork;
        std::vector<std::uint32_t> startupPixels;
        bool startupCancelled=false;
        {
            Budget startupBudget(30.);
            while(!gpuPreparation.ready()&&startupBudget.available()) {
                const auto begin=Clock::now();
                const bool ok=cpu.render(startupScene,160,90,90./160.,1,startupPixels,&startupBudget.revision,0);
                const auto end=Clock::now();
                if(!ok) {startupCancelled=true;break;}
                startupWork.emplace_back(begin,end);
            }
        }
        const bool preparationExceededWindow=!gpuPreparation.ready();
        if(preparationExceededWindow)
            out<<"BACKGROUND_CPU_WINDOW=ENDED JOINING_COMPILER=1 (CPU window limit does not terminate D3DCompile)\n"<<std::flush;
        const auto joinBegin=Clock::now();const bool haveGpu=gpuPreparation.wait();const double joinMs=ms(joinBegin);
        const double startupCpuEnd=processCpuMs();
        const auto compileBegin=gpuPreparation.startedAt(),compileEnd=gpuPreparation.finishedAt();
        unsigned beforeReady=0;double cpuWorkMs=0.,hostOverlapMs=0.;
        for(const auto& interval:startupWork) {
            cpuWorkMs+=ms(interval.first,interval.second);
            beforeReady+=interval.second<=compileEnd;
            const auto begin=std::max(interval.first,compileBegin),end=std::min(interval.second,compileEnd);
            if(end>begin)hostOverlapMs+=ms(begin,end);
        }
        out<<std::fixed<<std::setprecision(3)<<"BACKGROUND_GPU_READY="<<haveGpu
           <<" COMPILE_AND_DEVICE_PREP_MS="<<gpuPreparation.elapsedMs()
           <<" CPU_FRAMES="<<startupWork.size()<<" CPU_FRAMES_BEFORE_READY="<<beforeReady
           <<" CPU_WORKERS="<<cpu.cpuWorkerCount()<<" CPU_RENDER_WORK_MS="<<cpuWorkMs
           <<" HOST_INTERVAL_OVERLAP_MS="<<hostOverlapMs<<" JOIN_WAIT_MS="<<joinMs
           <<" STARTUP_HOST_MS="<<ms(startupBegin)
           <<" STARTUP_PROCESS_CPU_MS="<<(startupCpuBegin>=0.&&startupCpuEnd>=0.?startupCpuEnd-startupCpuBegin:-1.)
           <<" CPU_CANCELLED="<<startupCancelled<<" CPU_HASH="<<hexHash(startupCancelled?0:hashPixels(startupPixels))<<'\n'
           <<"BACKGROUND_OVERLAP_MEASURE=CPU_RENDER_HOST_INTERVALS_VS_PREPARATION_INTERVAL; process CPU time includes compiler and renderer.\n"
           <<"GPU_BACKEND_REUSED="<<haveGpu<<" DUPLICATE_RH_SHADER_COMPILES=0\n"<<std::flush;
        if(!haveGpu) {
            unavailableGpu.init(true,cpuThreads);unavailableGpu.diagnostic=gpuPreparation.error();
            out<<"BACKGROUND_GPU_ERROR="<<gpuPreparation.error()<<'\n';
        }
        amp::Renderer& gpu=haveGpu?*gpuPreparation.renderer():unavailableGpu;
        cpu.setThreadLimit(cpuThreads);gpu.setThreadLimit(cpuThreads);
        out<<"GPU_DEVICE="<<gpu.backend<<" GPU_AVAILABLE="<<haveGpu<<" GPU_INIT="<<gpu.diagnostic<<'\n';
        if(haveGpu) {
            // Compilation completed on the background device-only path. These
            // first GPU commands only warm model data and the two algebra flags.
            auto state=scenes.front().state;state.iterations=4;state.alpha=amp::C{};
            state.c=amp::C(.05,.03);state.cameras[state.plane]=amp::Camera{.1,.1,.02};
            Budget shaderBudget(30.);out<<"STAGE=GPU_DATA_AND_STABLE_WARMUP COMPILATION_ALREADY_DONE=1 TIME_LIMIT_S=30\n"<<std::flush;
            for(bool stable:{false,true}) {
                gpu.setGpuTuning(batchTuning(8192,stable));
                auto row=measured(gpu,{stable?"Stable_warmup":"Reference_warmup",state},8,4,1,
                    "data_warmup","GPU",batchName(8192,stable),0,shaderBudget,30.);
                csvSample(csv,row);out<<"GPU_DATA_WARMUP="<<row.status<<" STABLE="<<stable<<" PATH="<<row.path
                    <<" MS="<<row.total<<" DIAGNOSTIC="<<gpu.rhGpuDiagnostic()<<'\n'<<std::flush;
            }
        }
        Budget budget(options.maxSeconds);unsigned completedRows=0;bool incomplete=false;
        auto collect=[&](amp::Renderer& engine,const NamedState& named,int width,int height,
                         const std::string& suite,const std::string& backend,const std::string& policy,
                         bool warm=true,int repeatLimit=0,std::vector<std::uint32_t>* image=nullptr,
                         std::vector<std::uint8_t>* mask=nullptr) {
            std::vector<Sample> rows;if(!budget.available())return rows;
            if(warm) {
                auto row=measured(engine,named,width,height,cpuThreads,"warmup_"+suite,backend,policy,0,budget,options.rowSeconds);
                if(row.status!="OK") {
                    csvSample(csv,row);incomplete=true;
                    out<<"SKIP="<<named.name<<'/'<<policy<<" SIZE="<<width<<'x'<<height<<" REASON=WARMUP_"<<row.status
                       <<" WALL_MS="<<row.windowMs<<" DISPATCHES="<<row.gpu.dispatches<<" WAIT_MS="<<row.gpu.waitMs
                       <<" DRAIN_MS="<<row.drain<<'\n'<<std::flush;return rows;
                }
            }
            const int wanted=repeatLimit>0?repeatLimit:options.repeats;
            for(int repeat=1;repeat<=wanted&&budget.available();++repeat) {
                auto row=measured(engine,named,width,height,cpuThreads,suite,backend,policy,repeat,budget,
                    options.rowSeconds,0.,backend=="CPU"?100.:0.,image,mask);
                csvSample(csv,row);
                if(row.status!="OK") {
                    incomplete=true;out<<"SKIP="<<named.name<<'/'<<policy<<" SIZE="<<width<<'x'<<height
                       <<" REASON="<<row.status<<'\n'<<std::flush;break;
                }
                if(backend=="GPU"&&!row.usedGpu)
                    out<<"FALLBACK="<<named.name<<'/'<<policy<<" PATH="<<row.path
                       <<" CPU_WORKERS="<<row.effectiveWorkers<<" DIAGNOSTIC="<<engine.rhGpuDiagnostic()<<'\n';
                rows.push_back(row);
            }
            if(rows.size()!=std::size_t(wanted))incomplete=true;
            if(!rows.empty())++completedRows;return rows;
        };
        auto frameMedian=[](const std::vector<Sample>& rows) {
            std::vector<double> values;for(const auto& row:rows)values.push_back(row.total);return median(values);
        };
        auto hashesStable=[](const std::vector<Sample>& rows) {
            return !rows.empty()&&std::all_of(rows.begin(),rows.end(),[&](const Sample& row){return row.hash==rows.front().hash;});
        };
        struct Winner {
            bool valid=false;std::size_t scene=0,batch=0;int width=0,height=0;double cpuMs=0.,gpuMs=0.;
        } winner;
        auto runScale=[&](int multiplier,bool extraQueueCases) {
            const std::uint64_t wide=std::uint64_t(options.width)*multiplier,high=std::uint64_t(options.height)*multiplier;
            if(wide>4096||high>4096||wide*high>4000000ull) {
                out<<"SCALE="<<multiplier<<" SKIP=DIMENSION_LIMIT\n";return;
            }
            const int width=int(wide),height=int(high);
            for(std::size_t sceneIndex=0;sceneIndex<scenes.size()&&budget.available();++sceneIndex) {
                const auto& named=scenes[sceneIndex];std::vector<std::uint32_t> cpuImage,referenceImage,stableImage;
                std::vector<std::uint8_t> referenceMask,stableMask;
                auto cpuRows=collect(cpu,named,width,height,"scale","CPU","threads"+std::to_string(cpuThreads),true,0,&cpuImage);
                summary(out,cpuRows);const double cpuMs=frameMedian(cpuRows);const bool cpuHashValid=hashesStable(cpuRows);
                if(!haveGpu)continue;
                const std::vector<std::size_t> batches=multiplier==4?std::vector<std::size_t>{32768u}:
                    std::vector<std::size_t>{8192u,32768u};
                for(std::size_t batch:batches) {
                    if(!budget.available())break;
                    referenceImage.clear();referenceMask.clear();gpu.setGpuTuning(batchTuning(batch,false));
                    auto rows=collect(gpu,named,width,height,"queue","GPU",batchName(batch),true,0,&referenceImage,&referenceMask);
                    summary(out,rows,cpuMs,"CPU32");
                    if(!rows.empty())compareImages(out,named,width,height,batchName(batch)+"_vsCPU",cpuImage,referenceImage,referenceMask);
                    stableImage.clear();stableMask.clear();gpu.setGpuTuning(batchTuning(batch,true));
                    auto stableRows=collect(gpu,named,width,height,"algebra","GPU",batchName(batch,true),true,0,&stableImage,&stableMask);
                    summary(out,stableRows,cpuMs,"CPU32");
                    const auto quality=compareImages(out,named,width,height,batchName(batch,true)+"_vsCPU",cpuImage,stableImage,stableMask);
                    if(!referenceImage.empty()&&!stableImage.empty())
                        compareImages(out,named,width,height,batchName(batch,true)+"_vsReferenceGPU",referenceImage,stableImage,stableMask);
                    const double gpuMs=frameMedian(stableRows);
                    const bool actualGpu=!stableRows.empty()&&std::all_of(stableRows.begin(),stableRows.end(),[](const Sample& row){return row.usedGpu;});
                    const bool qualityMatch=quality.valid&&quality.full.maximum<=8&&cpuHashValid&&hashesStable(stableRows)
                        &&cpuRows.size()==std::size_t(options.repeats)&&stableRows.size()==std::size_t(options.repeats);
                    const bool eligible=qualityMatch&&actualGpu&&cpuMs>0.&&gpuMs>0.&&gpuMs<cpuMs;
                    out<<"HYBRID_CANDIDATE="<<named.name<<'/'<<batchName(batch,true)<<" SIZE="<<width<<'x'<<height
                       <<" MATCHED_CPU_TOL8="<<qualityMatch<<" FASTER_THAN_CPU="<<(cpuMs>0.&&gpuMs>0.&&gpuMs<cpuMs)
                       <<" ACTUAL_GPU="<<actualGpu<<" ELIGIBLE="<<eligible<<'\n';
                    if(eligible&&(!winner.valid||cpuMs/gpuMs>winner.cpuMs/winner.gpuMs))
                        winner={true,sceneIndex,batch,width,height,cpuMs,gpuMs};
                }
                if(extraQueueCases)for(std::size_t batch:{2048u,65536u}) {
                    if(!budget.available())break;
                    gpu.setGpuTuning(batchTuning(batch,false));referenceImage.clear();referenceMask.clear();
                    auto rows=collect(gpu,named,width,height,"queue_limits","GPU",batchName(batch),true,0,&referenceImage,&referenceMask);
                    summary(out,rows,cpuMs,"CPU32");
                    if(!rows.empty())compareImages(out,named,width,height,batchName(batch)+"_vsCPU",cpuImage,referenceImage,referenceMask);
                }
            }
        };
        runScale(1,true);runScale(2,false);
        if(options.full&&budget.available()&&budget.remaining()>30.)runScale(4,false);
        else out<<"OPTIONAL_1280_SCALE=SKIPPED BUDGET_OR_QUICK_MODE\n";
        if(haveGpu&&budget.available()) {
            // One large cached recolouring case isolates host queue overhead.
            const auto& named=scenes.front();const int width=std::min(1280,options.width*2),height=std::min(720,options.height*2);
            gpu.setGpuTuning(batchTuning(65536,false,true));
            auto warm=measured(gpu,named,width,height,cpuThreads,"warmup_cache","GPU","reference_batch65536",0,budget,options.rowSeconds);
            if(warm.status=="OK") {
                auto recolour=named;recolour.state.hue+=.17;recolour.state.rhExperiment.exposure=.8;
                recolour.state.rhExperiment.frequency+=.4;
                auto rows=collect(gpu,recolour,width,height,"cache","GPU","large_batch65536_recolour",false);
                summary(out,rows);
            } else {csvSample(csv,warm);incomplete=true;}
            if(budget.available()) {
                gpu.setGpuTuning(batchTuning(32768,true));
                cancellationProbe(gpu,named,std::max(640,options.width),std::max(360,options.height),cpuThreads,
                    "GPU","stable_batch32768",out,csv,options.rowSeconds,&budget);
            }
        }
        if(options.full&&winner.valid&&haveGpu&&budget.available()) {
            const auto& named=scenes[winner.scene];
            const double cpuScale=double(cpuThreads)/hybridCpu;
            const int cpuPercent=std::clamp(int(std::lround(100.*winner.gpuMs/(winner.cpuMs*cpuScale+winner.gpuMs))),10,90);
            gpu.setGpuTuning(batchTuning(winner.batch,true));
            out<<"HYBRID_POLICY=ONE_EXPERIMENT MATCHED_CPU_MAX_RGB_LE8=1 CPU_PERCENT="<<cpuPercent
               <<" FRACTION_SOURCE=MEASURED_COST_ESTIMATE (no GUI changes)\n";
            auto row=splitMeasured(cpu,gpu,named,winner.width,winner.height,cpuPercent,hybridCpu,1,budget,options.rowSeconds,hybridGpu);
            csvSample(csv,row);
            if(row.status=="OK") {summary(out,{row},winner.cpuMs,"CPU32");++completedRows;}else incomplete=true;
        } else out<<"HYBRID=SKIPPED REASON="<<(!options.full?"QUICK_MODE":!winner.valid?"NO_FASTER_GPU_WITH_FULL_CPU_TOL8_QUALITY":"BUDGET_OR_NO_GPU")<<'\n';
        if(haveGpu)gpu.drainGpu();
        out<<"ENGINE_PROBE="<<(!budget.available()||incomplete?"PARTIAL":"COMPLETE")
           <<" COMPLETED_ROWS="<<completedRows<<" CANCELLED_SAMPLES_SCORED=0\n"
           <<"REPORT="<<reportPath.string()<<"\nSAMPLES="<<csvPath.string()<<'\n'<<std::flush;
        csv.flush();reportFile.flush();
        if(!out||!csv||!reportFile)throw std::runtime_error("Engine probe output could not be completed");
        return 0;
    } catch(const std::exception& error) {
        out<<"ENGINE_PROBE=ERROR REASON="<<error.what()<<"\nREPORT="<<reportPath.string()<<'\n'<<std::flush;return 1;
    }
}
} // namespace engine_probe
