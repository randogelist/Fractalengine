#pragma once
#include "amp_gpu.hpp"
#include <ostream>

namespace engine_controls_tests {
inline bool run(std::ostream& out) {
    bool good=true;
    rhgpu::Engine engine;
    rhgpu::EngineTuning tuning;
    good=good&&rhgpu::validTuning(tuning);
    good=good&&tuning.waitPolicy==rhgpu::WaitPolicy::Responsive;
    tuning.schedule=rhgpu::Schedule::Legacy;
    tuning.waitPolicy=rhgpu::WaitPolicy::LegacySleep;
    tuning.metricCache=false;
    engine.setTuning(tuning);
    good=good&&engine.tuning().schedule==rhgpu::Schedule::Legacy&&!engine.tuning().metricCache&&
        engine.tuning().waitPolicy==rhgpu::WaitPolicy::LegacySleep;
    tuning.schedule=rhgpu::Schedule::Adaptive;
    tuning.waitPolicy=rhgpu::WaitPolicy::Responsive;
    tuning.minimumBatchPixels=2048;
    engine.setTuning(tuning);
    good=good&&engine.tuning().minimumBatchPixels==2048;
    tuning.targetBatchMs=0.;
    good=good&&!rhgpu::validTuning(tuning);
    bool rejected=false;
    try { engine.setTuning(tuning); }
    catch(const std::invalid_argument&) { rejected=true; }
    good=good&&rejected;
    tuning.targetBatchMs=4.;
    tuning.waitSliceMs=0.;
    good=good&&!rhgpu::validTuning(tuning);
    tuning.waitSliceMs=std::numeric_limits<double>::quiet_NaN();
    good=good&&!rhgpu::validTuning(tuning);
    tuning.waitSliceMs=.35;
    tuning.waitPolicy=static_cast<rhgpu::WaitPolicy>(-1);
    good=good&&!rhgpu::validTuning(tuning);
    tuning.waitPolicy=rhgpu::WaitPolicy::Responsive;
    good=good&&rhgpu::validTuning(tuning);
    tuning.batchPixels=31;
    good=good&&!rhgpu::validTuning(tuning);
    tuning.batchPixels=65537;
    good=good&&!rhgpu::validTuning(tuning);
    tuning.batchPixels=32768;
    tuning.captureRepairMask=true;
    tuning.stableCurvature=true;
    good=good&&rhgpu::validTuning(tuning);
    engine.setTuning(tuning);
    good=good&&engine.tuning().batchPixels==32768&&engine.tuning().captureRepairMask&&
        engine.tuning().stableCurvature&&engine.repairMask().empty();
    amp::Renderer renderer;
    renderer.setThreadLimit(100);
    good=good&&renderer.threadLimit()==64;
    renderer.setThreadLimit(0);
    good=good&&renderer.threadLimit()==0;
    out<<"ENGINE_CONTROLS="<<(good?"PASS":"FAIL")<<'\n';
    return good;
}
} // namespace engine_controls_tests
