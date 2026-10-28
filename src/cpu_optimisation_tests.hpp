#pragma once
#include "rh_master_math.hpp"
#include <ostream>

namespace cpu_optimisation_tests {
inline bool run(std::ostream& out){
    using C=rhmaster::C;
    using Diagnostics=rhmaster::OrbitDiagnostics;
    bool pass=true;
    auto check=[&](bool value,const char* label,int index=-1){
        if(!value){out<<"CPU_OPTIMISATION_FAIL="<<label<<" INDEX="<<index<<'\n';pass=false;}
    };
    // The reference evaluator forms each falling factorial in the original
    // inner loop. Cached derivatives must give the same direct/reverse jet.
    const std::array<rhfresh::detail::LC,6> samples{{
        {0.L,0.L},{.2L,-.3L},{1.L,0.L},{-1.L,0.L},{-1.0000001L,.00001L},{17.L,11.L}}};
    int jets=0;
    for(int degree:{4,8,16,32})for(int candidate:{0,3,5,6,7}){
        const auto data=rhfresh::detail::coefficientData(candidate,degree,.6);
        check(data&&data->valid,"COEFFICIENTS_VALID",candidate);
        if(!data||!data->valid)continue;
        const int order=rhfresh::detail::requiredThetaOrder(candidate);
        for(const auto y:samples){
            const bool reverse=std::abs(y)>1.L;
            const auto reference=rhfresh::detail::polynomialJets(data->theta,y,reverse);
            const auto cached=rhfresh::detail::polynomialJets(*data,y,reverse);
            for(int k=0;k<=order;++k){
                const long double tolerance=32.L*std::numeric_limits<long double>::epsilon()
                    *std::max(1.L,std::abs(reference[k]));
                check(std::abs(cached[k]-reference[k])<=tolerance,"CACHED_DERIVATIVE_JET",candidate);
                ++jets;
            }
            for(int k=order+1;k<=4;++k)
                check(cached[k]==rhfresh::detail::LC{},"UNUSED_DERIVATIVE_OMITTED",candidate);
        }
    }
    // Model instances retain distinct settings but reuse their polynomial
    // payload through level, scale, shift and strength changes.
    for(int candidate=0;candidate<rhfresh::CandidateCount;++candidate){
        auto settings=rhfresh::defaultSettings(candidate);
        rhfresh::Model before(settings);
        settings.level+=.125;settings.scale+=.01;settings.shift-=.02;settings.lambda+=C(.03,.01);
        rhfresh::Model after(settings);
        check(before.valid()&&after.valid()&&before.coefficientData()==after.coefficientData(),
              "IMMUTABLE_COEFFICIENT_REUSE",candidate);
        settings.radius+=.05;
        rhfresh::Model radiusEdit(settings);
        check(radiusEdit.valid()&&radiusEdit.coefficientData()!=after.coefficientData(),
              "RADIUS_INVALIDATES_COEFFICIENTS",candidate);
    }
    // Normalization is a colour observation, not a map update. Check axes,
    // signed directions, subnormal input and values with overflowing norm.
    const double tiny=std::numeric_limits<double>::denorm_min();
    const std::array<C,10> directions{{C(1.,0.),C(-1.,0.),C(0.,1.),C(0.,-1.),
        C(tiny,0.),C(tiny,-tiny),C(1e308,-1e308),C(1e-308,1e-309),C(-.7,.3),C(.2,-.8)}};
    for(const C z:directions){
        const C expected=std::polar(1.,std::arg(z));
        const C normalized=rhmaster::normalizedDirection(z);
        check(rhmaster::finitePoint(normalized)&&std::abs(normalized-expected)<2e-15,
              "ROBUST_FEEDBACK_DIRECTION");
        check(std::abs(std::abs(normalized)-1.)<5e-16,"UNIT_FEEDBACK_DIRECTION");
    }
    int orbits=0;
    for(int map=0;map<rhmaster::MapCount;++map)for(int coloring=0;coloring<rhmaster::ColouringCount;++coloring){
        rhmaster::Experiment experiment;experiment.map=map;experiment.coloring=coloring;
        experiment.twist=.2;experiment.blend=.37;
        rhmaster::Model model(1,16,.6,.9,.28,C(.15,.07),0.,experiment);
        std::vector<C> fullTrail,renderTrail;
        const auto full=model.orbit(C(.06,-.04),C(-.2,.1),96,&fullTrail);
        const auto render=model.orbit(C(.06,-.04),C(-.2,.1),96,&renderTrail,Diagnostics::Colour);
        const int index=map*rhmaster::ColouringCount+coloring;
        check(full.status==render.status&&full.endpoint==render.endpoint
              &&full.iterations==render.iterations&&full.smooth==render.smooth
              &&full.phase==render.phase&&fullTrail==renderTrail,"RENDER_PRESERVES_ORBIT",index);
        if(coloring==0)check(full.structure==render.structure,"ORIGINAL_STRUCTURE_PRESERVED",index);
        else check(render.structure==0.,"UNUSED_STRUCTURE_OMITTED",index);
        if(coloring==6)check(full.detectedPeriod==render.detectedPeriod
              &&full.periodResidual==render.periodResidual,"PERIOD_COLOUR_PRESERVED",index);
        else check(render.detectedPeriod==0,"UNUSED_PERIOD_OMITTED",index);
        check(full.circleTrap==render.circleTrap&&full.crossTrap==render.crossTrap
              &&full.winding==render.winding&&full.meanLogRadius==render.meanLogRadius
              &&full.feedbackMean==render.feedbackMean&&full.feedbackPhase==render.feedbackPhase
              &&full.feedbackPeak==render.feedbackPeak,"SELECTED_COLOUR_METRICS_PRESERVED",index);
        ++orbits;
    }
    rhmaster::Experiment periodExperiment;periodExperiment.coloring=3;
    rhmaster::Model periodModel(1,16,.6,.9,.28,C{},0.,periodExperiment);
    const auto report=periodModel.orbit(C{},C{},80);
    const auto render=periodModel.orbit(C{},C{},80,nullptr,Diagnostics::Colour);
    check(report.detectedPeriod==1&&report.periodRepeats==4&&report.structure>0.,
          "PUBLIC_REPORT_RETAINS_FULL_DIAGNOSTICS");
    check(render.detectedPeriod==0&&render.structure==0.,"RENDER_SKIPS_UNUSED_DIAGNOSTICS");
    rhmaster::Model lambdaModel(1,16,.6,.9,.28,C{},0.);
    const auto lambda=lambdaModel.orbitAtLambda(C{},C{},C{},80,nullptr,Diagnostics::LambdaColour);
    check(lambda.detectedPeriod==1&&lambda.periodRepeats==4&&lambda.structure>0.,
          "LAMBDA_ORIGINAL_RETAINS_PERIOD_COLOUR");
    out<<"CPU_OPTIMISATION_JETS="<<jets<<" ORBIT_CASES="<<orbits<<'\n';
    out<<"CPU_OPTIMISATION="<<(pass?"PASS":"FAIL")<<'\n';
    return pass;
}
} // namespace cpu_optimisation_tests
