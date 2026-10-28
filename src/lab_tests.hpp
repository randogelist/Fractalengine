#pragma once
#include "amp_core.hpp"
#include <chrono>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>

// Run on the user's Windows build. These fixtures were derived independently
// with Decimal arithmetic at 80 digits; no native executable was run to set them.
namespace lab_tests {
inline bool run(std::ostream& out) {
    using C=amp::C;
    bool pass=true;
    auto check=[&](bool ok,const char* label,int index=-1){
        if(!ok){out<<"LAB_TEST_FAIL="<<label<<" INDEX="<<index<<'\n';pass=false;}
    };
    // Independently optimized expressions may differ by a few rounding units.
    // Test numerical agreement, never bit identity across call contexts.
    constexpr double pointTolerance=64.*std::numeric_limits<double>::epsilon();
    auto nearPoint=[&](C a,C b){
        return rhmaster::finitePoint(a)&&rhmaster::finitePoint(b)
            &&std::abs(a-b)<=pointTolerance*(1.+std::abs(b));
    };
    auto nearColour=[](const rhmaster::Color& a,const rhmaster::Color& b){
        for(std::size_t k=0;k<a.size();++k)
            if(!std::isfinite(a[k])||!std::isfinite(b[k])||a[k]<0.||a[k]>1.
               ||b[k]<0.||b[k]>1.||std::abs(a[k]-b[k])>1e-12)return false;
        return true;
    };
    const C z0(.12,-.09),parameter(-.3,.1),strength(.2,.11);
    const std::array<C,6> expected{{
        C(.13978473000596947,.3858967865033045),
        C(.39179839746391665,.44868831437911505),
        C(.1923022091897494,.24723327538523737),
        C(-.477667966445858,-.1176915447916586),
        C(.14494009435920185,.41548000743725144),
        C(-.6135215956041903,-.1355717132746136)
    }};
    rhmaster::Experiment e;
    e.power=3;e.gain=.8;e.relaxation=.7;
    for(int map=0;map<static_cast<int>(expected.size());++map){
        e.map=map;
        const auto m=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
        const auto r=m->orbit(z0,parameter,2);
        check(r.iterations==2&&std::abs(r.endpoint-expected[map])<2e-10,
              "DECIMAL80_MAP_TWO_STEPS",map);
        if(map==0||map==3||map==4)
            check(std::abs(m->orbit(z0,parameter,1).endpoint-C(.1690987081260227,.19174496127568494))<2e-10,
                  "FIRST_STEP_SIGN_AND_ZERO_MEMORY",map);
    }
    const C display(.7,.4);
    for(int chart=0;chart<rhmaster::ChartCount;++chart){
        const auto physical=rhmaster::chartForward(display,chart);
        check(std::abs(rhmaster::chartInverse(physical,chart,display)-display)<3e-13,
              "CHART_ROUNDTRIP",chart);
    }
    check(!rhmaster::finitePoint(rhmaster::chartForward(C{},1))&&
          !rhmaster::finitePoint(rhmaster::chartForward(C{},3))&&
          !rhmaster::finitePoint(rhmaster::chartInverse(C{},2)),"CHART_POLES");
    const C logDisplay(.2,8.1);
    check(std::abs(rhmaster::chartInverse(rhmaster::chartForward(logDisplay,2),2,logDisplay)-logDisplay)<3e-13,
          "LOG_NEAREST_BRANCH");
    check(std::abs(rhmaster::chartInverse(C(1e10,0.),3,C(5e-11,0.))-C(5e-11,0.))<1e-22,
          "JOUKOWSKI_SMALL_SHEET");
    e={};
    const auto original=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
    const auto endpoint=original->orbit(z0,parameter,8).endpoint;
    for(int chart=0;chart<rhmaster::ChartCount;++chart){
        e.chart=chart;
        const auto m=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
        check(nearPoint(m->orbit(z0,parameter,8).endpoint,endpoint),"CHART_LEAVES_PHYSICAL_ORBIT",chart);
        // Isolate chart composition from sensitive long rational trajectories.
        // A smooth one-step, zero-gain endpoint colour detects a missing or
        // doubled chart without magnifying tiny differences in exp rounding.
        auto chartSettings=e;chartSettings.gain=0.;chartSettings.coloring=1;
        const auto chartModel=rhmaster::model(1,16,.6,.9,.28,strength,0.,chartSettings);
        auto cartesianSettings=chartSettings;cartesianSettings.chart=0;
        const auto cartesianModel=rhmaster::model(1,16,.6,.9,.28,strength,0.,cartesianSettings);
        const auto physical=rhmaster::chartForward(display,chart);
        const auto step=chartModel->orbit(physical,parameter,1);
        check(step.status==rhmaster::Status::Unresolved&&step.iterations==1
              &&nearPoint(step.endpoint,physical*physical+parameter),"CHART_PHYSICAL_BASELINE",chart);
        const auto a=chartModel->juliaColour(display,parameter,1,0,.2,1.);
        const auto b=cartesianModel->juliaColour(physical,parameter,1,0,.2,1.);
        check(nearColour(a,b),"CHART_COLOUR_PULLBACK",chart);
    }
    for(int color=0;color<rhmaster::ColouringCount;++color){
        e={};e.coloring=color;
        const auto m=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
        for(C point:std::array<C,3>{{z0,C{},C(18.,.2)}}){
            const auto c=m->juliaColour(point,parameter,64,2,.3,1.1);
            for(double v:c)check(std::isfinite(v)&&v>=0.&&v<=1.,"FINITE_COLOUR",color);
        }
    }
    e={};e.map=3;
    const auto alternating=rhmaster::model(1,16,.6,.9,.28,C{},0.,e)->orbit(C{},C{},80);
    check(alternating.detectedPeriod==2,"ALTERNATING_PHASE_PERIOD");
    e={};e.gain=0.;
    const auto bypass=rhmaster::model(1,16,.6,.9,.28,strength,0.,e)->orbit(z0,parameter,1);
    check(bypass.status==rhmaster::Status::Unresolved&&bypass.iterations==1
          &&nearPoint(bypass.endpoint,z0*z0+parameter),"ZERO_GAIN_BASELINE");
    e={};e.gain=std::numeric_limits<double>::quiet_NaN();
    check(!rhmaster::validExperiment(e),"REJECT_NAN");
    e={};e.map=rhmaster::MapCount;
    check(!rhmaster::validExperiment(e),"REJECT_MAP_BOUND");
    e={};e.relaxation=2.;
    check(!rhmaster::validExperiment(e),"REJECT_RELAXATION_BOUND");
    const auto cached=rhmaster::model(1,16,.6,.9,.28,strength,0.);
    e={};e.frequency=2.;
    const auto changed=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
    check(changed!=cached&&changed==rhmaster::model(1,16,.6,.9,.28,strength,0.,e),"EXPERIMENT_CACHE");

    const auto suffix=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    const auto path=std::filesystem::temp_directory_path()/("fraktalfaktur_lab_"+suffix+".view");
    try{
        auto s=amp::defaults(20,10);
        s.rhExperiment={};auto& v=s.rhExperiment;
        v.map=5;v.chart=3;v.coloring=7;v.power=6;v.gain=.7;v.relaxation=.6;
        v.memory=.3;v.trapRadius=1.2;v.trapAngle=-.4;v.exposure=1.4;v.frequency=2.3;
        s.rhLevel=3.125;
        amp::writeView(path,s);
        const auto restored=amp::readView(path);
        check(rhmaster::sameExperiment(restored.rhExperiment,s.rhExperiment)&&restored.rhLevel==s.rhLevel,
              "VIEW16_ROUNDTRIP");
        std::ifstream input(path);
        const std::string text((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
        input.close();
        check(text.rfind("FRAKTALFAKTUR_VIEW 16\n",0)==0,"VIEW16_HEADER");
        auto replaceFile=[&](const std::string& body){std::ofstream file(path);file<<body;};
        auto rejected=[&](const std::string& body){
            replaceFile(body);try{(void)amp::readView(path);return false;}catch(const std::exception&){return true;}
        };
        const std::array<std::string,13> keys{{"rh_map","rh_chart","rh_coloring","rh_power","rh_gain",
            "rh_relaxation","rh_memory","rh_trap_radius","rh_trap_angle","rh_exposure","rh_frequency","rh_twist","rh_blend"}};
        std::istringstream lines(text);std::string line,classic="FRAKTALFAKTUR_VIEW 14\n";
        std::getline(lines,line);
        while(std::getline(lines,line)){
            bool lab=false;for(const auto& key:keys)if(line.rfind(key+" ",0)==0)lab=true;
            if(!lab)classic+=line+'\n';
        }
        replaceFile(classic);
        const auto old=amp::readView(path);
        check(rhmaster::sameExperiment(old.rhExperiment,rhmaster::Experiment{})&&old.rhLevel==s.rhLevel,
              "VIEW14_DEFAULT_MIGRATION");
        check(rejected(classic+"rh_map 1\n"),"VIEW14_NEW_FIELD_REJECTED");
        check(rejected(text+"rh_gain 1\n"),"VIEW16_DUPLICATE_REJECTED");
        const auto pos=text.find("rh_gain "),end=text.find('\n',pos);
        std::string missing=text;missing.erase(pos,end-pos+1);
        check(rejected(missing),"VIEW16_MISSING_FIELD_REJECTED");
        std::string invalid=text;invalid.replace(pos,end-pos,"rh_gain 9");
        check(rejected(invalid),"VIEW16_RANGE_REJECTED");
    }catch(const std::exception& ex){out<<"LAB_TEST_ERROR="<<ex.what()<<'\n';pass=false;}
    std::error_code ignored;std::filesystem::remove(path,ignored);
    out<<"LAB_INTEGRATION="<<(pass?"PASS":"FAIL")<<'\n';return pass;
}
} // namespace lab_tests
