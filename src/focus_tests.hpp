#pragma once
#include "amp_core.hpp"
#include "focus_presets.hpp"
#include <chrono>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>

namespace focus_tests {
inline bool run(std::ostream& out) {
    using C=amp::C;
    bool pass=true;
    auto check=[&](bool ok,const char* label,int index=-1){
        if(!ok){out<<"FOCUS_TEST_FAIL="<<label<<" INDEX="<<index<<'\n';pass=false;}
    };
    const C initial(.12,-.09),parameter(-.3,.1),strength(.2,.11);
    // Independent Decimal80 reference; sin/cos(.4) are binary64 constants.
    const std::array<C,6> first{{
        C(-.4505640454510873,-.21532989644401251),C(-.5516085342627164,-.16675595277540636),
        C(-.5516085342627164,-.16675595277540636),C(-.33257037949746104,-.2379020883422114),
        C(.1483718624380169,.236150931287365),C(-.28604570345445524,-.22548395013859118)}};
    const std::array<C,6> second{{
        C(-.5126891407910501,.004547132111146501),C(-.6097849665151267,.04425079266799332),
        C(-.6840404917017755,-.14363108831603125),C(-.4324264139054102,-.14892460880524303),
        C(.18806314735343183,.29155308558132476),C(-.40379687613984444,-.13207371368976142)}};
    auto settings=[](){rhmaster::Experiment e;e.power=3;e.gain=.8;e.relaxation=.7;
        e.memory=.15;e.twist=.4;e.blend=.7;return e;};
    auto orbit=[&](rhmaster::Experiment e,int count){
        return rhmaster::model(1,16,.6,.9,.28,strength,0.,e)->orbit(initial,parameter,count);
    };
    for(int map=6;map<12;++map){
        auto e=settings();e.map=map;
        const auto a=orbit(e,1),b=orbit(e,2);
        check(a.status==rhmaster::Status::Unresolved&&a.iterations==1
              &&std::abs(a.endpoint-first[map-6])<2e-10,"DECIMAL80_FIRST",map);
        check(b.status==rhmaster::Status::Unresolved&&b.iterations==2
              &&std::abs(b.endpoint-second[map-6])<2e-10,"DECIMAL80_SECOND",map);
    }
    auto reduction=[&](int from,int to,double blend,double twist,double memory,const char* label){
        auto a=settings();a.map=from;a.blend=blend;a.twist=twist;a.memory=memory;
        auto b=a;b.map=to;
        const auto x=orbit(a,2),y=orbit(b,2);
        check(x.status==y.status&&x.iterations==y.iterations
              &&rhmaster::finitePoint(x.endpoint)&&rhmaster::finitePoint(y.endpoint)
              &&std::abs(x.endpoint-y.endpoint)<2e-12*(1.+std::abs(y.endpoint)),label,from);
    };
    reduction(8,5,.7,0.,0.,"ZERO_MEMORY_CAYLEY");
    reduction(9,5,.7,0.,.15,"ZERO_TWIST_CAYLEY");
    reduction(10,0,0.,.4,.15,"BLEND_ZERO_CLASSIC");
    reduction(10,2,1.,.4,.15,"BLEND_ONE_TRICORN");
    reduction(11,5,0.,0.,.15,"MIXED_CAYLEY_ZERO");
    reduction(11,6,1.,0.,.15,"MIXED_CAYLEY_CONJUGATE");
    auto phase=settings();phase.map=7;phase.gain=0.;
    const auto periodic=rhmaster::model(1,16,.6,.9,.28,C{},0.,phase)->orbit(C{},C{},80);
    check(periodic.detectedPeriod==2,"ALTERNATING_CAYLEY_PHASE");
    auto e=settings();e.map=9;
    const auto cached=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
    e.twist+=.1;const auto twisted=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
    check(cached!=twisted&&twisted->experiment.twist==e.twist,"TWIST_CACHE");
    e.blend+=.1;const auto mixed=rhmaster::model(1,16,.6,.9,.28,strength,0.,e);
    check(mixed!=twisted&&mixed->experiment.blend==e.blend,"BLEND_CACHE");
    e={};e.blend=1.01;check(!rhmaster::validExperiment(e),"INVALID_BLEND");
    e={};e.twist=4.;check(!rhmaster::validExperiment(e),"INVALID_TWIST");
    for(int bank=0;bank<focus_presets::BankCount;++bank)
        for(int index=0;index<focus_presets::PresetsPerBank;++index){
            auto s=amp::defaults(20,7);s.seed=C(.7,1.1);s.spectralD=.2;s.rhLevel=2.;s.spectralOrder=3;
            focus_presets::apply(s,bank,index);
            check(s.generation==7&&s.seed==C(.7,1.1)&&s.spectralD==.2&&s.rhLevel==2.
                  &&s.spectralOrder==3&&s.plane==1&&rhmaster::validExperiment(s.rhExperiment)
                  &&amp::scene(s).masterModel->validModel,"PRESET_PRESERVES_DRIVER",bank*4+index);
        }
    {
        auto s=amp::defaults(20,1);focus_presets::apply(s,1,0);
        check(s.rhExperiment.map==5&&s.c==C(-.56,.23)&&s.alpha==C(.16,0.)
              &&s.rhExperiment.coloring==5&&s.rhExperiment.frequency==3.,"FAVOURITE_CURRENTS");
        focus_presets::apply(s,2,0);
        check(s.rhExperiment.map==2&&s.c==C(-.12,.72)&&s.alpha==C(.18,0.)
              &&s.rhExperiment.coloring==3&&s.rhExperiment.frequency==2.5
              &&s.rhExperiment.trapAngle==.4,"FAVOURITE_TRICORN");
    }
    const auto suffix=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    const auto path=std::filesystem::temp_directory_path()/("fraktalfaktur_focus_"+suffix+".view");
    try{
        auto state=amp::defaults(20,1);state.rhExperiment=settings();state.rhExperiment.map=11;
        amp::writeView(path,state);
        const auto restored=amp::readView(path);
        check(rhmaster::sameExperiment(restored.rhExperiment,state.rhExperiment),"VIEW16_FOCUS_ROUNDTRIP");
        std::ifstream input(path);const std::string text((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());input.close();
        std::istringstream lines(text);std::string line,old="FRAKTALFAKTUR_VIEW 15\n";
        std::getline(lines,line);
        while(std::getline(lines,line)){
            if(line.rfind("rh_twist ",0)==0||line.rfind("rh_blend ",0)==0)continue;
            if(line.rfind("rh_map ",0)==0)line="rh_map 5";
            old+=line+'\n';
        }
        auto write=[&](const std::string& value){std::ofstream file(path);file<<value;};
        auto rejected=[&](const std::string& value){write(value);try{(void)amp::readView(path);return false;}
            catch(const std::exception&){return true;}};
        write(old);const auto migrated=amp::readView(path);
        check(migrated.rhExperiment.map==5&&migrated.rhExperiment.twist==0.
              &&migrated.rhExperiment.blend==.5,"VIEW15_DEFAULT_MIGRATION");
        check(rejected(old+"rh_twist 0.1\n"),"VIEW15_NEW_FIELD_REJECTED");
        check(rejected(text+"rh_blend 0.2\n"),"VIEW16_DUPLICATE_REJECTED");
        auto missing=text;const auto begin=missing.find("rh_twist "),end=missing.find('\n',begin);
        missing.erase(begin,end-begin+1);
        check(rejected(missing),"VIEW16_MISSING_TWIST_REJECTED");
    }catch(const std::exception& ex){out<<"FOCUS_TEST_ERROR="<<ex.what()<<'\n';pass=false;}
    std::error_code ignored;std::filesystem::remove(path,ignored);
    out<<"FOCUS_INTEGRATION="<<(pass?"PASS":"FAIL")<<'\n';return pass;
}
} // namespace focus_tests
