#pragma once
#include "amp_core.hpp"
#include "catalogue.hpp"
#include <chrono>
#include <fstream>
#include <iterator>
#include <ostream>

namespace fresh_classic_tests {
inline bool run(std::ostream& out) {
    bool pass=true;
    auto check=[&](bool ok,const char* label,int generation=0){
        if(!ok){out<<"CLASSIC_TEST_FAIL="<<label<<" GENERATION="<<generation<<'\n';pass=false;}
    };
    const auto entries=catalogue::builtins();
    check(entries.size()==rhfresh::CandidateCount+frk_extra::driver_count,"CATALOGUE_COUNT_29");
    for(int g=1;g<=rhfresh::CandidateCount;++g){
        const auto s=amp::defaults(20,g);
        const auto d=rhfresh::defaultSettings(g-1);
        check(s.generation==g&&s.seed.real()==d.radius&&s.rhLevel==d.level,"CANDIDATE_DEFAULTS",g);
        const auto scene=amp::scene(s);
        const auto model=scene.masterModel;
        check(model&&model->validModel&&model->generation==g,"MODEL_DISPATCH",g);
        if(!model||!model->validModel)continue;
        check(entries[g-1].id==rhfresh::registry()[g-1].id,"REGISTRY_ID",g);
        rhfresh::Model reference(d);
        const amp::C initial(.12,-.09),parameter(-.3,.1),strength(.2,.11);
        const auto step=model->orbitAtLambda(initial,parameter,strength,1);
        const auto expected=initial*initial+parameter-strength*reference.correction(initial);
        check(std::abs(step.endpoint-expected)<2e-12*(1.+std::abs(expected)),"ORIGINAL_LOOP_NEW_RESPONSE",g);
        const auto baseline=model->orbitAtLambda(initial,parameter,amp::C{},1);
        check(std::abs(baseline.endpoint-(initial*initial+parameter))<1e-14,"QUADRATIC_BASELINE",g);
        const auto same=rhmaster::model(g,16,d.radius,d.scale,d.shift,d.lambda,d.level);
        const auto sameAgain=rhmaster::model(g,16,d.radius,d.scale,d.shift,d.lambda,d.level);
        check(same==sameAgain,"MODEL_CACHE_REUSE",g);
        const auto levelEdit=rhmaster::model(g,16,d.radius,d.scale,d.shift,d.lambda,d.level+.125);
        check(levelEdit!=same&&levelEdit->level==d.level+.125,"MODEL_CACHE_LEVEL",g);
        const auto radiusEdit=rhmaster::model(g,16,d.radius+.05,d.scale,d.shift,d.lambda,d.level+.125);
        check(radiusEdit!=levelEdit&&radiusEdit->x==d.radius+.05,"MODEL_CACHE_RADIUS",g);
        out<<"CLASSIC_CANDIDATE="<<g<<" PASS="<<(model->validModel?"YES":"NO")<<'\n';
    }
    const auto suffix=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    const auto temp=std::filesystem::temp_directory_path()/("fraktalfaktur_classic_"+suffix+".view");
    try{
        auto s=amp::defaults(20,rhfresh::CandidateCount);
        s.seed.real(.75);s.rhLevel=3.125;s.spectralOrder=3;s.alpha=amp::C(.3,.15);
        amp::writeView(temp,s);
        const auto restored=amp::readView(temp);
        check(restored.family==20&&restored.generation==rhfresh::CandidateCount&&restored.rhLevel==s.rhLevel
              &&restored.seed==s.seed&&restored.alpha==s.alpha&&restored.spectralOrder==3,"CURRENT_VIEW_ROUNDTRIP");
        std::ifstream input(temp);std::string contents((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
        input.close();
        const auto pos=contents.find("FRAKTALFAKTUR_VIEW 16");
        check(pos==0,"VIEW16_HEADER");
        if(pos==0){
            contents.replace(0,std::string("FRAKTALFAKTUR_VIEW 16").size(),"FRAKTALFAKTUR_VIEW 9");
            {std::ofstream changed(temp);changed<<contents;}
            bool rejected=false;try{(void)amp::readView(temp);}catch(const std::exception&){rejected=true;}
            check(rejected,"OLD_ANGULAR_VIEW_REJECTED");
        }
    }catch(const std::exception& e){out<<"CLASSIC_TEST_ERROR="<<e.what()<<'\n';pass=false;}
    std::error_code ignored;std::filesystem::remove(temp,ignored);
    out<<"CLASSIC_INTEGRATION="<<(pass?"PASS":"FAIL")<<'\n';
    return pass;
}
} // namespace fresh_classic_tests
