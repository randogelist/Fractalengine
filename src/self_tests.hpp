#pragma once
#include "rh_candidates.hpp"
#include <array>
#include <cmath>
#include <ostream>

namespace rhfresh {
// Independently generated from an 80-decimal-digit polynomial reference.
// Reference input is the stored binary64 coefficient table, degree16.
inline bool runSelfTests(std::ostream& out) {
    const std::array<C,OriginalCandidateCount> expected{{
        C(-3.479070092220968,-.8023061629127474),
        C(.2437546429021348,.025820329116761243),
        C(.7870376428757812,.4312009567199926),
        C(.3786293916154225,.30098159625875204),
        C(.2629030448324114,.00914252184307315),
        C(-.34113148670729837,.7879404256130155),
        C(-.875116445205336,-.18952528882415984),
        C(-2.191423549686277,-.680625764656546),
        C(.21052907372328675,.016723448921424657),
        C(.08110231073894471,.004166849040393446)
    }};
    bool pass=true;
    auto check=[&](bool ok,const char* label,int index) {
        if(!ok) {out<<"SELF_TEST_FAIL="<<label<<" CANDIDATE="<<index+1<<'\n';pass=false;}
    };
    auto finite=[](C z){return std::isfinite(z.real())&&std::isfinite(z.imag());};
    const C sample(.31,.27);
    for(int i=0;i<CandidateCount;++i) {
        auto settings=defaultSettings(i);
        Model model(settings);
        check(model.valid(),"MODEL",i);
        if(!model.valid())continue;
        const C observed=model.correction(sample);
        if(i<OriginalCandidateCount)
            check(finite(observed)&&std::abs(observed-expected[i])<2e-7*(1.+std::abs(expected[i])),"REFERENCE_CORRECTION",i);
        else check(finite(observed),"NEW_BRANCH_FINITE_CORRECTION",i);
        const C reflected=model.correction(std::conj(sample));
        check(finite(reflected)&&std::abs(reflected-std::conj(observed))<1e-9*(1.+std::abs(observed)),"CONJUGATION",i);
        const C initial(.12,-.09),parameter(-.3,.1);
        const auto baseline=model.orbit(initial,parameter,C{},1);
        check(baseline.iterations==1&&std::abs(baseline.endpoint-(initial*initial+parameter))<1e-14,"ZERO_LAMBDA_BASELINE",i);
        const C strength(.2,.11);
        const auto one=model.orbit(initial,parameter,strength,1);
        const C next=initial*initial+parameter-strength*model.correction(initial);
        check(finite(one.endpoint)&&std::abs(one.endpoint-next)<1e-12*(1.+std::abs(next)),"COMMON_UPDATE",i);
        for(int degree: {8,16,24,32}) {
            settings.degree=degree;
            Model varied(settings);
            check(varied.valid()&&finite(varied.correction(sample)),"SUPPORTED_DEGREE",i);
        }
        out<<"SELF_TEST_CANDIDATE="<<registry()[i].id<<'\n';
    }
    // Challenge points probe cancellation near the critical boundary and outside the unit disk.
    const std::array<std::array<C,OriginalCandidateCount>,4> challenge{{
        {{C(-0.11354048685200728,0.072653401291636013),C(-2.4614968744554955,13.38084856635316),C(0.030248224109151886,-0.19579899278602478),C(-0.058369330288978921,0.10216199639552527),C(-0.078868211057335477,-0.064292560103699575),C(-0.14728038777701846,-0.013751701306380175),C(-0.0086359127701448636,0.069583576285098647),C(-0.010603853739611729,0.074148420120744435),C(-0.073456101664676421,-0.059392042625147121),C(-0.025716326498011168,-0.020667207402978019)}},
        {{C(-0.49134459462881669,0.31837791846389252),C(56.383206895223871,204.73360085857712),C(-1.2669949365588222,3.1966441816232489),C(-1.5895409806143448,-0.89636961785317026),C(0.17777809960163202,-0.087125858219074898),C(-0.14730727685059744,-0.013777572083657243),C(-0.020643993112640816,0.24784813520890331),C(-0.046834302263698543,0.31839896481969926),C(0.17859628848055401,-0.12385524050591357),C(0.059804572753646354,-0.032544595254263002)}},
        {{C(-0.49117671395407486,0.31850806554753458),C(56.385187636069766,204.72276137249818),C(-1.2626735941011562,3.2029887822443457),C(-1.5896011000929713,-0.89633551235831876),C(0.1778036437050399,-0.087154498057566226),C(-0.14730727679545547,-0.013777572004143278),C(-0.020856382305948114,0.24748395127075357),C(-0.046909275366608688,0.31780599806068921),C(0.1786244380325312,-0.12389595835200096),C(0.059812804562461477,-0.032557675062550286)}},
        {{C(-10.664867083754098,-3.821597682087754),C(0.31956505243847033,0.026696170466485227),C(6.0286435600823376,3.7978109481549271),C(3.3303118923215091,1.716286902561061),C(0.32614588404144895,0.025950742363910581),C(84833.28864692041,-132150.44196790538),C(-2.6958748594732351,-1.0119803383961528),C(-8.6529028898599023,-3.4776386875494527),C(0.29949508848933004,0.030893249508176866),C(0.10564369946010417,0.009262098101157348)}},
    }};
    const std::array<int,4> challengeDegree{{16,24,32,32}};
    const std::array<C,4> challengePoint{{C(-1.86,-.66),C(-1.86,-.66),C(-1.86,-.66),C(2.7,1.3)}};
    for(int row=0;row<4;++row) for(int i=0;i<OriginalCandidateCount;++i) {
        auto settings=defaultSettings(i);settings.degree=challengeDegree[row];Model model(settings);
        const C observed=model.correction(challengePoint[row]);const C target=challenge[row][i];
        check(model.valid()&&finite(observed)&&std::abs(observed-target)<2e-6*(1.+std::abs(target)),"BOUNDARY_CONDITIONING",i);
    }
    out<<"CANDIDATE_COUNT="<<CandidateCount<<'\n';
    out<<"SELF_TEST="<<(pass?"PASS":"FAIL")<<'\n';
    return pass;
}
} // namespace rhfresh
