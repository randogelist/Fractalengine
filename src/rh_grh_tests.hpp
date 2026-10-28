#pragma once
#include "garden_math.hpp"
#include "grh_observable_math.hpp"
#include <ostream>
#include <stdexcept>

namespace rh_grh_tests {
using Z = std::complex<long double>;
struct Jet { Z f{}, first{}, second{}; };
inline Jet add(Jet a, Jet b) {
    return {a.f+b.f,a.first+b.first,a.second+b.second};
}
inline Jet multiply(Jet a, Jet b) {
    return {a.f*b.f,a.first*b.f+a.f*b.first,
            a.second*b.f+2.L*a.first*b.first+a.f*b.second};
}
inline Jet scale(Jet a, Z c) { return {a.f*c,a.first*c,a.second*c}; }
inline std::pair<Jet,Jet> unscaled(Z z,int generation) {
    Jet p{z,1.L,0.L},q{1.L,0.L,0.L};
    for(int k=0;k<generation;k++) {
        Jet sum=add(p,q);
        p=scale(multiply(p,q),4.L);
        q=multiply(sum,sum);
    }
    return {p,q};
}
inline int run(std::ostream &out) {
    int checks=0;
    auto check=[&](bool ok,const char *name) {
        ++checks;
        if(!ok) throw std::runtime_error(std::string("RH/GRH garden test: ")+name);
    };
    for(int n=1;n<=4;n++) {
        auto data=garden::data(5,n);
        check(data->degree==(1<<n)&&data->records.size()==1,"Cayley data and upload record");
        for(int k=0;k<6;k++) {
            garden::C z(.15+1.4*std::sin(k*.73),.09+1.25*std::cos(k*.47));
            Z zz(z.real(),z.imag());
            auto normalized=cayley::preimage(z,n);
            auto raw=unscaled(zz,n);
            Z commonScale=raw.second.f/Z(normalized.q.f.real(),normalized.q.f.imag());
            for(int component=0;component<2;component++) {
                auto j=component ? normalized.q : normalized.p;
                auto value=[&](Z x) { auto v=unscaled(x,n); return component ? v.second.f : v.first.f; };
                long double h1=1e-6L,h2=1e-4L;
                Z first=(value(zz+h1)-value(zz-h1))/(2*h1);
                Z second=(value(zz+h2)-2.L*value(zz)+value(zz-h2))/(h2*h2);
                Z expectedFirst=commonScale*Z(j.first.real(),j.first.imag());
                Z expectedSecond=commonScale*Z(j.second.real(),j.second.imag());
                check(std::abs(first-expectedFirst)<3e-8L*(1+std::abs(expectedFirst)),
                      "Cayley first derivative versus finite difference");
                check(std::abs(second-expectedSecond)<3e-5L*(1+std::abs(expectedSecond)),
                      "Cayley second derivative versus finite difference");
            }
            garden::C rational=z;
            for(int j=0;j<n;j++) rational=4.*rational/((1.+rational)*(1.+rational));
            check(std::abs(normalized.p.f/normalized.q.f-rational)<3e-10*(1+std::abs(rational)),
                  "Homogeneous pair versus direct rational iteration");
            const garden::C target(.68,.24),relaxation(.81,-.17);
            Jet polynomial=add(raw.first,scale(raw.second,-Z(target.real(),target.imag())));
            for(int method=0;method<2;method++) {
                auto step=garden::correction(z,target,relaxation,*data,method);
                Z expected=Z(relaxation.real(),relaxation.imag())*(method ?
                    2.L*polynomial.f*polynomial.first/(2.L*polynomial.first*polynomial.first-polynomial.f*polynomial.second) :
                    polynomial.f/polynomial.first);
                check(step.valid&&std::abs(expected-Z(step.delta.real(),step.delta.imag()))<2e-10L*(1+std::abs(expected)),
                      "Relaxed Newton/Halley step versus unscaled polynomial jets");
            }
        }
        for(garden::C z:{garden::C(1e7,0),garden::C(-1e7,1e7)}) {
            auto v=cayley::preimage(z,n);
            bool finite=true;
            for(auto x:{v.p.f,v.p.first,v.p.second,v.q.f,v.q.first,v.q.second})
                finite=finite&&std::isfinite(x.real())&&std::isfinite(x.imag());
            check(finite,"Common normalization keeps large-coordinate jets finite");
        }
        for(int method=0;method<2;method++) {
            auto root=garden::orbit({1,0},{1,0},{1,0},*data,method,20);
            check(root.status==0&&root.iterations==0,"Exact target-one multiple root");
            auto zero=garden::orbit({0,0},{0,0},{1,0},*data,method,20);
            check(zero.status==0&&zero.iterations==0,"Exact target-zero root");
            if(n>1) {
                auto multiple=garden::orbit({-1,0},{0,0},{1,0},*data,method,20);
                check(multiple.status==0&&multiple.iterations==0,"Exact repeated target-zero root");
            }
        }
    }
    auto first=garden::data(5,1);
    for(garden::C a:{garden::C(0,0),garden::C(.68,.24),garden::C(1,0)}) {
        auto e=garden::evaluate({-1,0},a,{1,0},*first);
        check(std::abs(std::abs(e.f)/e.magnitude-1)<1e-14,"Cayley pole is not a false root");
    }
    for(int n=1;n<=4;n++) {
        auto data=garden::data(6,n);
        check(data->degree==0&&data->records.size()==1,"L-function data and upload record");
        for(garden::C z:{garden::C(-8.1,0),garden::C(12.1,0),garden::C(.5,80.1),garden::C(.5,-80.1)})
            check(garden::orbit(z,{0,0},{1,0},*data,0,20).status==3,"L-function window excludes unevaluated points");
    }
    auto zeta=garden::data(6,1);
    auto pole=garden::correction({1,0},{0,0},{1,0},*zeta,0);
    check(!pole.valid,"Zeta pole cannot be accepted as root");
    auto grh=garden::data(7,3);
    check(grh->degree==0&&grh->records.size()==1,"Primitive GRH Cayley data and upload record");
    {
        garden::C s0(.5,12.6749464170114);
        garden::C w=std::tanh(grh_cayley::log2*(s0-garden::C(.5,0)));
        garden::C mapped;
        check(grh_cayley::sheetPoint(w,3,mapped)&&std::abs(mapped-s0)<3e-12,
              "Cayley sheet inverse recovers primitive zero height");
        auto e=grh_cayley::evaluate(w,{0,0},3);
        check(e.valid&&e.residual<2e-11,"Primitive chi5 zero is a Cayley zero attractor");
        auto o=garden::orbit(w,{0,0},{1,0},*grh,0,8);
        check(o.status==0,"Primitive GRH Cayley accepts known numerical zero");
    }
    for (garden::C w:{garden::C(.23,.41),garden::C(-.31,.72),garden::C(.18,-.55)}) {
        auto n=grh_cayley::evaluate(w,{0,0},3);
        check(n.valid,"Primitive GRH Cayley normalized completed jet is finite");
        auto L=lfunction::evaluate(n.s,4);
        auto gamma=lfunction::detail::logGammaJet((n.s+1.0)*.5);
        garden::C A=std::exp(((n.s+1.0)*.5)*std::log(5.0/grh_cayley::pi)+gamma.logGamma);
        garden::C h1=.5*std::log(5.0/grh_cayley::pi)+.5*gamma.psi;
        garden::C h2=.25*gamma.trigamma;
        garden::C oneMinus=1.0-w*w;
        garden::C ds=1.0/(grh_cayley::log2*oneMinus);
        garden::C dds=2.0*w/(grh_cayley::log2*oneMinus*oneMinus);
        garden::C F=A*L.f;
        garden::C Fs=A*(L.first+h1*L.f);
        garden::C Fss=A*(L.second+2.0*h1*L.first+(h1*h1+h2)*L.f);
        garden::C Fw=Fs*ds;
        garden::C Fww=Fss*ds*ds+Fs*dds;
        auto newton=garden::correction(w,{0,0},{1,0},*grh,0);
        auto halley=garden::correction(w,{0,0},{1,0},*grh,1);
        garden::C expectedN=F/Fw;
        garden::C expectedH=2.0*F*Fw/(2.0*Fw*Fw-F*Fww);
        check(newton.valid&&std::abs(newton.delta-expectedN)<2e-10*(1+std::abs(expectedN)),
              "Primitive GRH Cayley Newton ratio equals direct completed-L jet");
        check(halley.valid&&std::abs(halley.delta-expectedH)<3e-10*(1+std::abs(expectedH)),
              "Primitive GRH Cayley Halley ratio equals direct completed-L jet");
    }
    auto pair=garden::data(8,1);
    check(pair->degree==0&&pair->records.size()==1,"Primitive character-pair data and upload record");
    for (garden::C s0:{garden::C(.5,6.18357819545086),garden::C(.5,-6.18357819545086)}) {
        auto e=primitive_pair::evaluate(s0,{0,0});
        check(e.valid&&e.residual<3e-10,"Primitive pair accepts chi/chibar numerical zero");
        auto o=garden::orbit(s0,{0,0},{1,0},*pair,0,8);
        check(o.status==0,"Primitive pair zero is a Newton attractor");
    }
    for (garden::C s0:{garden::C(.31,3.7),garden::C(.72,-5.1),garden::C(.48,10.2)}) {
        auto normalized=primitive_pair::evaluate(s0,{0,0});
        check(normalized.valid,"Primitive pair normalized completed jet is finite");
        auto L=lfunction::evaluate(s0,4);
        auto R=lfunction::evaluate(std::conj(s0),4);
        garden::C B0=std::conj(R.f),B1=std::conj(R.first),B2=std::conj(R.second);
        garden::C u0=L.f*B0;
        garden::C u1=L.first*B0+L.f*B1;
        garden::C u2=L.second*B0+2.0*L.first*B1+L.f*B2;
        auto gamma=lfunction::detail::logGammaJet((s0+1.0)*.5);
        garden::C A=std::exp(((s0+1.0)*.5)*std::log(5.0/primitive_pair::pi)+gamma.logGamma);
        garden::C h1=.5*std::log(5.0/primitive_pair::pi)+.5*gamma.psi;
        garden::C h2=.25*gamma.trigamma;
        garden::C F=A*A*u0;
        garden::C F1=A*A*(u1+2.0*h1*u0);
        garden::C F2=A*A*(u2+4.0*h1*u1+(4.0*h1*h1+2.0*h2)*u0);
        auto newton=garden::correction(s0,{0,0},{1,0},*pair,0);
        auto halley=garden::correction(s0,{0,0},{1,0},*pair,1);
        garden::C expectedN=F/F1;
        garden::C expectedH=2.0*F*F1/(2.0*F1*F1-F*F2);
        check(newton.valid&&std::abs(newton.delta-expectedN)<3e-10*(1+std::abs(expectedN)),
              "Primitive pair Newton ratio equals direct completed-product jet");
        check(halley.valid&&std::abs(halley.delta-expectedH)<4e-10*(1+std::abs(expectedH)),
              "Primitive pair Halley ratio equals direct completed-product jet");
    }
    {
        garden::C s(.43,7.25),target(.17,-.08);
        auto e=primitive_pair::evaluate(s,target);
        check(e.valid&&std::isfinite(e.residual),"Primitive pair normalized nonzero target remains finite");
    }
    // GRH phase-feedback garden: exact completed-L Cayley jets plus a bounded
    // vanishing phase kick. The deformation must preserve the zero attractors.
    auto feedbackData=garden::data(10,1);
    check(feedbackData->degree==0&&feedbackData->records.size()==1,"GRH phase-feedback data and upload record");
    {
        garden::C s0(.5,12.6749464170114);
        garden::C w=std::tanh(grh_cayley::log2*(s0-garden::C(.5,0)));
        auto step=grh_phase_feedback::correction(w,{.12,.8},{.9,0},0);
        check(step.valid&&step.residual<2e-11&&std::abs(step.delta)<1e-8,
              "Phase feedback vanishes at primitive completed-L zero");
        auto o=garden::orbit(w,{.12,.8},{.9,0},*feedbackData,0,8);
        check(o.status==0,"Phase-feedback garden preserves primitive zero attractor");
    }
    for (garden::C w:{garden::C(.23,.41),garden::C(-.31,.72),garden::C(.18,-.55)}) {
        auto e=grh_cayley::evaluate(w,{0,0},3);
        check(e.valid,"Phase-feedback sample has valid primitive completed-L Cayley jet");
        auto baseN=grh_phase_feedback::correction(w,{0,0},{.91,-.07},0);
        auto baseH=grh_phase_feedback::correction(w,{0,0},{.91,-.07},1);
        garden::C n=e.f/e.first;
        garden::C h=2.0*e.f*e.first/(2.0*e.first*e.first-e.f*e.second);
        check(baseN.valid&&std::abs(baseN.delta-garden::C(.91,-.07)*n)<3e-10*(1+std::abs(n)),
              "Zero-kick feedback Newton equals exact relaxed completed-L Newton");
        check(baseH.valid&&std::abs(baseH.delta-garden::C(.91,-.07)*h)<4e-10*(1+std::abs(h)),
              "Zero-kick feedback Halley equals exact relaxed completed-L Halley");
        auto phi=grh_phase_feedback::boundedPhase(e.f,.8);
        check(std::abs(phi)<=1.000000000001,"Bounded residual-phase kick has unit bound");
    }
    check(grh_phase_feedback::boundedPhase({0,0},1.2)==garden::C(0,0),
          "Residual-phase kick is exactly zero at f=0");

    // Literal completed-pair Cayley observable used by the A-tier GRH closure.
    for (auto xy : {std::pair<double,double>{4.7,1.3}, {5.6,2.4}, {6.3,3.2}, {7.0,4.0}}) {
        auto e=grh_observable::evaluate(xy.first,xy.second);
        check(e.valid,"Literal GRH Cayley observable is finite on interior sample");
        check(std::abs(e.zMinus-std::conj(e.zPlus))<2e-13*(1+e.a),
              "Cayley point minus is conjugate of point plus");
        check(e.symmetryError<2e-12,"Completed-pair logarithmic derivative respects conjugation symmetry");
        check(std::abs(e.observable.imag())<2e-11*(1+std::abs(e.observable.real())),
              "Literal two-point observable is real on the conjugate Cayley pair");
        double h=1e-6, rp=e.r+h, rm=e.r-h;
        auto logit=[](double r){return std::log(r/(1-r));};
        auto ep=grh_observable::evaluate(xy.first,logit(rp));
        auto em=grh_observable::evaluate(xy.first,logit(rm));
        double derivative=((ep.plus.logAbsTheta+ep.minus.logAbsTheta)-
                           (em.plus.logAbsTheta+em.minus.logAbsTheta))/(2*h);
        check(ep.valid&&em.valid&&std::abs(e.observable.real()+derivative)<2e-6*(1+std::abs(e.observable.real())),
              "Cayley observable equals minus r-derivative of symmetric completed-pair log product");
    }
    // Independent 50-decimal mpmath Dirichlet/Hurwitz references at s=.5+15i.
    struct Reference { int generation; garden::C f, first, second; };
    const Reference refs[] = {
        {1, {0.14710990704334914,0.70475224164321193}, {0.76952454010825733,-0.48055106252743174}, {-0.71995150809889807,0.29375560727409655}},
        {2, {0.28756914543790102,-1.4636518873819675}, {1.7501205413265943,1.8401380703884196}, {-3.5642343448113776,-2.7590858818961475}},
        {3, {2.2711259755793645,-0.94104385837824345}, {-2.2611759177772139,1.7880057689313071}, {4.4357636728384984,-3.6695169377589565}},
        {4, {0.078355004395224512,0.35236327703871967}, {1.9060819483422244,-0.88228514055329377}, {-4.7859159518838581,1.9735183228107167}}
    };
    for (const auto &r : refs) {
        auto e=lfunction::evaluate({.5,15},r.generation);
        check(e.valid && std::abs(e.f-r.f)<1e-10*(1+std::abs(r.f)), "L value versus independent high precision reference");
        check(e.valid && std::abs(e.first-r.first)<1e-10*(1+std::abs(r.first)), "L first derivative versus independent high precision reference");
        check(e.valid && std::abs(e.second-r.second)<1e-10*(1+std::abs(r.second)), "L second derivative versus independent high precision reference");
    }
    out<<"RH_GRH_GARDEN_CHECKS="<<checks<<" PASS\n";
    return checks;
}
} // namespace rh_grh_tests
