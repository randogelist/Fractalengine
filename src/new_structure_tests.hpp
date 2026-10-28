#pragma once
#include "rh_candidates.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

// Independent direct-power-polynomial oracle for the 1.8.0 generation.
// Does not call the shared-normalization evaluator being tested.
namespace new_structure_tests {
using Z=std::complex<long double>;
inline bool finite(Z z){return std::isfinite(z.real())&&std::isfinite(z.imag());}

inline std::array<Z,5> directJets(int degree,Z y){
    std::array<Z,5> jet{};
    for(int order=0;order<=4;++order){
        Z value{};
        for(int j=degree;j>=order;--j){
            long double a=static_cast<long double>(ff_rh::THETA_COEFFICIENTS[j]);
            for(int k=0;k<order;++k)a*=j-k;
            value=value*y+a;
        }
        jet[order]=value;
    }
    return jet;
}

inline Z target(int candidate,int degree,Z y,long double R){
    const auto t=directJets(degree,y);
    if(candidate==13){
        const auto a=directJets(degree,y+R),b=directJets(degree,y-R);
        return a[0]*b[0]-t[0]*t[0];
    }
    if(candidate==14)return t[0]*t[2]-t[1]*t[1];
    if(candidate==15){
        const auto a=directJets(degree,2.L*y),b=directJets(degree,.5L*y);
        return a[0]*b[0]-t[0]*t[0];
    }
    if(candidate==16){
        const auto b=directJets(degree,-y);
        return t[1]*b[0]-t[0]*b[1];
    }
    if(candidate==17)return 2.L*t[1]*t[3]-3.L*t[2]*t[2];
    return Z(std::numeric_limits<long double>::quiet_NaN(),0);
}
inline Z targetDerivative(int candidate,int degree,Z y,long double R){
    const auto t=directJets(degree,y);
    if(candidate==13){
        const auto a=directJets(degree,y+R),b=directJets(degree,y-R);
        return a[1]*b[0]+a[0]*b[1]-2.L*t[0]*t[1];
    }
    if(candidate==14)return t[0]*t[3]-t[1]*t[2];
    if(candidate==15){
        const auto a=directJets(degree,2.L*y),b=directJets(degree,.5L*y);
        return 2.L*a[1]*b[0]+.5L*a[0]*b[1]-2.L*t[0]*t[1];
    }
    if(candidate==16){
        const auto b=directJets(degree,-y);
        return t[2]*b[0]-2.L*t[1]*b[1]+t[0]*b[2];
    }
    if(candidate==17)return 2.L*t[1]*t[4]-4.L*t[2]*t[3];
    return Z(std::numeric_limits<long double>::quiet_NaN(),0);
}
inline bool run(std::ostream& out){
    int cases=0,failures=0;
    long double worst=0.L;
    auto check=[&](bool ok,const char* label,int index,int degree,Z y,long double error){
        ++cases;worst=std::max(worst,error);
        if(!ok){++failures;out<<"NEW_STRUCTURE_FAIL="<<label<<" CANDIDATE="<<index+1
            <<" N="<<degree<<" X="<<y.real()<<" Y="<<y.imag()
            <<" REL_ERROR="<<error<<'\n';}
    };
    const std::array<Z,6> points{{{.13L,.09L},{.31L,.27L},{-.42L,.13L},
                                 {.62L,-.2L},{1.7L,.4L},{3.2L,.2L}}};
    for(int candidate=13;candidate<=17;++candidate)
        for(int degree:{8,16,32})
            for(double R:{.12,.6,1.2}) {
                auto setting=rhfresh::defaultSettings(candidate);
                setting.degree=degree;setting.radius=R;setting.scale=.875;setting.shift=.125;
                const rhfresh::Model model(setting);
                for(Z z:points){
                    const Z y=static_cast<long double>(setting.scale)*z
                        +static_cast<long double>(setting.shift);
                    const Z f=target(candidate,degree,y,R);
                    const Z df=targetDerivative(candidate,degree,y,R);
                    const Z expected=f/df/static_cast<long double>(setting.scale);
                    const auto actual=model.correction({double(z.real()),double(z.imag())});
                    const Z observed(actual.real(),actual.imag());
                    const long double relative=std::abs(expected-observed)/(1.L+std::abs(expected));
                    check(model.valid()&&finite(expected)&&finite(observed)&&relative<3.e-9L,
                          "DIRECT_POWER_ORACLE",candidate,degree,z,relative);
                    // Independent symmetric difference verifies the derivative
                    // factor and the chain-rule multipliers (2 and 1/2).
                    constexpr long double h=2.e-6L;
                    const Z numerical=(target(candidate,degree,y+h,R)-target(candidate,degree,y-h,R))/(2.L*h);
                    const long double derivativeError=std::abs(numerical-df)/(1.L+std::abs(df));
                    check(finite(numerical)&&derivativeError<2.e-6L,
                          "COMPLEX_ANALYTIC_DERIVATIVE",candidate,degree,z,derivativeError);
                    const auto reflected=model.correction({double(z.real()),-double(z.imag())});
                    const Z conjugate(reflected.real(),reflected.imag());
                    const long double conjugateError=std::abs(conjugate-std::conj(observed))
                        /(1.L+std::abs(observed));
                    check(finite(conjugate)&&conjugateError<3.e-10L,
                          "CONJUGATE_SYMMETRY",candidate,degree,z,conjugateError);
                }
            }
    // Small-R finite shift converges quadratically to the Laguerre defect.
    for(int degree:{8,16,32}){
        const Z y(.34L,.19L);
        const long double R=.003L;
        const Z laguerre=target(14,degree,y,.6L);
        const Z echo=target(13,degree,y,R)/(R*R);
        const long double error=std::abs(echo-laguerre)/(1.L+std::abs(laguerre));
        check(error<5.e-4L,"ECHO_TO_LAGUERRE_LIMIT",13,degree,y,error);
    }
    // Antipodal current is odd, its derivative is even.
    for(int degree:{8,16,32}){
        const Z y(.24L,.4L);
        const Z f=target(16,degree,y,.6L),minus=target(16,degree,-y,.6L);
        const Z df=targetDerivative(16,degree,y,.6L),dMinus=targetDerivative(16,degree,-y,.6L);
        const long double error=(std::abs(f+minus)+std::abs(df-dMinus))
           /(1.L+std::abs(f)+std::abs(df));
        check(error<5.e-14L,"ANTIPODAL_ODD_CURRENT",16,degree,y,error);
    }
    // Projective invariance of S(T) under (aT+b)/(cT+d).
    for(int degree:{8,16,32}){
        const Z y(.31L,.27L);const auto t=directJets(degree,y);
        const Z a(1.1L,.2L),b(.2L,-.3L),c(.17L,.13L),d(.9L,-.1L);
        const Z D=a*d-b*c,p=c*t[0]+d;
        const Z g1=D*t[1]/(p*p);
        const Z g2=D*(t[2]/(p*p)-2.L*c*t[1]*t[1]/(p*p*p));
        const Z g3=D*(t[3]/(p*p)-6.L*c*t[1]*t[2]/(p*p*p)
            +6.L*c*c*t[1]*t[1]*t[1]/(p*p*p*p));
        const Z original=t[3]/t[1]-1.5L*std::pow(t[2]/t[1],2);
        const Z transformed=g3/g1-1.5L*std::pow(g2/g1,2);
        const long double error=std::abs(original-transformed)/(1.L+std::abs(original));
        check(error<2.e-11L,"SCHWARZIAN_MOBIUS_INVARIANT",17,degree,y,error);
    }
    out<<"NEW_STRUCTURE_TEST_CASES="<<cases<<" FAILURES="<<failures<<" WORST_REL_ERROR="<<worst<<'\n';
    out<<"NEW_STRUCTURE_TEST="<<(failures==0?"PASS":"FAIL")<<'\n';
    return failures==0;
}
} // namespace new_structure_tests
