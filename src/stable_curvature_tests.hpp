#pragma once
#include "rh_gpu_data.hpp"
#include <complex>
#include <ostream>

// Native algebra checks. These execute only when compiled by the user's
// C++ compiler. Hardware shader/repair/quality checks are a separate suite.
namespace stable_curvature_tests {
using LC=std::complex<long double>;
using Poly=std::vector<long double>;

inline LC evaluate(const Poly& coefficients,LC point,bool reversed=false) {
    LC value{};
    if(reversed)for(long double a:coefficients)value=value*point+a;
    else for(auto a=coefficients.rbegin();a!=coefficients.rend();++a)value=value*point+*a;
    return value;
}
inline Poly derivative(const Poly& coefficients) {
    if(coefficients.size()<2)return Poly{0.L};
    Poly result(coefficients.size()-1);
    for(std::size_t k=1;k<coefficients.size();++k)result[k-1]=k*coefficients[k];
    return result;
}
inline bool finite(LC value) {return std::isfinite(value.real())&&std::isfinite(value.imag());}
inline bool close(LC left,LC right,long double tolerance=2.e-10L) {
    return finite(left)&&finite(right)&&std::abs(left-right)<=tolerance*(1.L+std::max(std::abs(left),std::abs(right)));
}
inline long double majorant(const Poly& coefficients,LC point,bool reversed=false) {
    long double bound=0.L;
    const long double radius=std::abs(point);
    if(reversed)for(long double a:coefficients)bound=bound*radius+std::abs(a);
    else for(auto a=coefficients.rbegin();a!=coefficients.rend();++a)bound=bound*radius+std::abs(*a);
    return bound;
}
// Compare an algebraic identity against its absolute term majorant. This
// checks backward error even when the true result is severely cancelled.
// Use a conservative operation-count budget for coefficient construction,
// complex Horner evaluation and the independent derivative products below.
// This is a numerical regression bound, not an interval certificate.
inline bool identityClose(LC left,LC right,long double bound,int degree) {
    const long double roundoff=64.L*(3*degree+1)*std::numeric_limits<long double>::epsilon();
    const long double gamma=roundoff/(1.L-roundoff);
    return finite(left)&&finite(right)&&std::isfinite(bound)&&bound>=0.L
        &&std::abs(left-right)<=gamma*bound;
}
inline LC payloadCorrection(const rhgpu::StableCurvatureData& data,LC y,bool reverse) {
    const LC u=y/data.variableScale;
    const LC argument=reverse?LC(1.L,0.L)/u:u;
    LC result=data.quotientScale*evaluate(data.h,argument,reverse)*evaluate(data.theta,argument,reverse)
        /evaluate(data.j,argument,reverse);
    return reverse?u*result:result;
}
inline bool denominatorGuard(LC value,long double bound) {
    return !finite(value)||!std::isfinite(bound)||value==LC{}||
        std::abs(value)<3.e-5L*bound||std::abs(value)<1.e-30L;
}
struct GuardedEvaluation {
    bool finiteValues=false,rejected=true;
    long double thetaRatio=0.L,jRatio=0.L;
};
inline GuardedEvaluation shaderDomain(const rhgpu::StableCurvatureData& data,LC y,bool reverse) {
    const LC u=y/data.variableScale;
    if(reverse&&u==LC{})return {};
    const LC argument=reverse?LC(1.L,0.L)/u:u;
    const LC theta=evaluate(data.theta,argument,reverse),j=evaluate(data.j,argument,reverse);
    const long double thetaBound=majorant(data.theta,argument,reverse);
    const long double jBound=majorant(data.j,argument,reverse);
    GuardedEvaluation result;
    result.finiteValues=finite(theta)&&finite(j)&&std::isfinite(thetaBound)&&std::isfinite(jBound)
        &&thetaBound>0.L&&jBound>0.L;
    if(!result.finiteValues)return result;
    result.thetaRatio=std::abs(theta)/thetaBound;
    result.jRatio=std::abs(j)/jBound;
    // Match the production shader's theta and denominator guards. A small
    // H numerator remains valid: its regular roots are intended fixed points.
    result.rejected=result.thetaRatio<3.e-5L||denominatorGuard(j,jBound);
    return result;
}

inline bool run(std::ostream& out) {
    int cases=0,failures=0,acceptedPoints=0,guardedPoints=0;
    LC currentPoint{};
    auto check=[&](bool good,const char* name,int degree) {
        ++cases;
        if(!good){++failures;out<<"STABLE_CURVATURE_TEST_FAIL="<<name<<" DEGREE="<<degree
            <<" Y_RE="<<currentPoint.real()<<" Y_IM="<<currentPoint.imag()<<'\n';}
    };
    const std::array<LC,9> points{{LC{},LC(.31L,.27L),LC(-.4L,.6L),LC(2.7L,1.3L),
        LC(3.99L,.1L),LC(4.01L,.1L),LC(8.L,.7L),LC(0.L,4.2L),LC(-1.L,.5L)}};
    for(int degree:{8,16,32}) {
        currentPoint=LC{};
        const auto data=rhgpu::makeStableCurvatureData(degree);
        check(data.valid,"PAYLOAD_READY",degree);
        if(!data.valid)continue;
        check(data.theta.size()==std::size_t(degree+1)&&data.h.size()==std::size_t(2*degree-1)&&
            data.j.size()==std::size_t(3*degree-2),"POLYNOMIAL_DEGREES",degree);
        bool normal=true;
        for(const auto* polynomial:{&data.theta,&data.h,&data.j})for(long double a:*polynomial) {
            const float packed=static_cast<float>(a);
            normal=normal&&std::isfinite(packed)&&(a==0.L||std::abs(packed)>=std::numeric_limits<float>::min());
        }
        check(normal,"NO_DROPPED_FLOAT_COEFFICIENTS",degree);
        Poly p=data.theta;
        for(auto& a:p)a*=data.thetaScale;
        const auto p1=derivative(p),p2=derivative(p1),p3=derivative(p2);
        auto settings=rhfresh::defaultSettings(0);settings.degree=degree;settings.scale=1.;settings.shift=0.;
        const rhfresh::Model trusted(settings);
        int safeInDegree=0;
        for(const LC y:points) {
            currentPoint=y;
            const LC u=y/data.variableScale;
            const LC t=evaluate(p,u),t1=evaluate(p1,u),t2=evaluate(p2,u),t3=evaluate(p3,u);
            const LC h=t*t1+u*(t*t2-t1*t1);
            const LC hp=2.L*t*t2+u*(t*t3-t1*t2);
            const LC j=hp*t-2.L*h*t1;
            const long double tBound=majorant(p,u),t1Bound=majorant(p1,u);
            const long double t2Bound=majorant(p2,u),t3Bound=majorant(p3,u),radius=std::abs(u);
            const long double hBound=tBound*t1Bound+radius*(tBound*t2Bound+t1Bound*t1Bound);
            const long double hpBound=2.L*tBound*t2Bound+radius*(tBound*t3Bound+t1Bound*t2Bound);
            const long double jBound=hpBound*tBound+2.L*hBound*t1Bound;
            const LC expandedH=evaluate(data.h,u)*data.hScale;
            const LC expandedJ=evaluate(data.j,u)*data.jScale;
            // All nine points retain independent derivative-identity checks,
            // including the cancelled points which require CPU repair.
            check(identityClose(expandedH,h,hBound+majorant(data.h,u)*data.hScale,degree),
                "PAIR_H_BACKWARD_IDENTITY",degree);
            check(identityClose(expandedJ,j,jBound+majorant(data.j,u)*data.jScale,degree),
                "J_BACKWARD_IDENTITY",degree);
            const bool reverse=radius>1.L;
            const auto selected=shaderDomain(data,y,reverse);
            check(selected.finiteValues,"POINT_AND_MAJORANTS_FINITE",degree);
            if(selected.finiteValues&&!selected.rejected) {
                ++acceptedPoints;++safeInDegree;
                // Keep the original strict relative tolerances inside the
                // domain accepted by the renderer; do not loosen them.
                check(close(expandedH,h),"PAIR_H_VS_DERIVATIVES",degree);
                check(close(expandedJ,j),"J_VS_DERIVATIVES",degree);
                const LC direct=payloadCorrection(data,y,reverse);
                const rhfresh::C input(static_cast<double>(y.real()),static_cast<double>(y.imag()));
                const auto cpu=trusted.correction(input);
                check(close(direct,LC(cpu.real(),cpu.imag()),2.e-8L),"CORRECTION_VS_TRUSTED_CPU",degree);
                if(y!=LC{}&&!shaderDomain(data,y,!reverse).rejected)
                    check(close(direct,payloadCorrection(data,y,!reverse)),"ORDINARY_VS_REVERSED",degree);
            } else {
                ++guardedPoints;
                // Nonfinite/malformed payloads cannot masquerade as expected
                // precision repairs: these fixtures must show cancellation.
                check(selected.finiteValues&&(selected.thetaRatio<3.e-5L||selected.jRatio<3.e-5L),
                    "CANCELLED_POINT_REQUIRES_REPAIR",degree);
                out<<"STABLE_CURVATURE_GUARDED_POINT DEGREE="<<degree
                    <<" Y_RE="<<y.real()<<" Y_IM="<<y.imag()
                    <<" THETA_RATIO="<<selected.thetaRatio<<" J_RATIO="<<selected.jRatio<<'\n';
            }
            if(degree==32&&y==LC(0.L,4.2L))
                check(selected.finiteValues&&selected.rejected&&selected.jRatio<3.e-5L,
                    "DEGREE32_IMAGINARY_CANCELLATION_GUARDED",degree);
        }
        check(safeInDegree>=6,"SAFE_POINT_COVERAGE",degree);
        // The actual affine model uses y=b*z+h and divides the correction by b.
        settings.scale=.9;settings.shift=.28;
        const rhfresh::Model affine(settings);
        const rhfresh::C z(.31,.27);
        const LC y=static_cast<long double>(settings.scale)*LC(z.real(),z.imag())+
            static_cast<long double>(settings.shift);
        currentPoint=y;
        const auto cpu=affine.correction(z);
        check(close(payloadCorrection(data,y,std::abs(y/data.variableScale)>1.L)/static_cast<long double>(settings.scale),
            LC(cpu.real(),cpu.imag()),2.e-8L),"AFFINE_SCALE_AND_SHIFT",degree);
        auto state=amp::defaults(20,1);state.spectralOrder=degree==8?1:degree==16?2:3;
        currentPoint=LC{};
        const auto gpu=rhgpu::makeData(state);
        const auto metadata=gpu.texels[std::size_t(rhgpu::CurvatureMetaRow)*rhgpu::CoefficientWidth+rhgpu::CurvatureMetaColumn];
        check(gpu.valid&&metadata[0]==1.f&&metadata[1]==4.f&&metadata[2]>0.f,
            "PACKED_READY_AND_SCALE",degree);
    }
    check(!rhgpu::makeStableCurvatureData(0).valid&&!rhgpu::makeStableCurvatureData(33).valid,
        "INVALID_DEGREE_DECLINED",0);
    // Synthetic T=1+y: H=1,J=-2. Exercise the same factored evaluator at
    // its simple theta root and far outside the ordinary-Horner domain.
    rhgpu::StableCurvatureData linear;
    linear.variableScale=1.L;linear.theta={1.L,1.L};linear.h={1.L};linear.j={-2.L};
    check(close(payloadCorrection(linear,LC(-1.L),false),LC{}),"SIMPLE_THETA_ROOT_CONTINUATION",1);
    const LC distantPoint(4.L,2.L);
    check(close(payloadCorrection(linear,distantPoint,true),-(1.L+distantPoint)/2.L),"REVERSED_LINEAR_FINITE_VALUE",1);
    // Synthetic T=1+y^2: H=4y,J=4-12y^2. A Newton pole must be repaired,
    // while the exact numerator root y=0 and simple theta roots remain finite.
    rhgpu::StableCurvatureData quadratic;
    quadratic.variableScale=1.L;quadratic.theta={1.L,0.L,1.L};
    quadratic.h={0.L,4.L};quadratic.j={4.L,0.L,-12.L};
    const long double pole=1.L/std::sqrt(3.L);
    check(denominatorGuard(evaluate(quadratic.j,LC(pole)),4.L+12.L*pole*pole),"NEWTON_POLE_GUARDED",2);
    check(close(payloadCorrection(quadratic,LC{},false),LC{}),"REGULAR_NUMERATOR_ROOT",2);
    check(!denominatorGuard(LC(4.L),4.L),"REGULAR_DENOMINATOR_ACCEPTED",2);
    check(denominatorGuard(LC{},0.L)&&denominatorGuard(LC(1.e-35L),1.e-35L),"EXACT_AND_TINY_POLES_GUARDED",2);
    const LC imaginary(0.L,1.L);
    check(close(payloadCorrection(quadratic,imaginary,false),LC{}),
        "QUADRATIC_THETA_ROOT_CONTINUATION",2);
    const Poly multipleTheta{1.L,2.L,1.L},multipleJ{-4.L,-12.L,-12.L,-4.L};
    check(evaluate(multipleTheta,LC(-1.L))==LC{}&&
        denominatorGuard(evaluate(multipleJ,LC(-1.L)),32.L),"MULTIPLE_THETA_ROOT_DELEGATED",2);
    check(guardedPoints>0,"CANCELLED_POINT_COVERAGE",32);
    out<<"STABLE_CURVATURE_ALGEBRA="<<(failures?"FAIL":"PASS")<<" CASES="<<cases
        <<" ACCEPTED_POINTS="<<acceptedPoints<<" GUARDED_POINTS="<<guardedPoints<<'\n';
    return failures==0;
}
} // namespace stable_curvature_tests
