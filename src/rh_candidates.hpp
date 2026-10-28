#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include "theta_coefficients.hpp"

// The coefficient table defines T(y)=Theta_N(200y)/Theta(0).  These are
// finite, rounded Taylor models, not RH certificates.  Every correction is
// the Newton quotient of the stated F, including its quotient derivative.
// The displayed iteration is z^2+c-lambda*F(y)/(b*F_y(y)), y=b*z+h.
// No square-root sampling, conjugation, denominator epsilon, or clipping
// changes that rational map.  A radius-16 stop is a display cutoff only.
namespace rhfresh {
using C = std::complex<double>;
constexpr int OriginalCandidateCount = 10;
constexpr int CandidateCount = 18;
constexpr int NewGenerationStart = 13;

struct Entry {
    const char* id;
    const char* name;
    const char* description;
    const char* formula;
};

inline const std::array<Entry, CandidateCount>& registry() {
    static const std::array<Entry, CandidateCount> entries{{
        {"curvature", "Normalized boundary curvature", "Derivative of the logarithmic radial count.",
         "g=y*T'/T; F=q=g'; correction=H*T/(H'*T-2*H*T')"},
        {"equipotential", "Positive equipotentials", "Angular levels of the exact polynomial lift, at fixed R.",
         "F=E_R'(y)-level; E_R(x)=T(Rw)T(R/w), 2x=w+1/w"},
        {"logslope", "Logarithmic equipotentials", "Angular slope divided by the potential itself.",
         "F=E_R'(y)/E_R(y)-level; correction=(E_R'-level*E_R)*E_R/(E_R''*E_R-E_R'^2)"},
        {"integer", "Integer count levels", "Unit levels of the logarithmic radial count.",
         "F=y*T'(y)/T(y)-level; correction=(y*T'-level*T)*T/H"},
        {"angular_fold", "Angular fold", "Stationary points of the angular positivity profile.",
         "F=E_R''(y); correction=E_R''/E_R'''"},
        {"boundary_source", "Critical-boundary sources", "Boundary angular sources continued in complex radius.",
         "F=B(y)-level; B(y)=2*y*H(-y); H=T*T'+y*T*T''-y*T'^2"},
        {"curvature_jet", "Curvature jet", "One additional derivative of radial logarithmic curvature.",
         "F=q'=g''; J=H'*T-2*H*T'; correction=J*T/(J'*T-3*J*T')"},
        {"log_curvature", "Dilation curvature", "Derivative of y times radial curvature.",
         "F=(y*q)'; K=H*T+y*(H'*T-2*H*T'); correction=K*T/(K'*T-3*K*T')"},
        {"radial_fold", "Radial fold", "Radial stationary points of the angular profile, at fixed R.",
         "F=partial_R E_R'(y); correction=(partial_R E_R')/(partial_R E_R'')"},
        {"equipotential_inflection", "Equipotential inflection", "Zero signed-curvature numerator of U(R,x)=constant, U=E_R'.",
         "F=U_RR*U_x^2-2*U_Rx*U_R*U_x+U_xx*U_R^2; correction=F/F_x"},
        {"log_angular_fold", "Logarithmic angular fold", "Stationary logarithmic angular slope of the finite Theta lift.",
         "D=E*E''-E'^2; F=D/E^2-level; correction=(D-level*E^2)*E/(D'*E-2*D*E')"},
        {"dilation_gap", "Dyadic dilation gap", "Dyadic change of the finite-model logarithmic radial count.",
         "C(y)=y*T'/T; F=C(2y)-C(y)-level; correction=F/(2*C'(2y)-C'(y))"},
        {"equipotential_flow", "Equipotential flow tilt", "Radial-to-angular slope ratio of U=E_R', continued to complex coordinates.",
         "F=U_R/U_x-level; correction=(U_R-level*U_x)*U_x/(U_Rx*U_x-U_R*U_xx)"},
        {"spectral_echo", "Spectral echo determinant", "Symmetric finite-shift correlation defect of the Theta polynomial.",
         "F=T(y+R)*T(y-R)-T(y)^2; F'=T'(y+R)*T(y-R)+T(y+R)*T'(y-R)-2*T(y)*T'(y)"},
        {"laguerre_gram", "Laguerre differential determinant", "Infinitesimal spectral echo: logarithmic second-derivative determinant.",
         "F=T*T''-T'^2; F'=T*T'''-T'*T''"},
        {"dyadic_holonomy", "Dyadic scale holonomy", "Multiplicative failure of homogeneity across three dyadic radii.",
         "F=T(2y)*T(y/2)-T(y)^2; F'=2*T'(2y)*T(y/2)+(1/2)*T(2y)*T'(y/2)-2*T(y)*T'(y)"},
        {"antipodal_wronskian", "Antipodal Wronskian current", "Odd reflection current between positive and negative spectral sheets.",
         "F=T'(y)*T(-y)-T(y)*T'(-y); F'=T''(y)*T(-y)-2*T'(y)*T'(-y)+T(y)*T''(-y)"},
        {"projective_schwarzian", "Projective Schwarzian curvature", "Möbius-invariant projective derivative curvature, numerator-zero basins.",
         "F=2*T'*T'''-3*T''^2; F'=2*T'*T''''-4*T''*T'''"}
    }};
    return entries;
}

struct Settings {
    int candidate = 0;
    int degree = 16;
    double radius = 1.;
    double scale = .9;
    double shift = .28;
    double level = 1.;
    C lambda{.45, 0.};
    C parameter{-.56, .23};
    int iterations = 160;
};

inline Settings defaultSettings(int candidate) {
    Settings s;
    s.candidate = candidate;
    s.radius = .6; // Physical radial coordinate is r=200R.
    s.level = 0.;
    // Literal anchor levels in the normalized degree-16 model. They remain
    // fixed when degree or radius is edited; there is no hidden rescaling.
    if (candidate == 1) s.level = 8.3393668337974311;  // E_0.6'(0)
    if (candidate == 2) s.level = 5.1651956222738728;  // E_0.6'(0)/E_0.6(0)
    if (candidate == 3) s.level = 1.;                // Integer count level
    if (candidate == 5) s.level = .21520193813860753; // B(0.28)
    if (candidate == 12) s.level = 1.;              // Equipotential flow level
    return s;
}

enum class Status { Cutoff, Unresolved, Numerical, Invalid };
struct Result {
    C endpoint{};
    int iterations = 0;
    double smooth = 0., structure = 0., phase = 0.;
    Status status = Status::Unresolved;
};

namespace detail {
using Poly = std::vector<long double>; // Increasing powers or Chebyshev indices.
using LC = std::complex<long double>;

inline void trim(Poly& a) {
    while (a.size() > 1 && a.back() == 0.L) a.pop_back();
    if (a.empty()) a.push_back(0.L);
}
inline bool finite(C z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
inline bool finite(LC z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
inline C nonfinite() {
    const double q = std::numeric_limits<double>::quiet_NaN();
    return C(q,q);
}
inline LC badLong() {
    const long double q = std::numeric_limits<long double>::quiet_NaN();
    return LC(q,q);
}
inline LC quotient(LC numerator, LC denominator) {
    if (!finite(numerator) || !finite(denominator)) return badLong();
    // At a multiple zero or a pole of any stated Newton target, F/F'
    // extends to zero. All common zeros in these factored formulas have
    // that form; a nonzero numerator over exact zero is a numerical pole.
    if (denominator == LC(0.L,0.L))
        return numerator == LC(0.L,0.L) ? LC(0.L,0.L) : badLong();
    return numerator/denominator;
}

// E_R(x)=sum a_n^2 R^(2n)+sum_{k>=1} amplitude_k*T_k(x).
// Ordinary radial derivatives use (p)_j R^(p-j), not Euler derivatives.
inline Poly angularAmplitudes(int degree, long double radius, int radialOrder) {
    std::array<long double,65> powers{};
    powers[0] = 1.L;
    for (int j = 1; j <= 2*degree; ++j) powers[j] = powers[j-1]*radius;
    auto weight = [&](int exponent) {
        if (exponent < radialOrder) return 0.L;
        long double v = powers[exponent-radialOrder];
        for (int j = 0; j < radialOrder; ++j) v *= exponent-j;
        return v;
    };
    Poly amplitude(degree+1,0.L);
    for (int n = 0; n <= degree; ++n) {
        const long double a = ff_rh::THETA_COEFFICIENTS[n];
        amplitude[0] += a*a*weight(2*n);
    }
    for (int k = 1; k <= degree; ++k)
        for (int n = 0; n+k <= degree; ++n)
            amplitude[k] += 2.L*static_cast<long double>(ff_rh::THETA_COEFFICIENTS[n])
                *static_cast<long double>(ff_rh::THETA_COEFFICIENTS[n+k])*weight(2*n+k);
    return amplitude;
}
inline Poly chebyshevDerivative(const Poly& a) {
    const int degree = static_cast<int>(a.size())-1;
    if (degree < 1) return Poly{0.L};
    Poly d(degree,0.L);
    d[degree-1] = 2.L*degree*a[degree];
    if (degree >= 2) d[degree-2] = 2.L*(degree-1)*a[degree-1];
    for (int k = degree-3; k >= 0; --k) d[k] = d[k+2]+2.L*(k+1)*a[k+1];
    d[0] *= .5L;
    trim(d);
    return d;
}
inline LC clenshaw(const Poly& coefficients, LC y) {
    LC next(0.L,0.L), after(0.L,0.L);
    for (int k = static_cast<int>(coefficients.size())-1; k >= 1; --k) {
        const LC current = 2.L*y*next-after+coefficients[k];
        after = next;
        next = current;
    }
    return y*next-after+coefficients[0];
}

// For |y|<=1, return derivatives through order four by long-double Horner.
// Outside that disk, return jets[k]=T^(k)(y)/y^(N-k), each computed by
// reversed Horner. The caller balances the powers before division, so no
// y^N or expanded products need be formed. Rounded coefficients are the
// exact binary64 table entries promoted to long double.
inline std::array<LC,5> polynomialJets(const Poly& coefficients, LC y, bool reversed,
                                     int maximumOrder=4) {
    const int degree = static_cast<int>(coefficients.size())-1;
    const LC argument = reversed ? LC(1.L,0.L)/y : y;
    std::array<LC,5> jet{};
    for (int order = 0; order <= maximumOrder && order <= 4 && order <= degree; ++order) {
        LC value(0.L,0.L);
        if (reversed) {
            for (int j = order; j <= degree; ++j) {
                long double a = coefficients[j];
                for (int k = 0; k < order; ++k) a *= j-k;
                value = value*argument+a;
            }
        } else {
            for (int j = degree; j >= order; --j) {
                long double a = coefficients[j];
                for (int k = 0; k < order; ++k) a *= j-k;
                value = value*argument+a;
            }
        }
        jet[order] = value;
    }
    return jet;
}

inline int requiredThetaOrder(int candidate) {
    if(candidate==13||candidate==15)return 1;
    if(candidate==16)return 2;
    if(candidate==14)return 3;
    if(candidate==17)return 4;
    return candidate==3||candidate==11?2:(candidate==0||candidate==5?3:4);
}

// Only the degree/radius/candidate affect this immutable payload. Model
// instances for palettes, exposure, gain, level or chart edits reuse it.
// Derivative coefficients keep the exact multiplication order of the
// original evaluator; Horner no longer forms factorial weights per orbit.
struct CoefficientData {
    int candidate=0,degree=0;
    double radius=0.;
    Poly theta;
    std::array<Poly,5> thetaDerivatives{};
    // Candidate 11 also needs derivatives of B(u)=u^N*T(1/u).
    // Its count offset C-N stays resolvable far beyond |y|=1.
    std::array<Poly,3> inverseThetaDerivatives{};
    std::array<std::array<Poly,5>,3> angular{};
    bool valid=false;

    CoefficientData(int selected,int count,double r):candidate(selected),degree(count),radius(r) {
        if(selected<0||selected>=CandidateCount||count<4||count>32||!std::isfinite(r)||r<=0.)return;
        theta.resize(count+1);
        for(int j=0;j<=count;++j)theta[j]=ff_rh::THETA_COEFFICIENTS[j];
        const bool usesTheta=selected==0||selected==3||selected==5||selected==6||selected==7||selected==11||selected>=13;
        if(usesTheta){
            for(int order=0;order<=requiredThetaOrder(selected);++order){
                auto& values=thetaDerivatives[order];
                values.resize(count-order+1);
                for(int j=order;j<=count;++j){
                    long double a=theta[j];
                    for(int k=0;k<order;++k)a*=j-k;
                    values[j-order]=a;
                }
                if(selected==11){
                    auto& inverse=inverseThetaDerivatives[order];
                    inverse.resize(count-order+1);
                    for(int k=order;k<=count;++k){
                        long double a=theta[count-k];
                        for(int j=0;j<order;++j)a*=k-j;
                        inverse[k-order]=a;
                    }
                }
            }
        }else{
            const std::array<int,3> orders=selected==9?std::array<int,3>{{4,3,2}}
                :selected==12?std::array<int,3>{{3,2,-1}}
                :selected==8?std::array<int,3>{{-1,2,-1}}
                :std::array<int,3>{{selected==4||selected==10?3:2,-1,-1}};
            for(int radial=0;radial<=2;++radial){
                if(orders[radial]<0)continue;
                angular[radial][0]=angularAmplitudes(count,r,radial);
                for(int k=1;k<=orders[radial];++k)
                    angular[radial][k]=chebyshevDerivative(angular[radial][k-1]);
                for(int k=0;k<=orders[radial];++k)
                    for(long double a:angular[radial][k])if(!std::isfinite(a))return;
            }
        }
        valid=true;
    }
};

inline std::shared_ptr<const CoefficientData> coefficientData(int candidate,int degree,double radius) {
    if(candidate<0||candidate>=CandidateCount||degree<4||degree>32||!std::isfinite(radius)||radius<=0.)
        return {};
    static std::array<std::shared_ptr<const CoefficientData>,CandidateCount> cache{};
    auto& slot=cache[static_cast<std::size_t>(candidate)];
    auto data=std::atomic_load_explicit(&slot,std::memory_order_acquire);
    if(data&&data->degree==degree&&data->radius==radius)return data;
    data=std::make_shared<const CoefficientData>(candidate,degree,radius);
    // Concurrent callers may build the same small payload once; each sees
    // a complete immutable snapshot, and there is no mutex in evaluation.
    std::atomic_store_explicit(&slot,data,std::memory_order_release);
    return data;
}

inline std::array<LC,5> polynomialJets(const CoefficientData& data,LC y,bool reversed) {
    const LC argument=reversed?LC(1.L,0.L)/y:y;
    std::array<LC,5> jet{};
    for(int order=0;order<=requiredThetaOrder(data.candidate);++order){
        const auto& coefficients=data.thetaDerivatives[order];
        LC value{};
        if(reversed){
            for(long double a:coefficients)value=value*argument+a;
        }else{
            for(auto a=coefficients.rbegin();a!=coefficients.rend();++a)value=value*argument+*a;
        }
        jet[order]=value;
    }
    return jet;
}

inline std::array<LC,3> dilationJets(const CoefficientData& data,LC y,bool reversed) {
    std::array<LC,3> jet{};
    if(!reversed){
        const auto ordinary=polynomialJets(data,y,false);
        for(int order=0;order<=2;++order)jet[order]=ordinary[order];
        return jet;
    }
    const LC inverse=LC(1.L,0.L)/y;
    for(int order=0;order<=2;++order){
        const auto& coefficients=data.inverseThetaDerivatives[order];
        LC value{};
        for(auto a=coefficients.rbegin();a!=coefficients.rend();++a)value=value*inverse+*a;
        jet[order]=value;
    }
    return jet;
}
// Evaluate T^(k)(w)/M^N in a shared real scale M>=max(1, |w|).
// Horner differentiates the same scaled polynomial. All multi-argument
// products share M^(2N), which cancels exactly in their Newton quotient.
// This avoids forming y^N or destructive independent normalizations.
inline std::array<LC,5> boundedThetaJets(const CoefficientData& data,LC w,
                                         long double M,int maximumOrder) {
    std::array<LC,5> jet{};
    if (!finite(w)||!std::isfinite(M)||M<1.L||maximumOrder<0||maximumOrder>4)
        return {{badLong(),badLong(),badLong(),badLong(),badLong()}};
    const LC x=w/M;
    const long double inverse=1.L/M;
    long double inversePower=1.L;
    for(int j=data.degree;j>=0;--j) {
        for(int k=maximumOrder;k>=1;--k)
            jet[k]=x*jet[k]+static_cast<long double>(k)*jet[k-1];
        jet[0]=x*jet[0]+data.theta[static_cast<std::size_t>(j)]*inversePower;
        inversePower*=inverse;
    }
    long double derivativeScale=inverse;
    for(int k=1;k<=maximumOrder;++k) {
        jet[k]*=derivativeScale;
        derivativeScale*=inverse;
    }
    return jet;
}
} // namespace detail

class Model {
    Settings settings_;
    std::shared_ptr<const detail::CoefficientData> coefficients_;
    bool valid_ = false;

    detail::LC evaluateCorrection(detail::LC y) const {
        using namespace detail;
        const int candidate = settings_.candidate;
        const long double level = settings_.level;
        const bool reversed = std::abs(y) > 1.L;
        if (candidate == 0 || candidate == 3 || candidate == 6 || candidate == 7) {
            const auto t = polynomialJets(*coefficients_,y,reversed);
            if (!finite(t[0])) return badLong();
            // g and its indicated derivatives have poles at every nonzero
            // zero of T. Their Newton correction has the exact value zero.
            if (t[0] == LC(0.L,0.L)) return LC(0.L,0.L);
            const LC a = t[1]/t[0], b = t[2]/t[0];
            const LC l2 = b-a*a;
            const LC q = reversed?a+l2:a+y*l2;
            if(candidate==3)return reversed?y*quotient(a-level,q):quotient(y*a-level,q);
            const LC c = t[3]/t[0];
            const LC l3 = c-3.L*a*b+2.L*a*a*a;
            const LC qp = reversed?2.L*l2+l3:2.L*l2+y*l3;
            if(candidate==0)return reversed?y*quotient(q,qp):quotient(q,qp);
            const LC d = t[4]/t[0];
            const LC l4 = d-4.L*a*c-3.L*b*b+12.L*a*a*b-6.L*a*a*a*a;
            if (reversed) {
                const LC qpp = 3.L*l3+l4;
                if (candidate == 6) return y*quotient(qp,qpp);
                if (candidate == 7) return y*quotient(q+qp,2.L*qp+qpp);
            } else {
                const LC qpp = 3.L*l3+y*l4;
                if (candidate == 6) return quotient(qp,qpp);
                if (candidate == 7) return quotient(q+y*qp,2.L*qp+y*qpp);
            }
            return badLong();
        }
        if (candidate == 5) {
            const LC v = -y;
            const auto t = polynomialJets(*coefficients_,v,reversed);
            // H=T*T'+v*T*T''-v*T'^2 and
            // H'=2*T*T''+v*(T*T'''-T'*T''). Keep products factored;
            // expanding their high-degree coefficients is ill-conditioned.
            if (reversed) {
                const LC h = t[0]*t[1]+t[0]*t[2]-t[1]*t[1];
                const LC hp = 2.L*t[0]*t[2]+t[0]*t[3]-t[1]*t[2];
                LC inversePower(1.L,0.L);
                const LC inverse = LC(1.L,0.L)/v;
                for (int j = 0; j < 2*settings_.degree; ++j) inversePower *= inverse;
                return v*quotient(-2.L*h-level*inversePower,2.L*(h+hp));
            }
            const LC h = t[0]*t[1]+v*(t[0]*t[2]-t[1]*t[1]);
            const LC hp = 2.L*t[0]*t[2]+v*(t[0]*t[3]-t[1]*t[2]);
            return quotient(2.L*y*h-level,2.L*h-2.L*y*hp);
        }
        if (candidate == 11) {
            // Each evaluation chooses its own normal/reversed jet. This
            // matters when |y|<=1 but |2y|>1. Reversed jets divide out the
            // otherwise enormous common powers before forming C and C'.
            const LC doubled=2.L*y;
            if(!finite(y)||!finite(doubled))return badLong();
            const auto t=dilationJets(*coefficients_,y,reversed);
            const bool doubleReversed=std::abs(doubled)>1.L;
            const auto u=dilationJets(*coefficients_,doubled,doubleReversed);
            if(!finite(t[0])||!finite(u[0]))return badLong();
            const bool firstPole=t[0]==LC{},secondPole=u[0]==LC{};
            // Joint poles can cancel; no continuation is asserted there.
            if(firstPole&&secondPole)return badLong();
            if(firstPole||secondPole)return LC{};
            const auto countJet=[](const std::array<LC,3>& jet,LC w,bool reverse){
                const LC a=jet[1]/jet[0],b=jet[2]/jet[0];
                const LC l2=b-a*a;
                if(reverse){
                    const LC inverse=LC(1.L,0.L)/w;
                    // B(u)=u^N*T(1/u): C=N-u*B'/B. Keep only the
                    // offset from N, so dyadic subtraction cannot erase
                    // the gap by rounding both large counts to N.
                    return std::array<LC,2>{{-inverse*a,inverse*inverse*(a+inverse*l2)}};
                }
                return std::array<LC,2>{{w*a,a+w*l2}};
            };
            const auto atY=countJet(t,y,reversed),atDouble=countJet(u,doubled,doubleReversed);
            const int baseDifference=(doubleReversed?settings_.degree:0)-(reversed?settings_.degree:0);
            return quotient(atDouble[0]-atY[0]+static_cast<long double>(baseDifference)-level,
                            2.L*atDouble[1]-atY[1]);
        }
        if (candidate >= 13 && candidate <= 17) {
            // Structural generation 1.8: finite exact polynomial targets.
            // A shared positive scale M removes the common M^(2N) in F/F'.
            LC a=y, b=y, c=y;
            int order=0;
            if(candidate==13){a=y+static_cast<long double>(settings_.radius);b=y-static_cast<long double>(settings_.radius);order=1;}
            if(candidate==14){order=3;}
            if(candidate==15){a=2.L*y;b=.5L*y;order=1;}
            if(candidate==16){b=-y;order=2;}
            if(candidate==17){order=4;}
            long double M=std::max(1.L,std::abs(a));
            if(candidate==13||candidate==15||candidate==16)M=std::max(M,std::abs(b));
            if(candidate==13||candidate==15)M=std::max(M,std::abs(c));
            if(!std::isfinite(M))return badLong();
            const auto A=boundedThetaJets(*coefficients_,a,M,order);
            if(!finite(A[0])||!finite(A[order]))return badLong();
            LC F{},Fp{};
            if(candidate==14) {
                F=A[0]*A[2]-A[1]*A[1];
                Fp=A[0]*A[3]-A[1]*A[2];
            } else if(candidate==17) {
                F=2.L*A[1]*A[3]-3.L*A[2]*A[2];
                Fp=2.L*A[1]*A[4]-4.L*A[2]*A[3];
            } else {
                const auto B=boundedThetaJets(*coefficients_,b,M,order);
                if(!finite(B[0])||!finite(B[order]))return badLong();
                if(candidate==16) {
                    F=A[1]*B[0]-A[0]*B[1];
                    Fp=A[2]*B[0]-2.L*A[1]*B[1]+A[0]*B[2];
                } else {
                    const auto C=boundedThetaJets(*coefficients_,c,M,order);
                    if(!finite(C[0])||!finite(C[order]))return badLong();
                    F=A[0]*B[0]-C[0]*C[0];
                    Fp=(candidate==13?1.L:2.L)*A[1]*B[0]
                        +(candidate==13?1.L:.5L)*A[0]*B[1]-2.L*C[0]*C[1];
                }
            }
            return quotient(F,Fp);
        }
        const auto e = [&](int radialOrder,int angularOrder) {
            return clenshaw(coefficients_->angular[radialOrder][angularOrder],y);
        };
        if (candidate == 1) return quotient(e(0,1)-level,e(0,2));
        if (candidate == 2) {
            const LC E = e(0,0), U = e(0,1), Ux = e(0,2);
            return quotient((U-level*E)*E,Ux*E-U*U);
        }
        if (candidate == 4) return quotient(e(0,2),e(0,3));
        if (candidate == 8) return quotient(e(1,1),e(1,2));
        if (candidate == 9) {
            const LC Ux = e(0,2), Uxx = e(0,3), Uxxx = e(0,4);
            const LC Ur = e(1,1), Urx = e(1,2), Urxx = e(1,3);
            const LC Urr = e(2,1), Urrx = e(2,2);
            const LC F = Urr*Ux*Ux-2.L*Urx*Ur*Ux+Uxx*Ur*Ur;
            const LC Fx = Urrx*Ux*Ux+2.L*Urr*Ux*Uxx
                -2.L*Urxx*Ur*Ux-2.L*Urx*Urx*Ux+Uxxx*Ur*Ur;
            return quotient(F,Fx);
        }
        if (candidate == 10) {
            const LC E=e(0,0),U=e(0,1),Ux=e(0,2),Uxx=e(0,3);
            if(!finite(E)||!finite(U)||!finite(Ux)||!finite(Uxx))return badLong();
            // A zero of E is a double pole of (E'/E)' regardless of its
            // multiplicity, so the Newton correction extends to zero.
            if(E==LC{})return LC{};
            const LC D=E*Ux-U*U,Dp=E*Uxx-U*Ux;
            return quotient((D-level*E*E)*E,Dp*E-2.L*D*U);
        }
        if (candidate == 12) {
            const LC Ux=e(0,2),Uxx=e(0,3),Ur=e(1,1),Urx=e(1,2);
            if(!finite(Ux)||!finite(Uxx)||!finite(Ur)||!finite(Urx))return badLong();
            // A common zero may remove the pole of U_R/U_x, with a
            // nonzero quotient. Leave that indeterminate point numerical.
            if(Ur==LC{}&&Ux==LC{})return badLong();
            return quotient((Ur-level*Ux)*Ux,Urx*Ux-Ur*Uxx);
        }
        return badLong(); // New registry entries require their own implementation.
    }

public:
    explicit Model(const Settings& s) : settings_(s) {
        if (s.candidate < 0 || s.candidate >= CandidateCount || s.degree < 4 || s.degree > 32
            || !std::isfinite(s.radius) || s.radius <= 0.
            || !std::isfinite(s.scale) || s.scale <= 0.
            || !std::isfinite(s.shift) || !std::isfinite(s.level)
            || !detail::finite(s.lambda) || !detail::finite(s.parameter)) return;
        coefficients_=detail::coefficientData(s.candidate,s.degree,s.radius);
        valid_ = coefficients_&&coefficients_->valid;
    }
    bool valid() const { return valid_; }
    const Settings& settings() const { return settings_; }
    const std::shared_ptr<const detail::CoefficientData>& coefficientData()const{return coefficients_;}
    C correction(C z) const {
        if (!valid_ || !detail::finite(z)) return detail::nonfinite();
        const detail::LC y = static_cast<long double>(settings_.scale)
            *detail::LC(z.real(),z.imag())+static_cast<long double>(settings_.shift);
        const detail::LC q = evaluateCorrection(y)/static_cast<long double>(settings_.scale);
        if (!detail::finite(q)) return detail::nonfinite();
        const C result(static_cast<double>(q.real()),static_cast<double>(q.imag()));
        return detail::finite(result) ? result : detail::nonfinite();
    }
    Result orbit(C initial, C parameter, C strength, int budget,
                 std::vector<C>* trail = nullptr) const {
        Result result;
        result.endpoint = initial;
        if (trail) trail->clear();
        if (!valid_ || !detail::finite(parameter) || !detail::finite(strength)) {
            result.status = Status::Invalid; return result;
        }
        if (!detail::finite(initial)) { result.status = Status::Numerical; return result; }
        budget = std::clamp(budget,0,4096);
        C z = initial;
        double minimumTrap = std::numeric_limits<double>::infinity();
        for (int n = 0; n <= budget; ++n) {
            result.endpoint = z;
            result.iterations = n;
            result.smooth = static_cast<double>(n);
            if (trail) trail->push_back(z);
            minimumTrap = std::min(minimumTrap,std::abs(z-C(-.21,.37)));
            result.structure = std::exp(-2.*minimumTrap);
            result.phase = .5+std::arg(z)/6.2831853071795864769;
            const double component = std::max(std::abs(z.real()),std::abs(z.imag()));
            if (component > 16. || std::norm(z) > 256.) {
                result.status = Status::Cutoff;
                const double m = std::max(std::abs(z.real()),std::abs(z.imag()));
                const C scaled = z/m;
                const double logRadius = std::log(m)+.5*std::log(std::norm(scaled));
                result.smooth = std::max(0.,n+1.-std::log(std::max(1.,logRadius))/std::log(2.));
                return result;
            }
            if (n == budget) return result;
            C next = z*z+parameter;
            // EXACT bypass: zero strength does not evaluate a correction,
            // including at any pole of the rational carrier.
            if (strength != C(0.,0.)) {
                const C q = correction(z);
                if (!detail::finite(q)) { result.status = Status::Numerical; return result; }
                next -= strength*q;
            }
            if (!detail::finite(next)) { result.status = Status::Numerical; return result; }
            z = next;
        }
        return result;
    }
};

inline uint32_t colour(const Result& r, int palette = 0, double hue = 0.) {
    if (r.status == Status::Invalid) return 0x862330;
    if (r.status == Status::Numerical) return 0xe052c5;
    constexpr double tau = 6.2831853071795864769;
    if (!std::isfinite(hue)) hue = 0.;
    const double smooth = std::isfinite(r.smooth) ? r.smooth : 0.;
    const double structure = std::clamp(r.structure,0.,1.);
    const double phase = std::isfinite(r.phase) ? r.phase : 0.;
    double h, saturation, value;
    if (r.status == Status::Unresolved) {
        h = hue+.61+.14*phase+.08*structure;
        saturation = .55;
        value = .018+.115*structure;
    } else {
        const double stripe = .5+.5*std::cos(.57*smooth+2.4*structure);
        h = palette%3 == 1 ? hue+.66-.023*smooth
            : palette%3 == 2 ? hue+.65*structure+.18*phase
            : hue+.57+.30*phase+.012*smooth;
        saturation = .72+.11*std::sin(tau*structure);
        value = .29+.69*stripe;
    }
    h -= std::floor(h);
    saturation = std::clamp(saturation,0.,1.);
    value = std::clamp(value,0.,1.);
    const double sector = h*6.;
    const int k = static_cast<int>(sector);
    const double f = sector-k, p = value*(1.-saturation);
    const double q = value*(1.-saturation*f), t = value*(1.-saturation*(1.-f));
    double red,green,blue;
    switch (k%6) {
    case 0: red=value;green=t;blue=p;break;
    case 1: red=q;green=value;blue=p;break;
    case 2: red=p;green=value;blue=t;break;
    case 3: red=p;green=q;blue=value;break;
    case 4: red=t;green=p;blue=value;break;
    default:red=value;green=p;blue=q;break;
    }
    return (static_cast<uint32_t>(255.*red+.5)<<16)
        | (static_cast<uint32_t>(255.*green+.5)<<8)
        | static_cast<uint32_t>(255.*blue+.5);
}
} // namespace rhfresh
