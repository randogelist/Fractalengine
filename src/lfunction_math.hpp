#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <utility>

// Numerical values and analytic s-derivatives of actual L-functions.
// Generation 1: zeta; 2: chi mod 3; 3: chi mod 4 (Dirichlet beta);
// 4: primitive chi mod 5 with chi(2)=i, hence [1,i,-i,-1].
// Sources: https://dlmf.nist.gov/25.11.E7 (Euler--Maclaurin),
// https://dlmf.nist.gov/25.15.E3 and https://dlmf.nist.gov/25.15.E5,
// https://dlmf.nist.gov/5.11.E1 (log Gamma asymptotics).
// Operational domain: -8 <= Re(s) <= 12, |Im(s)| <= 80. The zeta pole
// |s-1| < 1e-8 is excluded. valid is a numerical-domain/finite-value flag,
// not an interval certificate. No RH/GRH hypothesis or preloaded zeros.
// A basin renderer MUST use an absolute |f| residual, not a residual
// relative to a cancelled Dirichlet sum, to accept a zero.
namespace lfunction {
using C = std::complex<double>;
struct Evaluation { C f{}, first{}, second{}; bool valid = false; };
namespace detail {
constexpr double pi = 3.141592653589793238462643383279502884;
struct J { C v{}, d{}, dd{}; };
inline J add(const J& a, const J& b) { return {a.v+b.v,a.d+b.d,a.dd+b.dd}; }
inline J scale(const J& a, C b) { return {a.v*b,a.d*b,a.dd*b}; }
inline J mul(const J& a, const J& b) {
    return {a.v*b.v, a.d*b.v+a.v*b.d,
            a.dd*b.v+2.0*a.d*b.d+a.v*b.dd};
}
inline J reciprocal(C z) {
    C r=1.0/z;
    return {r,-r*r,2.0*r*r*r};
}
inline J exponential(C s, double coefficient) {
    C v=std::exp(s*coefficient);
    return {v,coefficient*v,coefficient*coefficient*v};
}
inline bool finite(C z) { return std::isfinite(z.real())&&std::isfinite(z.imag()); }
inline int modulus(int generation) { return generation==1 ? 1 : generation+1; }
inline C character(int n, int generation, bool conjugate=false) {
    if (generation==1) return C(1,0);
    int r=n%modulus(generation);
    if (r==0) return C(0,0);
    if (generation==2) return C(r==1 ? 1:-1,0);
    if (generation==3) return C(r==1 ? 1:(r==3 ? -1:0),0);
    if (r==1) return C(1,0);
    if (r==4) return C(-1,0);
    return C(0,(r==2 ? 1:-1)*(conjugate ? -1:1));
}
inline const std::array<double,256>& logTable() {
    static const std::array<double,256> t=[] {
        std::array<double,256> a{};
        for (int n=1;n<256;++n) a[n]=std::log(double(n));
        return a;
    }();
    return t;
}
struct Sum {
    J s{}, correction{};
    void append(J a) {
        C y=a.v-correction.v,t=s.v+y; correction.v=(t-s.v)-y; s.v=t;
        y=a.d-correction.d; t=s.d+y; correction.d=(t-s.d)-y; s.d=t;
        y=a.dd-correction.dd; t=s.dd+y; correction.dd=(t-s.dd)-y; s.dd=t;
    }
};
// (x^(1-s)-X^(1-s))/(q*(s-1)), with the removable singularity
// evaluated as an entire exprel series near s=1. The subtracted common
// term cancels because every nonprincipal character here has sum zero.
inline J integralDifference(C s, double lx, double lX, double q) {
    C z=s-1.0;
    J R=exponential(z,-lX);
    if (std::abs(z)>0.5) {
        J E=exponential(z,-lx);
        return scale(mul(add(E,scale(R,-1.0)),reciprocal(z)),1.0/q);
    }
    double h=lx-lX;
    J u{-z*h,C(-h,0),C(0,0)};
    // Horner series sum_{k=0}^{14} u^k/(k+1)!; |u| < .031 here.
    double c=1.0;
    for (int j=2;j<=15;++j) c/=j;
    J e{C(c,0),{}, {}};
    for (int k=13;k>=0;--k) {
        c*=k+2;
        e=add(mul(e,u),J{C(c,0),{}, {}});
    }
    return scale(mul(R,e),-h/q);
}
inline J eulerMaclaurin(C s,int generation,bool conjugate=false) {
    static constexpr double bc[16]={
        8.33333333333333333333e-2,-1.38888888888888888889e-3,
        3.30687830687830687831e-5,-8.26719576719576719577e-7,
        2.08767569878680989792e-8,-5.28419013868749318485e-10,
        1.33825365306846788328e-11,-3.38968029632258286683e-13,
        8.58606205627784456414e-15,-2.17486869855806187304e-16,
        5.50900282836022951520e-18,-1.39544646858125233407e-19,
        3.53470703962946747169e-21,-8.95351742703754685040e-23,
        2.26795245233768306031e-24,-5.74479066887220244526e-26};
    const int q=modulus(generation);
    const int N=16+int(std::ceil(std::abs(s.imag())/3.0));
    const auto& logs=logTable();
    Sum total;
    // Prefix is evaluated in physical integers: no subtraction of enormous
    // independent Hurwitz values. Reflection keeps Re(s) >= -0.25 here.
    int last=q==1 ? N : q*N-1;
    for (int n=1;n<=last;++n) {
        C chi=character(n,generation,conjugate);
        if (chi!=C(0,0)) total.append(scale(exponential(s,-logs[n]),chi));
    }
    const int residues=q==1 ? 1 : q-1;
    for (int r=1;r<=residues;++r) {
        C chi=character(r,generation,conjugate);
        if (chi==C(0,0)) continue;
        double x=q*N+r, lx=logs[int(x)], ratio=double(q)/x;
        J base=exponential(s,-lx);
        J tail=scale(base,0.5);
        if (q==1) tail=add(tail,scale(mul(base,reciprocal(s-1.0)),x));
        else tail=add(tail,integralDifference(s,lx,logs[q*N],double(q)));
        J rising{C(1,0),{}, {}};
        double rpower=1.0;
        for (int j=0;j<31;++j) {
            rising=mul(rising,J{s+double(j),C(1,0),C(0,0)});
            rpower*=ratio;
            if ((j&1)==0) {
                int k=j/2;
                tail=add(tail,scale(mul(base,rising),bc[k]*rpower));
            }
        }
        total.append(scale(tail,chi));
    }
    return total.s;
}
// Fused quartic-character pair evaluation. Both primitive characters modulo 5
// share the same powers, Euler--Maclaurin tail, and derivative basis; only the
// character coefficient is conjugated. This is algebraically identical to two
// separate eulerMaclaurin calls but avoids duplicating the expensive complex
// exponentials and rising-factorial construction.
inline std::pair<J,J> eulerMaclaurinQuarticPair(C s) {
    static constexpr double bc[16]={
        8.33333333333333333333e-2,-1.38888888888888888889e-3,
        3.30687830687830687831e-5,-8.26719576719576719577e-7,
        2.08767569878680989792e-8,-5.28419013868749318485e-10,
        1.33825365306846788328e-11,-3.38968029632258286683e-13,
        8.58606205627784456414e-15,-2.17486869855806187304e-16,
        5.50900282836022951520e-18,-1.39544646858125233407e-19,
        3.53470703962946757169e-21,-8.95351742703754685040e-23,
        2.26795245233768306031e-24,-5.74479066887220244526e-26};
    constexpr int q=5;
    const int N=16+int(std::ceil(std::abs(s.imag())/3.0));
    const auto& logs=logTable();
    Sum chi,bar;
    int last=q*N-1;
    for (int n=1;n<=last;++n) {
        C c=character(n,4,false);
        if (c==C(0,0)) continue;
        J basis=exponential(s,-logs[n]);
        chi.append(scale(basis,c));
        bar.append(scale(basis,std::conj(c)));
    }
    for (int r=1;r<q;++r) {
        C c=character(r,4,false);
        if (c==C(0,0)) continue;
        double x=q*N+r, lx=logs[int(x)], ratio=double(q)/x;
        J base=exponential(s,-lx);
        J tail=scale(base,0.5);
        tail=add(tail,integralDifference(s,lx,logs[q*N],double(q)));
        J rising{C(1,0),{}, {}};
        double rpower=1.0;
        for (int j=0;j<31;++j) {
            rising=mul(rising,J{s+double(j),C(1,0),C(0,0)});
            rpower*=ratio;
            if ((j&1)==0) {
                int k=j/2;
                tail=add(tail,scale(mul(base,rising),bc[k]*rpower));
            }
        }
        chi.append(scale(tail,c));
        bar.append(scale(tail,std::conj(c)));
    }
    return {chi.s,bar.s};
}
struct GammaJet { C logGamma,psi,trigamma; };
inline GammaJet logGammaJet(C z) {
    // Recurrence into Re(z)>=12, then differentiated Stirling expansion.
    static constexpr double c[10]={1.0/12.0,-1.0/360.0,1.0/1260.0,
        -1.0/1680.0,1.0/1188.0,-691.0/360360.0,1.0/156.0,
        -3617.0/122400.0,43867.0/244188.0,-174611.0/125400.0};
    C correction=0,psiCorrection=0,triCorrection=0;
    while(z.real()<12.0) {
        correction+=std::log(z); psiCorrection+=1.0/z;
        triCorrection+=1.0/(z*z); z+=1.0;
    }
    C inv=1.0/z, inv2=inv*inv, power=inv;
    C lg=(z-0.5)*std::log(z)-z+0.5*std::log(2.0*pi);
    C psi=std::log(z)-0.5*inv, tri=inv+0.5*inv2;
    for (int k=0;k<10;++k) {
        int p=2*k+1;
        lg+=c[k]*power;
        psi-=double(p)*c[k]*power*inv;
        tri+=double(p*(p+1))*c[k]*power*inv2;
        power*=inv2;
    }
    return {lg-correction,psi-psiCorrection,tri+triCorrection};
}
inline C rootNumber(int generation) {
    if (generation!=4) return C(1,0);
    static const C w=[] {
        C tau=0;
        for(int r=1;r<5;++r) tau+=character(r,4)*std::exp(C(0,2.0*pi*r/5.0));
        return tau/C(0,std::sqrt(5.0));
    }();
    return w;
}
inline std::pair<J,J> reflectedQuarticPair(C s) {
    constexpr int q=5;
    GammaJet gamma=logGammaJet(1.0-s);
    double logScale=std::log(2.0*pi/double(q));
    C common=(std::sqrt(double(q))/pi)*std::exp(gamma.logGamma+s*logScale);
    C w=rootNumber(4);
    C dlog=logScale-gamma.psi;
    J Achi{w*common,w*common*dlog,w*common*(dlog*dlog+gamma.trigamma)};
    C wb=std::conj(w);
    J Abar{wb*common,wb*common*dlog,wb*common*(dlog*dlog+gamma.trigamma)};
    C arg=(pi/2.0)*s;
    C c=std::cos(arg),dc=-(pi/2.0)*std::sin(arg);
    J T{c,dc,-(pi*pi/4.0)*c};
    auto pair=eulerMaclaurinQuarticPair(1.0-s);
    pair.first.d=-pair.first.d;
    pair.second.d=-pair.second.d;
    // Functional equation swaps chi and chibar on the reflected Dirichlet side.
    return {mul(mul(Achi,T),pair.second),mul(mul(Abar,T),pair.first)};
}
inline J reflected(C s,int generation) {
    int q=modulus(generation);
    GammaJet gamma=logGammaJet(1.0-s);
    double logScale=std::log(2.0*pi/double(q));
    C a=rootNumber(generation)*(std::sqrt(double(q))/pi)*
        std::exp(gamma.logGamma+s*logScale);
    C dlog=logScale-gamma.psi;
    J A{a,a*dlog,a*(dlog*dlog+gamma.trigamma)};
    C arg=(pi/2.0)*s;
    C t=generation==1 ? std::sin(arg) : std::cos(arg);
    C dt=(pi/2.0)*(generation==1 ? std::cos(arg) : -std::sin(arg));
    J T{t,dt,-(pi*pi/4.0)*t};
    J L=eulerMaclaurin(1.0-s,generation,generation==4);
    L.d=-L.d;
    return mul(mul(A,T),L);
}
} // namespace detail
inline std::pair<Evaluation,Evaluation> evaluateQuarticPair(C s) {
    if (!detail::finite(s)||s.real() < -8.0||s.real()>12.0||std::abs(s.imag())>80.0)
        return {};
    std::pair<detail::J,detail::J> j = s.real() < -0.25
        ? detail::reflectedQuarticPair(s)
        : detail::eulerMaclaurinQuarticPair(s);
    bool oka=detail::finite(j.first.v)&&detail::finite(j.first.d)&&detail::finite(j.first.dd);
    bool okb=detail::finite(j.second.v)&&detail::finite(j.second.d)&&detail::finite(j.second.dd);
    return {{j.first.v,j.first.d,j.first.dd,oka},{j.second.v,j.second.d,j.second.dd,okb}};
}
inline Evaluation evaluate(C s,int generation) {
    if (generation<1||generation>4||!detail::finite(s)||s.real() < -8.0||
        s.real()>12.0||std::abs(s.imag())>80.0||
        (generation==1&&std::abs(s-1.0)<1e-8)) return {};
    detail::J j=s.real() < -0.25 ? detail::reflected(s,generation) :
                                               detail::eulerMaclaurin(s,generation);
    bool ok=detail::finite(j.v)&&detail::finite(j.d)&&detail::finite(j.dd);
    return {j.v,j.d,j.dd,ok};
}
} // namespace lfunction
