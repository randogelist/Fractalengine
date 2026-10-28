#pragma once
#include "lfunction_math.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <vector>

// Finite completed-L carrier for the project-linked GRH rational map.
//
// chi is the odd primitive quartic character modulo 5, chi(2)=i.
// With theta_chi(t)=sum_{n>=1} n chi(n) exp(-pi n^2 t/5), the split
// Mellin formula at the centre is
//
// Lambda_chi(1/2+w) = integral_0^infinity exp(3u/4) *
//   [theta_chi(exp(u)) exp(wu/2)
//    + epsilon_chi theta_chibar(exp(u)) exp(-wu/2)] du.
//
// Its coefficient integrals avoid Cauchy-circle cancellation and gamma
// factors at removable trivial zeros.  The paired functional equation gives
// P(1/2+w)=P(1/2-w), where P=Lambda_chi Lambda_chibar, so the target
// H(w)=P(1/2+w)P(1/2-w) is even.  We form its finite Taylor polynomial.
//
// Native runtime construction uses a fixed finite theta sum and quadrature,
// NOT an interval enclosure.  The result defines a finite rational model
// with numerically constructed coefficients; its roots are not certified L-function
// zeros.  Its Julia dynamics must not be labelled the exact infinite H map.
// Increasing |w| does not establish escape for a logarithmic-derivative map.
namespace project_grh {
using C = std::complex<double>;
constexpr double pi = 3.141592653589793238462643383279502884;

inline bool finite(C z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

struct Jet {
    C f{}, first{}, second{};
    bool valid = false;
};

struct Carrier {
    // H_M(w), up to a nonzero constant, in powers of t=w/radius.
    // 'radius' is solely a coefficient scale, not a certified domain.
    std::vector<C> coefficients;
    double radius = 12.0;
    int degree = 0;
    bool valid = false;
    std::string status;

    Jet evaluate(C w) const {
        if (!valid || coefficients.empty() || !finite(w)) return {};
        const C t=w/radius;
        C f=coefficients.back(), d=0, dd=0;
        for (int n=int(coefficients.size())-2;n>=0;--n) {
            dd=dd*t+2.0*d;
            d=d*t+f;
            f=f*t+coefficients[std::size_t(n)];
        }
        d/=radius;
        dd/=(radius*radius);
        return {f,d,dd,finite(f)&&finite(d)&&finite(dd)};
    }
};

namespace detail {
using LC = std::complex<long double>;
constexpr long double piL = 3.141592653589793238462643383279502884L;

// Long-double accumulators reduce summation loss in coefficient construction.
// No claimed precision beyond the unvalidated finite quadrature is implied.
struct Sum {
    LC value{}, correction{};
    void add(LC x) {
        const LC y=x-correction, next=value+y;
        correction=(next-value)-y;
        value=next;
    }
};

inline LC character(int n) {
    switch (n%5) {
    case 1: return LC(1,0);
    case 2: return LC(0,1);
    case 3: return LC(0,-1);
    case 4: return LC(-1,0);
    default: return LC(0,0);
    }
}

inline LC rootNumber() {
    Sum gauss;
    for (int r=1;r<5;++r)
        gauss.add(character(r)*std::exp(LC(0,2.0L*piL*r/5.0L)));
    return gauss.value/LC(0,std::sqrt(5.0L));
}

inline LC theta(long double t) {
    Sum s;
    // At t>=1 the first omitted absolute term is already exponentially
    // small.  This remains a finite numerical truncation, not a certificate.
    for (int n=1;n<=16;++n) {
        const LC chi=character(n);
        if (chi!=LC(0,0))
            s.add(chi*(static_cast<long double>(n)*
                std::exp(-piL*static_cast<long double>(n*n)*t/5.0L)));
    }
    return s.value;
}
} // namespace detail

inline Carrier makeCarrier(int requestedDegree=24, double radius=12.0) {
    Carrier out;
    out.radius=radius;
    // Even degrees only: odd coefficients vanish by the paired equation.
    if (requestedDegree<4 || requestedDegree>64 || (requestedDegree&1) ||
        !std::isfinite(radius) || radius<1.0 || radius>32.0) {
        out.status="invalid GRH carrier degree/scale";
        return out;
    }
    const int degree=requestedDegree;
    out.degree=degree;
    constexpr int panels=2048;
    const long double upper=std::log(128.0L);
    const long double h=upper/static_cast<long double>(panels);
    const detail::LC epsilon=detail::rootNumber();
    std::vector<detail::Sum> integral(std::size_t(degree+1));
    for (int j=0;j<=panels;++j) {
        const long double u=h*static_cast<long double>(j);
        const long double t=std::exp(u);
        const detail::LC theta=detail::theta(t);
        const detail::LC reflected=epsilon*std::conj(theta);
        const long double simpson=(j==0 || j==panels)?1.0L:((j&1)?4.0L:2.0L);
        const long double weight=simpson*std::exp(0.75L*u)*(h/3.0L);
        long double power=1.0L;
        const long double argument=static_cast<long double>(radius)*u/2.0L;
        for (int n=0;n<=degree;++n) {
            integral[std::size_t(n)].add(weight*power*
                (theta+((n&1)?-reflected:reflected)));
            power*=argument/static_cast<long double>(n+1);
        }
    }

    // Lambda_chibar has conjugated Taylor coefficients at the real centre.
    // P's odd coefficients are analytically zero; imposing that symmetry
    // removes numerical antisymmetric roundoff without moving a parameter.
    std::vector<long double> pair(std::size_t(degree+1),0.0L);
    for (int n=0;n<=degree;n+=2) {
        detail::Sum coefficient;
        for (int j=0;j<=n;++j)
            coefficient.add(integral[std::size_t(j)].value*
                            std::conj(integral[std::size_t(n-j)].value));
        pair[std::size_t(n)]=coefficient.value.real();
    }
    std::vector<long double> result(std::size_t(degree+1),0.0L);
    long double scale=0;
    for (int n=0;n<=degree;n+=2) {
        long double value=0, correction=0;
        for (int j=0;j<=n;j+=2) {
            const long double y=pair[std::size_t(j)]*pair[std::size_t(n-j)]-correction;
            const long double next=value+y;
            correction=(next-value)-y;
            value=next;
        }
        result[std::size_t(n)]=value;
        scale=std::max(scale,std::abs(value));
    }
    if (!(scale>0) || !std::isfinite(scale)) {
        out.status="GRH carrier coefficient construction failed";
        return out;
    }
    out.coefficients.resize(std::size_t(degree+1));
    for (int n=0;n<=degree;++n) {
        const double value=static_cast<double>(result[std::size_t(n)]/scale);
        if (!std::isfinite(value)) {
            out.coefficients.clear();
            out.status="nonfinite GRH carrier coefficient";
            return out;
        }
        out.coefficients[std::size_t(n)]=C(value,0);
    }
    if (out.coefficients.front()==C(0,0) || out.coefficients.back()==C(0,0)) {
        out.coefficients.clear();
        out.status="degenerate GRH carrier coefficients";
        return out;
    }
    out.valid=true;
    out.status="chi5 completed-pair H: even Taylor degree "+std::to_string(degree)+
        "; finite theta/Mellin quadrature (16 terms, 2048 panels, t<=128)";
    return out;
}
} // namespace project_grh
