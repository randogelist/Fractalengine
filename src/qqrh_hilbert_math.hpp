#pragma once
#include "theta_coefficients.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

// Finite numerical realization of the recovered QQRH Hardy/shift construction.
// Theta(z)=sum a_n z^n; the existing table stores b_n=200^n*a_n/a_0.
// Here r=200*R and v_n=b_n*R^n / sqrt(sum b_j^2*R^(2j)), 0<=n<=12.
// Thus c_k=<B^k v,v>. These finite vectors and their Gram matrices do not
// establish positivity of the full infinite QQRH angular kernel or RH.
namespace qqrh_hilbert {
using C = std::complex<double>;
inline constexpr int degree=12;

inline double bounded_radius(double value) {
    return std::clamp(std::isnan(value) ? .6 : value,.1,1.25);
}

struct Model {
    std::array<double,degree+1> v{},c{};
    double radius;
    double angular_norm=0.;

    explicit Model(double R) : radius(bounded_radius(R)) {
        double power=1.,norm2=0.;
        for(int n=0;n<=degree;++n) {
            v[n]=ff_rh::THETA_COEFFICIENTS[n]*power;
            norm2+=v[n]*v[n];
            power*=radius;
        }
        const double inverse_norm=1./std::sqrt(norm2);
        for(double& value:v)value*=inverse_norm;
        for(int k=0;k<=degree;++k)
            for(int n=0;n+k<=degree;++n)c[k]+=v[n]*v[n+k];
        for(int k=1;k<=degree;++k)angular_norm+=double(k*k)*c[k];
    }

    // Holomorphic overlap H(z)=<w(z),v> for w(z)=(v_n*z^n)_n.
    // Its zeros are orthogonality points, not zeros of the full Theta.
    C kernel(C z,C* derivative=nullptr) const {
        C value(v[degree]*v[degree],0.),d(0.,0.);
        for(int n=degree-1;n>=0;--n) {
            d=d*z+value;
            value=value*z+v[n]*v[n];
        }
        if(derivative)*derivative=d;
        return value;
    }

    // Pullback of the Fubini-Study metric, with ds^2=g(z)|dz|^2:
    // g=(||w||^2 ||w'||^2-|<w',w>|^2)/||w||^4.
    // A nonnegative pair sum avoids subtracting almost equal moments.
    double projective_metric(C z) const {
        if(!std::isfinite(z.real()) || !std::isfinite(z.imag()))return 0.;
        const double magnitude=std::abs(z);
        if(magnitude<=1.) {
            const double t=magnitude*magnitude;
            std::array<double,2*degree> powers{};
            powers[0]=1.;
            for(int k=1;k<2*degree;++k)powers[k]=powers[k-1]*t;
            double norm2=0.,numerator=0.;
            for(int n=0;n<=degree;++n)norm2+=v[n]*v[n]*powers[n];
            for(int i=0;i<degree;++i)for(int j=i+1;j<=degree;++j)
                numerator+=double((j-i)*(j-i))*v[i]*v[i]*v[j]*v[j]*powers[i+j-1];
            return numerator/(norm2*norm2);
        }
        const double scale=std::max(std::abs(z.real()),std::abs(z.imag()));
        const double inverse=(1./scale)/std::abs(z/scale);
        const double u=inverse*inverse;
        std::array<double,degree+1> weights{};
        double power=1.,total=0.;
        for(int k=0;k<=degree;++k) {
            weights[k]=v[degree-k]*v[degree-k]*power;
            total+=weights[k];power*=u;
        }
        for(double& weight:weights)weight/=total;
        double variance=0.;
        for(int i=0;i<degree;++i)for(int j=i+1;j<=degree;++j)
            variance+=double((j-i)*(j-i))*weights[i]*weights[j];
        return u*variance;
    }

    // Angle between the projective rays [v] and [w(z)], in [0,pi/2].
    double projective_distance(C z) const {
        return std::acos(std::sqrt(coherence(z)));
    }

    // A(z)=sum_{k>=1} c_k z^k/c_1; A(0)=0 and A'(0)=1.
    C correlation(C z,C* derivative=nullptr) const {
        C value(c[degree]/c[1],0.),d(0.,0.);
        for(int k=degree-1;k>=0;--k) {
            d=d*z+value;
            value=value*z+(k ? c[k]/c[1] : 0.);
        }
        if(derivative)*derivative=d;
        return value;
    }

    // F(x)=sum k*c_k*U_(k-1)(x)/sum k^2*c_k, so F(1)=1.
    // For real x in (-1,1), d_x |sum v_n exp(i*n*acos x)|^2
    // equals 2*angular_norm*F(x). The polynomial is continued to complex x
    // only to define the separate Newton visualization.
    C angular(C z,C* derivative=nullptr) const {
        C previous(0.,0.),current(1.,0.);
        C dprevious(0.,0.),dcurrent(0.,0.);
        C value(0.,0.),d(0.,0.);
        for(int k=1;k<=degree;++k) {
            value+=double(k)*c[k]*current;
            d+=double(k)*c[k]*dcurrent;
            const C next=2.*z*current-previous;
            const C dnext=2.*current+2.*z*dcurrent-dprevious;
            previous=current;current=next;
            dprevious=dcurrent;dcurrent=dnext;
        }
        if(derivative)*derivative=d/angular_norm;
        return value/angular_norm;
    }

    // Squared normalized overlap between v and (v_n*z^n)_n.
    // The outside-unit-circle branch removes |z|^degree before evaluation.
    double coherence(C z) const {
        if(!std::isfinite(z.real()) || !std::isfinite(z.imag()))return 0.;
        const double magnitude=std::abs(z);
        C numerator(0.,0.);
        double denominator=0.;
        if(magnitude<=1.) {
            numerator=C(v[degree]*v[degree],0.);
            denominator=v[degree]*v[degree];
            for(int n=degree-1;n>=0;--n) {
                numerator=numerator*z+v[n]*v[n];
                denominator=denominator*magnitude*magnitude+v[n]*v[n];
            }
        } else {
            // magnitude may overflow only when the complex modulus exceeds
            // DBL_MAX; first scale its two finite components in that case.
            const double scale=std::max(std::abs(z.real()),std::abs(z.imag()));
            const C reduced=z/scale;
            const double reduced_magnitude=std::abs(reduced);
            const C unit=reduced/reduced_magnitude;
            const double inverse_magnitude=(1./scale)/reduced_magnitude;
            C unit_power(1.,0.);
            numerator=C(v[0]*v[0],0.);
            denominator=v[0]*v[0];
            for(int n=1;n<=degree;++n) {
                unit_power*=unit;
                numerator=numerator*inverse_magnitude+v[n]*v[n]*unit_power;
                denominator=denominator*inverse_magnitude*inverse_magnitude+v[n]*v[n];
            }
        }
        const double result=std::norm(numerator)/(denominator*c[0]);
        return std::isfinite(result) ? std::clamp(result,0.,1.) : 0.;
    }
};

inline const Model& cached(double lambda) {
    const double radius=bounded_radius(lambda);
    thread_local Model model(radius);
    if(model.radius!=radius)model=Model(radius);
    return model;
}
} // namespace qqrh_hilbert
