#pragma once
#include "qqrh_hilbert_math.hpp"
#include <iostream>
#include <limits>

namespace qqrh_hilbert_tests {
inline bool run(std::ostream& out) {
    using namespace qqrh_hilbert;
    int checks=0,failures=0;
    auto check=[&](bool passed,const char* name) {
        ++checks;
        if(!passed){++failures;out<<"QQRH_HILBERT_CHECK=FAIL NAME="<<name<<'\n';}
    };
    auto close=[](C a,C b,double tolerance=2e-8) {
        return std::abs(a-b)<=tolerance*(1.+std::abs(b));
    };
    for(double radius:{.1,.35,.6,1.,1.25}) {
        Model model(radius);
        check(close(model.c[0],1.,1e-14),"unit_vector");
        C derivative;
        check(close(model.correlation(0.,&derivative),0.,1e-14) &&
              close(derivative,1.,1e-14),"normalized_correlation_at_origin");
        check(close(model.angular(1.),1.,1e-14),"angular_normalization");
        check(close(model.coherence(1.),1.,1e-14),"identical_vectors_overlap");
        check(close(model.coherence(0.),model.v[0]*model.v[0],1e-14),"origin_overlap");
        check(close(model.kernel(1.),1.,1e-14),"kernel_identical_vector_inner_product");
        check(close(model.projective_metric(0.),std::pow(model.v[1]/model.v[0],2),1e-14),
              "projective_metric_origin_limit");
        check(model.projective_distance(1.)<2e-8,"projective_distance_same_ray");
        for(C z:{C(.17,.23),C(-.8,.31),C(1.3,-.2)}) {
            const double h=2e-6;
            model.correlation(z,&derivative);
            check(close(derivative,(model.correlation(z+h)-model.correlation(z-h))/(2.*h)),
                  "correlation_derivative_finite_difference");
            model.angular(z,&derivative);
            check(close(derivative,(model.angular(z+h)-model.angular(z-h))/(2.*h)),
                  "angular_derivative_finite_difference");
            C numerator(0.,0.),power(1.,0.);
            double denominator=0.;
            for(int n=0;n<=degree;++n) {
                numerator+=model.v[n]*model.v[n]*power;
                denominator+=model.v[n]*model.v[n]*std::norm(power);
                power*=z;
            }
            const double independent=std::norm(numerator)/(denominator*model.c[0]);
            check(close(model.coherence(z),independent,2e-13),"overlap_direct_sum");
            check(close(model.kernel(z),numerator,2e-13),"kernel_direct_inner_product");
            model.kernel(z,&derivative);
            check(close(derivative,(model.kernel(z+h)-model.kernel(z-h))/(2.*h)),
                  "kernel_derivative_finite_difference");
            // Independently form the two-vector Gram determinant from w,w'.
            C derivative_inner(0.,0.);
            double derivative_norm=0.;
            C previous_power(1.,0.);
            for(int n=1;n<=degree;++n) {
                const C tangent=double(n)*model.v[n]*previous_power;
                const C component=model.v[n]*previous_power*z;
                derivative_norm+=std::norm(tangent);
                derivative_inner+=tangent*std::conj(component);
                previous_power*=z;
            }
            const double gram=(denominator*derivative_norm-std::norm(derivative_inner))
                              /(denominator*denominator);
            check(close(model.projective_metric(z),gram,2e-12),"projective_metric_Gram_determinant");
            check(close(model.projective_metric(z),model.projective_metric(std::abs(z)),2e-13),
                  "projective_metric_rotation_invariance");
            const double angle=model.projective_distance(z);
            check(close(std::pow(std::cos(angle),2),independent,2e-13),"projective_distance_overlap_identity");
        }
        for(C z:{C(0.,0.),C(.4,.5),C(-2.,.7),C(1e100,-2e100),C(1e308,1e308)}) {
            double overlap=model.coherence(z);
            check(std::isfinite(overlap) && overlap>=0. && overlap<=1.,"Cauchy_Schwarz_bounds");
            check(close(overlap,model.coherence(std::conj(z)),1e-13),"overlap_conjugation_symmetry");
            const double metric=model.projective_metric(z);
            check(std::isfinite(metric) && metric>=0.,"projective_metric_finite_nonnegative");
            const double angle=model.projective_distance(z);
            check(std::isfinite(angle) && angle>=0. && angle<=std::acos(-1.)/2.,
                  "projective_distance_bounds");
        }
        for(double tiny:{1e-8,1e-80,1e-200})
            check(close(model.projective_metric(tiny),model.projective_metric(0.),2e-12),
                  "projective_metric_tiny_argument_limit");
        check(close(model.projective_metric(1.-1e-8),model.projective_metric(1.+1e-8),2e-7),
              "projective_metric_branch_continuity");
        for(double x:{-.9,-.3,.2,.8}) {
            // Evaluate d_x |T_N(R exp(i theta))|^2 from polynomial jets,
            // independently of the correlation/Chebyshev implementation.
            const C phase(x,std::sqrt(1.-x*x)),z=radius*phase;
            C value(ff_rh::THETA_COEFFICIENTS[degree],0.),jet(0.,0.);
            for(int n=degree-1;n>=0;--n) {
                jet=jet*z+value;
                value=value*z+ff_rh::THETA_COEFFICIENTS[n];
            }
            const C dzdx=radius*C(1.,-x/std::sqrt(1.-x*x));
            const double derivative_raw=2.*std::real(jet*dzdx*std::conj(value));
            double norm2=0.,power=1.;
            for(int n=0;n<=degree;++n) {
                const double term=ff_rh::THETA_COEFFICIENTS[n]*power;
                norm2+=term*term;power*=radius;
            }
            const double expected_psi=derivative_raw/(2.*radius*norm2);
            const double actual_psi=model.angular_norm*model.angular(x).real()/radius;
            check(close(actual_psi,expected_psi,4e-11),"angular_identity_polynomial_jet");
        }
        for(int p=0;p<4;++p)for(int q=0;q<4;++q) {
            double gram=0.,shifted_gram=0.;
            for(int n=0;n+std::max(p,q)<=degree;++n)gram+=model.v[n+p]*model.v[n+q];
            for(int n=0;n+std::max(p,q)+1<=degree;++n)shifted_gram+=model.v[n+p+1]*model.v[n+q+1];
            check(close(gram-shifted_gram,model.v[p]*model.v[q],2e-14),"Gram_rank_one_defect");
        }
        for(int m:{1,4,degree+1}) {
            double lost=0.,tail=0.;
            for(int n=0;n<m;++n)lost+=model.v[n]*model.v[n];
            for(int n=m;n<=degree;++n)tail+=model.v[n]*model.v[n];
            check(close(lost+tail,1.,2e-14),"shift_energy_telescoping");
        }
    }
    check(Model(-3.).radius==.1 && Model(8.).radius==1.25,"radius_limits");
    check(std::isfinite(Model(std::numeric_limits<double>::quiet_NaN()).radius),"NaN_radius_is_finite");
    check(cached(.31).radius==.31 && cached(.87).radius==.87,"thread_cache_rekeys");
    {
        Model model(.65);
        C derivative;
        const C value=model.angular(-1.,&derivative);
        check(std::abs(value)<1e-4 && std::abs(value/derivative)>.09,
              "small_normalized_residual_is_not_root_convergence");
        const C seeds[]={C(-.5,0.),C(-.4,.3),C(.2,.6)};
        const C reference_roots[]={C(-1.0928407284,0.),C(-1.0538875,.1780159),
                                   C(-.8915953,.4349980)};
        for(int sample=0;sample<3;++sample) {
            C z=seeds[sample];
            bool converged=false;
            for(int iteration=0;iteration<64;++iteration) {
                const C f=model.angular(z,&derivative);
                if(!std::isfinite(std::abs(f)) || std::abs(derivative)<1e-13)break;
                const C correction=f/derivative;
                if(std::abs(f)<1e-6 && std::abs(correction)<1e-6*(1.+std::abs(z))) {
                    converged=true;break;
                }
                z-=correction;
            }
            check(converged && close(z,reference_roots[sample],2e-6),
                  "angular_Newton_reference_root_with_correction_check");
        }
    }
    out<<"QQRH_HILBERT_TEST="<<(failures ? "FAIL" : "PASS")
       <<" CHECKS="<<checks<<" FAILURES="<<failures<<" N="<<degree<<'\n';
    return failures==0;
}
} // namespace qqrh_hilbert_tests
