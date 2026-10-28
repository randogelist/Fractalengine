#pragma once
// FRAKTALFAKTUR_RH_LAB_GPU v1.8 additive mathematical candidates, C++17.
// Standalone extension API: no modification or assumption about the base renderer.
// Driver IDs are stable 18..31 (base v1.8 has 18 existing drivers, IDs 0..17).
// All systems here are exploratory complex maps, not proofs of RH or Collatz.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include "qqrh_hilbert_math.hpp"

namespace frk_extra {
using C = std::complex<double>;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr int first_driver = 18;
constexpr int driver_count = 14;
constexpr int last_driver = first_driver + driver_count - 1;
static constexpr std::array<const char*, driver_count> names = {{
  "19 Cayley commutator", "20 Tricorn Cayley", "21 GP echo memory",
  "22 Spectral Newton web", "23 Halley quintic web",
  "24 Entire Collatz parity", "25 Dirichlet eta Newton",
  "26 Schwarzian pole web", "27 QQRH Hilbert overlap",
  "28 QQRH Hardy shift fractal", "29 QQRH angular kernel Newton",
  "30 QQRH projective metric", "31 QQRH overlap kernel Newton",
  "32 QQRH kernel feedback fractal"
}};
static constexpr std::array<const char*,driver_count> algebra = {{
  "M(z)=A B A^-1 B^-1(z); z_next=M(z)^2+c",
  "z_next=K(conj(z))^2+c+0.14*lambda*K(conj(z))",
  "z_next=z^2+c+lambda*(1+0.22i)*(z-z_previous)",
  "z_next=z-1/sum_j 1/(z-r_j); r_j=-(gamma_j/gamma_1)^2",
  "z_next=z-2*P(z)*P'(z)/(2*P'(z)^2-P(z)*P''(z)), P=z^5-1+0.30*lambda*z^2",
  "z_next=z/2+(z+1/2)*sin(pi*z/2)^2+lambda*sin(pi*z)/8",
  "z_next=z-eta9(z)/eta9'(z), eta9=sum_n=1..9 (-1)^(n-1)n^(-z)",
  "z_next=z^2+c+lambda*S(z^3-3z), S=(-4z^2-2)/(z^2-1)^2",
  "coherence=|<v,v(z)>|^2/(||v||^2*||v(z)||^2); v_n=b_n*R^n, n=0..12",
  "c_k=<B^k v,v>; A(z)=sum_(k=1..12) c_k*z^k/c_1; z_next=z^2+c-.35*A/A'",
  "F(x)=sum_(k=1..12) k*c_k*U_(k-1)(x)/sum k^2*c_k; z_next=z-F/F'",
  "g(z)=partial_z partial_conj(z) log K(z,z); K=sum v_n^2*|z|^(2n); color=g/(1+g)",
  "H(z)=<w(z),v>=sum v_n^2*z^n; z_next=z-H/H'; kernel-orthogonality basins",
  "H(z)=sum v_n^2*z^n; z_next=z^2+c-.35*H/H'; finite Hilbert kernel feedback"
}};
inline bool finite(C a) { return std::isfinite(a.real()) && std::isfinite(a.imag()); }
inline C divide(C a, C b) {
  const double n=std::norm(b);
  if(!(n>1e-26) || !std::isfinite(n)) return {1e12,1e12};
  return a/b;
}
struct Mat { C a,b,c,d; };
inline Mat multiply(const Mat &x,const Mat &y) {
  return {x.a*y.a+x.b*y.c, x.a*y.b+x.b*y.d,
          x.c*y.a+x.d*y.c, x.c*y.b+x.d*y.d};
}
inline Mat inverse_projective(Mat m) { return {m.d,-m.b,-m.c,m.a}; }
inline Mat cayley_commutator(double lambda) {
  C a=C(.48,.13)*(1.+.20*lambda);
  C b=C(-.26,.43)*(1.+.15*lambda);
  Mat A{1.,a,std::conj(a),1.}, B{1.,b,std::conj(b),1.};
  Mat R=multiply(multiply(multiply(A,B),inverse_projective(A)),inverse_projective(B));
  double s=std::max({std::abs(R.a),std::abs(R.b),std::abs(R.c),std::abs(R.d),1e-12});
  R.a/=s;R.b/=s;R.c/=s;R.d/=s;
  return R;
}
inline C apply(Mat M,C z) { return divide(M.a*z+M.b,M.c*z+M.d); }
inline C cayley(C z,C a) { return divide(z-a,1.-std::conj(a)*z); }
inline C collatz(C z, double lambda=0.) {
  C s=std::sin((pi*.5)*z);
  return z*.5 + (z+.5)*s*s + (lambda*.125)*std::sin(pi*z);
}
inline C schwarzian(C z) {
  const C z2=z*z;
  return divide(-4.*z2-2.,(z2-1.)*(z2-1.));
}
inline C eta9(C z,C *deriv=nullptr) {
  C s(0.,0.), ds(0.,0.);
  for(int j=1;j<=9;++j){
    const double logj=std::log(double(j)), sign=(j&1)?1.:-1.;
    C t=sign*std::exp(-z*logj);
    s+=t; ds-=logj*t;
  }
  if(deriv) *deriv=ds;
  return s;
}
inline double spectral_root(int j) {
  constexpr double g[5]={14.134725141734694,21.022039638771555,
                 25.010857580145688,30.424876125859514,32.93506158773919};
  const double q=g[std::max(0,std::min(4,j))]/g[0];
  return -q*q;
}
inline C spectral_newton(C z) {
  C ld(0.,0.);
  for(int j=0;j<5;++j) ld+=divide(C(1.,0.),z-C(spectral_root(j),0.));
  return z-divide(C(1.,0.),ld);
}
inline C halley(C z,double lambda) {
  const C z2=z*z,z3=z2*z,z4=z2*z2;
  const C P=z4*z-1.+(.30*lambda)*z2;
  const C d=5.*z4+(.60*lambda)*z;
  const C dd=20.*z3+(.60*lambda);
  return z-divide(2.*P*d,2.*d*d-P*dd);
}
struct State {
  C z, previous, c;
  int n=0;
};
struct Params {
  double lambda=.3;
  C julia_c={-.16,.19};
  bool julia=false;
  int iterations=160;
  double bailout=64.;
};
inline bool is_field_driver(int driver){return driver==26 || driver==29;}
inline bool seed_at_pixel(int driver){return (driver>=21 && driver<=24) || is_field_driver(driver) || driver==28 || driver==30;}
inline State start(int driver,C pixel,const Params &p) {
  const bool seed=seed_at_pixel(driver);
  C z=(seed || p.julia) ? pixel : C(0.,0.);
  C c=seed ? C(0.,0.) : (p.julia ? p.julia_c : pixel);
  return {z,C(0.,0.),c,0};
}
inline State step(int driver,State s,double lambda) {
  const C z=s.z;
  switch(driver) {
    case 18: { C t=apply(cayley_commutator(lambda),z);s.z=t*t+s.c; break; }
    case 19: { C a=C(.28,.22)*(1.+.17*lambda);
               C t=cayley(std::conj(z),a);
               s.z=t*t+s.c+(.14*lambda)*t; break; }
    case 20: {s.z=z*z+s.c+C(lambda,.22*lambda)*(z-s.previous);
              s.previous=z;break;}
    case 21: s.z=spectral_newton(z);break;
    case 22: s.z=halley(z,lambda);break;
    case 23: s.z=collatz(z,lambda);break;
    case 24: {C d;C f=eta9(z,&d);s.z=z-divide(f,d);break;}
    case 25: s.z=z*z+s.c+lambda*schwarzian(z);break;
    case 26: break; // Direct field: evaluate() has no orbit.
    case 27: {C d;const C a=qqrh_hilbert::cached(lambda).correlation(z,&d);
              s.z=z*z+s.c-.35*divide(a,d);break;}
    case 28: {C d;const C f=qqrh_hilbert::cached(lambda).angular(z,&d);
              s.z=z-divide(f,d);break;}
    case 29: break;
    case 30: {C d;const C h=qqrh_hilbert::cached(lambda).kernel(z,&d);
              s.z=z-divide(h,d);break;}
    case 31: {C d;const C h=qqrh_hilbert::cached(lambda).kernel(z,&d);
              s.z=z*z+s.c-.35*divide(h,d);break;}
    default: s.z=z*z+s.c;break;
  }
  ++s.n;
  return s;
}
struct Result {
  int iterations=0, basin=-1;
  bool escaped=false,converged=false;
  C final_z=0.;
  double trap=1e12;
  bool field=false;
  double value=0.;
};
inline bool is_root_driver(int driver){return driver==21 || driver==22 || driver==24 || driver==28 || driver==30;}
inline Result evaluate(int driver,C pixel,const Params &p) {
  Result o;
  if(is_field_driver(driver)) {
    o.field=true;
    if(driver==29) {const double g=qqrh_hilbert::cached(p.lambda).projective_metric(pixel);o.value=g/(1.+g);}
    else o.value=qqrh_hilbert::cached(p.lambda).coherence(pixel);
    o.final_z=pixel;return o;
  }
  State s=start(driver,pixel,p);
  const int limit=std::max(1,std::min(4096,p.iterations));
  const double bail=(driver==23)?72.:(is_root_driver(driver)?2e4:p.bailout);
  for(int n=0;n<limit;++n){
    if(!finite(s.z) || std::norm(s.z)>bail*bail ||
       (driver==23 && std::abs(s.z.imag())>16.) ||
       (driver==24 && (s.z.real() < -24. || s.z.real() > 28.))){
      o.iterations=n;o.escaped=true;break;
    }
    o.trap=std::min(o.trap,std::abs(s.z));
    if(driver==21){
      int closest=-1;double best=1e90;
      for(int j=0;j<5;++j){
        double d=std::abs(s.z-C(spectral_root(j),0.));
        if(d<best){closest=j;best=d;}
      }
      if(best<1e-6){o.iterations=n;o.basin=closest;o.converged=true;break;}
    } else if(driver==22){
      C f=std::pow(s.z,5)-1.+(.30*p.lambda)*s.z*s.z;
      if(std::abs(f)<1e-6){o.iterations=n;o.converged=true;
        o.basin=int(std::floor((std::arg(s.z)+pi)/(2*pi)*5.));break;}
    } else if(driver==24){
      C f=eta9(s.z);
      if(std::abs(f)<1e-6){o.iterations=n;o.converged=true;
        o.basin=int(std::floor((std::arg(s.z)+pi)/(2*pi)*8.));break;}
    } else if(driver==28 || driver==30){
      C d;const auto& model=qqrh_hilbert::cached(p.lambda);
      const C f=driver==28 ? model.angular(s.z,&d) : model.kernel(s.z,&d);
      if(finite(f) && finite(d) && std::norm(d)>1e-26 && std::abs(f)<1e-6 &&
         std::abs(f/d)<1e-6*(1.+std::abs(s.z))){o.iterations=n;o.converged=true;
        o.basin=driver==30 ? -2 : int(std::floor((std::arg(s.z)+pi)/(2*pi)*12.));break;}
    }
    s=step(driver,s,p.lambda);
    o.iterations=n+1;
  }
  o.final_z=s.z;
  return o;
}
inline uint32_t color(Result o, int palette=0) {
  if(o.field) {
    const double q=std::clamp(o.value,0.,1.);
    double red=.035+.90*q,green=.08+.65*std::sqrt(q),blue=.16+.75*(1.-q);
    const double stripe=.78+.22*std::cos(12.*pi*q);
    red*=stripe;green*=stripe;blue*=stripe;
    if(palette==1){const double t=red;red=green;green=blue;blue=t;}
    if(palette==2)std::swap(red,blue);
    auto byte=[](double t){return uint32_t(std::lround(255.*std::clamp(t,0.,1.)));};
    return 0xff000000u|(byte(red)<<16)|(byte(green)<<8)|byte(blue);
  }
  if(!o.converged&&!o.escaped) return 0xff061321;
  const double phase=finite(o.final_z)?std::arg(o.final_z):0.;
  const double mag=finite(o.final_z)?std::abs(o.final_z):0.;
  double t=o.converged?(.173*o.basin+.009*o.iterations+.05*phase):
    (.023*o.iterations+.16*phase/(2*pi)+.012*std::log1p(mag));
  // Continuous endpoint colour for overlap roots: avoid an artificial branch
  // cut on the negative real axis where several of these roots lie.
  if(o.converged && o.basin==-2)t=.47*std::log1p(mag)+.10*std::sin(phase)+.006*o.iterations;
  if(palette==1)t+=.19;
  if(palette==2)t=-t;
  const double light=o.converged ?(.69+.25*std::exp(-o.iterations/25.)):
                   (.42+.52*(1.-std::exp(-o.iterations/30.)));
  auto channel=[&](double off){
    double v=light*(.54+.46*std::cos(2*pi*t+off));
    return uint32_t(std::lround(255*std::max(0.,std::min(1.,v))));
  };
  return 0xff000000u | (channel(.1)<<16) |(channel(2.1)<<8)|channel(4.2);
}
} // namespace frk_extra
