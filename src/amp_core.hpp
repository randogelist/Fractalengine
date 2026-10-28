#pragma once
#include "garden_math.hpp"
#include "grh_observable_math.hpp"
#include "spectral_math.hpp"
#include "quantum_math.hpp"
#include "spherical_math.hpp"
#include "spin_math.hpp"
#include "project_math.hpp"
#include "riemann_witness_math.hpp"
#include "moment_cubic_math.hpp"
#include "moment_quotient_math.hpp"
#include "moment_triple_math.hpp"
#include "rh_master_math.hpp"
#include "bkqr_lambda_math.hpp"
#include "qqrh_direct_math.hpp"
#include "frk_extra_catalogue.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace amp {
using C = std::complex<double>;
constexpr double pi = 3.14159265358979323846;
struct Q {
    C x{0, 0}, y{1, 0};
};
inline Q normal(Q q) {
    double s = std::sqrt(std::norm(q.x) + std::norm(q.y));
    if (!(s > 0) || !std::isfinite(s))
        throw std::runtime_error("Invalid homogeneous point");
    return {q.x / s, q.y / s};
}
inline Q point(C z) { return normal({z, C(1, 0)}); }
inline Q infinity() { return {C(1, 0), C(0, 0)}; }
inline double chordal(Q a, Q b) { return std::abs(a.x * b.y - b.x * a.y); }
inline bool affine(Q q, C &z) {
    if (std::norm(q.y) < 1e-28)
        return false;
    z = q.x / q.y;
    return true;
}
inline Q advance(Q q, C c, C alpha, double *ly = nullptr) {
    if (std::norm(alpha) == 0) {
        if (ly)
            *ly = std::log(1e-24);
        return point(c);
    }
    C xx = q.x * q.x, xy = q.x * q.y, yy = q.y * q.y;
    C d = xx - 6.0 * xy + 8.0 * yy;
    C n = 0.6 * alpha * (8.0 * xy - 3.0 * xx) + c * d;
    double den = std::norm(n) + std::norm(d);
    if (ly)
        *ly = std::log(std::max(1e-24, 0.6 * std::abs(alpha) * std::abs(10.0 * xx - 48.0 * xy + 64.0 * yy) /
                                           std::max(den, 1e-30)));
    return normal({n, d});
}
struct OrbitInfo {
    int period = 0;
    double stretch = 0;
    Q last;
    std::vector<Q> cycle;
};
inline OrbitInfo analyse(C c, C alpha, Q seed, int budget = 768) {
    Q q = seed;
    std::vector<Q> tail;
    double sum = 0;
    for (int i = 0; i < budget; i++) {
        double ly;
        q = advance(q, c, alpha, &ly);
        if (i >= budget - 128)
            sum += ly;
        if (i >= budget - 80)
            tail.push_back(q);
    }
    OrbitInfo a;
    a.last = q;
    a.stretch = sum / 128;
    for (int p = 1; p <= 16; p++) {
        bool good = true;
        for (int j = 0; j < std::max(4, 2 * p); j++)
            if (chordal(tail[tail.size() - 1 - j], tail[tail.size() - 1 - j - p]) > 1e-9) {
                good = false;
                break;
            }
        if (good) {
            a.period = p;
            a.cycle.assign(tail.end() - p, tail.end());
            break;
        }
    }
    return a;
}
struct Camera {
    double x = 2.1, y = -1.3, span = 8.2;
};
struct State {
    int family = 0, generation = 1, method = 0;
    C c{1.2, .36}, alpha{1, 0}, seed{2.4, .8};
    int plane = 0, orbitMode = 0, palette = 0, iterations = 256, preset = 0;
    double hue = 0, bands = 1, spectralD = .25, spectralWeight = .35;
    int spectralOrder = 1;
    int bkqrQuality = 0, bkqrPrecision = 0;
    double qqrhSlope = 11;
    double qqrhScale = 2;
    double rhLevel = 0;
    rhmaster::Experiment rhExperiment{};
    quantum::Settings quantum;
    spherical::Settings sphere;
    spin::Settings spin;
    projects::Settings project;
    riemann::Settings riemann;
    std::array<Camera, 3> cameras{{{2.1, -1.3, 8.2}, {1.3, .3, 5.2}, {1, 0, 5.2}}};
};
inline bkqrlambda::Settings bkqrSettings(const State &s) {
    bkqrlambda::Settings q;
    if(s.family==22 || s.family==23) {
        q.c=s.c;q.host=bkqrlambda::JULIA;q.power=2;
        q.branch=s.orbitMode;q.coupling=s.spectralOrder-1;
        q.quality=s.bkqrQuality;q.precision=s.bkqrPrecision;
        q.iterations=s.iterations;q.metric=s.palette;
        if(s.family==23) {q.qqrhChannel=s.generation;q.qqrhSlope=s.qqrhSlope;}
        else q.rhChannel=s.generation;
        q.rhN=8 << std::clamp(s.method,0,2);
        q.rhX=s.seed.real();q.rhB=s.seed.imag();q.rhH=s.spectralD;q.eta=s.alpha;
        return q;
    }
    q.c=s.c;q.host=s.method;q.branch=s.orbitMode;q.coupling=s.spectralOrder-1;
    q.power=int(std::lround(s.alpha.real()));q.quality=int(std::lround(s.seed.real()));
    q.precision=int(std::lround(s.seed.imag()));q.iterations=s.iterations;q.metric=s.palette;
    return q;
}
inline const char* qqrhTitle(int generation) {
    static const char* names[]={"QQRH", "QQRH / Curvature", "QQRH / Dyadic",
        "QQRH / Laplacian 1", "QQRH / Laplacian 2", "QQRH / Moving boundary",
        "QQRH / Scale-ratio curvature"};
    return names[std::clamp(generation,0,6)];
}
inline std::string qqrhCarrierFormula(int generation) {
    std::string result="T(y)=Theta_N(200*y)/Theta(0); y=b*z+h; ";
    if(generation==1)
        result+="J(y)=T(-y); F(y)=-d_y(y*J'(y)/J(y)); Q=c+eta*(F(y)-F(h)); ";
    else {
        result+="rho+/-=x+/-i*sqrt(1-x^2); P_G(y,x)=G(y*rho+)*G(y*rho-); ";
        if(generation==2)
            result+="F(y)=P_T(2*y,x)/P_T(y,x); Q=c+eta*(F(y)-F(h)); ";
        else {
            result+="Psi_G(y,x)=d_x(P_G(y,x))/(2*y), with polynomial extension; ";
            if(generation==3) result+="F(y)=Psi_(T')(y,x); ";
            else if(generation==4) result+="F(y)=Psi_(T'')(y,x); ";
            else result+="a=(1-m^2)/(1+m^2); F(y)=(Psi_T(y,x)-Psi_T(y,a))/(y*(x-a)), with polynomial extension; ";
            result+="Q=z^2+c-eta*F(y)/(b*F'(y)); ";
        }
    }
    return result;
}
inline std::string qqrhDirectFormula(int generation) {
    if(generation==6)
        return "T(y)=Theta_N(200*y)/Theta(0); J(y)=T(-y); C_J(y)=-d_y(y*J'(y)/J(y)); K_kappa(y)=kappa*C_J(kappa*y)-C_J(y); y=b*z+h; z_next=c+eta*(K_kappa(y)-K_kappa(h)); kappa=1 gives the constant host c; parameter portrait uses z0=0";
    return qqrhCarrierFormula(generation)+"z_next=Q; Julia initial state is the pixel; parameter portrait uses z0=0";
}
inline std::string qqrhFormula(int generation) {
    return qqrhCarrierFormula(generation)+"z_next=z+(1+sigma*H(lambda))*(Q-z); inner eta is independent of outer pixel lambda";
}
inline std::string familyName(int family) {
    return family == 1   ? "GP / BKQR"
           : family == 2 ? "Percolation"
           : family == 4 ? "Restart / spectral"
           : family == 5 ? "Cayley / RH kernel"
           : family == 6 ? "RH / Dirichlet L"
           : family == 7 ? "Primitive GRH Cayley"
           : family == 8 ? "Primitive Character Pair"
           : family == 9 ? "GRH Cayley Observable"
           : family == 10 ? "GRH Phase Feedback"
           : family == 12 ? "Quantum Measurement Histories"
           : family == 13 ? "Spherical Fractals"
           : family == 14 ? "Dissipative Spin Chaos"
           : family == 15 ? "Five Project Fractals"
           : family == 25 ? "FRAKTALFAKTUR / new complex fractals"
           : family == 20 ? "RH / finite-model candidate drivers"
           : family == 21 ? "BKQR / original complex-lambda phase map"
           : family == 22 ? "BKQR lambda / RH Master / Sub_1 / Sub_2"
           : family == 23 ? "BKQR lambda / QQRH five carriers"
           : family == 24 ? "QQRH / direct Julia and parameter portraits"
           : family == 19 ? "Riemann Moment Dynamics / Three Ways"
           : family == 18 ? "Riemann Moment Quotient Dynamics"
           : family == 17 ? "Riemann Moment Julia Dynamics"
           : family == 16 ? "Riemann Witness Fractals"
           : family == 11 ? "Spectral Jet Determinant"
           : family == 3 ? "CROWN"
                         : "AMP";
}
inline std::string lfunctionName(int g) {
    static const char* names[] = {"", "Riemann zeta", "L / chi modulo 3", "Dirichlet beta / modulo 4", "L / quartic chi modulo 5"};
    return g >= 1 && g <= 4 ? names[g] : "Dirichlet L";
}
inline std::string formula(const State &s) {
    if(s.family==25) {
        if(s.generation>=9) return frk_extra::algebra[std::clamp(s.generation,1,frk_extra::driver_count)-1];
        static const char* expressions[]={
          "z_next=M_commutator(z)^2+c, where M=[A,B], A and B are Cayley/Mobius maps",
          "z_next=Cayley(conj(z),a)^2+c+.14*lambda*Cayley(conj(z),a)",
          "z_next=z^2+c+(lambda+.22*i*lambda)*(z-z_previous) (four-real-dimensional state)",
          "P(z)=product_{j=1..5}(z+gamma_j^2/gamma_1^2), z_next=z-P/P'",
          "P(z)=z^5-1+.30*lambda*z^2, z_next=z-2*P*P'/(2*P'^2-P*P'')",
          "f(z)=z/2+(z+1/2)*sin(pi*z/2)^2+lambda*sin(pi*z)/8; f(n)=Collatz(n) for integer n, lambda=0",
          "eta_9(z)=sum_{n=1..9}(-1)^(n-1)*n^(-z), z_next=z-eta_9(z)/eta_9'(z)",
          "z_next=z^2+c+lambda*S(z), S(z)=(-4*z^2-2)/(z^2-1)^2"
        };
        return expressions[std::clamp(s.generation,1,8)-1];
    }
    if(s.family==24) return qqrhDirectFormula(s.generation);
    if(s.family==23) return qqrhFormula(s.generation);
    if(s.family==22)
        return "Q=z^2+c+eta*RH_response_j(z); z_next=z+(1+sigma*H(lambda))*(Q-z); inner eta is independent of outer pixel lambda";
    if (s.family==21)
        return "z_next=z+(1+sigma*h(lambda))*(F(z)-z); lambda=pixel; original BKQR v0.4.7 host/branch/coupling ensemble and global percentile metrics";
    if (s.family==20) return rhmaster::formula(s.generation);
    if (s.family==19) return momenttriple::formula(s.generation);
    if (s.family==18)
        return "z_next=lambda*(z-b-1/z)+c; lambda=1/sqrt(h), 0<h<1; critical starts +/-i; pole 0 -> infinity";
    if (s.family==17)
        return "z_next=z^3-A*z+c; A=1+h+b^2/3 from five fixed RH moments; critical starts +/-sqrt(A/3)";
    if (s.family == 16) {
        static const char *carriers[]={"",
            "P(z)=product_j(z-q_j/q_max)",
            "P(z)=product_j(kappa(z)-u_j); u_j=q_j/(1+sqrt(1-q_j))^2",
            "P(z)=product_j(q_j*(1+kappa(z))^2-4*kappa(z))",
            "P(z)=(z^2-b*z-1)^2+h*z^2",
            "P(z)=z^4-(2+b^2+h)*z^2+1+h",
            "P(z)=z*(z^2-b*z-1-h); Jacobi terminal diagonal=0",
            "P(z)=z^4*A((z+1/z)/2); A(t)=(t^2-b*t-1)^2+h*t^2; polynomial continuation at z=0"};
        const int g=std::clamp(s.generation,1,7);
        return std::string(carriers[g])
            +(g<=3 ? "; kappa(z)=27*z*(1-z)^2/4; "
                    : "; b=E[t^3], h=E[t^4]-b^2-1, t=(q-mu)/sigma; ")
            +(s.method==0 ? "z_next=z-alpha*P/P'" : "z_next=z-alpha*2*P*P'/(2*P'^2-P*P'')")
            +"; finite rounded actual-xi moments";
    }

    if (s.family == 15) {
        switch (s.generation) {
        case 1: return "w_next=a+b*H_M'(w)/H_M(w); H=completedPair(1/2+w)*completedPair(1/2-w); finite even Taylor carrier";
        case 2: return "D(z)=sum_t [beta_(2*(1024+offset+t))-beta_(2*(1024+offset+t)+1)-1024*c(1024+offset+t)]*z^t; z_next=z-D/D'";
        case 3: return "P_t(z)=sum_j (j+1)^(-slope)*exp(-2*nu_t*4^j)*z^j; z_next=z-P_t/P_t'; exact periodic NSE shear energies";
        case 4: return "F(z)=z+1/4-(z/2+1/4)*cos(pi*z); z_next=F(F(z)); exact integer shortcut blocks, entire extension";
        default: return "Q(z)=[p/(2-p)]*z^5-z^4-z^2-3*z-15; z=2*d-1; z_next=z-Q/Q'; q=2 GP-4 critical-dimension basins";
        }
    }
    if (s.family == 14)
        return "q=R_y(beta)*r; u=R_z(k*q.z)*q; r_next=(sech(tau)*u.x,sech(tau)*u.y,u.z-tanh(tau))/(1-u.z*tanh(tau)); finite-time tangent diagnostics";
    if (s.family == 13)
        return s.generation == 3
            ? "f(z)=z^3+lambda/z^3; six equally weighted inverse branches on CP1; lambda=10^logMagnitude*exp(i*phase)"
            : "I_n,a(r)=((1-a^2)*r+2*(a-n.r)*n)/(1+a^2-2*a*n.r); tetrahedral circle reflections; no immediate repeated generator";
    if (s.family == 12)
        return "r=R_kick*r; p_i=(1+e^2+2e*n_i.r)/(N*(1+e^2)); r_next=((1-e^2)*r+2e*(1+e*n_i.r)*n_i)/(1+e^2+2e*n_i.r); finite history density";
    if (s.family == 11 && s.generation >= 2)
        return "x_next=x^2+c+eta*(R_s*R_(s+2)/R_(s+1)^2-1); R_k=sum_i w_i/(x-t_i)^k; t=d*(-1,.35,1), w=((1-m)/2,m,(1-m)/2); q=1/2+t";
    if (s.family == 11)
        return "x_next=x^2+c+eta*d^2*(x^2-d^2)/(x^2+d^2)^2; x=z-1/2, q_+-=1/2+-d; finite numerical portrait";
    if (s.family == 5)
        return "P'=4PQ; Q'=(P+Q)^2; P0=z, Q0=1; f=P[depth]-target*Q[depth]";
    if (s.family == 6)
        return lfunctionName(s.generation) + ": f(s)=L(s,chi)-target; analytic Newton / Halley";
    if (s.family == 7)
        return "s_k(w)=1/2+(atanh(w)+i*pi*k)/log(2); f=completed L(s_k,chi_5), common completion factor removed";
    if (s.family == 8)
        return "P(s)=Lambda(s,chi_5)*Lambda(s,chibar_5); common A(s)^2 factor removed from Newton/Halley jets";
    if (s.family == 9)
        return "C_chi(a,r)=4a[z_-R(z_-)-z_+R(z_+)]/(z_--z_+), R=Theta'/Theta, Theta(z)=completedPair(1/2+sqrt(z))";
    if (s.family == 10)
        return "w_next=w-rho*N_f(w)+kappa*Phi(f;tau)/(1+0.4|w|^2), f=completed primitive L on Cayley sheet 3";
    if (s.family == 1)
        return "P[n+1] = 1 + a*z*P[n](z^2) + b*z^4*P[n](z^4); G = (1-z)*P[n]";
    if (s.family == 4)
        return "D(z) = z^L - feedback * sum[h_l/M * z^(L-l)]; spectral factor of Q + feedback*u*b^T";
    if (s.family == 3)
        return "G(z) = product over valuation words of (z - binary_embedding(q))";
    if (s.family == 2)
        return "G(p) = C[L](p) - target; exact finite strip connectivity";
    return "R(z) = c + (3/5)*alpha*z*(8-3*z)/((2-z)*(4-z))";
}
inline State defaults(int family, int generation = 1) {
    State s;
    s.family = family;
    s.generation = generation;
    if(family==25) {
        const auto& q=frk_extra::extra_catalogue[std::clamp(generation,1,frk_extra::driver_count)-1];
        s.method=0;s.plane=1;s.palette=0;s.iterations=q.suggested_iterations;
        s.c={-.16,.19};s.alpha={q.default_lambda,0};s.seed={0,0};
        // Original camera spans are horizontal; the addon presets use vertical spans.
        const Camera camera{q.center_real,q.center_imag,q.vertical_span * (16.0/9.0)};
        s.cameras={camera,camera,camera};
        return s;
    }
    if (family == 23 || family == 24) {
        static const C constants[]={C(-.25,.2),C(-.30,.15),C(-.56,.23),C(-.68,.27),C(-.56,.23),C(.7,.2)};
        static const double strengths[]={.15,.025,.45,.45,.45,.08};
        const int index=std::clamp(generation,1,family==24 ? 6 : 5)-1;
        s.method=1;s.orbitMode=0;s.spectralOrder=1;s.palette=0;s.plane=family==24 ? 1 : 2;s.iterations=family==24 ? 160 : 96;
        s.c=constants[index];s.alpha={strengths[index],0};s.seed={0,.9};s.spectralD=.28;
        s.bkqrQuality=0;s.bkqrPrecision=0;s.qqrhSlope=generation==5 ? 2.0 : 11.0;
        s.cameras={{{0,0,4},{0,0,4},{0,0,4}}};
        if(family==24 && generation==1)s.cameras[1]={1.3,0,1.6};
        if(family==24 && generation==6)s.cameras[1]={.7,0,1.8};
        return s;
    }
    if (family == 22) {
        static const C constants[]={C(-.56,.23),C(-.68,.27),C(-.30,.15)};
        s.method=1;s.orbitMode=0;s.spectralOrder=1;s.palette=0;s.plane=2;s.iterations=96;
        s.c=constants[std::clamp(generation,1,3)-1];s.alpha={.45,0};s.seed={0,.9};s.spectralD=.28;
        s.bkqrQuality=0;s.bkqrPrecision=0;
        s.cameras={{{0,0,4},{0,0,4},{0,0,4}}};
        return s;
    }
    if (family == 21) {
        s.method=0;s.orbitMode=0;s.spectralOrder=1;s.palette=0;s.plane=2;s.iterations=96;
        s.c={-.8,.156};s.alpha={2,0};s.seed={0,0};
        s.cameras={{{0,0,4},{0,0,4},{0,0,4}}};
        return s;
    }
    if (family == 20) {
        const auto q=rhfresh::defaultSettings(std::clamp(generation,1,rhfresh::CandidateCount)-1);
        s.method=0;s.plane=1;s.palette=0;s.iterations=q.iterations;
        s.c=q.parameter;s.alpha=q.lambda;s.seed={q.radius,q.scale};
        s.spectralD=q.shift;s.spectralOrder=q.degree==8 ? 1 : q.degree==32 ? 3 : 2;
        s.rhLevel=q.level;
        s.rhExperiment=rhmaster::Experiment{};
        s.cameras={{{0,0,4},{0,0,4},{0,0,4}}};
        return s;
    }
    if (family == 19) {
        s.method=0;s.plane=0;s.palette=0;s.c={0,0};s.alpha={1,0};
        s.seed={.2,.3};s.iterations=256;
        s.cameras={{{0,0,6},{0,0,5},{0,0,5}}};
        return s;
    }
    if (family == 18) {
        s.method=0; s.plane=0; s.palette=0; s.c={0,0}; s.alpha={1,0};
        s.seed={.2,.3}; s.iterations=256;
        s.cameras={{{1,0,8},{0,0,6},{0,0,6}}};
        return s;
    }
    if (family == 17) {
        s.method=0; s.plane=0; s.palette=0; s.c={0,0}; s.alpha={1,0};
        s.seed={.2,.3}; s.iterations=256;
        s.cameras={{{0,0,4},{0,0,4},{0,0,4}}};
        return s;
    }
    if (family == 16) {
        s.method=0; s.plane=0; s.palette=0;
        s.c={0,0}; s.alpha={1,0}; s.seed={.4,.25}; s.iterations=224;
        Camera v=generation==1 ? Camera{.5,0,2.6}
            : generation==2 ? Camera{.65,0,1.8} : generation==3 ? Camera{3,0,58}
            : generation==4 ? Camera{.25,0,5} : generation==5 ? Camera{0,0,5.2}
            : generation==6 ? Camera{.3,0,5.2} : Camera{0,0,9};
        s.cameras={v,v,v};
        return s;
    }
    if (family == 15) {
        s.method=0; s.plane=0; s.palette=0;
        s.c={.2,.15}; s.alpha={6,.8}; s.seed={.3,.2};
        s.iterations=(generation==1 || generation==4)?192:160;
        Camera v=generation==1 ? Camera{0,0,28} : generation==2 ? Camera{0,0,3.5}
            : generation==3 ? Camera{0,0,6} : generation==4 ? Camera{1.5,0,6} : Camera{1,0,14};
        s.cameras={v,v,v};
        return s;
    }
    if (family == 1) {
        s.c = {1, 0};
        s.alpha = {1, 0};
        s.seed = {.35, .31};
        s.iterations = 160;
        s.cameras = {{{0, 0, 5.4}, {0, 0, 5.4}, {0, 0, 5.4}}};
    }
    if (family == 2) {
        s.c = {.5, 0};
        s.alpha = {1, 0};
        s.seed = {.37, .23};
        s.iterations = 128;
        s.cameras = {{{.5, 0, 4.8}, {.5, 0, 4.8}, {1, 0, 4.8}}};
    }
    if (family == 3) {
        s.c = {.45, .35};
        s.alpha = {1, 0};
        s.seed = {.17, .23};
        s.iterations = 128;
        s.cameras = {{{0, 0, 5.2}, {0, 0, 2.4}, {1, 0, 4.8}}};
    }
    if (family == 4) {
        s.c = {.85, .35};
        s.alpha = {1, 0};
        s.seed = {.27, .19};
        s.iterations = 128;
        s.cameras = {{{0, 0, 3.6}, {0, 0, 4.8}, {1, 0, 4.8}}};
    }
    if (family == 5) {
        s.c = {.55, .22};
        s.alpha = {1, 0};
        s.seed = {.37, .26};
        s.iterations = 128;
        s.cameras = {{{0, 0, 5.4}, {.5, 0, 3.6}, {1, 0, 4.8}}};
    }
    if (family == 6) {
        static const double heights[] = {0, 18, 11, 9, 10};
        s.c = {0, 0};
        s.alpha = {1, 0};
        s.seed = {.7, heights[generation]};
        s.iterations = 64;
        s.palette = 1;
        s.cameras = {{{.5, heights[generation], 15}, {0, 0, 4.8}, {1, 0, 3.6}}};
    }
    if (family == 7) {
        s.c = {0, 0};
        s.alpha = {1, 0};
        s.seed = {.22, .35};
        s.iterations = 64;
        s.palette = 0;
        s.cameras = {{{0, .20, 4.4}, {0, 0, 4.8}, {1, 0, 3.6}}};
    }
    if (family == 8) {
        s.c = {0, 0};
        s.alpha = {1, 0};
        s.seed = {.5, 2.0};
        s.iterations = 64;
        s.palette = 0;
        s.cameras = {{{.5, 1.5, 14.0}, {0, 0, 4.8}, {1, 0, 3.6}}};
    }
    if (family == 9) {
        s.c = {0, 0};
        s.alpha = {1, 0};
        s.seed = {0, 0};
        s.iterations = 32;
        s.palette = 0;
        s.cameras = {{{6.25, 4.1, 5.8}, {6.25, 4.1, 5.8}, {6.25, 4.1, 5.8}}};
    }
    if (family == 10) {
        s.c = {.10, .80};       // kick kappa, phase twist tau
        s.alpha = {.90, 0};     // complex Newton/Halley relaxation rho
        s.seed = {.20, .20};
        s.iterations = 64;
        s.palette = 0;
        s.bands = 1.15;
        s.cameras = {{{0, 0, 3.6}, {0, 0, 3.6}, {0, 0, 3.6}}};
    }
    if (family == 11) {
        s.c = {-.30, .55};
        s.alpha = {.05, 0}; // eta is the spectral feedback strength.
        s.seed = {0, 0};
        s.spectralD = .25;
        s.iterations = 256;
        s.palette = 0;
        s.cameras = {{{0, 0, 6.0}, {0, 0, 4.8}, {0, 0, 3.6}}};
        if (generation >= 2) {
            s.alpha = {.08, .025};
            s.spectralD = .32;
            s.spectralWeight = .35;
            s.spectralOrder = 1;
        }
    }
    if (family == 12 || family == 13 || family == 14) {
        s.c = {0, 0}; s.alpha = {1, 0}; s.seed = {0, 0};
        s.quantum = quantum::Settings{};
        s.plane = 0; s.method = 1; s.palette = 0;
        s.cameras = {{{0, 0, 2.4}, {0, 0, 2.4}, {0, 0, 2.4}}};
        if (family == 13 && generation == 2) s.sphere.gap = .015;
        if (family == 14 && generation == 2) s.palette = 1;
    }
    return s;
}
struct Scene {
    State state;
    bkqrlambda::Settings bkqr;
    std::shared_ptr<const garden::Data> gardenData;
    std::shared_ptr<const projects::Model> projectModel;
    std::shared_ptr<const riemann::Model> riemannModel;
    std::shared_ptr<const momentcubic::Model> cubicModel;
    std::shared_ptr<const momentquotient::Model> quotientModel;
    std::shared_ptr<const momenttriple::Model> tripleModel;
    std::shared_ptr<const rhmaster::Model> masterModel;
    std::shared_ptr<const qqrhdirect::Model> qqrhDirectModel;
    std::vector<C> crownRoots;
    OrbitInfo a, b;
    bool same = false;
};
inline Scene scene(const State &s) {
    Scene r;
    r.state = s;
    if(s.family==25) {
        auto data=std::make_shared<garden::Data>();
        data->family=25;data->generation=s.generation;data->degree=2;
        data->records.push_back({0,0,0,0});r.gardenData=data;
        return r;
    }
    if (s.family==24) {
        const int degree=8 << std::clamp(s.method,0,2);
        r.qqrhDirectModel=qqrhdirect::model(s.generation,degree,s.seed.real(),s.seed.imag(),s.spectralD,s.alpha,s.qqrhSlope,s.qqrhScale);
        auto data=std::make_shared<garden::Data>();
        data->family=24;data->generation=s.generation;data->degree=r.qqrhDirectModel->degree;
        data->records.push_back({0,0,0,0});r.gardenData=data;return r;
    }
    if (s.family==21 || s.family==22 || s.family==23) {
        r.bkqr=bkqrSettings(s);
        auto data=std::make_shared<garden::Data>();
        data->family=s.family;data->generation=s.generation;data->degree=r.bkqr.power;
        data->records.push_back({0,0,0,0});r.gardenData=data;
        return r; // Full ensemble scan is performed only by the render worker.
    }
    if (s.family==20) {
        const int degree=8 << (std::clamp(s.spectralOrder,1,3)-1);
        r.masterModel=rhmaster::model(s.generation,degree,s.seed.real(),s.seed.imag(),s.spectralD,s.alpha,s.rhLevel,s.rhExperiment);
        auto data=std::make_shared<garden::Data>();
        data->family=20;data->generation=s.generation;data->degree=r.masterModel->degree;
        data->records.push_back({0,0,0,0});r.gardenData=data;return r;
    }
    if (s.family==19) {
        r.tripleModel=momenttriple::model(s.generation);
        auto data=std::make_shared<garden::Data>();
        data->family=19;data->generation=s.generation;data->degree=r.tripleModel->degree;
        data->records.push_back({0,0,0,0});r.gardenData=data;return r;
    }
    if (s.family==18) {
        r.quotientModel=momentquotient::model();
        auto data=std::make_shared<garden::Data>();
        data->family=18; data->generation=1; data->degree=2;
        data->records.push_back({0,0,0,0}); r.gardenData=data;
        return r;
    }
    if (s.family==17) {
        r.cubicModel=momentcubic::model();
        auto data=std::make_shared<garden::Data>();
        data->family=17; data->generation=1; data->degree=3;
        data->records.push_back({0,0,0,0}); r.gardenData=data;
        return r;
    }
    if (s.family == 16) {
        r.riemannModel=riemann::model(s.generation,s.riemann,s.alpha,s.method);
        auto data=std::make_shared<garden::Data>();
        data->family=16; data->generation=s.generation; data->degree=r.riemannModel->degree;
        data->records.push_back({0,0,0,0}); r.gardenData=data;
        return r;
    }
    if (s.family == 15) {
        r.projectModel=projects::model(s.generation,s.project,s.c,s.alpha);
        auto data=std::make_shared<garden::Data>();
        data->family=15; data->generation=s.generation; data->degree=r.projectModel->degree;
        data->records.push_back({0,0,0,0}); r.gardenData=data;
        return r;
    }
    if (s.family == 12 || s.family == 13 || s.family == 14) {
        // Histories use an accumulated spherical measure, never a per-pixel
        // garden iteration. Keep shared metadata safe without garden::data.
        auto data = std::make_shared<garden::Data>();
        data->family = s.family; data->generation = s.generation;
        data->degree = s.family == 13 && s.generation == 3 ? 6 : 0;
        data->records.push_back({0, 0, 0, 0});
        r.gardenData = data;
        return r;
    }
    if (s.family == 11) {
        // Keep the common catalogue/GPU record contract without asking the
        // Newton-garden factory to interpret a different dynamical family.
        auto data = std::make_shared<garden::Data>();
        data->family = 11;
        data->generation = s.generation;
        // Degree is generic: exceptional algebraic common factors can lower it.
        data->degree = s.spectralD == 0 || (s.plane != 2 && s.alpha == C(0, 0)) ? 2
                     : s.generation >= 2 ? 4 * s.spectralOrder + 6 : 6;
        data->records.push_back({0, 0, 0, 0});
        r.gardenData = data;
        return r;
    }
    if (s.family) {
        r.gardenData = garden::data(s.family, s.generation);
        if (s.family == 3 && s.plane != 1)
            r.crownRoots = crown::roots(r.gardenData->cylinders, s.c);
        return r;
    }
    r.a = analyse(s.c, s.alpha, point({2.4, .8}));
    r.b = analyse(s.c, s.alpha, point({2.4, -.8}));
    r.same =
        r.a.period && r.a.period == r.b.period &&
        std::any_of(r.a.cycle.begin(), r.a.cycle.end(), [&](Q q) { return chordal(q, r.b.last) < 1e-7; });
    return r;
}
inline void preset(State &s, int i) {
    if(s.family==25) {
        const int generation=s.generation,iterations=s.iterations;
        const C c=s.c;const double lambda=s.alpha.real();
        const double hue=s.hue,bands=s.bands;
        s=defaults(25,generation);s.iterations=iterations;s.hue=hue;s.bands=bands;
        if(generation>=9) {
            const int driver=generation+17;
            if(frk_extra::is_field_driver(driver)) {
                if(i==1)s.alpha={.25,0};
                if(i==2)s.alpha={1.10,0};
                if(i==3){s.alpha={.65,0};for(auto& camera:s.cameras)camera={1.,0.,1.8};s.palette=1;}
            } else if(frk_extra::is_root_driver(driver)) {
                if(i==1)s.alpha={.35,0};
                if(i==2)s.alpha={1.10,0};
                if(i==3){for(auto& camera:s.cameras)camera={-.8,.15,1.8};s.palette=1;}
            } else {
                if(i==1)s.plane=0;
                if(i==2){s.c=c;s.alpha={lambda,0};s.plane=1;}
                if(i==3){s.c={-.58,.36};s.alpha={.85,0};s.palette=1;s.plane=1;}
            }
            s.preset=i;return;
        }
        if(i==1)s.plane=frk_extra::seed_at_pixel(generation+17) ? 1 : 0;
        if(i==2){s.c=c;s.alpha={lambda,0};s.plane=1;}
        if(i==3){s.palette=1;s.plane=1;}
        s.preset=i;return;
    }
    if (s.family==24) {
        // Preserve the analytic host while selecting its view and truncation.
        s.preset=i;s.plane=i==1 ? 0 : 1;
        if(i==2)s.method=0;
        if(i==3)s.method=2;
        return;
    }
    if (s.family==22 || s.family==23) {
        // Change wrapper policy while preserving the chosen RH host and view.
        s.preset=i;s.orbitMode=i==1 ? 1 : 0;
        s.spectralOrder=i==2 ? 2 : i==3 ? 3 : 1;
        return;
    }
    if (s.family==21) {
        const int iterations=s.iterations;
        s=defaults(21,1);s.iterations=iterations;s.preset=i;
        if(i==1)s.orbitMode=1; // Legacy SUBTRACT branch.
        if(i==2)s.spectralOrder=2; // SELFDUAL coupling.
        if(i==3)s.spectralOrder=3; // COMPACT coupling.
        return;
    }
    if (s.family==20) {
        const int g=s.generation,iterations=s.iterations;
        const double hue=s.hue,bands=s.bands;
        s=defaults(20,g);s.iterations=iterations;s.hue=hue;s.bands=bands;s.preset=i;
        if(i==1)s.plane=0;
        if(i==2)s.spectralOrder=1;
        if(i==3)s.spectralOrder=3;
        return;
    }
    if (s.family==19) {
        const int g=s.generation,iterations=s.iterations;
        const double hue=s.hue,bands=s.bands;
        s=defaults(19,g);s.iterations=iterations;s.hue=hue;s.bands=bands;s.preset=i;
        if(i>0)s.plane=1;
        if(i==2)s.c=momenttriple::model(g)->criticalFixedParameter();
        if(i==3)s.c={.2,.55};
        return;
    }
    if (s.family==18) {
        const int iterations=s.iterations;
        const double hue=s.hue, bands=s.bands;
        s=defaults(18,1); s.iterations=iterations; s.hue=hue; s.bands=bands; s.preset=i;
        if (i>0) s.plane=1;
        if (i==1 || i==2) {
            const auto m=momentquotient::model();
            s.c={m->lambda*m->b,i==2 ? 1-2*m->lambda : 0};
        }
        if (i==3) s.c={.4,.6};
        return;
    }
    if (s.family==17) {
        const int iterations=s.iterations;
        const double hue=s.hue, bands=s.bands;
        s=defaults(17,1); s.iterations=iterations; s.hue=hue; s.bands=bands; s.preset=i;
        if (i>0) s.plane=1;
        if (i==2) {
            const auto m=momentcubic::model();
            s.c={m->critical*(1+2*m->A/3),0};
            s.cameras[1].span=5;
        }
        if (i==3) s.c={.15,.35};
        return;
    }
    if (s.family == 16) {
        const int g=s.generation, iterations=s.iterations;
        const double hue=s.hue, bands=s.bands;
        s=defaults(16,g); s.iterations=iterations; s.hue=hue; s.bands=bands; s.preset=i;
        if (g==1) {
            if (i==1) s.cameras[0]={.5,0,1.5};
            if (i==2) s.riemann.order=1;
            if (i==3) s.alpha={.85,.2};
        } else if (g==2) {
            if (i==1) s.cameras[0]={.0005,0,.002};
            if (i==2) s.cameras[0]={1,0,.15};
            if (i==3) s.method=1;
        } else if (g==3) {
            if (i==1) s.cameras[0]={.65,0,2.4};
            if (i==2) s.cameras[0]={10,0,32};
            if (i==3) s.alpha={.9,.2};
        } else {
            if (i==1) s.cameras[0].span=g==7 ? 3.2 : 2.6;
            if (i==2) s.method=1;
            if (i==3) s.alpha={.85,.2};
        }
        s.cameras[1]=s.cameras[2]=s.cameras[0];
        return;
    }
    if (s.family == 15) {
        const int g=s.generation,iterations=s.iterations;
        const double hue=s.hue,bands=s.bands;
        s=defaults(15,g); s.iterations=iterations; s.hue=hue; s.bands=bands; s.preset=i;
        if (g==1) {
            if (i==1) {s.c={-.25,.4};s.alpha={-8,3};}
            if (i==2) {s.project.grhDegree=32;s.project.grhRadius=16;s.cameras[0]={0,0,38};}
            if (i==3) {s.project.grhDegree=40;s.cameras[0]={0,5.5,8};}
        } else if (g==2) {
            if (i==1) s.project.piTerms=32;
            if (i==2) s.project.piTerms=128;
            if (i==3) s.project.piOffset=384;
        } else if (g==3) {
            if (i==1) s.project.nseTime=0;
            if (i==2) s.project.nseTime=.002;
            if (i==3) s.project.nseSlope=2;
        } else if (g==4) {
            if (i==1) s.cameras[0]={3,0,14};
            if (i==2) s.cameras[0]={-.75,0,3};
            if (i==3) s.cameras[0]={1.5,0,2.5};
        } else {
            if (i==1) s.project.pcdProbability=4972./19293.;
            if (i==2) s.project.pcdProbability=29620./175861.;
            if (i==3) s.cameras[0]={-1,1.5,6};
        }
        s.cameras[1]=s.cameras[2]=s.cameras[0];
        return;
    }
    if (s.family) {
        int family = s.family, generation = s.generation, it = s.iterations;
        double hue = s.hue, bands = s.bands;
        const int samplePower = s.quantum.samplePower;
        const uint32_t randomSeed = s.quantum.randomSeed;
        const auto spinSettings = s.spin;
        s = defaults(family, generation);
        s.iterations = it;
        s.hue = hue;
        s.bands = bands;
        s.preset = i;
        if (family == 1) {
            if (i == 1)
                s.alpha = {-.7, 0};
            if (i == 2) {
                s.c = {.65, .35};
                s.alpha = {.8, -.2};
            }
            if (i == 3)
                s.c = {0, 0};
        } else if (family == 5) {
            if (i == 1) s.c = {.5, 0};
            if (i == 2) { s.c = {-.4, .45}; s.method = 1; }
            if (i == 3) s.alpha = {.88, .22};
        } else if (family == 6) {
            if (i == 1) { s.method = 1; s.palette = 0; }
            if (i == 2) { s.c = {.3, .2}; s.palette = 0; }
            if (i == 3) { s.alpha = {.86, .18}; s.palette = 1; }
        } else if (family == 7) {
            if (i == 1) { s.method = 1; s.palette = 0; }
            if (i == 2) { s.cameras[0] = {0, 0, 1.8}; s.palette = 0; }
            if (i == 3) { s.alpha = {.86, .18}; s.method = 1; s.palette = 1; }
        } else if (family == 8) {
            if (i == 1) { s.method = 1; s.palette = 0; }
            if (i == 2) { s.cameras[0] = {.5, 0.8, 7.2}; s.palette = 0; }
            if (i == 3) { s.alpha = {.84, .20}; s.method = 1; s.palette = 1; }
        } else if (family == 9) {
            s.plane = 0; s.method = 0;
            if (i == 0) { s.cameras[0] = {6.25, 4.1, 5.8}; s.palette = 0; }
            if (i == 1) { s.cameras[0] = {5.1, 2.1, 8.0}; s.palette = 0; }
            if (i == 2) { s.cameras[0] = {7.7, 5.0, 5.2}; s.palette = 1; }
            if (i == 3) { s.cameras[0] = {6.5, 7.2, 4.2}; s.palette = 2; }
            s.cameras[1] = s.cameras[2] = s.cameras[0];
        } else if (family == 10) {
            s.plane = 0; s.method = 0;
            if (i == 0) { s.c={.10,.80}; s.alpha={.90,0}; s.cameras[0]={0,0,3.6}; }
            if (i == 1) { s.c={.05,.60}; s.alpha={.92,0}; s.cameras[0]={-.28,-.18,1.55}; s.bands=1.5; }
            if (i == 2) { s.c={.16,1.00}; s.alpha={.88,0}; s.cameras[0]={.20,.18,2.25}; s.bands=2.0; }
            if (i == 3) { s.c={.08,.80}; s.alpha={.94,0}; s.method=1; s.cameras[0]={0,0,3.6}; s.bands=1.25; }
            s.cameras[1]=s.cameras[2]=s.cameras[0];
        } else if (family == 14) {
            s.quantum.samplePower = samplePower;
            s.quantum.randomSeed = randomSeed;
            s.spin.burnPower = spinSettings.burnPower;
            s.spin.windowPower = spinSettings.windowPower;
            s.spin.gridPower = spinSettings.gridPower;
            if (generation == 1) {
                if (i == 1) s.spin.damping = .15;
                if (i == 2) s.spin.twist = 8;
                if (i == 3) s.spin.damping = 0;
            } else if (generation == 2) {
                if (i == 1) s.spin.twist = 8;
                if (i == 2) { s.spin.twist = 0; s.spin.tilt = 0; }
                if (i == 3) s.spin.damping = 0;
            } else {
                if (i == 1 || i == 2) {
                    s.spin.kMin = 2; s.spin.kMax = 7;
                    s.spin.tauMin = .05; s.spin.tauMax = .4;
                    s.spin.tilt = i == 1 ? 1.5 : 1.75;
                    s.spin.twist = 4.5; s.spin.damping = .18;
                }
                if (i == 3) {
                    s.spin.kMin = 6; s.spin.kMax = 14;
                    s.spin.tauMin = .1; s.spin.tauMax = 1;
                    s.spin.twist = 8;
                }
            }
        } else if (family == 13) {
            s.quantum.samplePower = samplePower;
            s.quantum.randomSeed = randomSeed;
            if (generation == 1) {
                if (i == 1) s.sphere.boost = .45;
                if (i == 2) s.plane = 2;
                if (i == 3) s.plane = 1;
            } else if (generation == 2) {
                static const double gaps[] = {.015, .06, .16, 0};
                s.sphere.gap = gaps[std::clamp(i, 0, 3)];
            } else {
                if (i == 1) s.sphere.lambdaLog = -3;
                if (i == 2) s.sphere.lambdaPhase = 1.15;
                if (i == 3) { s.sphere.lambdaLog = -1.2; s.sphere.lambdaPhase = 0; }
            }
        } else if (family == 12) {
            s.quantum.geometry = std::clamp(i, 0, 3);
            if (i == 1) s.quantum.strength = .58;
            if (i == 2) s.quantum.strength = .75;
            if (i == 3) { s.quantum.strength = .78; s.quantum.kick = .23; }
        } else if (family == 11) {
            // Presets are finite numerical exploration settings, not certified basins.
            s.plane = 0;
            s.method = 0;
            if (generation >= 2) {
                if (i == 1) { s.c = {-.4, .2}; s.alpha = {.16, .06}; s.spectralD = .28; s.spectralWeight = .60; }
                if (i == 2) { s.c = {-.25, .6}; s.alpha = {.06, .025}; s.spectralD = .32; s.spectralWeight = .35; s.spectralOrder = 2; }
                if (i == 3) { s.c = {-.3, .55}; s.alpha = {.025, .015}; s.spectralD = .30; s.spectralWeight = .45; s.spectralOrder = 3; }
            } else {
                if (i == 1) { s.c = {-.4, .2}; s.alpha = {.2, .02}; s.spectralD = .25; }
                if (i == 2) { s.c = {-.6, .24}; s.alpha = {.25, .08}; s.spectralD = .30; }
                if (i == 3) { s.c = {-.3, .55}; s.alpha = {0, 0}; s.spectralD = .25; }
            }
        } else if (family == 4) {
            if (i == 1)
                s.c = {1, 0};
            if (i == 2)
                s.c = {-.65, .55};
            if (i == 3) {
                s.c = {.85, .35};
                s.alpha = {.86, .18};
            }
        } else if (family == 3) {
            if (i == 1)
                s.c = {.30, .55};
            if (i == 2)
                s.c = {-.48, .34};
            if (i == 3) {
                s.c = {.45, .35};
                s.alpha = {.86, .16};
            }
        } else {
            if (i == 1)
                s.c = {.25, 0};
            if (i == 2)
                s.c = {.75, 0};
            if (i == 3)
                s.c = {.5, .2};
        }
        return;
    }
    static const C values[4] = {{1.2, .36}, {1.5, .54}, {1.7, .96}, {.6, .06}};
    static const Camera views[4] = {{2.1, -1.3, 8.2}, {2.2, -.8, 6.8}, {1.85, -.6, 5.4}, {2.1, -.8, 8.2}};
    s.c = values[std::clamp(i, 0, 3)];
    s.alpha = {1, 0};
    s.plane = 0;
    s.preset = i;
    s.palette = 0;
    s.orbitMode = 0;
    s.cameras[0] = views[i];
    s.cameras[1] = {s.c.real(), s.c.imag(), 3.8};
    s.cameras[2] = {1, 0, 5.2};
}
inline C world(const Camera &v, double x, double y, int w, int h) {
    return {v.x + (x / double(w) - .5) * v.span, v.y + (.5 - y / double(h)) * v.span * h / double(w)};
}
using Color = std::array<double, 3>;
inline Color ramp(double t) {
    return {.51 + .46 * std::cos(2 * pi * t), .51 + .46 * std::cos(2 * pi * (t + .105)),
            .51 + .46 * std::cos(2 * pi * (t + .19))};
}
inline Color scale(Color c, double a) {
    for (auto &v : c)
        v *= a;
    return c;
}
inline std::array<Color, 3> grhObservableColors(C p, const State &s) {
    auto e = grh_observable::evaluate(p.real(), p.imag());
    if (!e.valid)
        return {Color{.035,.04,.055}, Color{.035,.04,.055}, Color{.035,.04,.055}};
    double contrast = s.palette == 0 ? 1.0 : s.palette == 1 ? 1.35 : 1.75;
    double phaseTheta = e.plus.phaseTheta / (2.0 * pi);
    double thetaLight = .24 + .63 * (.5 + .5 * std::tanh((e.plus.logAbsTheta + 1.0) * .22 * contrast));
    auto left = scale(ramp(s.hue + s.bands * phaseTheta), thetaLight);

    double rm = std::log1p(std::abs(e.plus.R));
    double rp = std::arg(e.plus.R) / (2.0 * pi);
    auto middle = scale(ramp(s.hue + .17 + s.bands * rp), .18 + .77 * (1.0 - std::exp(-.36 * contrast * rm)));

    double value = e.observable.real();
    double u = std::asinh(value) * contrast;
    double signedLevel = std::tanh(.23 * u);
    double contour = .035 * s.bands * u;
    auto right = scale(ramp(s.hue + (value >= 0 ? .04 : .54) + .13 * signedLevel + contour),
                       .20 + .75 * (1.0 - std::exp(-.28 * std::abs(u))));
    if (e.nearZero) {
        for (auto *c : {&left, &middle, &right}) {
            (*c)[0] = std::min(1.0, (*c)[0] + .35);
            (*c)[1] = std::min(1.0, (*c)[1] + .35);
            (*c)[2] = std::min(1.0, (*c)[2] + .35);
        }
    }
    return {left, middle, right};
}
inline Color stretchColor(double ly, const State &s) {
    return scale(ramp(s.hue + .48 + s.bands * .20 * (ly * .65) / (1 + std::abs(ly * .65))),
                 .32 + .5 * (1 - std::exp(-std::abs(ly) * 1.8)));
}
struct Inspect {
    double ly = 0, trap = 1;
    int period = 0;
};
inline Inspect inspect(Q q, C c, C alpha, int n) {
    Inspect a;
    double sum = 0;
    int count = 0;
    Q poleA = point({2, 0}), poleB = point({4, 0});
    for (int i = 0; i < n; i++) {
        double ly;
        q = advance(q, c, alpha, &ly);
        if (i > n * .72) {
            sum += ly;
            count++;
        }
        a.trap = std::min(a.trap, std::min(chordal(q, poleA), chordal(q, poleB)));
    }
    a.ly = sum / std::max(count, 1);
    Q base = q;
    for (int j = 1; j <= 16; j++) {
        q = advance(q, c, alpha);
        if (a.period == 0 && chordal(q, base) < .00002)
            a.period = j;
    }
    return a;
}
inline Color trapColor(double d, const State &s) {
    double t = -std::log(std::max(1e-12, d));
    return scale(ramp(s.hue + .12 + t * .095 * s.bands), .36 + .48 * (1 - std::exp(-t * .22)));
}
inline Color renderColor(C p, const Scene &r) {
    const State &s = r.state;
    if(s.family==25) {
        frk_extra::Params q;
        q.lambda=s.alpha.real();q.julia_c=s.c;q.julia=(s.plane==1);
        q.iterations=s.iterations;
        const auto result=frk_extra::evaluate(s.generation+17,p,q);
        const uint32_t rgb=frk_extra::color(result,s.palette);
        return {double((rgb>>16)&255)/255.0,
                double((rgb>>8)&255)/255.0,double(rgb&255)/255.0};
    }
    if (s.family==24)
        return s.plane==0 ? r.qqrhDirectModel->atlasColour(p,s.iterations,s.palette,s.hue,s.bands)
                         : r.qqrhDirectModel->juliaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
    if (s.family==21 || s.family==22 || s.family==23) return {1,0,1}; // Full-grid metric path; never a per-pixel orbit.
    if (s.family==20) {
        if (s.plane==2)
            return r.masterModel->lambdaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
        return s.plane==0 ? r.masterModel->atlasColour(p,s.iterations,s.palette,s.hue,s.bands)
                         : r.masterModel->juliaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
    }
    if (s.family==19)
        return s.plane==0 ? r.tripleModel->atlasColour(p,s.iterations,s.palette,s.hue,s.bands)
                         : r.tripleModel->juliaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
    if (s.family==18)
        return s.plane==0 ? r.quotientModel->atlasColour(p,s.iterations,s.palette,s.hue,s.bands)
                         : r.quotientModel->juliaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
    if (s.family==17)
        return s.plane==0 ? r.cubicModel->atlasColour(p,s.iterations,s.palette,s.hue,s.bands)
                         : r.cubicModel->juliaColour(p,s.c,s.iterations,s.palette,s.hue,s.bands);
    if (s.family == 16)
        return projects::colour(r.riemannModel->orbit(p,s.iterations),s.palette,s.hue,s.bands);
    if (s.family == 15)
        return projects::colour(r.projectModel->orbit(p,s.iterations),s.palette,s.hue,s.bands);
    if (s.family == 12 || s.family == 13 || s.family == 14)
        return {1, 0, 1}; // Diagnostic only: these families use cached sphere clouds.
    if (s.family == 11) {
        const C z = s.plane == 0 ? p : s.seed;
        const C c = s.plane == 1 ? p : s.c;
        const C eta = s.plane == 2 ? p : s.alpha;
        const auto orbit = s.generation >= 2
            ? spectral::orbitGen2(z, c, eta, s.spectralD, s.spectralWeight, s.spectralOrder, s.iterations)
            : spectral::orbit(z, c, eta, s.spectralD, s.iterations);
        return spectral::colour(orbit, s.palette, s.hue, s.bands);
    }
    if (s.family == 9)
        return grhObservableColors(p, s)[std::clamp(s.palette, 0, 2)];
    if (s.family) {
        C z = s.plane == 0 ? p : s.seed, a = s.plane == 1 ? p : s.c, b = s.plane == 2 ? p : s.alpha;
        auto o = garden::orbit(z, a, b, *r.gardenData, s.method, s.iterations, nullptr,
                               r.crownRoots.empty() ? nullptr : &r.crownRoots);
        return garden::colour(o, s.palette, s.hue, s.bands, s.family);
    }
    if (s.plane != 0) {
        C c = s.plane == 1 ? p : s.c, alpha = s.plane == 2 ? p : s.alpha;
        Inspect a, b;
        if (s.orbitMode == 3)
            a = b = inspect(point(s.seed), c, alpha, s.iterations);
        else {
            if (s.orbitMode != 2)
                a = inspect(point({2.4, .8}), c, alpha, s.iterations);
            if (s.orbitMode != 1)
                b = inspect(point({2.4, -.8}), c, alpha, s.iterations);
            if (s.orbitMode == 1)
                b = a;
            if (s.orbitMode == 2)
                a = b;
        }
        double ly = std::max(a.ly, b.ly);
        if (s.palette == 1)
            return stretchColor(ly, s);
        if (s.palette == 2)
            return trapColor(std::min(a.trap, b.trap), s);
        int period = std::max(a.period, b.period);
        double t = s.hue + .08 + s.bands * (period * .069 + std::log1p(std::abs(ly)) * .10);
        double light = .32 + .47 * (1 - std::exp(-std::abs(ly) * 1.3));
        if (period == 0) {
            t = s.hue + .54 + .10 * std::sin(ly * 4 * s.bands);
            light = .43 + .4 * std::clamp(ly, 0., 1.);
        }
        auto col = scale(ramp(t), light);
        double slow = std::exp(-std::abs(ly) * 14);
        col[0] += .12 * slow;
        col[1] += .17 * slow;
        col[2] += .17 * slow;
        return col;
    }
    Q q = point(p);
    double ly = 0, sum = 0, distance = 1, trap = 1;
    int count = 0, hit = 0, at = 0;
    Q poleA = point({2, 0}), poleB = point({4, 0});
    for (int i = 0; i < s.iterations; i++) {
        q = advance(q, s.c, s.alpha, &ly);
        if (i > 8) {
            sum += ly;
            count++;
        }
        if (s.palette == 2)
            trap = std::min(trap, std::min(chordal(q, poleA), chordal(q, poleB)));
        if (s.palette == 0) {
            double da = r.a.period ? chordal(q, r.a.last) : 2,
                   db = r.b.period && !r.same ? chordal(q, r.b.last) : 2;
            double d = std::min(da, db);
            if (d < .00004) {
                hit = da <= db ? 1 : 2;
                at = i;
                distance = std::max(d, 1e-12);
                break;
            }
        }
    }
    if (s.palette == 2)
        return trapColor(trap, s);
    if (hit) {
        double smooth =
            at + std::clamp(std::log(distance / .00004) / std::min(sum / std::max(count, 1), -.08), -3., 3.);
        return scale(ramp(s.hue + .085 + std::log1p(std::max(smooth, 0.)) * .105 * s.bands + (hit - 1) * .40),
                     .4 + .33 * (1 - std::exp(-std::max(smooth, 0.) * .055)));
    }
    return stretchColor(sum / std::max(count, 1), s);
}
inline uint32_t pack(Color c) {
    uint32_t u[3];
    for (int i = 0; i < 3; i++)
        u[i] = uint32_t(std::clamp(std::pow(std::max(c[i], 0.), .88) * 255 + .5, 0., 255.));
    return u[2] | (u[1] << 8) | (u[0] << 16);
}
inline uint32_t pixel(int x, int y, int w, int h, const Scene &r, int samples = 1) {
    Color col{};
    for (int sy = 0; sy < samples; sy++)
        for (int sx = 0; sx < samples; sx++) {
            auto c = renderColor(
                world(r.state.cameras[r.state.plane], x + (sx + .5) / samples, y + (sy + .5) / samples, w, h),
                r);
            for (int k = 0; k < 3; k++)
                col[k] += c[k] / double(samples * samples);
        }
    return pack(col);
}
inline void writeView(const std::filesystem::path &p, const State &s) {
    std::ofstream f(p);
    if (!f)
        throw std::runtime_error("Cannot write view");
    f << std::setprecision(17) << "FRAKTALFAKTUR_VIEW " << (s.family==25 ? 17 : s.family==20 ? 16 : s.family==24 ? (s.generation==6 ? 13 : 12) : s.family==23 ? 11 : s.family==22 ? 10 : 9) << '\n';
    f << "family " << s.family << "\ngeneration " << s.generation << "\nmethod " << s.method << "\n";
    f << "c " << s.c.real() << ' ' << s.c.imag() << "\nalpha " << s.alpha.real() << ' ' << s.alpha.imag()
      << "\nseed " << s.seed.real() << ' ' << s.seed.imag();
    f << "\nplane " << s.plane << "\ncritical " << s.orbitMode << "\ncolour " << s.palette << "\niterations "
      << s.iterations << "\nhue " << s.hue << "\nbands " << s.bands
      << "\nspectral_d " << s.spectralD << "\nspectral_weight " << s.spectralWeight
      << "\nspectral_order " << s.spectralOrder << '\n';
    const auto& q = s.quantum;
    f << "quantum_strength " << q.strength << "\nquantum_kick " << q.kick
      << "\nquantum_axis_azimuth " << q.axisAzimuth << "\nquantum_axis_elevation " << q.axisElevation
      << "\nquantum_yaw " << q.yaw << "\nquantum_pitch " << q.pitch
      << "\nquantum_zoom " << q.zoom << "\nquantum_exposure " << q.exposure
      << "\nquantum_geometry " << q.geometry << "\nquantum_sample_power " << q.samplePower
      << "\nquantum_random_seed " << q.randomSeed << '\n';
    f << "sphere_gap " << s.sphere.gap << "\nsphere_lambda_log " << s.sphere.lambdaLog
      << "\nsphere_lambda_phase " << s.sphere.lambdaPhase << "\nsphere_boost " << s.sphere.boost << '\n';
    f << "spin_twist " << s.spin.twist << "\nspin_tilt " << s.spin.tilt
      << "\nspin_damping " << s.spin.damping << "\nspin_burn_power " << s.spin.burnPower
      << "\nspin_window_power " << s.spin.windowPower << "\nspin_grid_power " << s.spin.gridPower
      << "\nspin_k_min " << s.spin.kMin << "\nspin_k_max " << s.spin.kMax
      << "\nspin_tau_min " << s.spin.tauMin << "\nspin_tau_max " << s.spin.tauMax << '\n';
    f << "project_pi_terms " << s.project.piTerms << "\nproject_pi_offset " << s.project.piOffset
      << "\nproject_nse_modes " << s.project.nseModes << "\nproject_grh_degree " << s.project.grhDegree
      << "\nproject_nse_time " << s.project.nseTime << "\nproject_nse_slope " << s.project.nseSlope
      << "\nproject_pcd_probability " << s.project.pcdProbability << "\nproject_grh_radius " << s.project.grhRadius << '\n';
    f << "rh_witness_order " << s.riemann.order << "\nrh_radau_anchor " << s.riemann.anchor << '\n';
    if(s.family==20) {
        const auto& e=s.rhExperiment;
        f << "rh_level " << s.rhLevel << '\n'
          << "rh_map " << e.map << "\nrh_chart " << e.chart
          << "\nrh_coloring " << e.coloring << "\nrh_power " << e.power
          << "\nrh_gain " << e.gain << "\nrh_relaxation " << e.relaxation
          << "\nrh_memory " << e.memory << "\nrh_trap_radius " << e.trapRadius
          << "\nrh_trap_angle " << e.trapAngle << "\nrh_exposure " << e.exposure
          << "\nrh_frequency " << e.frequency
          << "\nrh_twist " << e.twist << "\nrh_blend " << e.blend << '\n';
    }
    if(s.family==22 || s.family==23)
        f << "bkqr_quality " << s.bkqrQuality << "\nbkqr_precision " << s.bkqrPrecision << '\n';
    if(s.family==23 || s.family==24)
        f << "qqrh_slope " << s.qqrhSlope << '\n';
    if(s.family==24 && s.generation==6)
        f << "qqrh_scale " << s.qqrhScale << '\n';
    for (int i = 0; i < 3; i++)
        f << "camera" << i << ' ' << s.cameras[i].x << ' ' << s.cameras[i].y << ' ' << s.cameras[i].span
          << '\n';
    if (!f)
        throw std::runtime_error("View write failed");
}
inline State readView(const std::filesystem::path &p) {
    std::ifstream f(p);
    if (!f)
        throw std::runtime_error("Cannot read view");
    std::string magic;
    int version = 0;
    if (!(f >> magic >> version))
        throw std::runtime_error("Malformed native view header");
    bool legacy = magic == "AMP_NATIVE_VIEW" && version == 1;
    if (!legacy && !(magic == "FRAKTALFAKTUR_VIEW" && (version >= 2 && version <= 17)))
        throw std::runtime_error("Not a supported native view");
    State s;
    std::string key;
    int found = 0;
    std::set<std::string> fields;
    const std::array<const char*,11> experimentFields{{
        "rh_map","rh_chart","rh_coloring","rh_power","rh_gain","rh_relaxation",
        "rh_memory","rh_trap_radius","rh_trap_angle","rh_exposure","rh_frequency"
    }};
    const std::array<const char*,2> focusFields{{"rh_twist","rh_blend"}};
    while (f >> key) {
        if (legacy && (key == "family" || key == "generation" || key == "method"))
            throw std::runtime_error("Version-1 AMP view contains a version-2 field");
        if (version < 3 && key == "spectral_d")
            throw std::runtime_error("Older view contains a version-3 field");
        if (version < 4 && (key == "spectral_weight" || key == "spectral_order"))
            throw std::runtime_error("Older view contains a version-4 field");
        if (version < 5 && key.rfind("quantum_", 0) == 0)
            throw std::runtime_error("Older view contains a version-5 field");
        if (version < 6 && key.rfind("sphere_", 0) == 0)
            throw std::runtime_error("Older view contains a version-6 field");
        if (version < 7 && key.rfind("spin_", 0) == 0)
            throw std::runtime_error("Older view contains a version-7 field");
        if (version < 8 && key.rfind("project_", 0) == 0)
            throw std::runtime_error("Older view contains a version-8 field");
        if (version < 9 && key.rfind("rh_", 0) == 0)
            throw std::runtime_error("Older view contains a version-9 field");
        if (version < 10 && key.rfind("bkqr_", 0) == 0)
            throw std::runtime_error("Older view contains a version-10 field");
        if (version < 11 && key.rfind("qqrh_", 0) == 0)
            throw std::runtime_error("Older view contains a version-11 field");
        if (version < 13 && key == "qqrh_scale")
            throw std::runtime_error("Older view contains a version-13 field");
        if (version < 14 && key == "rh_level")
            throw std::runtime_error("Older view contains a version-14 field");
        if (version < 15 && std::find(experimentFields.begin(),experimentFields.end(),key)!=experimentFields.end())
            throw std::runtime_error("Older view contains a version-15 field");
        if (version < 16 && std::find(focusFields.begin(),focusFields.end(),key)!=focusFields.end())
            throw std::runtime_error("Older view contains a version-16 field");
        if (!fields.insert(key).second)
            throw std::runtime_error("Duplicate view field");
        if (key == "c" || key == "alpha" || key == "seed") {
            double a = 0, b = 0;
            f >> a >> b;
            if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a) > 1e6 || std::abs(b) > 1e6)
                throw std::runtime_error("Invalid complex parameter");
            (key == "c" ? s.c : key == "alpha" ? s.alpha : s.seed) = {a, b};
        } else if (key == "bkqr_quality")
            f >> s.bkqrQuality;
        else if (key == "bkqr_precision")
            f >> s.bkqrPrecision;
        else if (key == "qqrh_slope")
            f >> s.qqrhSlope;
        else if (key == "qqrh_scale")
            f >> s.qqrhScale;
        else if (key == "rh_level")
            f >> s.rhLevel;
        else if (key == "rh_map")
            f >> s.rhExperiment.map;
        else if (key == "rh_chart")
            f >> s.rhExperiment.chart;
        else if (key == "rh_coloring")
            f >> s.rhExperiment.coloring;
        else if (key == "rh_power")
            f >> s.rhExperiment.power;
        else if (key == "rh_gain")
            f >> s.rhExperiment.gain;
        else if (key == "rh_relaxation")
            f >> s.rhExperiment.relaxation;
        else if (key == "rh_memory")
            f >> s.rhExperiment.memory;
        else if (key == "rh_trap_radius")
            f >> s.rhExperiment.trapRadius;
        else if (key == "rh_trap_angle")
            f >> s.rhExperiment.trapAngle;
        else if (key == "rh_exposure")
            f >> s.rhExperiment.exposure;
        else if (key == "rh_frequency")
            f >> s.rhExperiment.frequency;
        else if (key == "rh_twist")
            f >> s.rhExperiment.twist;
        else if (key == "rh_blend")
            f >> s.rhExperiment.blend;
        else if (key == "family")
            f >> s.family;
        else if (key == "generation")
            f >> s.generation;
        else if (key == "method")
            f >> s.method;
        else if (key == "plane")
            f >> s.plane;
        else if (key == "critical")
            f >> s.orbitMode;
        else if (key == "colour")
            f >> s.palette;
        else if (key == "iterations")
            f >> s.iterations;
        else if (key == "hue")
            f >> s.hue;
        else if (key == "bands")
            f >> s.bands;
        else if (key == "spectral_d")
            f >> s.spectralD;
        else if (key == "spectral_weight")
            f >> s.spectralWeight;
        else if (key == "spectral_order")
            f >> s.spectralOrder;
        else if (key == "sphere_gap") f >> s.sphere.gap;
        else if (key == "sphere_lambda_log") f >> s.sphere.lambdaLog;
        else if (key == "sphere_lambda_phase") f >> s.sphere.lambdaPhase;
        else if (key == "sphere_boost") f >> s.sphere.boost;
        else if (key == "spin_twist") f >> s.spin.twist;
        else if (key == "spin_tilt") f >> s.spin.tilt;
        else if (key == "spin_damping") f >> s.spin.damping;
        else if (key == "spin_burn_power") f >> s.spin.burnPower;
        else if (key == "spin_window_power") f >> s.spin.windowPower;
        else if (key == "spin_grid_power") f >> s.spin.gridPower;
        else if (key == "spin_k_min") f >> s.spin.kMin;
        else if (key == "spin_k_max") f >> s.spin.kMax;
        else if (key == "spin_tau_min") f >> s.spin.tauMin;
        else if (key == "spin_tau_max") f >> s.spin.tauMax;
        else if (key == "project_pi_terms") f >> s.project.piTerms;
        else if (key == "project_pi_offset") f >> s.project.piOffset;
        else if (key == "project_nse_modes") f >> s.project.nseModes;
        else if (key == "project_grh_degree") f >> s.project.grhDegree;
        else if (key == "project_nse_time") f >> s.project.nseTime;
        else if (key == "project_nse_slope") f >> s.project.nseSlope;
        else if (key == "project_pcd_probability") f >> s.project.pcdProbability;
        else if (key == "project_grh_radius") f >> s.project.grhRadius;
        else if (key == "rh_witness_order") f >> s.riemann.order;
        else if (key == "rh_radau_anchor") f >> s.riemann.anchor;
        else if (key == "quantum_strength") f >> s.quantum.strength;
        else if (key == "quantum_kick") f >> s.quantum.kick;
        else if (key == "quantum_axis_azimuth") f >> s.quantum.axisAzimuth;
        else if (key == "quantum_axis_elevation") f >> s.quantum.axisElevation;
        else if (key == "quantum_yaw") f >> s.quantum.yaw;
        else if (key == "quantum_pitch") f >> s.quantum.pitch;
        else if (key == "quantum_zoom") f >> s.quantum.zoom;
        else if (key == "quantum_exposure") f >> s.quantum.exposure;
        else if (key == "quantum_geometry") f >> s.quantum.geometry;
        else if (key == "quantum_sample_power") f >> s.quantum.samplePower;
        else if (key == "quantum_random_seed") {
            std::string value; f >> value;
            if (value.empty() || value.size() > 10 || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid quantum random seed");
            const auto seed = std::stoull(value);
            if (seed > UINT32_MAX) throw std::runtime_error("Quantum random seed out of range");
            s.quantum.randomSeed = uint32_t(seed);
        }
        else if (key == "camera0" || key == "camera1" || key == "camera2") {
            auto &v = s.cameras[key.back() - '0'];
            f >> v.x >> v.y >> v.span;
        } else
            throw std::runtime_error("Unknown view field: " + key);
        if (!f)
            throw std::runtime_error("Malformed view");
        found++;
    }
    const bool experimentComplete=std::all_of(experimentFields.begin(),experimentFields.end(),
        [&](const char* field){return fields.count(field)!=0;});
    const bool focusComplete=std::all_of(focusFields.begin(),focusFields.end(),
        [&](const char* field){return fields.count(field)!=0;});
    if (found != (legacy ? 12 : version == 2 ? 15 : version == 3 ? 16 : version == 4 ? 18 : version == 5 ? 29 : version == 6 ? 33 : version == 7 ? 43 : version == 8 ? 51 : version == 9 ? 53 : version == 10 ? 55 : version == 11 ? 56 : version == 12 ? 54 : version == 13 ? 55 : version == 14 ? 54 : version == 15 ? 65 : version == 17 ? 53 : 67) || s.family < 0 || s.family > 25 || s.generation < 1 ||
        (version == 10 && (s.family!=22 || !fields.count("bkqr_quality") || !fields.count("bkqr_precision"))) ||
        (version == 11 && (s.family!=23 || !fields.count("bkqr_quality") || !fields.count("bkqr_precision") || !fields.count("qqrh_slope"))) ||
        (version == 12 && (s.family!=24 || s.generation>5 || !fields.count("qqrh_slope") || fields.count("bkqr_quality") || fields.count("bkqr_precision"))) ||
        (version == 13 && (s.family!=24 || s.generation!=6 || !fields.count("qqrh_slope") || !fields.count("qqrh_scale") || fields.count("bkqr_quality") || fields.count("bkqr_precision"))) ||
        (version >= 14 && version <= 16 && (s.family!=20 || !fields.count("rh_level") || fields.count("qqrh_slope") || fields.count("qqrh_scale") || fields.count("bkqr_quality") || fields.count("bkqr_precision"))) ||
        (version >= 15 && version <= 16 && !experimentComplete) || (version == 16 && !focusComplete) ||
        (version == 17 && (s.family!=25 || fields.count("rh_level") || fields.count("qqrh_slope") || fields.count("bkqr_quality") || fields.count("rh_map") || fields.count("rh_twist"))) ||
        (s.family == 13 && version < 6) || (s.family == 14 && version < 7) ||
        (s.family == 15 && (version < 8 || s.plane != 0 || s.method != 0 || (s.generation == 1 && s.alpha == C(0,0)))) ||
        (s.family == 24 && (version != (s.generation==6 ? 13 : 12) || s.plane>1 ||
            std::abs(s.c)>1e6 || std::abs(s.alpha)>1e6 ||
            s.seed.real() < -1 || s.seed.real() > 1 || s.seed.imag() < .16 || s.seed.imag() > 2.4 ||
            s.spectralD < .025 || s.spectralD > .49 ||
            !std::isfinite(s.qqrhSlope) || s.qqrhSlope<0 || s.qqrhSlope>64 ||
            (s.generation==6 && (!std::isfinite(s.qqrhScale) || s.qqrhScale<.25 || s.qqrhScale>4)))) ||
        (s.family == 23 && (version != 11 || s.plane!=2 || s.orbitMode>2 ||
            std::abs(s.c)>1e6 || std::abs(s.alpha)>1e6 ||
            s.seed.real() < -1 || s.seed.real() > 1 || s.seed.imag() < .16 || s.seed.imag() > 2.4 ||
            s.spectralD < .025 || s.spectralD > .49 ||
            s.bkqrQuality<0 || s.bkqrQuality>2 || s.bkqrPrecision<0 || s.bkqrPrecision>4 ||
            !std::isfinite(s.qqrhSlope) || s.qqrhSlope<0 || s.qqrhSlope>64)) ||
        (s.family == 22 && (version != 10 || s.plane!=2 || s.orbitMode>2 ||
            std::abs(s.c)>1e6 || std::abs(s.alpha)>1e6 ||
            s.seed.real() < -1 || s.seed.real() > 1 || s.seed.imag() < .16 || s.seed.imag() > 2.4 ||
            s.spectralD < .025 || s.spectralD > .49 ||
            s.bkqrQuality<0 || s.bkqrQuality>2 || s.bkqrPrecision<0 || s.bkqrPrecision>4)) ||
        (s.family == 21 && (version < 9 || s.generation!=1 || s.plane!=2 || std::abs(s.c)>1e6 ||
            s.alpha.imag()!=0 || s.alpha.real()<2 || s.alpha.real()>8 || s.alpha.real()!=std::floor(s.alpha.real()) ||
            s.seed.real()<0 || s.seed.real()>2 || s.seed.real()!=std::floor(s.seed.real()) ||
            s.seed.imag()<0 || s.seed.imag()>4 || s.seed.imag()!=std::floor(s.seed.imag()) || s.orbitMode>2)) ||
        (s.family == 25 && (version!=17 || s.plane>1 || s.method!=0 ||
            std::abs(s.c)>1e6 || !std::isfinite(s.alpha.real()) ||
            s.alpha.real() < -2. || s.alpha.real() > 2. || s.alpha.imag()!=0 ||
            (s.generation>=9 && (s.alpha.real()<.1 || s.alpha.real()>1.25)) ||
            (s.generation>=9 && frk_extra::seed_at_pixel(s.generation+17) && s.plane!=1))) ||
        (s.family == 20 && ((version != 14 && version != 15 && version != 16) || s.plane>2 || s.method!=0 || std::abs(s.c)>1e6 ||
            std::abs(s.alpha)>1e6 ||
            s.seed.real() < .05 || s.seed.real() > 3 || s.seed.imag() < .16 || s.seed.imag() > 2.4 ||
            s.spectralD < .025 || s.spectralD > .49 || !std::isfinite(s.rhLevel) || s.rhLevel < -12 || s.rhLevel > 12 ||
            !rhmaster::validExperiment(s.rhExperiment))) ||
        (s.family == 19 && (version < 9 || s.plane>1 || s.method!=0 || std::abs(s.c)>1e6)) ||
        (s.family == 18 && (version < 9 || s.plane>1 || s.method!=0 || std::abs(s.c)>1e6)) ||
        (s.family == 17 && (version < 9 || s.plane>1 || s.method!=0 || std::abs(s.c)>1e6)) ||
        (s.family == 16 && (version < 9 || s.plane != 0 || std::abs(s.alpha)==0 || std::abs(s.alpha)>2 || (s.generation>=4 && s.riemann.order!=2))) ||
        s.generation > (s.family == 25 ? frk_extra::driver_count : s.family == 24 ? 6 : s.family == 23 ? 5 : s.family == 22 ? 3 : s.family == 21 ? 1 : s.family == 20 ? rhfresh::CandidateCount : s.family == 19 ? 3 : s.family == 18 ? 1 : s.family == 17 ? 1 : s.family == 16 ? 7 : s.family == 15 ? 5 : s.family == 14 ? 3 : s.family == 13 ? 3 : s.family == 12 ? 1 : s.family == 11 ? 2 : (s.family == 1 || s.family >= 3) ? 4
                        : s.family == 2                  ? 6
                                                         : 1) ||
        s.method < 0 || s.method > ((s.family==22 || s.family==23 || s.family==24) ? 2 : s.family==21 ? 4 : 1) || s.plane < 0 || s.plane > 2 || s.orbitMode < 0 || s.orbitMode > 3 ||
        s.palette < 0 || s.palette > ((s.family==21 || s.family==22 || s.family==23) ? 8 : 2) || s.iterations < 32 || s.iterations > ((s.family==21 || s.family==22 || s.family==23) ? 4096 : 1024) || !std::isfinite(s.hue) ||
        s.hue < 0 || s.hue > 1 || !std::isfinite(s.bands) || s.bands < .2 || s.bands > 3 ||
        !std::isfinite(s.spectralD) || s.spectralD < 0 || s.spectralD > .49 ||
        !std::isfinite(s.spectralWeight) || s.spectralWeight < .01 || s.spectralWeight > .98 ||
        s.spectralOrder < 1 || s.spectralOrder > 3 || !quantum::valid(s.quantum) || !spherical::valid(s.sphere) || !spin::valid(s.spin) || !projects::valid(s.project) || !riemann::valid(s.riemann))
        throw std::runtime_error("View values out of range");
    for (auto v : s.cameras)
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.span) || std::abs(v.x) > 1e6 ||
            std::abs(v.y) > 1e6 || v.span < 1e-10 || v.span > 1e5)
            throw std::runtime_error("Invalid camera");
    if(s.family!=20 && s.family!=25) throw std::runtime_error("This catalogue accepts only registered fractal candidate views.");
    s.preset = -1;
    return s;
}
inline void writeBmp(const std::filesystem::path &p, const std::vector<uint32_t> &data, int w, int h) {
    std::ofstream f(p, std::ios::binary);
    if (!f)
        throw std::runtime_error("Cannot write BMP");
    auto u16 = [&](uint16_t v) {
        for (int i = 0; i < 2; i++)
            f.put(char(v >> (8 * i)));
    };
    auto u32 = [&](uint32_t v) {
        for (int i = 0; i < 4; i++)
            f.put(char(v >> (8 * i)));
    };
    f.write("BM", 2);
    u32(54 + w * h * 4);
    u32(0);
    u32(54);
    u32(40);
    u32(w);
    u32(uint32_t(-h));
    u16(1);
    u16(32);
    u32(0);
    u32(w * h * 4);
    u32(2835);
    u32(2835);
    u32(0);
    u32(0);
    f.write(reinterpret_cast<const char *>(data.data()), std::streamsize(data.size() * 4));
    if (!f)
        throw std::runtime_error("BMP write failed");
}
inline int selfTest(std::ostream &out) {
    int checks = 0;
    auto check = [&](bool ok, const char *name) {
        checks++;
        if (!ok)
            throw std::runtime_error(std::string("Self-test: ") + name);
    };
    C c{1.2, .36}, alpha{1, .17};
    for (int i = 0; i < 250; i++) {
        C z{std::sin(i * .77) * 5.3, std::cos(i * .33) * 4.7};
        Q q = advance(point(z), c, alpha);
        C expected = c + alpha * .6 * (z / (2.0 - z) + 2.0 * z / (4.0 - z));
        check(chordal(q, point(expected)) < 2e-13, "affine/homogeneous agreement");
        double ly;
        advance(point(z), c, alpha, &ly);
        C derivative = alpha * .6 * (2.0 / ((2.0 - z) * (2.0 - z)) + 8.0 / ((4.0 - z) * (4.0 - z)));
        double sph = std::abs(derivative) * (1 + std::norm(z)) / (1 + std::norm(expected));
        check(std::abs(std::exp(ly) - sph) < 1e-10 * std::max(1., sph), "spherical derivative");
    }
    for (double pole : {2., 4.})
        check(chordal(advance(point({pole, 0}), c, alpha), infinity()) < 1e-14, "pole maps to infinity");
    check(chordal(advance(infinity(), c, alpha), point(c - 1.8 * alpha)) < 1e-14, "infinity returns");
    for (Q q : {point({2, 0}), point({4, 0}), infinity()})
        check(chordal(advance(q, c, {0, 0}), point(c)) < 1e-14, "zero multiplier is constant");
    for (double sign : {-1., 1.}) {
        double ly;
        advance(point({2.4, .8 * sign}), c, alpha, &ly);
        check(std::exp(ly) < 1e-13, "critical derivative zero");
    }
    State s;
    auto r = scene(s);
    check(r.a.period == 1 && r.b.period == 4, "fourfold critical periods");
    preset(s, 1);
    r = scene(s);
    check(r.a.period == 1 && r.b.period == 7, "seventh critical periods");
    preset(s, 2);
    r = scene(s);
    check(r.a.period == 1 && r.b.period == 3, "triplet critical periods");
    {
        const C x{.71, .19};
        const double d = .25, t = d * d;
        const C rm = 1.0 / (x - d), rp = 1.0 / (x + d);
        const C r1 = .5 * (rm + rp), r2 = .5 * (rm * rm + rp * rp);
        const C r3 = .5 * (rm * rm * rm + rp * rp * rp);
        const C determinant = r1 * r3 / (r2 * r2) - 1.0;
        const C kernel = t * (x * x - t) / ((x * x + t) * (x * x + t));
        check(std::abs(determinant - kernel) < 2e-13, "spectral jet determinant identity");
        const auto pole = spectral::orbit({0, d}, {0, 0}, {.1, 0}, d, 32);
        check(pole.status == spectral::Status::Escaped && pole.poleHit, "spectral pole maps to infinity");
        const auto control = spectral::orbit({0, d}, {0, 0}, {0, 0}, d, 1);
        check(control.status == spectral::Status::Unresolved && !control.poleHit &&
              std::abs(control.last - C(-t, 0)) < 1e-14, "eta zero removes spectral poles");
        const auto collapsed = spectral::orbit({.5, 0}, {0, 0}, {.1, .2}, 0, 1);
        check(collapsed.status == spectral::Status::Unresolved &&
              std::abs(collapsed.last - C(.25, 0)) < 1e-14, "coincident spectral atoms are quadratic");
        const auto atom = spectral::orbit({d, 0}, {0, .1}, {.3, .2}, d, 1);
        check(atom.status == spectral::Status::Unresolved &&
              std::abs(atom.last - C(t, .1)) < 1e-14, "spectral atom location is removable");
        const auto escaped = spectral::orbit({4, 0}, {0, 0}, {.1, 0}, d, 32);
        check(escaped.status == spectral::Status::Escaped && escaped.steps == 0,
              "spectral proven escape radius");
        const auto invalid = spectral::orbit({std::numeric_limits<double>::quiet_NaN(), 0},
                                             {0, 0}, {.1, 0}, d, 32);
        check(invalid.status == spectral::Status::Numerical, "spectral NaN is not escape");
        std::vector<C> trail;
        const auto unresolved = spectral::orbit({0, 0}, {0, 0}, {0, 0}, d, 32, &trail);
        check(unresolved.status == spectral::Status::Unresolved && trail.size() == 33 &&
              std::all_of(trail.begin(), trail.end(), [](C z) { return spectral::finite(z); }),
              "spectral finite trail and unresolved budget");
        for (int palette = 0; palette < 3; ++palette) {
            const auto color = spectral::colour(escaped, palette, .25, 1.2);
            check(std::all_of(color.begin(), color.end(), [](double v) {
                return std::isfinite(v) && v >= 0 && v <= 1;
            }), "spectral colours remain finite");
        }
    }
    {
        const double d = .32, m = .35, w = (1.0 - m) * .5;
        const std::array<double, 3> atoms{{-d, .35 * d, d}}, weights{{w, m, w}};
        const C x{.71, .29}, c2{-.3, .55}, eta{.08, .025};
        for (int order = 1; order <= 3; ++order) {
            std::array<C, 3> resolvent{};
            for (int j = 0; j < 3; ++j) {
                const C inverse = 1.0 / (x - atoms[j]);
                C power{1, 0};
                for (int k = 0; k < order; ++k) power *= inverse;
                for (int k = 0; k < 3; ++k) {
                    resolvent[k] += weights[j] * power;
                    power *= inverse;
                }
            }
            const auto kernel = spectral::gen2Kernel(x, d, m, order);
            const C expected = resolvent[0] * resolvent[2] / (resolvent[1] * resolvent[1]) - 1.0;
            const C actual = (kernel.numerator / kernel.denominator) / kernel.denominator;
            check(std::abs(actual - expected) < 2e-12, "generation-2 raw/cancelled jet agreement");
            for (double atom : atoms) {
                const auto step = spectral::orbitGen2({atom, 0}, c2, eta, d, m, order, 1);
                check(step.status == spectral::Status::Unresolved &&
                      std::abs(step.last - (C(atom * atom, 0) + c2)) < 2e-12,
                      "generation-2 spectral atom singularity removed");
            }
            const auto edge = spectral::gen2Kernel(C(8 * d, 0), d, m, order);
            check(std::abs((edge.numerator / edge.denominator) / edge.denominator) < 1,
                  "generation-2 jet satisfies exterior estimate");
        }
        const C start{.11, .27};
        const auto noFeedback = spectral::orbitGen2(start, c2, {0, 0}, d, m, 3, 1);
        const auto collapsed = spectral::orbitGen2(start, c2, eta, 0, m, 3, 1);
        check(noFeedback.status == spectral::Status::Unresolved &&
              std::abs(noFeedback.last - (start * start + c2)) < 1e-14,
              "generation-2 eta-zero quadratic bypass");
        check(collapsed.status == spectral::Status::Unresolved &&
              std::abs(collapsed.last - (start * start + c2)) < 1e-14,
              "generation-2 coincident-atom quadratic bypass");
        const double radius = spectral::escapeRadiusGen2(c2, eta, d);
        const auto escaped = spectral::orbitGen2({radius + 1, 0}, c2, eta, d, m, 3, 32);
        check(escaped.status == spectral::Status::Escaped && escaped.steps == 0,
              "generation-2 proved escape radius");
    }
    {
        const quantum::Vec3 r = quantum::normalized({.23, -.51, .72});
        for (int geometry = 0; geometry < 4; ++geometry) {
            const auto& ns = quantum::directions(geometry);
            const int sizes[] = {4, 6, 12, 20};
            check(int(ns.size()) == sizes[geometry], "quantum Platonic direction count");
            quantum::Vec3 sum{0,0,0};
            double second[3][3]{};
            for (const auto n : ns) {
                check(std::abs(quantum::length(n)-1) < 2e-14, "quantum unit measurement direction");
                sum = sum+n;
                const double v[] = {n.x,n.y,n.z};
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) second[a][b]+=v[a]*v[b];
            }
            check(quantum::length(sum) < 2e-14, "quantum zero direction sum gives Kraus completeness");
            for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                check(std::abs(second[a][b]-(a==b?ns.size()/3.0:0)) < 3e-14,
                      "quantum isotropic second moment");
            for (double e : {0.0, .5, 1.0}) {
                double probabilitySum=0;
                quantum::Vec3 mean{0,0,0};
                for (const auto n : ns) {
                    const double p=quantum::probability(r,n,e,int(ns.size()));
                    check(p>=0 && p<=1, "quantum valid Born outcome probability");
                    probabilitySum+=p;
                    if (p>0) {
                        const auto next=quantum::measure(r,n,e);
                        check(std::abs(quantum::length(next)-1) < 2e-14, "quantum update preserves pure-state norm");
                        mean=mean+next*p;
                    }
                }
                check(std::abs(probabilitySum-1) < 2e-14, "quantum Born probabilities sum to one");
                check(quantum::length(mean-r*quantum::meanFactor(e)) < 3e-14,
                      "quantum exact averaged channel contraction");
            }
            const auto n=ns.front();
            check(quantum::probability(-n,n,1,int(ns.size())) < 2e-16,
                  "quantum antipodal projective outcome has zero probability");
            check(quantum::length(quantum::measure(r,n,0)-r) == 0,
                  "quantum zero measurement strength leaves state unchanged");
        }
        const auto axis=quantum::kickAxis(quantum::Settings{});
        const auto rotated=quantum::rotate(r,axis,.23);
        check(std::abs(quantum::length(rotated)-1) < 2e-14, "quantum unitary kick preserves Bloch norm");
        check(quantum::length(quantum::rotate(rotated,axis,-.23)-r) < 2e-14,
              "quantum inverse unitary kick");
    }
    checks += spherical::selfTest(out);
    checks += spin::selfTest(out);
    checks += projects::selfTest(out);
    checks += project_pi::selfTest(out);
    out << "PASS: " << checks << " independent mathematical checks\n";
    return checks;
}
} // namespace amp
