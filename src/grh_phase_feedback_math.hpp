#pragma once
#include "grh_cayley_math.hpp"
#include <cmath>
#include <complex>

// Exploratory nonlinear dynamics built from the exact primitive completed-L
// Cayley jets.  This is NOT an additional GRH equivalence.  It preserves every
// zero of the completed primitive L-function as a fixed point while feeding a
// bounded residual-phase term back into Newton/Halley motion.
//
//   f(w) = Lambda(s_3(w),chi_5) / A(s_3(w))       (common factor removed)
//   N_f  = f/f'  or  2ff'/(2f'^2-ff'')
//   Phi(f;tau) = tanh(beta |f|) f/|f| exp(i tau arg f)
//
//   w_{n+1} = w_n - rho N_f(w_n)
//                   + kappa Phi(f(w_n);tau)/(1+eta |w_n|^2).
//
// Phi -> 0 as f -> 0, so the GRH zero attractors are not displaced.
namespace grh_phase_feedback {
using C = std::complex<double>;
constexpr int sheet = 3;
constexpr double beta = 0.8;
constexpr double eta = 0.4;

struct Step {
    C delta{};          // orbit update uses w <- w - delta
    double residual = 1;
    double phase = 0;
    bool valid = false;
};

inline C boundedPhase(C f, double twist) {
    double a = std::abs(f);
    if (!(a > 1e-300) || !std::isfinite(a))
        return C(0,0);
    C unit = f / a;
    return std::tanh(beta * a) * unit * std::polar(1.0, twist * std::arg(f));
}

inline Step correction(C w, C feedback, C relaxation, int method) {
    Step out;
    auto e = grh_cayley::evaluate(w, C(0,0), sheet);
    if (!e.valid)
        return out;
    out.residual = e.residual;
    out.phase = std::arg(e.f);
    if (out.residual < 3e-6) {
        out.valid = true;
        return out;
    }
    C base;
    if (method == 1) {
        C den = 2.0*e.first*e.first - e.f*e.second;
        if (std::norm(den) < 1e-26)
            return out;
        base = 2.0*e.f*e.first/den;
    } else {
        if (std::norm(e.first) < 1e-26)
            return out;
        base = e.f/e.first;
    }
    // feedback.real = kick kappa; feedback.imag = phase twist tau.
    C kick = feedback.real() * boundedPhase(e.f, feedback.imag()) /
             (1.0 + eta * std::norm(w));
    out.delta = relaxation * base - kick;
    out.valid = std::isfinite(out.delta.real()) && std::isfinite(out.delta.imag());
    return out;
}
} // namespace grh_phase_feedback
