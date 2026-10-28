#pragma once
#include "lfunction_math.hpp"
#include <cmath>
#include <complex>

// Primitive GRH Cayley coordinate used by renderer family 7.
// Character: primitive quartic chi modulo 5, chi(2)=i (lfunction generation 4).
//
//   w = tanh(log(2) * (s - 1/2))
//   s_k(w) = 1/2 + (atanh(w) + i*pi*k) / log(2)
//
// tanh has period i*pi, so k labels the vertical Cayley sheet.  On every
// sheet Re(s)=1/2 maps exactly to Re(w)=0.  At target=0 the analytic function
// is the completed primitive L-function up to a nowhere-zero common factor:
//
//   Lambda(s,chi) = A(s) L(s,chi),
//   A(s) = (5/pi)^((s+1)/2) Gamma((s+1)/2).
//
// Newton/Halley is evaluated after dividing f,f',f'' by the same A(s).  This
// leaves both corrections exactly unchanged while avoiding the exponential
// Gamma decay from being mistaken for a small zero residual.
namespace grh_cayley {
using C = std::complex<double>;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double log2 = 0.693147180559945309417232121458176568;
constexpr int lGeneration = 4;

struct Evaluation {
    C f{}, first{}, second{}, s{};
    double residual = 0;
    bool valid = false;
};

inline bool sheetPoint(C w, int sheet, C &s, C *ds = nullptr, C *dds = nullptr) {
    if (sheet < 1 || sheet > 4 || !std::isfinite(w.real()) || !std::isfinite(w.imag()))
        return false;
    C oneMinus = 1.0 - w * w;
    if (std::abs(oneMinus) < 1e-12)
        return false;
    C a = std::atanh(w);
    s = C(.5, 0) + (a + C(0, pi * sheet)) / log2;
    if (!std::isfinite(s.real()) || !std::isfinite(s.imag()) || s.real() < -8.0 || s.real() > 12.0 ||
        std::abs(s.imag()) > 80.0)
        return false;
    if (ds)
        *ds = 1.0 / (log2 * oneMinus);
    if (dds)
        *dds = 2.0 * w / (log2 * oneMinus * oneMinus);
    return true;
}

inline Evaluation evaluate(C w, C target, int sheet) {
    Evaluation out;
    C ds, dds;
    if (!sheetPoint(w, sheet, out.s, &ds, &dds))
        return out;
    auto L = lfunction::evaluate(out.s, lGeneration);
    if (!L.valid)
        return out;

    // All primitive nonprincipal characters currently shipped here are odd;
    // for chi modulo 5, chi(-1)=-1, hence a_chi=1.
    C gammaArg = (out.s + 1.0) * .5;
    auto gamma = lfunction::detail::logGammaJet(gammaArg);
    C h1 = .5 * std::log(5.0 / pi) + .5 * gamma.psi;
    C h2 = .25 * gamma.trigamma;

    // f/A = L-target and its s-derivatives after the common completion factor
    // A has been divided out.  Newton and Halley ratios are unchanged.
    C u = L.f - target;
    C fs1 = L.first + h1 * u;
    C fs2 = L.second + 2.0 * h1 * L.first + (h1 * h1 + h2) * u;

    out.f = u;
    out.first = fs1 * ds;
    out.second = fs2 * ds * ds + fs1 * dds;
    out.residual = std::abs(u);
    out.valid = std::isfinite(out.f.real()) && std::isfinite(out.f.imag()) &&
                std::isfinite(out.first.real()) && std::isfinite(out.first.imag()) &&
                std::isfinite(out.second.real()) && std::isfinite(out.second.imag());
    return out;
}
} // namespace grh_cayley
