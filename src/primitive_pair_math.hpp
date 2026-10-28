#pragma once
#include "lfunction_math.hpp"
#include <cmath>
#include <complex>

// Primitive character-pair completed-L dynamics.
//
// Character: primitive quartic chi modulo 5 with chi(2)=i, paired with
// its complex conjugate character.  With
//
//   A(s) = (5/pi)^((s+1)/2) Gamma((s+1)/2),
//   Lambda_chi(s)    = A(s) L(s,chi),
//   Lambda_chibar(s) = A(s) L(s,chibar),
//
// the zero-target object is the completed pair product
//
//   P(s) = Lambda_chi(s) Lambda_chibar(s).
//
// For numerical stability all jets are divided by the common nowhere-zero
// factor A(s)^2. Newton and Halley corrections are therefore unchanged.
// A nonzero UI target is a normalized exploratory deformation
// A(s)^2 (L(s,chi)L(s,chibar)-target), not P(s)-target.
namespace primitive_pair {
using C = std::complex<double>;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr int lGeneration = 4;

struct Evaluation {
    C f{}, first{}, second{};
    double residual = 0;
    bool valid = false;
};

inline Evaluation evaluate(C s, C target = C(0,0)) {
    Evaluation out;
    auto pair = lfunction::evaluateQuarticPair(s);
    const auto &L = pair.first;
    const auto &B = pair.second;
    if (!L.valid || !B.valid)
        return out;

    C B0 = B.f;
    C B1 = B.first;
    C B2 = B.second;

    C u0 = L.f * B0 - target;
    C u1 = L.first * B0 + L.f * B1;
    C u2 = L.second * B0 + 2.0 * L.first * B1 + L.f * B2;

    auto gamma = lfunction::detail::logGammaJet((s + 1.0) * .5);
    C h1 = .5 * std::log(5.0 / pi) + .5 * gamma.psi;
    C h2 = .25 * gamma.trigamma;
    C k1 = 2.0 * h1; // (log A^2)'
    C k2 = 2.0 * h2; // (log A^2)''

    out.f = u0;
    out.first = u1 + k1 * u0;
    out.second = u2 + 2.0 * k1 * u1 + (k1 * k1 + k2) * u0;
    out.residual = std::abs(u0);
    out.valid = lfunction::detail::finite(out.f) && lfunction::detail::finite(out.first) &&
                lfunction::detail::finite(out.second);
    return out;
}
} // namespace primitive_pair
