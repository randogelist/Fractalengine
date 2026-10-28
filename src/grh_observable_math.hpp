#pragma once
#include "primitive_pair_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

// Literal completed-pair Cayley observable used by the all-character GRH
// closure, specialized here to the primitive quartic character mod 5.
//
//   P_chi(s) = Lambda(s,chi) Lambda(s,chibar)
//   Theta_chi(z) = P_chi(1/2 + sqrt(z))       (principal sqrt / slit plane)
//   R_chi(z) = Theta'_chi(z) / Theta_chi(z)
//
//   z_+(a,r) = a - 2ar + 2 i a sqrt(r(1-r))
//   z_-(a,r) = conjugate(z_+)
//
//   C_chi(a,r) = 4a ( z_- R(z_-) - z_+ R(z_+) ) / (z_- - z_+).
//
// The native diagnostic parameterization is
//      x = log2(a),  y = logit(r),
// so the entire theorem domain a>0, 0<r<1 fits in an unconstrained plane.
namespace grh_observable {
using C = std::complex<double>;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double log2c = 0.693147180559945309417232121458176568;

struct ThetaPoint {
    C z{}, sqrtz{}, s{};
    C normalizedTheta{}; // P/A^2 = L_chi L_chibar; same zero set.
    C R{};               // literal Theta'/Theta.
    double logAbsTheta = 0;
    double phaseTheta = 0;
    bool valid = false;
    bool nearZero = false;
};

struct Evaluation {
    double log2a = 0, logitR = 0, a = 1, r = .5;
    C zPlus{}, zMinus{};
    ThetaPoint plus, minus;
    C observable{};
    double symmetryError = 0;
    bool valid = false;
    bool nearZero = false;
};

inline double logistic(double y) {
    if (y >= 0) {
        double e = std::exp(-std::min(y, 700.0));
        return 1.0 / (1.0 + e);
    }
    double e = std::exp(std::max(y, -700.0));
    return e / (1.0 + e);
}
inline bool finite(C z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }

inline ThetaPoint thetaPoint(C z) {
    ThetaPoint out;
    out.z = z;
    out.sqrtz = std::sqrt(z); // principal square root: slit-plane pullback.
    if (std::abs(out.sqrtz) < 1e-14 || !finite(out.sqrtz))
        return out;
    out.s = C(.5, 0) + out.sqrtz;
    auto pair = primitive_pair::evaluate(out.s, C(0, 0));
    if (!pair.valid || !finite(pair.f) || !finite(pair.first))
        return out;
    out.normalizedTheta = pair.f;
    double fAbs = std::abs(pair.f);
    out.nearZero = fAbs < 1e-9;
    if (!(fAbs > 1e-300))
        return out;

    // primitive_pair::evaluate returns P/A^2 and P'/A^2. Their ratio is
    // exactly the completed-pair logarithmic derivative P'/P.
    C logDerivativeS = pair.first / pair.f;
    out.R = logDerivativeS / (2.0 * out.sqrtz);

    // Literal completed-pair phase/magnitude without exponentiating A(s)^2.
    auto gamma = lfunction::detail::logGammaJet((out.s + 1.0) * .5);
    C logA = ((out.s + 1.0) * .5) * std::log(5.0 / pi) + gamma.logGamma;
    C logTheta = 2.0 * logA + std::log(pair.f);
    out.logAbsTheta = logTheta.real();
    out.phaseTheta = std::remainder(logTheta.imag(), 2.0 * pi);
    out.valid = finite(out.R) && std::isfinite(out.logAbsTheta) && std::isfinite(out.phaseTheta);
    return out;
}

inline Evaluation evaluate(double log2a, double logitR) {
    Evaluation out;
    out.log2a = log2a;
    out.logitR = logitR;
    // Keep the analytic L evaluator inside its validated numerical window.
    if (!std::isfinite(log2a) || !std::isfinite(logitR) || log2a < -14.0 || log2a > 10.0 ||
        logitR < -12.0 || logitR > 14.0)
        return out;
    out.a = std::exp(log2a * log2c);
    out.r = logistic(logitR);
    double rr = out.r * (1.0 - out.r);
    if (!(out.a > 0) || !(rr > 1e-14))
        return out;
    double imag = 2.0 * out.a * std::sqrt(rr);
    double real = out.a * (1.0 - 2.0 * out.r);
    out.zPlus = {real, imag};
    out.zMinus = {real, -imag};
    out.plus = thetaPoint(out.zPlus);
    out.minus = thetaPoint(out.zMinus);
    if (!out.plus.valid || !out.minus.valid)
        return out;
    C denominator = out.zMinus - out.zPlus;
    if (std::abs(denominator) < 1e-14 * std::max(1.0, out.a))
        return out;
    out.observable = 4.0 * out.a * (out.zMinus * out.minus.R - out.zPlus * out.plus.R) / denominator;
    if (!finite(out.observable))
        return out;
    out.symmetryError = std::abs(out.minus.R - std::conj(out.plus.R)) /
                        std::max(1.0, std::abs(out.plus.R));
    out.nearZero = out.plus.nearZero || out.minus.nearZero;
    out.valid = true;
    return out;
}

} // namespace grh_observable
