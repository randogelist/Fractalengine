#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

// Finite numerical portraits of two- and three-atom spectral jet determinants.
// With q_+- = b +- d and x=z-b, D_1(x)=d^2(x^2-d^2)/(x^2+d^2)^2.
// The iterated map is F(x)=x^2+c+eta*D_1(x).  This is new dynamics,
// not a claim of a formal theorem about its Julia set.
namespace spectral {
using C = std::complex<double>;
constexpr double tau = 6.283185307179586476925286766559;
enum class Status { Unresolved, Escaped, Numerical };
struct Result {
    Status status = Status::Unresolved;
    int steps = 0;
    double smooth = 0, poleTrap = 1, phase = 0;
    C last{0, 0}; // Last representable finite point; infinity is not stored.
    bool poleHit = false;
};
inline bool finite(C z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
inline double escapeRadius(C c, C eta, double d) {
    // If r>=2d, |eta D_1(x)| <= 5|eta|/9.  For r>R below,
    // |F(x)| >= r^2/2 > r, and the subsequent orbit tends to infinity.
    return std::max({2.0, 2.0 * d, std::sqrt(2.0 * (std::abs(c) + 5.0 * std::abs(eta) / 9.0))});
}
inline Result orbit(C start, C c, C eta, double d, int n, std::vector<C> *trail = nullptr) {
    Result o;
    o.last = start;
    if (!finite(start) || !finite(c) || !finite(eta) || !std::isfinite(d) || d < 0 || n < 0) {
        o.status = Status::Numerical;
        return o;
    }
    const double t = d * d;
    const bool quadratic = d == 0 || eta == C(0, 0);
    const double radius = escapeRadius(c, quadratic ? C(0, 0) : eta, d);
    if (!std::isfinite(t) || !std::isfinite(radius) || (!quadratic && t == 0)) {
        o.status = Status::Numerical;
        return o;
    }
    C z = start, phaseSum{0, 0};
    auto finish = [&](Status status) {
        o.status = status;
        o.phase = std::abs(phaseSum) > 0 ? std::arg(phaseSum) / tau : 0;
        return o;
    };
    for (int i = 0; i <= n; ++i) {
        if (!finite(z)) return finish(Status::Numerical);
        o.last = z;
        o.steps = i;
        if (trail && trail->size() < 192) trail->push_back(z);
        const double r = std::abs(z);
        if (!std::isfinite(r)) return finish(Status::Numerical);
        if (r > 0) phaseSum += z / r;
        o.poleTrap = std::min(o.poleTrap,
            std::min(std::abs(z - C(0, d)), std::abs(z + C(0, d))) / (1.0 + r));
        if (r > radius) {
            o.smooth = std::max(0.0, i + 1.0 - std::log(std::log(std::max(r, 2.0))) / std::log(2.0));
            return finish(Status::Escaped);
        }
        if (i == n) return finish(Status::Unresolved);
        const C w = z * z, base = w + c;
        if (!finite(w) || !finite(base)) return finish(Status::Numerical);
        if (quadratic) {
            z = base; // The eta=0 and d=0 degeneracies have no poles.
            continue;
        }
        const C a = w + t;
        const double denominator = std::abs(a);
        if (a == C(0, 0)) {
            o.poleHit = true;
            o.steps = i + 1;
            o.smooth = o.steps;
            return finish(Status::Escaped); // Genuine double pole maps to infinity.
        }
        // Near a pole, only certify escape when the reverse triangle
        // inequality proves |F(x)|>R.  A small denominator alone is insufficient.
        const double numerator = std::abs(eta) * t * std::abs(w - t);
        if (denominator <= 1e-8 * (1.0 + std::abs(w) + t) &&
            numerator > (radius + std::abs(base)) * std::norm(a) * (1.0 + 1e-12)) {
            o.steps = i + 1;
            o.smooth = o.steps;
            return finish(Status::Escaped);
        }
        const C next = base + eta * t * ((w - t) / a) / a;
        if (!finite(next)) {
            const double logTerm = std::log(std::abs(eta)) + std::log(t) +
                                   std::log(std::abs(w - t)) - 2.0 * std::log(denominator);
            if (logTerm > std::log(radius + std::abs(base)) + 1e-12) {
                o.steps = i + 1;
                o.smooth = o.steps;
                return finish(Status::Escaped);
            }
            return finish(Status::Numerical); // Never colour a NaN as escape.
        }
        z = next;
    }
    return finish(Status::Unresolved);
}
// Generation 2: three positive atoms q_i=1/2+t_i, with
// t_i=d*(-1,.35,1), w_i=((1-m)/2,m,(1-m)/2), and jet order s=1..3.
// For B_s=sum_i w_i product_{j!=i}(x-t_j)^(s+1),
// N_s=sum_{i<j}w_i*w_j*(t_i-t_j)^2*(x-t_i)^s*(x-t_j)^s*(x-t_k)^(2s+2),
// D_s=N_s/B_s^2.  This cancels the removable atom poles before evaluation.
// The scaled coordinates below cancel the same common power of L from N and B^2.
struct Gen2Kernel {
    C numerator{0, 0}, denominator{1, 0};
    double denominatorScale = 1;
};
inline Gen2Kernel gen2Kernel(C x, double d, double m, int s) {
    if (d == 0) return {};
    const double length = std::max(std::abs(x), d);
    const double beta = d / length;
    if (beta == 0 || beta * beta == 0) return {C(0, 0), C(0, 0), 0};
    const C v = x / length;
    const std::array<C, 3> u{{v + beta, v - .35 * beta, v - beta}};
    std::array<C, 3> power{{u[0], u[1], u[2]}}, next;
    for (int k = 1; k < s; ++k)
        for (int j = 0; j < 3; ++j) power[j] *= u[j];
    for (int j = 0; j < 3; ++j) next[j] = power[j] * u[j];
    const double w = (1.0 - m) * .5;
    const C b0 = w * next[1] * next[2];
    const C b1 = m * next[0] * next[2];
    const C b2 = w * next[0] * next[1];
    Gen2Kernel out;
    out.denominator = b0 + b1 + b2;
    out.denominatorScale = std::abs(b0) + std::abs(b1) + std::abs(b2);
    out.numerator = beta * beta *
        (w * m * (1.35 * 1.35) * power[0] * power[1] * next[2] * next[2] +
         w * w * 4.0 * power[0] * power[2] * next[1] * next[1] +
         m * w * (.65 * .65) * power[1] * power[2] * next[0] * next[0]);
    return out;
}
inline double escapeRadiusGen2(C c, C eta, double d) {
    // For r>=8d, s<=3 and positive unit-sum weights:
    // |D_s| <= 16777216/24423364 < 1.  Indeed the scaled R_(s+1)
    // has modulus at least 2-(8/7)^4, and the pair numerator uses
    // sum_{i<j}w_i*w_j*(t_i-t_j)^2 = Var(t) <= d^2.
    // Beyond R, |F(x)| >= r^2-|c|-|eta| >= r^2/2 > r.
    return std::max({2.0, 8.0 * d, std::sqrt(2.0 * (std::abs(c) + std::abs(eta)))});
}
inline Result orbitGen2(C start, C c, C eta, double d, double m, int s, int n,
                        std::vector<C> *trail = nullptr) {
    Result o;
    o.last = start;
    if (!finite(start) || !finite(c) || !finite(eta) || !std::isfinite(d) ||
        d < 0 || d > .49 || !std::isfinite(m) || m < .01 || m > .98 ||
        s < 1 || s > 3 || n < 0) {
        o.status = Status::Numerical;
        return o;
    }
    const bool quadratic = d == 0 || eta == C(0, 0);
    const double radius = quadratic ? escapeRadius(c, C(0, 0), 0)
                                    : escapeRadiusGen2(c, eta, d);
    if (!std::isfinite(radius)) {
        o.status = Status::Numerical;
        return o;
    }
    C z = start, phaseSum{0, 0};
    auto finish = [&](Status status) {
        o.status = status;
        o.phase = std::abs(phaseSum) > 0 ? std::arg(phaseSum) / tau : 0;
        return o;
    };
    for (int i = 0; i <= n; ++i) {
        if (!finite(z)) return finish(Status::Numerical);
        o.last = z;
        o.steps = i;
        if (trail && trail->size() < 192) trail->push_back(z);
        const double r = std::abs(z);
        if (!std::isfinite(r)) return finish(Status::Numerical);
        if (r > 0) phaseSum += z / r;
        if (r > radius) {
            o.smooth = std::max(0.0, i + 1.0 - std::log(std::log(std::max(r, 2.0))) / std::log(2.0));
            return finish(Status::Escaped);
        }
        // The denominator trap measures cancellation, not distance to atoms.
        const auto kernel = quadratic ? Gen2Kernel{} : gen2Kernel(z, d, m, s);
        if (!finite(kernel.numerator) || !finite(kernel.denominator) ||
            !std::isfinite(kernel.denominatorScale) || !(kernel.denominatorScale > 0))
            return finish(Status::Numerical);
        const double denominator = std::abs(kernel.denominator);
        o.poleTrap = std::min(o.poleTrap,
            std::clamp(denominator / kernel.denominatorScale, 0.0, 1.0));
        if (i == n) return finish(Status::Unresolved);
        const C base = z * z + c;
        if (!finite(base)) return finish(Status::Numerical);
        if (quadratic) {
            z = base;
            continue;
        }
        if (kernel.denominator == C(0, 0)) {
            if (kernel.numerator == C(0, 0)) return finish(Status::Numerical);
            o.poleHit = true;
            o.steps = i + 1;
            o.smooth = o.steps;
            return finish(Status::Escaped);
        }
        const double numerator = std::abs(eta) * std::abs(kernel.numerator);
        if (denominator <= 1e-8 * kernel.denominatorScale &&
            numerator > (radius + std::abs(base)) * std::norm(kernel.denominator) * (1.0 + 1e-12)) {
            o.steps = i + 1;
            o.smooth = o.steps;
            return finish(Status::Escaped);
        }
        const C next = base + eta * (kernel.numerator / kernel.denominator) / kernel.denominator;
        if (!finite(next)) {
            const double logTerm = std::log(std::abs(eta)) + std::log(std::abs(kernel.numerator)) -
                                   2.0 * std::log(denominator);
            if (logTerm > std::log(radius + std::abs(base)) + 1e-12) {
                o.steps = i + 1;
                o.smooth = o.steps;
                return finish(Status::Escaped);
            }
            return finish(Status::Numerical);
        }
        z = next;
    }
    return finish(Status::Unresolved);
}
inline std::array<double, 3> colour(const Result &o, int palette, double hue, double bands) {
    if (o.status == Status::Numerical) return {.35, .02, .20};
    const bool escaped = o.status == Status::Escaped;
    if (palette == 0 && !escaped) return {0, 0, 0};
    const double depth = -std::log(std::max(o.poleTrap, 1e-12));
    double h = hue + .08 + .075 * bands * o.smooth;
    double light = .45 + .50 * (1.0 - std::exp(-.045 * o.smooth));
    if (palette == 1) {
        h = hue + .18 + .12 * bands * depth;
        light = .20 + .70 * (1.0 - std::exp(-.22 * depth));
    } else if (palette == 2) {
        h = hue + bands * o.phase + .025 * bands * std::log1p(o.steps);
        light = escaped ? .40 + .50 * (1.0 - std::exp(-.06 * o.smooth))
                        : .26 + .38 * (1.0 - std::exp(-.18 * depth));
    }
    return {light * (.52 + .46 * std::cos(tau * h)),
            light * (.52 + .46 * std::cos(tau * (h + .33))),
            light * (.52 + .46 * std::cos(tau * (h + .67)))};
}
} // namespace spectral
