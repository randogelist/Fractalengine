#pragma once
#include "theta_coefficients.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

// Five finite QQRH carriers derived from the same rounded Taylor polynomial
// T(y)=Theta_N(200*y)/Theta(0).  These are visualization constructions, not
// assertions of RH, positivity, or completeness of a numerical zero search.
namespace qqrh {
using C = std::complex<double>;

class Channel {
    struct Polynomial {
        std::array<double, 65> coefficient{};
        int degree = 0;
        bool valid = false;

        void assign(const std::array<long double, 65>& source, int top, bool rescale = true) {
            degree = top;
            while (degree > 0 && source[degree] == 0.0L) --degree;
            long double largest = 0.0L;
            for (int n = 0; n <= degree; ++n)
                largest = std::max(largest, std::abs(source[n]));
            if (!(largest > 0.0L) || !std::isfinite(largest)) return;
            for (int n = 0; n <= degree; ++n)
                coefficient[n] = static_cast<double>(source[n] / (rescale ? largest : 1.0L));
            valid = true;
        }

        // Ordinary or reversed polynomial, together with its first two jets.
        // Reversed polynomial R(s)=s^degree*P(1/s).
        void jets(C argument, bool reversed, C& value, C& first, C& second) const {
            value = C(coefficient[reversed ? 0 : degree], 0.0);
            first = second = C{};
            for (int k = degree - 1; k >= 0; --k) {
                second = second * argument + 2.0 * first;
                first = first * argument + value;
                const int index = reversed ? degree - k : k;
                value = value * argument + coefficient[index];
            }
        }

        C evaluate(C argument, bool reversed = false) const {
            C value(coefficient[reversed ? 0 : degree], 0.0);
            for (int k = degree - 1; k >= 0; --k)
                value = value * argument + coefficient[reversed ? degree - k : k];
            return value;
        }
    };

    int mode_ = 0;
    double scale_ = 1.0, shift_ = 0.0, angularAnchor_ = 0.0;
    Polynomial polynomial_;
    bool valid_ = false, centerValid_ = false;
    C centerValue_{};
    C rhoPlus_{1.0, 0.0}, rhoMinus_{1.0, 0.0};

    static bool finite(C value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    }

    static C nonfinite() {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return C(nan, nan);
    }

    C curvature(C y) const {
        C value, first, second;
        if (std::abs(y) <= 1.0) {
            polynomial_.jets(y, false, value, first, second);
            if (value == C{}) return nonfinite();
            const C logarithmic = first / value;
            // C_H=-d_y(y H'/H), with H(y)=T(-y).
            return -logarithmic - y * (second / value - logarithmic * logarithmic);
        }
        // For H(y)=y^n R(1/y), the leading monomial contributes zero:
        // C_H=-s^2*(R'/R+s*(R''/R-(R'/R)^2)), s=1/y.
        const C inverse = C(1.0, 0.0) / y;
        polynomial_.jets(inverse, true, value, first, second);
        if (value == C{}) return nonfinite();
        const C logarithmic = first / value;
        return -inverse * inverse * (logarithmic + inverse *
            (second / value - logarithmic * logarithmic));
    }

    C doublingRatio(C y) const {
        if (std::abs(y) <= .5) {
            const C denominator = polynomial_.evaluate(y);
            if (denominator == C{}) return nonfinite();
            return polynomial_.evaluate(2.0 * y) / denominator;
        }
        // T(2y)/T(y)=2^degree R(1/(2y))/R(1/y).  Evaluate the
        // scaled numerator before division to avoid tiny-ratio underflow.
        const C inverse = C(1.0, 0.0) / y;
        const C denominator = polynomial_.evaluate(inverse, true);
        if (denominator == C{}) return nonfinite();
        const C numerator = std::ldexp(1.0, polynomial_.degree) *
                            polynomial_.evaluate(.5 * inverse, true);
        return numerator / denominator;
    }

    C dyadic(C y) const {
        // Preserve the exact product definition.  Expanding P and rounding
        // its coefficients would lose the square structure at x=+/-1 and
        // can create large cancellation errors near its double zeros.
        const C plus = doublingRatio(y * rhoPlus_);
        const C minus = doublingRatio(y * rhoMinus_);
        if (!finite(plus) || !finite(minus)) return nonfinite();
        return plus * minus;
    }

    C newtonDisplacement(C y) const {
        C value, first, second;
        if (std::abs(y) <= 1.0) {
            polynomial_.jets(y, false, value, first, second);
            if (first == C{})
                return value == C{} ? C{} : nonfinite();
            return -(value / first) / scale_;
        }
        const C inverse = C(1.0, 0.0) / y;
        polynomial_.jets(inverse, true, value, first, second);
        const C denominator = static_cast<double>(polynomial_.degree) * value - inverse * first;
        if (denominator == C{})
            return value == C{} ? C{} : nonfinite();
        return -(y / scale_) * (value / denominator);
    }

public:
    // 1 curvature; 2 dyadic ratio; 3/4 Laplacian derivative levels 1/2;
    // 5 angular divided difference anchored at a=(1-m^2)/(1+m^2).
    Channel(int mode, int N, double x, double b, double h, double slope = 11.0)
        : mode_(mode), scale_(b), shift_(h) {
        if (mode < 1 || mode > 5 || (N != 8 && N != 16 && N != 32) ||
            !std::isfinite(x) || x < -1.0 || x > 1.0 ||
            !std::isfinite(b) || !(b > 0.0) || !std::isfinite(h) ||
            !std::isfinite(slope) || slope < 0.0 || slope > 64.0) return;

        std::array<long double, 65> carrier{};
        std::array<long double, 33> a{};
        for (int n = 0; n <= N; ++n)
            a[n] = static_cast<long double>(ff_rh::THETA_COEFFICIENTS[n]);
        const long double xx = static_cast<long double>(x);
        const long double m = static_cast<long double>(slope);
        const long double anchor = (1.0L - m * m) / (1.0L + m * m);
        angularAnchor_ = static_cast<double>(anchor);

        if (mode == 1) {
            for (int n = 0; n <= N; ++n) carrier[n] = (n & 1) ? -a[n] : a[n];
            polynomial_.assign(carrier, N, false);
        } else if (mode == 2) {
            for (int n = 0; n <= N; ++n) carrier[n] = a[n];
            polynomial_.assign(carrier, N, false);
            const double imaginary = std::sqrt(std::max(0.0, 1.0 - x * x));
            rhoPlus_ = C(x, imaginary); rhoMinus_ = C(x, -imaginary);
        } else {
            const int derivativeLevel = mode == 3 ? 1 : mode == 4 ? 2 : 0;
            const int order = N - derivativeLevel;
            std::array<long double, 33> coefficient{};
            for (int n = 0; n <= order; ++n) {
                long double factor = 1.0L;
                for (int j = 1; j <= derivativeLevel; ++j) factor *= n + j;
                coefficient[n] = factor * a[n + derivativeLevel];
            }
            std::array<long double, 33> U{}, atAnchor{}, V{};
            U[0] = atAnchor[0] = 1.0L;
            U[1] = 2.0L * xx; atAnchor[1] = 2.0L * anchor;
            V[0] = 0.0L; V[1] = 2.0L;
            for (int k = 1; k < order - 1; ++k) {
                U[k + 1] = 2.0L * xx * U[k] - U[k - 1];
                atAnchor[k + 1] = 2.0L * anchor * atAnchor[k] - atAnchor[k - 1];
                // (U_k(x)-U_k(a))/(x-a), including x=a, without division.
                V[k + 1] = 2.0L * xx * V[k] - V[k - 1] + 2.0L * atAnchor[k];
            }
            for (int k = mode == 5 ? 2 : 1; k <= order; ++k) {
                const long double angular = static_cast<long double>(k) *
                                             (mode == 5 ? V[k - 1] : U[k - 1]);
                for (int n = 0; n + k <= order; ++n) {
                    const int exponent = 2 * n + k - (mode == 5 ? 2 : 1);
                    carrier[exponent] += angular * coefficient[n] * coefficient[n + k];
                }
            }
            polynomial_.assign(carrier, mode == 5 ? 2 * N - 4 : 2 * order - 2);
        }
        valid_ = polynomial_.valid && polynomial_.degree > 0;
        if (!valid_) return;
        if (direct()) {
            centerValue_ = mode == 1 ? curvature(C(h, 0.0)) : dyadic(C(h, 0.0));
            centerValid_ = finite(centerValue_);
        } else {
            centerValid_ = true;
        }
    }

    // Structural validity deliberately does not reject a pole at the direct
    // centering point.  eta=0 and a zero forced outer BKQR factor must still
    // bypass that inactive host; active evaluation reports nonfinite instead.
    bool valid() const { return valid_; }
    bool centerValid() const { return centerValid_; }
    bool direct() const { return mode_ == 1 || mode_ == 2; }
    // For the factorized dyadic channel this is deg P=2N, not deg T=N.
    int degree() const { return mode_ == 2 ? 2 * polynomial_.degree : polynomial_.degree; }
    double angularAnchor() const { return angularAnchor_; }

    C response(C z) const {
        if (!valid_ || !finite(z)) return nonfinite();
        if (direct()) {
            if (!centerValid_) return nonfinite();
            if (z == C{}) return C{};
        }
        const C y = scale_ * z + shift_;
        if (!finite(y)) return nonfinite();
        const C result = mode_ == 1 ? curvature(y) - centerValue_
            : mode_ == 2 ? dyadic(y) - centerValue_ : newtonDisplacement(y);
        return finite(result) ? result : nonfinite();
    }
};

} // namespace qqrh
