#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

#include "theta_coefficients.hpp"

// Finite Taylor realization of the three RH-branch channels.
// The coefficient table is for T(y) = Theta_N(200*y)/Theta(0).
// The common constant normalization has no effect on F/F'.
namespace ff_rh {

class Channel {
public:
    using Complex = std::complex<double>;

    // mode: 0 = Master, 1 = Sub_1, 2 = Sub_2.
    // N is the Taylor degree, x is the angular coordinate in [-1,1].
    // The complex radius used by the renderer is y = b*z + h.
    Channel(int mode = 0, int N = 16, double x = 0.0,
            double b = 1.0, double h = 0.0)
        : scale_(b), shift_(h) {
        if (mode < 0 || mode > 2 || N < 2 || N > 32 ||
            !std::isfinite(x) || x < -1.0 || x > 1.0 ||
            !std::isfinite(b) || b == 0.0 || !std::isfinite(h)) {
            return;
        }

        // U_m is Chebyshev of the second kind.  V_m is its polynomial
        // divided difference (U_m(x)-U_m(-1))/(x+1).  Constructing V
        // directly avoids cancellation and works exactly at x=-1.
        std::array<long double, 32> U{};
        std::array<long double, 32> V{};
        const long double xx = static_cast<long double>(x);
        U[0] = 1.0L;
        U[1] = 2.0L * xx;
        V[0] = 0.0L;
        V[1] = 2.0L;
        for (int m = 1; m < N - 1; ++m) {
            U[m + 1] = 2.0L * xx * U[m] - U[m - 1];
            const long double boundary = (m & 1) ? -(m + 1.0L)
                                                       : (m + 1.0L);
            V[m + 1] = 2.0L * xx * V[m] - V[m - 1]
                       + 2.0L * boundary;
        }

        std::array<long double, 63> accum{};
        // c_k(y) = sum_n a_n*a_(n+k)*y^(2n+k).
        // A = sum_k k*c_k*U_(k-1), B = A(y,-1).
        // Remove universal factors before forming Newton corrections:
        // Master F=A/y, Sub_1 F=(A-B)/(y^2*(x+1)), Sub_2 F=B/y.
        for (int k = 1; k <= N; ++k) {
            if (mode == 1 && k == 1) continue;
            long double angular = 0.0L;
            if (mode == 0) {
                angular = static_cast<long double>(k) * U[k - 1];
            } else if (mode == 1) {
                angular = static_cast<long double>(k) * V[k - 1];
            } else {
                angular = static_cast<long double>(k) * k;
                if ((k & 1) == 0) angular = -angular;
            }
            for (int n = 0; n + k <= N; ++n) {
                const int exponent = 2 * n + k - ((mode == 1) ? 2 : 1);
                accum[exponent] += angular
                    * static_cast<long double>(THETA_COEFFICIENTS[n])
                    * static_cast<long double>(THETA_COEFFICIENTS[n + k]);
            }
        }

        degree_ = (mode == 1) ? 2 * N - 4 : 2 * N - 2;
        while (degree_ > 0 && accum[degree_] == 0.0L) --degree_;
        long double largest = 0.0L;
        for (int j = 0; j <= degree_; ++j)
            largest = std::max(largest, std::abs(accum[j]));
        if (!(largest > 0.0L) || !std::isfinite(largest)) return;
        for (int j = 0; j <= degree_; ++j)
            coefficients_[j] = static_cast<double>(accum[j] / largest);
        valid_ = degree_ > 0;
    }

    bool valid() const { return valid_; }
    int degree() const { return degree_; }

    // -F(y)/(b*F_y'(y)), the Newton displacement in the z coordinate.
    // The caller supplies its existing outer map z^p+c+lambda*response(z).
    // No clipping or denominator regularization changes the rational map.
    Complex response(Complex z) const {
        if (!valid_) return nonfinite();
        const Complex y = scale_ * z + shift_;
        if (!finite(y)) return nonfinite();

        Complex value;
        Complex derivative;
        Complex ratio;
        if (std::abs(y) <= 1.0) {
            value = Complex(coefficients_[degree_], 0.0);
            derivative = Complex(0.0, 0.0);
            for (int j = degree_ - 1; j >= 0; --j) {
                derivative = derivative * y + value;
                value = value * y + coefficients_[j];
            }
            if (derivative == Complex(0.0, 0.0))
                return value == Complex(0.0, 0.0) ? Complex(0.0, 0.0)
                                                  : nonfinite();
            ratio = value / derivative;
        } else {
            // Reversed Horner: F/F' = y*R(1/y)/(d*R(1/y)-(1/y)*R'(1/y)).
            // Avoids forming y^degree for radii outside the unit disk.
            const Complex inverse = Complex(1.0, 0.0) / y;
            value = Complex(coefficients_[0], 0.0);
            derivative = Complex(0.0, 0.0);
            for (int j = 1; j <= degree_; ++j) {
                derivative = derivative * inverse + value;
                value = value * inverse + coefficients_[j];
            }
            const Complex denominator = static_cast<double>(degree_) * value
                                        - inverse * derivative;
            if (denominator == Complex(0.0, 0.0))
                return value == Complex(0.0, 0.0) ? Complex(0.0, 0.0)
                                                  : nonfinite();
            ratio = y * (value / denominator);
        }
        return -ratio / scale_;
    }

private:
    static bool finite(Complex z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }

    static Complex nonfinite() {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Complex(nan, nan);
    }

    std::array<double, 63> coefficients_{};
    int degree_ = 0;
    double scale_ = 1.0;
    double shift_ = 0.0;
    bool valid_ = false;
};

} // namespace ff_rh
