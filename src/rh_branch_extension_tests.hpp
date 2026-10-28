#pragma once

#include "rh_master_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

// Independent native algebra checks for the three branch extensions. This
// oracle builds E in ordinary powers from its defining Laurent product;
// it does not call the production Chebyshev/Clenshaw, jet, or quotient code.
namespace rh_branch_extension_tests {
namespace oracle {
using Z = std::complex<long double>;
using Polynomial = std::vector<long double>;

inline bool finite(Z z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

inline Polynomial angular(int degree, long double radius, int radialOrder) {
    // S_k(x)=w^k+w^(-k), 2x=w+w^(-1). Pairing each (i,j)
    // with (j,i) in T(Rw)T(R/w) removes w without a square root.
    std::vector<Polynomial> symmetric(std::size_t(degree) + 1);
    symmetric[0] = Polynomial{2.L};
    symmetric[1] = Polynomial{0.L, 2.L};
    for (int k = 1; k < degree; ++k) {
        auto& next = symmetric[std::size_t(k + 1)];
        next.assign(std::size_t(k) + 2, 0.L);
        for (std::size_t j = 0; j < symmetric[std::size_t(k)].size(); ++j)
            next[j + 1] += 2.L * symmetric[std::size_t(k)][j];
        for (std::size_t j = 0; j < symmetric[std::size_t(k - 1)].size(); ++j)
            next[j] -= symmetric[std::size_t(k - 1)][j];
    }
    std::vector<long double> powers(std::size_t(2 * degree) + 1, 1.L);
    for (std::size_t k = 1; k < powers.size(); ++k) powers[k] = powers[k - 1] * radius;
    Polynomial result(std::size_t(degree) + 1, 0.L);
    for (int i = 0; i <= degree; ++i)
        for (int j = 0; j <= i; ++j) {
            const int exponent = i + j;
            if (exponent < radialOrder) continue;
            long double amplitude = static_cast<long double>(ff_rh::THETA_COEFFICIENTS[i])
                * static_cast<long double>(ff_rh::THETA_COEFFICIENTS[j])
                * powers[std::size_t(exponent - radialOrder)];
            for (int derivative = 0; derivative < radialOrder; ++derivative)
                amplitude *= exponent - derivative;
            if (i == j) result[0] += amplitude;
            else {
                const auto& pair = symmetric[std::size_t(i - j)];
                for (std::size_t k = 0; k < pair.size(); ++k) result[k] += amplitude * pair[k];
            }
        }
    return result;
}

inline Z derivative(const Polynomial& polynomial, Z y, int order) {
    Z value{};
    for (int k = int(polynomial.size()) - 1; k >= order; --k) {
        long double coefficient = polynomial[std::size_t(k)];
        for (int j = 0; j < order; ++j) coefficient *= k - j;
        value = value * y + coefficient;
    }
    return value;
}

inline std::array<Z, 4> theta(int degree, Z y) {
    std::array<Z, 4> derivatives{};
    // Evaluate independently differentiated ordinary polynomials, rather
    // than the production reversed jet recurrence.
    for (int order = 0; order < 4; ++order)
        for (int k = degree; k >= order; --k) {
            long double coefficient = ff_rh::THETA_COEFFICIENTS[k];
            for (int j = 0; j < order; ++j) coefficient *= k - j;
            derivatives[std::size_t(order)] = derivatives[std::size_t(order)] * y + coefficient;
        }
    return derivatives;
}

inline std::array<Z, 2> count(int degree, Z y) {
    const auto t = theta(degree, y);
    const Z logarithmic = t[1] / t[0];
    return {{y * logarithmic,
        logarithmic + y * (t[2] / t[0] - logarithmic * logarithmic)}};
}

inline Z correction(int candidate, int degree, const Polynomial& e,
                    const Polynomial& er, Z y, long double level) {
    if (candidate == 11) {
        const auto one = count(degree, y), two = count(degree, 2.L * y);
        return (two[0] - one[0] - level) / (2.L * two[1] - one[1]);
    }
    if (candidate == 0) {
        const auto t = theta(degree, y);
        const Z a = t[1] / t[0], b = t[2] / t[0], c = t[3] / t[0];
        const Z second = b - a * a;
        const Z third = c - 3.L * a * b + 2.L * a * a * a;
        return (a + y * second) / (2.L * second + y * third);
    }
    const Z E = derivative(e, y, 0), U = derivative(e, y, 1);
    const Z Ux = derivative(e, y, 2), Uxx = derivative(e, y, 3);
    if (candidate == 1) return (U - level) / Ux;
    if (candidate == 2) return (U / E - level) / (Ux / E - (U / E) * (U / E));
    if (candidate == 10) {
        // F=(log E)''-level, F'=(log E)'''.
        const Z slope = U / E;
        const Z second = Ux / E - slope * slope;
        const Z third = Uxx / E - 3.L * slope * Ux / E + 2.L * slope * slope * slope;
        return (second - level) / third;
    }
    if (candidate == 12) {
        const Z Ur = derivative(er, y, 1), Urx = derivative(er, y, 2);
        // The target is Ur/Ux. Differentiating that ratio here fixes its
        // orientation independently of the production factored quotient.
        const Z ratio = Ur / Ux;
        const Z ratioDerivative = Urx / Ux - ratio * Uxx / Ux;
        return (ratio - level) / ratioDerivative;
    }
    const long double nan = std::numeric_limits<long double>::quiet_NaN();
    return Z(nan, nan);
}
} // namespace oracle

inline bool run(std::ostream& out) {
    try {
        int cases = 0, failures = 0;
        long double worstRelative = 0.L;
        constexpr long double tolerance = 2.e-9L;
        const std::array<int, 6> candidates{{0, 1, 2, 10, 11, 12}};
        const std::array<rhfresh::C, 5> points{{
            {.31, .27}, {-.4, .19}, {.55, .19}, {1.25, .23}, {4.2, .13}
        }};
        auto check = [&](bool good, const std::string& label, int candidate, int degree,
                         double radius, double level, rhfresh::C point, long double error,
                         long double allowed = 2.e-9L) {
            ++cases;
            worstRelative = std::max(worstRelative, error);
            if (!good) {
                ++failures;
                out << "RH_BRANCH_EXTENSION_TEST_FAIL=" << label
                    << " CANDIDATE=" << candidate + 1 << " DEGREE=" << degree
                    << " RADIUS=" << radius << " LEVEL=" << level
                    << " Z_REAL=" << point.real() << " Z_IMAG=" << point.imag()
                    << " RELATIVE_ERROR=" << error << " TOLERANCE=" << allowed << '\n';
            }
        };
        for (int degree : {8, 16, 32})
            for (double radius : {.35, .6, .9}) {
                const auto e = oracle::angular(degree, radius, 0);
                const auto er = oracle::angular(degree, radius, 1);
                for (int candidate : candidates)
                    for (double level : {0., .125, -.25}) {
                        auto settings = rhfresh::defaultSettings(candidate);
                        settings.degree = degree;
                        settings.radius = radius;
                        settings.level = level;
                        settings.scale = .875;
                        settings.shift = .125;
                        const rhfresh::Model model(settings);
                        const rhmaster::CandidateChannel channel(candidate + 1, degree, radius,
                            settings.scale, settings.shift, settings.lambda, level);
                        for (auto point : points) {
                            const oracle::Z z(point.real(), point.imag());
                            const oracle::Z y = static_cast<long double>(settings.scale) * z
                                + static_cast<long double>(settings.shift);
                            const oracle::Z expected = oracle::correction(candidate, degree, e, er, y, level)
                                / static_cast<long double>(settings.scale);
                            const auto actual = model.correction(point);
                            const oracle::Z got(actual.real(), actual.imag());
                            const long double error = std::abs(got - expected) / (1.L + std::abs(expected));
                            check(model.valid() && channel.valid() && oracle::finite(expected)
                                && oracle::finite(got) && error <= tolerance,
                                "INDEPENDENT_POWER_ORACLE", candidate, degree, radius, level, point, error);

                            const auto reflected = model.correction(std::conj(point));
                            const long double conjugateError = std::abs(
                                oracle::Z(reflected.real(), reflected.imag()) - std::conj(got))
                                / (1.L + std::abs(got));
                            check(oracle::finite(oracle::Z(reflected.real(), reflected.imag()))
                                && conjugateError <= tolerance,
                                "REAL_COEFFICIENT_CONJUGATION", candidate, degree, radius, level,
                                point, conjugateError);

                            const auto response = channel.response(point);
                            const long double responseError = std::abs(
                                oracle::Z(response.real(), response.imag()) + got) / (1.L + std::abs(got));
                            check(oracle::finite(oracle::Z(response.real(), response.imag()))
                                && responseError <= tolerance,
                                "NEGATIVE_RESPONSE_AND_AFFINE_SCALE", candidate, degree, radius, level,
                                point, responseError);
                        }
                    }
            }
        // A dilation gap can use different representations at y and 2y.
        // These samples cross both switch points without using the
        // production inverse-polynomial identity in the oracle.
        const std::array<rhfresh::C, 9> transitionPoints{{
            {0., 0.}, {.49, 0.}, {.5, 0.}, {.500001, 0.}, {.75, 0.},
            {.999999, 0.}, {1., 0.}, {1.000001, 0.}, {1., .2}
        }};
        for (int degree : {8, 16, 32})
            for (double level : {0., .125, -.25}) {
                auto settings = rhfresh::defaultSettings(11);
                settings.degree = degree;
                settings.level = level;
                settings.scale = 1.;
                settings.shift = 0.;
                const rhfresh::Model model(settings);
                for (auto point : transitionPoints) {
                    const oracle::Z y(point.real(), point.imag());
                    const auto expected = oracle::correction(11, degree, {}, {}, y, level);
                    const auto result = model.correction(point);
                    const oracle::Z actual(result.real(), result.imag());
                    const long double error = std::abs(actual - expected) / (1.L + std::abs(expected));
                    check(model.valid() && oracle::finite(expected) && oracle::finite(actual)
                        && error <= tolerance, "DILATION_MIXED_REPRESENTATION_BOUNDARIES",
                        11, degree, settings.radius, level, point, error);
                }
                // If T(y)=a_N*y^N*(1+A/y+O(y^-2)), then
                // C(2y)-C(y)=A/(2y)+O(y^-2). Therefore the zero-level
                // Newton quotient is -y+O(1), and a nonzero level gives
                // (2*level/A)*y^2+O(y). These independent asymptotics catch
                // loss of the small difference between the two counts.
                const long double A = static_cast<long double>(ff_rh::THETA_COEFFICIENTS[degree - 1])
                    / static_cast<long double>(ff_rh::THETA_COEFFICIENTS[degree]);
                for (double coordinate : {1.e8, 1.e30, 1.e100}) {
                    const long double y = coordinate;
                    const long double expected = level == 0. ? -y
                        : (2.L * static_cast<long double>(level) / A) * y * y;
                    const auto result = model.correction(rhfresh::C(coordinate, 0.));
                    const oracle::Z actual(result.real(), result.imag());
                    const long double error = std::abs(actual - expected) / (1.L + std::abs(expected));
                    const long double asymptoticTolerance = coordinate == 1.e8 ? 1.e-5L : tolerance;
                    check(model.valid() && oracle::finite(actual) && error <= asymptoticTolerance,
                        "DILATION_LARGE_ARGUMENT_ASYMPTOTIC", 11, degree, settings.radius,
                        level, rhfresh::C(coordinate, 0.), error, asymptoticTolerance);
                }
            }
        out << "RH_BRANCH_EXTENSION_TEST_CASES=" << cases << " FAILURES=" << failures
            << " MAX_RELATIVE_ERROR=" << worstRelative << " TOLERANCE=" << tolerance << '\n';
        out << "RH_BRANCH_EXTENSION_TEST=" << (failures == 0 ? "PASS" : "FAIL") << '\n';
        return failures == 0;
    } catch (const std::exception& error) {
        out << "RH_BRANCH_EXTENSION_TEST=FAIL STAGE=EXCEPTION REASON=" << error.what() << '\n';
        return false;
    } catch (...) {
        out << "RH_BRANCH_EXTENSION_TEST=FAIL STAGE=EXCEPTION REASON=Unknown exception\n";
        return false;
    }
}
} // namespace rh_branch_extension_tests
