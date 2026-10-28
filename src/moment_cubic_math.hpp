#pragma once
#include "riemann_witness_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// A direct polynomial iteration derived from the five fixed RH point moments.
// Let b be standardized skewness and h=E[X^4]-b^2-1 for E[X]=0,
// E[X^2]=1. The existing Jacobi closure has characteristic polynomial
// P(w)=w^3-b*w^2-(1+h)*w. Substitution w=z+b/3 gives
// P(z+b/3)=z^3-A*z-2*b^3/27-b*(1+h)/3, A=1+h+b^2/3.
// Keeping A and allowing the constant to vary defines the NEW dynamical family
//                  f_c(z)=z^3-A*z+c.
// This is direct iteration, not Newton/Halley iteration and not a claim of
// dynamical conjugacy to the unshifted characteristic-polynomial map.
// Its two critical points are +/-sqrt(A/3), since f'_c(z)=3*z^2-A.
//
// For r > R=max(2,sqrt(2*(A+1)),cbrt(2*|c|)),
// |f_c(z)| >= r^3-A*r-|c| >= r*(r^2/2-A) > r.
// Iteration then diverges: the last multiplicative lower bound increases
// with r and already exceeds one. Runtime inflates R by 1+1e-12 to avoid
// classifying a rounded equality as a strict inequality. This is a numerical
// application of the analytic escape criterion, NOT interval certification.
// Finite-budget nonescape remains unresolved; it never certifies boundedness,
// connectedness, attracting cycles, RH, or any infinite moment condition.
namespace momentcubic {
using C = std::complex<double>;
using Color = std::array<double, 3>;

enum class Status { Escaped, Unresolved, NumericalLimit, InvalidModel };
inline const char *statusName(Status status) {
    switch (status) {
    case Status::Escaped: return "numerical polynomial escape criterion reached";
    case Status::Unresolved: return "unresolved within iteration budget";
    case Status::NumericalLimit: return "floating-point evaluation limit (not escape)";
    default: return "invalid model or parameter (finite |c| <= 1e6 required)";
    }
}

struct Result {
    Status status = Status::Unresolved;
    C endpoint{};
    int iterations = 0;
    double smooth = 0, structure = 0, phase = 0;
};

class Model {
    static bool finite(C z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }
    static bool validParameter(C c) {
        return finite(c) && std::abs(c) <= 1e6;
    }
    static double safeHue(double hue) { return std::isfinite(hue) ? hue : 0; }
    static double safeBands(double bands) {
        return std::isfinite(bands) ? std::clamp(bands, 0.0, 100.0) : 1;
    }
    static Color errorColour(Status status) {
        return status == Status::InvalidModel ? Color{{.42, .035, .055}}
                                              : Color{{.83, .12, .72}};
    }
    static Color unresolvedColour(double structure, double hue) {
        // Very dark finite-budget nonescape. No basin/period is inferred.
        return projects::hsv(hue + .61, .48, .019 + .017 * structure);
    }

public:
    std::shared_ptr<const riemann::Model> source;
    double A = 0, critical = 0;
    bool validModel = false;
    std::string description, statusText;

    Model() : source(riemann::model(6, riemann::Settings{})) {
        description = "Riemann point-moment centered cubic: direct iteration of z^3-A*z+c";
        if (!source || !source->validModel) {
            statusText = "The fixed H2 point-moment source is unavailable";
            return;
        }
        const double b = source->standardizedSkew, h = source->standardizedH;
        A = 1 + h + b * b / 3;
        if (!std::isfinite(b) || !(h > 0) || !std::isfinite(h) ||
            !(A > 0) || !std::isfinite(A)) {
            statusText = "The point-moment cubic coefficient is invalid";
            return;
        }
        critical = std::sqrt(A / 3);
        if (!(critical > 0) || !std::isfinite(critical)) {
            statusText = "The point-moment cubic critical points are invalid";
            return;
        }
        validModel = true;
        statusText = "Fixed rounded H2 point moments; direct cubic iteration with two critical orbits. "
            "Numerical escape colours; finite-budget nonescape is unresolved, not an RH proof.";
    }

    double escapeRadius(C parameter) const {
        if (!validModel || !validParameter(parameter))
            return std::numeric_limits<double>::quiet_NaN();
        return (1 + 1e-12) * std::max({2.0, std::sqrt(2 * (A + 1)),
                                     std::cbrt(2 * std::abs(parameter))});
    }

    Result orbit(C initial, C parameter, int budget,
                 std::vector<C> *trail = nullptr) const {
        Result result;
        result.endpoint = initial;
        if (trail) trail->clear();
        if (!validModel || !validParameter(parameter)) {
            result.status = Status::InvalidModel;
            return result;
        }
        const double radius = escapeRadius(parameter);
        if (!std::isfinite(radius)) {
            result.status = Status::NumericalLimit;
            return result;
        }
        budget = std::clamp(budget, 32, 1024);
        C z = initial;
        bool trailOpen = finite(z) && std::abs(z) < 1e12;
        if (trail && trailOpen) trail->push_back(z);
        double structureSum = 0, weightSum = 0, sampleWeight = 1;
        // n counts completed updates. Including n==budget ensures that an
        // escape caused by the final permitted update is not lost.
        for (int n = 0; n <= budget; ++n) {
            result.endpoint = z;
            result.iterations = n;
            result.smooth = static_cast<double>(n);
            if (!finite(z)) { result.status = Status::NumericalLimit; break; }
            const double r = std::abs(z);
            if (!std::isfinite(r)) { result.status = Status::NumericalLimit; break; }
            // A bounded, weighted orbit-shape observable; it is only colour
            // data and never contributes to dynamical classification.
            const double shape = .55 * (r / (1 + r)) +
                                 .45 * (std::abs(z.imag()) / (1 + r));
            structureSum += sampleWeight * shape;
            weightSum += sampleWeight;
            sampleWeight *= .97;
            if (r > radius) {
                result.status = Status::Escaped;
                // Degree-three potential smoothing relative to the same R.
                // A colour coordinate, not an extra orbit or proof step.
                result.smooth = std::max(0.0, n + 1.0 -
                    std::log(std::log(r) / std::log(radius)) / std::log(3.0));
                break;
            }
            if (n == budget) break;
            // The escape check precedes cubing, keeping evaluated iterates
            // small for the supported parameter range |c|<=1e6.
            z = z * z * z - A * z + parameter;
            if (!finite(z) || std::abs(z) >= 1e12) trailOpen = false;
            if (trail && trailOpen) trail->push_back(z);
        }
        result.structure = weightSum > 0 ? structureSum / weightSum : 0;
        if (finite(result.endpoint))
            result.phase = .5 + std::arg(result.endpoint) / (2 * projects::pi);
        return result;
    }

    Color juliaColour(C initial, C parameter, int budget, int palette,
                      double hue, double bands) const {
        const Result result = orbit(initial, parameter, budget);
        hue = safeHue(hue); bands = safeBands(bands);
        if (result.status == Status::InvalidModel || result.status == Status::NumericalLimit)
            return errorColour(result.status);
        if (result.status == Status::Unresolved)
            return unresolvedColour(result.structure, hue);
        const double stripes = .5 + .5 * std::cos(bands *
            (.57 * result.smooth + 2.4 * result.structure));
        if (palette == 1)
            return projects::hsv(hue + .66 - .023 * result.smooth, .76,
                                 .24 + .74 * stripes);
        if (palette == 2)
            return projects::hsv(hue + .65 * result.structure + .18 * result.phase,
                                 .73, .26 + .72 * stripes);
        return projects::hsv(hue + .57 + .30 * result.phase + .012 * result.smooth,
                             .74, .26 + .72 * stripes);
    }

    Color atlasColour(C parameter, int budget, int palette,
                      double hue, double bands) const {
        // A cubic has two critical points. Inspect both rather than using a
        // single preferred start as a surrogate for its parameter portrait.
        const Result plus = orbit(C(critical, 0), parameter, budget);
        const Result minus = orbit(C(-critical, 0), parameter, budget);
        hue = safeHue(hue); bands = safeBands(bands);
        if (plus.status == Status::InvalidModel || minus.status == Status::InvalidModel)
            return errorColour(Status::InvalidModel);
        if (plus.status == Status::NumericalLimit || minus.status == Status::NumericalLimit)
            return errorColour(Status::NumericalLimit);
        const bool plusEscaped = plus.status == Status::Escaped;
        const bool minusEscaped = minus.status == Status::Escaped;
        const double structure = (plus.structure + minus.structure) / 2;
        if (!plusEscaped && !minusEscaped)
            return unresolvedColour(structure, hue);
        const bool bothEscaped = plusEscaped && minusEscaped;
        // Orange: +critical escaped; cyan: -critical escaped; violet: both.
        // "Only" here means only observed within this budget. The other
        // critical orbit remains unresolved, never certified bounded.
        const double fateHue = bothEscaped ? .80 : plusEscaped ? .065 : .54;
        const double time = bothEscaped ? (plus.smooth + minus.smooth) / 2
                                       : plusEscaped ? plus.smooth : minus.smooth;
        const double difference = bothEscaped ? std::abs(plus.smooth - minus.smooth) : 0;
        const double stripes = .5 + .5 * std::cos(bands *
            (.57 * time + .13 * difference + 2.4 * structure));
        if (palette == 1)
            return projects::hsv(hue + fateHue + .07 * std::sin(.24 * time),
                                 .73, .21 + .77 * stripes);
        if (palette == 2)
            return projects::hsv(hue + fateHue + .13 * (structure - .5),
                                 .62 + .30 * structure, .26 + .72 * stripes);
        return projects::hsv(hue + fateHue, .74, .28 + .70 * stripes);
    }
};

inline std::shared_ptr<const Model> model() {
    static const std::shared_ptr<const Model> value = std::make_shared<const Model>();
    return value;
}
} // namespace momentcubic
