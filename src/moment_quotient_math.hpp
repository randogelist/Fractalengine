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

// The five fixed RH point moments give E[X]=0, E[X^2]=1, standardized
// skewness b=E[X^3], and h=E[X^4]-b^2-1. The associated Jacobi recurrence
// begins p_1(z)=z, p_2(z)=z^2-b*z-1. Its quotient p_2/p_1 is therefore
// z-b-1/z. This NEW degree-two rational family directly iterates
//                  f_c(z)=lambda*(z-b-1/z)+c,
// with lambda=1/sqrt(h). It is not a Newton/Halley iteration. We require
// 0<h<1 rather than silently substituting some other multiplier or map.
// Its critical points are +/-i, since f'_c(z)=lambda*(1+1/z^2).
// Infinity is a fixed point with local multiplier 1/lambda in w=1/z.
//
// Set d=c-lambda*b. If
// R=max(1,4*|d|/(lambda-1),sqrt(4*lambda/(lambda-1))) and r>R, then
// |f_c(z)| >= lambda*r-|d|-lambda/r > (lambda+1)*r/2 > r.
// This proves escape to infinity once the strict bound holds. The runtime
// inflates R by 1+1e-12 and evaluates that bound in floating point; its
// diagnostics are NOT interval certificates. Nonescape at a finite budget
// is unresolved. No polynomial connectedness criterion is applied here.
//
// To include the pole without unstable affine division, iterate on P^1:
// [Z:W] -> [lambda*Z^2+d*Z*W-lambda*W^2 : Z*W].
// A common positive coefficient scale and a common coordinate normalization
// do not change the projective point. In particular, 0 -> infinity -> infinity;
// encountering the pole is a valid orbit event, not a numerical failure.
namespace momentquotient {
using C = std::complex<double>;
using Color = std::array<double, 3>;

enum class Status { Escaped, Unresolved, NumericalLimit, InvalidModel };
inline const char *statusName(Status status) {
    switch (status) {
    case Status::Escaped: return "numerical rational escape criterion reached";
    case Status::Unresolved: return "unresolved within iteration budget";
    case Status::NumericalLimit: return "floating-point evaluation limit (not escape)";
    default: return "invalid model or parameter (finite |c| <= 1e6 required)";
    }
}

struct Result {
    Status status = Status::Unresolved;
    C endpoint{}, endpointZ{}, endpointW{1, 0};
    bool atInfinity = false, endpointRepresentable = true;
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
        return projects::hsv(hue + .61, .48, .019 + .017 * structure);
    }
    static void storeEndpoint(Result &result, C Z, C W) {
        result.endpointZ = Z;
        result.endpointW = W;
        result.atInfinity = W == C{};
        result.endpointRepresentable = false;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result.endpoint = C(nan, nan);
        if (result.atInfinity || !finite(Z) || !finite(W)) return;
        // Only the display endpoint uses affine division. It cannot change
        // classification or feed back into the homogeneous orbit. Wider
        // intermediates retain representable components near a small W.
        using L = std::complex<long double>;
        const L point = L(Z.real(), Z.imag()) / L(W.real(), W.imag());
        const long double limit = std::numeric_limits<double>::max();
        if (!std::isfinite(point.real()) || !std::isfinite(point.imag()) ||
            std::abs(point.real()) > limit || std::abs(point.imag()) > limit) return;
        result.endpoint = C(static_cast<double>(point.real()),
                            static_cast<double>(point.imag()));
        result.endpointRepresentable = finite(result.endpoint);
    }

public:
    std::shared_ptr<const riemann::Model> source;
    double b = 0, h = 0, lambda = 0, infinityMultiplier = 0;
    bool validModel = false;
    std::string description, statusText;

    Model() : source(riemann::model(6, riemann::Settings{})) {
        description = "Riemann point-moment Jacobi quotient: direct rational iteration of lambda*(z-b-1/z)+c";
        if (!source || !source->validModel) {
            statusText = "The fixed H2 point-moment source is unavailable";
            return;
        }
        b = source->standardizedSkew;
        h = source->standardizedH;
        if (!std::isfinite(b) || !std::isfinite(h) || !(h > 0 && h < 1)) {
            statusText = "The point-moment quotient requires finite b and 0<h<1";
            return;
        }
        lambda = 1 / std::sqrt(h);
        infinityMultiplier = 1 / lambda;
        if (!std::isfinite(lambda) || !(lambda > 1) ||
            !std::isfinite(lambda * b) || !(infinityMultiplier > 0)) {
            statusText = "The point-moment quotient multiplier is invalid";
            return;
        }
        validModel = true;
        statusText = "Fixed rounded H2 point moments; direct rational iteration with critical points +i and -i. "
            "Infinity attracts; finite-budget nonescape is unresolved, not an RH proof.";
    }

    double escapeRadius(C parameter) const {
        if (!validModel || !validParameter(parameter))
            return std::numeric_limits<double>::quiet_NaN();
        const C d = parameter - lambda * b;
        return (1 + 1e-12) * std::max({1.0, 4 * std::abs(d) / (lambda - 1),
                                     std::sqrt(4 * lambda / (lambda - 1))});
    }

    Result orbit(C initial, C parameter, int budget,
                 std::vector<C> *trail = nullptr) const {
        Result result;
        result.endpoint = initial;
        result.endpointRepresentable = finite(initial);
        if (trail) trail->clear();
        if (!validModel || !validParameter(parameter)) {
            result.status = Status::InvalidModel;
            return result;
        }
        if (!finite(initial)) {
            result.status = Status::NumericalLimit;
            return result;
        }
        const double radius = escapeRadius(parameter);
        const C d = parameter - lambda * b;
        const double coefficientScale = std::max({1.0, lambda, std::abs(d)});
        if (!std::isfinite(radius) || !(radius > 0) || !finite(d) ||
            !std::isfinite(coefficientScale)) {
            result.status = Status::NumericalLimit;
            return result;
        }
        const double scaledLambda = lambda / coefficientScale;
        const C scaledD = d / coefficientScale;
        budget = std::clamp(budget, 32, 1024);
        // Componentwise scaling also accepts finite inputs whose Euclidean
        // magnitude would overflow. It does not impose an affine cutoff.
        const double initialScale = std::max({1.0, std::abs(initial.real()),
                                              std::abs(initial.imag())});
        C Z = initial / initialScale, W(1 / initialScale, 0);
        bool trailOpen = true;
        double structureSum = 0, weightSum = 0, sampleWeight = 1;
        for (int n = 0; n <= budget; ++n) {
            result.endpointZ = Z;
            result.endpointW = W;
            result.atInfinity = W == C{};
            result.iterations = n;
            result.smooth = static_cast<double>(n);
            if (!finite(Z) || !finite(W) || (Z == C{} && W == C{})) {
                result.status = Status::NumericalLimit;
                break;
            }
            if (trail && trailOpen) {
                storeEndpoint(result, Z, W);
                if (result.atInfinity || !result.endpointRepresentable ||
                    std::abs(result.endpoint) >= 1e12) trailOpen = false;
                if (trailOpen) trail->push_back(result.endpoint);
            }

            const double zMagnitude = std::abs(Z), wMagnitude = std::abs(W);
            // A bounded projective shape observable for colouring only.
            const double mass = zMagnitude * zMagnitude + wMagnitude * wMagnitude;
            const double shape = .55 * zMagnitude / (zMagnitude + wMagnitude) +
                .45 * 2 * std::abs((Z * std::conj(W)).imag()) / mass;
            structureSum += sampleWeight * std::clamp(shape, 0.0, 1.0);
            weightSum += sampleWeight;
            sampleWeight *= .97;

            // The same affine escape inequality is tested without Z/W.
            // Exact infinity is a valid escaped endpoint, including 0's image.
            if (result.atInfinity || zMagnitude / radius > wMagnitude) {
                result.status = Status::Escaped;
                if (!result.atInfinity) {
                    // Infinity has linear expansion lambda in the z chart,
                    // not polynomial degree-two log-log potential smoothing.
                    const double excess = std::log(zMagnitude) - std::log(wMagnitude) -
                                          std::log(radius);
                    result.smooth = std::max(0.0, n + 1.0 - excess / std::log(lambda));
                }
                break;
            }
            // Includes the endpoint produced by the final permitted update.
            if (n == budget) break;
            const C product = Z * W;
            const C nextZ = scaledLambda * (Z * Z - W * W) + scaledD * product;
            const C nextW = product / coefficientScale;
            const double coordinateScale = std::max(std::abs(nextZ), std::abs(nextW));
            if (!finite(nextZ) || !finite(nextW) || !std::isfinite(coordinateScale) ||
                !(coordinateScale > 0)) {
                result.status = Status::NumericalLimit;
                break;
            }
            const C normalizedZ = nextZ / coordinateScale;
            const C normalizedW = nextW / coordinateScale;
            // Loss of a nonzero denominator to floating-point underflow must
            // not masquerade as reaching the exact pole/infinity. With exact
            // arithmetic nextW vanishes precisely when Z or W vanishes.
            if (Z != C{} && W != C{} && normalizedW == C{}) {
                result.status = Status::NumericalLimit;
                break;
            }
            Z = normalizedZ;
            W = normalizedW;
        }
        storeEndpoint(result, result.endpointZ, result.endpointW);
        result.structure = weightSum > 0 ? structureSum / weightSum : 0;
        if (!result.atInfinity && finite(result.endpointZ) && finite(result.endpointW) &&
            result.endpointZ != C{} && result.endpointW != C{}) {
            result.phase = .5 + (std::arg(result.endpointZ) - std::arg(result.endpointW)) /
                                 (2 * projects::pi);
            result.phase -= std::floor(result.phase);
        }
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
        const Result plus = orbit(C(0, 1), parameter, budget);
        const Result minus = orbit(C(0, -1), parameter, budget);
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
        // Orange: +i escaped; cyan: -i escaped; violet: both escaped.
        // Any other critical orbit is unresolved within this finite budget.
        // These are orbit diagnostics, not connectedness or cycle labels.
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
} // namespace momentquotient
