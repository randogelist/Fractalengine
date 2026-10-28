#pragma once
#include "qqrh_math.hpp"
#include "rh_master_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

// Direct iteration of six finite-Taylor QQRH hosts.  There is no outer
// BKQR step, sign branch, or lambda coupling in this adapter.  Parameter
// portraits follow the chosen seed z=0, not all critical points.  |z|>16 is
// only a display cutoff, not a certified escape radius for these maps.
namespace qqrhdirect {
using C = std::complex<double>;
using Color = std::array<double, 3>;
using Status = rhmaster::Status;
using Result = rhmaster::Result;
using rhmaster::statusName;

inline const char* title(int generation) {
    switch (generation) {
    case 1: return "QQRH direct / radial curvature";
    case 2: return "QQRH direct / dyadic ratio";
    case 3: return "QQRH direct / Laplacian level 1";
    case 4: return "QQRH direct / Laplacian level 2";
    case 5: return "QQRH direct / moving boundary";
    case 6: return "QQRH / scale-ratio curvature / direct";
    default: return "QQRH direct / invalid channel";
    }
}

inline const char* formula(int generation) {
    switch (generation) {
    case 1:
        return "y=b*z+h; H(y)=T(-y); C_H=-d_y[y*H'/H]; z_next=c+eta*(C_H(y)-C_H(h))";
    case 2:
        return "y=b*z+h; P(y,x)=T(y*rho+)*T(y*rho-); D=P(2y,x)/P(y,x); z_next=c+eta*(D(y,x)-D(h,x))";
    case 3:
        return "y=b*z+h; F(y)=Psi_(T')(y,x); z_next=z^2+c-eta*F(y)/(b*F'(y))";
    case 4:
        return "y=b*z+h; F(y)=Psi_(T'')(y,x); z_next=z^2+c-eta*F(y)/(b*F'(y))";
    case 5:
        return "y=b*z+h; a=(1-m^2)/(1+m^2); F(y)=[Psi_T(y,x)-Psi_T(y,a)]/[y*(x-a)], removable factors cancelled; z_next=z^2+c-eta*F(y)/(b*F'(y))";
    case 6:
        return "y=b*z+h; H(y)=T(-y); C_H=-d_y[y*H'/H]; K_kappa(y)=kappa*C_H(kappa*y)-C_H(y); z_next=c+eta*(K_kappa(y)-K_kappa(h)); kappa=1 gives exactly c";
    default: return "Invalid direct QQRH channel";
    }
}

class Model {
    const qqrh::Channel channel_;
    // Only generation 6 needs a second immutable carrier.  Composing two
    // existing curvature channels preserves their stable polynomial jets.
    const std::unique_ptr<const qqrh::Channel> scaledChannel_;

    static bool finite(C z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }

    static double logRadius(C z) {
        const double scale = std::max(std::abs(z.real()), std::abs(z.imag()));
        if (scale == 0.0) return -std::numeric_limits<double>::infinity();
        const double re = z.real() / scale, im = z.imag() / scale;
        return std::log(scale) + 0.5 * std::log(re * re + im * im);
    }

    static Color errorColour(Status status) {
        return status == Status::InvalidModel ? Color{{.42, .035, .055}}
                                              : Color{{.83, .12, .72}};
    }

    // Match the existing RH adapter's orbit display, including dark unresolved
    // points.  Finite-tail near-period detections remain available as metadata.
    static Color colour(const Result& result, int palette, double hue, double bands) {
        if (result.status == Status::InvalidModel ||
            result.status == Status::NumericalLimit)
            return errorColour(result.status);
        hue = std::isfinite(hue) ? hue : 0.0;
        bands = std::isfinite(bands) ? std::clamp(bands, 0.0, 100.0) : 1.0;
        if (result.status == Status::Unresolved)
            return projects::hsv(hue + .61, .48, .019 + .017 * result.structure);
        const double stripes = .5 + .5 * std::cos(
            bands * (.57 * result.smooth + 2.4 * result.structure));
        const double colorHue = palette == 1 ? hue + .66 - .023 * result.smooth
            : palette == 2 ? hue + .65 * result.structure + .18 * result.phase
            : hue + .57 + .30 * result.phase + .012 * result.smooth;
        return projects::hsv(colorHue, .74, .26 + .72 * stripes);
    }

    static void detectTailPeriod(Result& result, const std::array<C, 80>& tail) {
        if (result.status != Status::Unresolved || result.iterations < 32) return;
        const int available = std::min(result.iterations + 1, 80);
        const int last = result.iterations % 80;
        constexpr double tolerance = 1e-7;
        // Compare four entire consecutive blocks, requiring five blocks of
        // stored states.  This is recurrence within a finite budget, not a
        // proof of exact periodicity or stability of later iterates.
        for (int period = 1; period <= 16 && 5 * period <= available; ++period) {
            double residual = 0.0;
            bool repeated = true;
            for (int lag = 0; lag < 4 * period; ++lag) {
                const C recent = tail[(last + 80 - lag) % 80];
                const C previous = tail[(last + 80 - lag - period) % 80];
                const double error = std::abs(recent - previous) /
                    (1.0 + std::max(std::abs(recent), std::abs(previous)));
                if (!std::isfinite(error) || error > tolerance) {
                    repeated = false;
                    break;
                }
                residual = std::max(residual, error);
            }
            if (repeated) {
                result.detectedPeriod = period;
                result.periodRepeats = 4;
                result.periodResidual = residual;
                return;
            }
        }
    }

public:
    int generation = 1;
    int N = 16;
    double x = 0.0;
    double b = 1.0;
    double h = 0.0;
    C eta{.45, 0.0};
    double slope = 11.0;
    double scale = 2.0;
    bool validModel = false;
    int degree = 0; // Carrier polynomial degree, not the full rational-map degree.
    std::string description;
    std::string statusText;

    Model(int g, int truncation, double angular, double radialScale,
          double shift, C strength, double movingSlope = 11.0,
          double curvatureScale = 2.0)
        : channel_(g == 6 ? 1 : g, truncation, angular, radialScale, shift, movingSlope),
          scaledChannel_(g == 6 ? std::make_unique<const qqrh::Channel>(
              1, truncation, angular, curvatureScale * radialScale,
              curvatureScale * shift, movingSlope) : nullptr),
          generation(g), N(truncation), x(angular), b(radialScale), h(shift),
          eta(strength), slope(movingSlope), scale(curvatureScale),
          degree(g == 6 ? 2 * channel_.degree() : channel_.degree()),
          description(title(g)) {
        // A direct channel's centering point may be a pole.  Keep structural
        // validity separate so eta=0 can still bypass the unused carrier.
        validModel = g >= 1 && g <= 6 && channel_.valid() && finite(eta) &&
            (g != 6 || (std::isfinite(scale) && scale >= .25 && scale <= 4.0));
        if (g == 6) {
            // degree records deg[H(kappa*y)H(y)]=2N, not the full map's
            // degree (generically 4N).  Common poles may cancel; kappa=1
            // cancels the entire carrier before either response is evaluated.
            if (!validModel) {
                statusText = "Require N=8/16/32, -1<=x<=1, finite b>0, h, complex eta, 0<=m<=64 and .25<=kappa<=4.";
            } else if (eta == C{} || scale == 1.0) {
                statusText = "Exact constant map z_next=c: eta=0 or kappa=1 bypasses both curvature evaluations. |z|>16 remains a display cutoff only.";
            } else if (!channel_.centerValid() || !scaledChannel_->valid() ||
                       !scaledChannel_->centerValid()) {
                statusText = "An active curvature carrier has a nonfinite scaled parameter or centering value; evaluation reports a numerical limit. eta=0 or kappa=1 bypasses both carriers.";
            } else {
                statusText = "Rounded finite Taylor coefficients; scale-ratio logarithmic curvature with fixed complex eta. Product carrier degree 2N; generic rational-map degree 4N. |z|>16 is a display cutoff only; finite-tail detections do not certify stability or RH.";
            }
        } else if (!validModel) {
            statusText = "Require generation 1..5, N=8/16/32, -1<=x<=1, finite b>0, h, complex eta and 0<=m<=64.";
        } else if (eta != C{} && !channel_.centerValid()) {
            statusText = "The active direct carrier has a nonfinite centering value at h; evaluation reports a numerical limit. eta=0 bypasses the carrier.";
        } else {
            statusText = "Rounded finite Taylor coefficients; direct Q iteration with fixed complex eta. |z|>16 is a display cutoff only. Seed-0 parameter slices; finite-tail near-periods do not certify stability or RH.";
        }
    }

    Result orbit(C initial, C parameter, int budget,
                 std::vector<C>* trail = nullptr) const {
        Result result;
        result.endpoint = initial;
        if (trail) trail->clear();
        if (!validModel || !finite(parameter)) {
            result.status = Status::InvalidModel;
            return result;
        }
        if (!finite(initial)) {
            result.status = Status::NumericalLimit;
            return result;
        }

        budget = std::clamp(budget, 0, 1024);
        const bool constantCarrier = generation == 6 && scale == 1.0;
        C z = initial;
        std::array<C, 80> tail{};
        double minimumTrap = std::numeric_limits<double>::infinity();
        for (int n = 0; n <= budget; ++n) {
            result.endpoint = z;
            result.iterations = n;
            result.smooth = static_cast<double>(n);
            tail[n % 80] = z;
            const double largestComponent = std::max(std::abs(z.real()), std::abs(z.imag()));
            if (trail && largestComponent < 1e12)
                trail->push_back(z);
            minimumTrap = std::min(minimumTrap, std::abs(z - C(-.21, .37)));

            // Only finite endpoints reach here.  Avoid squaring large
            // components, and keep smoothing representable at any magnitude.
            if (largestComponent > 16.0 || std::norm(z) > 256.0) {
                const double logR = logRadius(z);
                result.status = Status::Escaped;
                // Inherited cosmetic smoothing, not a certified escape rate.
                result.smooth = std::max(0.0, n + 1.0 -
                    std::log(logR / std::log(16.0)) / std::log(2.0));
                break;
            }
            if (n == budget) break;

            C correction{};
            if (eta != C{} && !constantCarrier) {
                if (generation == 6) {
                    // Curvature of R_kappa(y)=H(kappa*y)/H(y):
                    // -d_y[y*R'/R]=kappa*C_H(kappa*y)-C_H(y).
                    // Each existing response already subtracts its own center.
                    const C scaled = scaledChannel_->response(z);
                    const C original = channel_.response(z);
                    correction = scale * scaled - original;
                } else {
                    correction = channel_.response(z);
                }
                if (!finite(correction)) {
                    result.status = Status::NumericalLimit;
                    break;
                }
            }
            // eta=0 gives exactly c for curvature/dyadic and exactly z^2+c
            // for the other three hosts, even at an unused carrier pole.
            const C baseline = channel_.direct() ? parameter : z * z + parameter;
            const C feedback = eta.imag() == 0.0
                ? eta.real() * correction : eta * correction;
            const C next = eta == C{} || constantCarrier ? baseline : baseline + feedback;
            if (!finite(next)) {
                result.status = Status::NumericalLimit;
                break;
            }
            z = next;
        }
        result.structure = std::exp(-2.5 * minimumTrap);
        result.phase = .5 + std::arg(result.endpoint) / (2.0 * projects::pi);
        result.phase -= std::floor(result.phase);
        detectTailPeriod(result, tail);
        return result;
    }

    Color juliaColour(C initial, C parameter, int budget, int palette,
                      double hue, double bands) const {
        return colour(orbit(initial, parameter, budget), palette, hue, bands);
    }

    Color atlasColour(C parameter, int budget, int palette,
                      double hue, double bands) const {
        return colour(orbit(C{}, parameter, budget), palette, hue, bands);
    }
};

// One immutable cached model per generation.  Model construction is shared
// by a frame; neither carrier coefficients nor centering jets rebuild per pixel.
inline std::shared_ptr<const Model> model(int generation, int N = 16,
        double x = 0.0, double b = 1.0, double h = 0.0,
        C eta = C(.45, 0.0), double slope = 11.0, double scale = 2.0) {
    struct Entry {
        int N = 0;
        double x = 0.0, b = 0.0, h = 0.0;
        C eta{};
        double slope = 0.0;
        double scale = 2.0;
        std::shared_ptr<const Model> value;
    };
    static std::array<Entry, 6> cache{};
    static std::mutex guard;
    if (generation < 1 || generation > 6)
        return std::make_shared<const Model>(generation, N, x, b, h, eta, slope, scale);
    std::lock_guard<std::mutex> lock(guard);
    Entry& entry = cache[static_cast<std::size_t>(generation - 1)];
    if (entry.value && entry.N == N && entry.x == x && entry.b == b &&
        entry.h == h && entry.eta == eta && entry.slope == slope &&
        (generation != 6 || entry.scale == scale))
        return entry.value;
    auto value = std::make_shared<const Model>(generation, N, x, b, h, eta, slope, scale);
    entry = {N, x, b, h, eta, slope, scale, value};
    return value;
}

} // namespace qqrhdirect
