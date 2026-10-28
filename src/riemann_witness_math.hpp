#pragma once
#include "project_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

// Finite, rounded-input moment witnesses and their complex root iterations.
// These floating-point calculations are not interval certificates or RH proofs.
namespace riemann {
using C = std::complex<double>;

struct Settings { int order = 2; double anchor = .0001; };
inline bool valid(const Settings &s) {
    return s.order >= 1 && s.order <= 2 && std::isfinite(s.anchor) &&
           s.anchor >= 1e-6 && s.anchor <= 1e-3;
}
inline bool sameSettings(const Settings &a, const Settings &b) {
    return a.order == b.order && a.anchor == b.anchor;
}

class Model {
    Settings settings;
    C relaxation{1, 0};
    int method = 0;
    std::vector<double> factorTargets;
    std::vector<C> roots;

    struct Jet {
        C value{1, 0}, first{}, second{}, newton{};
        double residual = 1;
        bool valid = false, exactRoot = false, critical = false;
    };

    // A common scale for all three jet components leaves both iteration
    // quotients unchanged. No expanded high-degree coefficients are formed.
    static bool normalise(C &p, C &d, C &dd) {
        if (!projects::finite(p) || !projects::finite(d) || !projects::finite(dd))
            return false;
        const double scale = std::max({std::abs(p), std::abs(d), std::abs(dd)});
        if (!(scale > 0) || !std::isfinite(scale)) return false;
        p /= scale; d /= scale; dd /= scale;
        return projects::finite(p) && projects::finite(d) && projects::finite(dd);
    }

    Jet evaluateRootProduct(C z) const {
        Jet out;
        if (!projects::finite(z) || std::abs(z) > 1e90 || roots.empty()) return out;
        out.residual = std::numeric_limits<double>::infinity();
        for (C root : roots) {
            C f = z - root, df{1, 0}, ddf{};
            const double residual = std::abs(f) / (1 + std::abs(z) + std::abs(root));
            if (!std::isfinite(residual)) return out;
            out.residual = std::min(out.residual, residual);
            // Only an individual zero factor is an exact root. A product
            // rounded to zero without one is a numerical failure instead.
            out.exactRoot = out.exactRoot || f == C{};
            if (!normalise(f, df, ddf)) return out;
            const C p = out.value * f;
            const C d = out.first * f + out.value * df;
            const C dd = out.second * f + 2.0 * out.first * df + out.value * ddf;
            if (p == C{} && !out.exactRoot) return out;
            out.value = p; out.first = d; out.second = dd;
            if (!normalise(out.value, out.first, out.second)) return out;
            if (out.value == C{} && !out.exactRoot) return out;
        }
        out.critical = out.first == C{};
        if (!out.critical) out.newton = out.value / out.first;
        out.valid = std::isfinite(out.residual);
        return out;
    }

    Jet evaluate(C z) const {
        if (generation >= 4) return evaluateRootProduct(z);
        Jet out;
        // A conservative evaluation limit, never a mathematical escape test.
        if (!projects::finite(z) || std::abs(z) > 1e90) return out;
        C carrier = z, first{1, 0}, second{};
        if (generation != 1) {
            const C oneMinus = C(1) - z;
            carrier = 6.75 * z * oneMinus * oneMinus;
            first = 6.75 * oneMinus * (C(1) - 3.0 * z);
            second = 6.75 * (6.0 * z - C(4));
        }
        if (!projects::finite(carrier) || !projects::finite(first) ||
            !projects::finite(second)) return out;
        out.residual = std::numeric_limits<double>::infinity();
        for (double target : factorTargets) {
            C f = carrier - C(target), df = first, ddf = second;
            const double denominator = 1 + std::abs(carrier) + std::abs(target);
            const double residual = std::abs(f) / denominator;
            if (!std::isfinite(residual)) return out;
            out.residual = std::min(out.residual, residual);
            // Product underflow must never be mistaken for an exact root.
            out.exactRoot = out.exactRoot || f == C{};
            if (!normalise(f, df, ddf)) return out;
            const C p = out.value * f;
            const C d = out.first * f + out.value * df;
            const C dd = out.second * f + 2.0 * out.first * df + out.value * ddf;
            out.value = p; out.first = d; out.second = dd;
            if (!normalise(out.value, out.first, out.second)) return out;
        }
        out.critical = out.first == C{};
        if (!out.critical) out.newton = out.value / out.first;
        out.valid = std::isfinite(out.residual);
        return out;
    }

    static long double kappa(long double z) {
        return 6.75L * z * (1 - z) * (1 - z);
    }
    static long double bisect(long double target, long double lo,
                              long double hi, bool increasing) {
        // Construction-time scalar bisection, not a pixel/orbit computation.
        for (int n = 0; n < 160; ++n) {
            const long double mid = lo + (hi - lo) / 2;
            if (mid == lo || mid == hi) break;
            if ((kappa(mid) < target) == increasing) lo = mid;
            else hi = mid;
        }
        return lo + (hi - lo) / 2;
    }
    bool appendCubicRoots(double target) {
        const long double t = target;
        if (!(t > 0) || !std::isfinite(t) || t == 1) return false;
        if (t < 1) {
            roots.emplace_back(static_cast<double>(bisect(t, 0, 1.L / 3, true)), 0);
            roots.emplace_back(static_cast<double>(bisect(t, 1.L / 3, 1, false)), 0);
            roots.emplace_back(static_cast<double>(bisect(t, 1, 4.L / 3, true)), 0);
        } else {
            const long double hi = std::max(2.L, std::cbrt(4 * t / 27) + 2);
            if (!std::isfinite(hi) || kappa(hi) < t) return false;
            const long double realRoot = bisect(t, 4.L / 3, hi, true);
            // Divide z^3-2z^2+z-4t/27 by its real-root factor. The remaining
            // roots are 1-r/2 +/- i*sqrt(r(3r-4))/2; this avoids Cardano loss.
            const long double pairReal = 1 - realRoot / 2;
            const long double pairImag = std::sqrt(realRoot * (3 * realRoot - 4)) / 2;
            if (!(pairImag > 0) || !std::isfinite(pairImag)) return false;
            roots.emplace_back(static_cast<double>(realRoot), 0);
            roots.emplace_back(static_cast<double>(pairReal), static_cast<double>(pairImag));
            roots.emplace_back(static_cast<double>(pairReal), -static_cast<double>(pairImag));
        }
        return true;
    }
    double rootPhase(C endpoint) const {
        if (roots.empty() || !projects::finite(endpoint)) return 0;
        std::size_t nearest = 0;
        double distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < roots.size(); ++i) {
            const double candidate = std::abs(endpoint - roots[i]);
            if (candidate < distance) { distance = candidate; nearest = i; }
        }
        return (static_cast<double>(nearest) + .25) / static_cast<double>(roots.size());
    }

    bool constructAtoms() {
        // Oct. 11 witness inputs. With D=r*d/dr and
        // Q_n(X)=X^2*product_{j=1..n}(4j^2-X^2)/(2n+1)!,
        // c_n=[Q_n(D) log xi(1/2+r)]_{r=1}.
        // c0..c2: conductor-1 entry of GRH_Witness_Search_v1.zip's
        // certificates.json. c3,c4: rh_h2_finite_certificate.py assertion
        // boxes. Only these rounded decimal midpoints are used here.
        constexpr std::array<long double, 5> c{{
            .09212431120995805L,
            .00058563952247175L,
            .00000897679033715L,
            .00000016329845283513060L,
            .000000003132575952496475L
        }};
        const std::size_t count = static_cast<std::size_t>(settings.order + 1);
        std::vector<long double> q(count), omega(count);
        if (settings.order == 1) {
            const long double mu = c[1] / c[0];
            const long double variance = c[2] / c[0] - mu * mu;
            if (!(mu > 0) || !(variance > 0) || !std::isfinite(variance)) return false;
            q[0] = mu / 2;
            q[1] = mu + 2 * variance / mu;
            omega[0] = c[0] * (q[1] - mu) / (q[1] - q[0]);
            omega[1] = c[0] * (mu - q[0]) / (q[1] - q[0]);
        } else {
            const long double a = settings.anchor;
            std::array<long double, 4> y{};
            for (std::size_t k = 0; k < y.size(); ++k) y[k] = c[k + 1] - a * c[k];
            const long double delta = y[0] * y[2] - y[1] * y[1];
            if (!(delta > 0) || !std::isfinite(delta)) return false;
            const long double A = (y[0] * y[3] - y[1] * y[2]) / delta;
            const long double B = (y[1] * y[3] - y[2] * y[2]) / delta;
            const long double discriminant = A * A - 4 * B;
            if (!(discriminant > 0) || !std::isfinite(discriminant)) return false;
            q[0] = a;
            q[1] = (A - std::sqrt(discriminant)) / 2;
            q[2] = (A + std::sqrt(discriminant)) / 2;
            for (std::size_t j = 0; j < q.size(); ++j) {
                const std::size_t u = (j + 1) % 3, v = (j + 2) % 3;
                const long double denominator = (q[j] - q[u]) * (q[j] - q[v]);
                if (denominator == 0 || !std::isfinite(denominator)) return false;
                omega[j] = (c[2] - (q[u] + q[v]) * c[1] + q[u] * q[v] * c[0]) /
                           denominator;
            }
        }
        for (std::size_t j = 0; j < count; ++j) {
            if (!(q[j] > 0 && q[j] < 1) || !(omega[j] > 0) ||
                !std::isfinite(q[j]) || !std::isfinite(omega[j])) return false;
            for (std::size_t k = 0; k < j; ++k) if (q[j] == q[k]) return false;
            const double node = static_cast<double>(q[j]);
            const double ordinary = static_cast<double>(omega[j]);
            const double weight = static_cast<double>(omega[j] / q[j]);
            if (!(node > 0 && node < 1) || !(ordinary > 0) || !(weight > 0) ||
                !std::isfinite(node) || !std::isfinite(ordinary) || !std::isfinite(weight))
                return false;
            for (double other : nodes) if (node == other) return false;
            nodes.push_back(node); ordinaryWeights.push_back(ordinary); weights.push_back(weight);
        }
        const std::size_t momentCount = static_cast<std::size_t>(2 * settings.order + 1);
        maxMomentResidual = 0;
        for (std::size_t n = 0; n < momentCount; ++n) {
            moments.push_back(static_cast<double>(c[n]));
            long double reconstructed = 0;
            for (std::size_t j = 0; j < count; ++j) {
                long double power = nodes[j];
                for (std::size_t k = 0; k < n; ++k) power *= nodes[j];
                reconstructed += static_cast<long double>(weights[j]) * power;
            }
            const long double reference = moments[n];
            const long double residual = std::abs(reconstructed - reference) / std::abs(reference);
            if (!std::isfinite(residual)) return false;
            maxMomentResidual = std::max(maxMomentResidual, static_cast<double>(residual));
        }
        // A floating-point consistency guard, not an interval enclosure.
        return maxMomentResidual <= 1e-8;
    }

    using LC = std::complex<long double>;
    static bool finiteLong(LC z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }
    static bool quadraticRoots(LC linear, LC constant, std::array<C, 2> &pair) {
        if (!finiteLong(linear) || !finiteLong(constant)) return false;
        const LC discriminant = linear * linear - 4.L * constant;
        if (!finiteLong(discriminant)) return false;
        const LC radical = std::sqrt(discriminant);
        const LC plus = -linear + radical, minus = -linear - radical;
        // Choose the numerator with less cancellation and recover its
        // partner from the product of the roots, rather than subtract again.
        LC first = (std::abs(plus) >= std::abs(minus) ? plus : minus) / 2.L;
        LC second;
        if (first == LC{}) {
            if (constant != LC{}) return false;
            second = -linear;
        } else second = constant / first;
        if (!finiteLong(first) || !finiteLong(second)) return false;
        pair = {{C(static_cast<double>(first.real()), static_cast<double>(first.imag())),
                 C(static_cast<double>(second.real()), static_cast<double>(second.imag()))}};
        for (C root : pair) {
            if (!projects::finite(root)) return false;
            const LC z(root.real(), root.imag());
            const long double denominator = 1 + std::norm(z) +
                std::abs(linear) * std::abs(z) + std::abs(constant);
            if (!(denominator > 0) || !std::isfinite(denominator)) return false;
            const long double residual = std::abs(z * z + linear * z + constant) / denominator;
            if (!std::isfinite(residual) || residual > 1e-12L) return false;
        }
        return true;
    }

    bool constructStandardizedRoots() {
        if (settings.order != 2 || moments.size() != 5) return false;
        const long double mass = moments[0];
        if (!(mass > 0) || !std::isfinite(mass)) return false;
        const long double mu = static_cast<long double>(moments[1]) / mass;
        const long double m2 = static_cast<long double>(moments[2]) / mass;
        const long double m3 = static_cast<long double>(moments[3]) / mass;
        const long double m4 = static_cast<long double>(moments[4]) / mass;
        const long double variance = m2 - mu * mu;
        if (!(variance > 0) || !std::isfinite(variance)) return false;
        const long double sigma = std::sqrt(variance);
        const long double central3 = m3 - 3 * mu * m2 + 2 * mu * mu * mu;
        const long double central4 = m4 - 4 * mu * m3 + 6 * mu * mu * m2 -
                                    3 * mu * mu * mu * mu;
        const long double b = central3 / (variance * sigma);
        const long double h = central4 / (variance * variance) - b * b - 1;
        if (!std::isfinite(mu) || !(sigma > 0) || !std::isfinite(sigma) ||
            !std::isfinite(b) || !(h > 0) || !std::isfinite(h)) return false;
        standardizedMean = static_cast<double>(mu);
        standardizedScale = static_cast<double>(sigma);
        standardizedSkew = static_cast<double>(b);
        standardizedH = static_cast<double>(h);
        if (!std::isfinite(standardizedMean) || !(standardizedScale > 0) ||
            !std::isfinite(standardizedScale) || !std::isfinite(standardizedSkew) ||
            !(standardizedH > 0) || !std::isfinite(standardizedH)) return false;

        // All coefficients below depend only on the five fixed point moments,
        // not on the arbitrary Radau anchor used for the witness diagnostics.
        std::array<C, 2> pair{};
        if (generation == 4 || generation == 7) {
            // P4=(z^2-b*z-1)^2+h*z^2. Build one complex quadratic
            // and its exact floating-point conjugate to preserve symmetry.
            if (!quadraticRoots(LC(-b, std::sqrt(h)), LC(-1, 0), pair)) return false;
            const std::array<C, 4> pencil{{pair[0], pair[1],
                                          std::conj(pair[0]), std::conj(pair[1])}};
            if (generation == 4) {
                roots.assign(pencil.begin(), pencil.end());
                description = "Riemann finite point moments: orthogonal-polynomial pencil root basins";
            } else {
                // Product of (z^2-2*r*z+1) is 16*z^4*P4((z+1/z)/2).
                // This is a degree-eight polynomial with nonzero constant;
                // its Newton/Halley iteration has no artificial root at zero.
                for (C r : pencil) {
                    if (!quadraticRoots(-2.L * LC(r.real(), r.imag()), LC(1, 0), pair))
                        return false;
                    if (pair[0] == C{} || pair[1] == C{}) return false;
                    roots.push_back(pair[0]); roots.push_back(pair[1]);
                }
                description = "Riemann finite point moments: polynomial Joukowski lift of the orthogonal pencil";
            }
        } else if (generation == 5) {
            // P5=z^4-(2+b^2+h)*z^2+1+h. Its quadratic in y=z^2
            // has two positive roots; expand the discriminant without loss
            // from subtracting nearly equal positive coefficients.
            const long double sum = 2 + b * b + h;
            const long double discriminant = b * b * b * b + 2 * b * b * (h + 2) + h * h;
            if (!(discriminant > 0) || !std::isfinite(discriminant)) return false;
            const long double large = (sum + std::sqrt(discriminant)) / 2;
            const long double small = (1 + h) / large;
            if (!(small > 0 && large > small) || !std::isfinite(small) ||
                !std::isfinite(large)) return false;
            for (long double y : {small, large}) {
                const long double denominator = 1 + y * y + sum * y + 1 + h;
                if (!(denominator > 0) || !std::isfinite(denominator)) return false;
                const long double residual = std::abs(y * y - sum * y + 1 + h) / denominator;
                if (!std::isfinite(residual) || residual > 1e-12L) return false;
                if (!quadraticRoots(LC{}, LC(-y, 0), pair)) return false;
                roots.push_back(pair[0]); roots.push_back(pair[1]);
            }
            description = "Riemann finite point moments: mirrored Christoffel-kernel polynomial root basins";
        } else if (generation == 6) {
            // The symmetric tridiagonal Jacobi closure has diagonal (0,b,0)
            // and off-diagonal (1,sqrt(h)); this is its characteristic polynomial.
            if (!quadraticRoots(LC(-b, 0), LC(-1 - h, 0), pair)) return false;
            roots = {C{}, pair[0], pair[1]};
            description = "Riemann finite point moments: Jacobi-closure characteristic-polynomial root basins";
        } else return false;
        for (std::size_t j = 0; j < roots.size(); ++j) {
            if (!projects::finite(roots[j])) return false;
            for (std::size_t k = 0; k < j; ++k) {
                const double separation = std::abs(roots[j] - roots[k]);
                const double scale = 1 + std::max(std::abs(roots[j]), std::abs(roots[k]));
                if (!(separation > 64 * std::numeric_limits<double>::epsilon() * scale))
                    return false;
            }
        }
        degree = generation == 6 ? 3 : generation == 7 ? 8 : 4;
        return true;
    }

public:
    int generation = 1, degree = 0;
    bool validModel = false;
    std::string description, statusText, coefficientOrigin;
    std::vector<double> moments, nodes, weights, ordinaryWeights;
    double maxMomentResidual = 0;
    double standardizedMean = 0, standardizedScale = 0;
    double standardizedSkew = 0, standardizedH = 0;

    Model(int gen, const Settings &s, C alpha = {1, 0}, int iterationMethod = 0)
        : settings(s), relaxation(alpha), method(iterationMethod), generation(gen) {
        coefficientOrigin = "Fixed rounded decimal midpoints: c0..c2 from GRH_Witness_Search_v1.zip "
            "certificates.json, conductor 1; c3,c4 from rh_h2_finite_certificate.py assertion boxes. "
            "Runtime uses point values, not an interval certificate.";
        if (!valid(s) || gen < 1 || gen > 7 || method < 0 || method > 1 ||
            !projects::finite(alpha) || !(std::abs(alpha) > 0) || std::abs(alpha) > 2) {
            statusText = "Invalid Riemann witness settings, relaxation, or iteration method";
            return;
        }
        if (gen >= 4 && s.order != 2) {
            statusText = "Standardized moment fractals require H2 and all five fixed point moments";
            return;
        }
        if (!constructAtoms()) {
            statusText = "Rounded-input witness construction failed positivity or moment consistency";
            return;
        }
        if (gen >= 4) {
            if (!constructStandardizedRoots()) {
                statusText = "Standardized point-moment polynomial failed positivity or root consistency";
                return;
            }
        } else if (gen == 1) {
            const double maximum = *std::max_element(nodes.begin(), nodes.end());
            for (double q : nodes) factorTargets.push_back(q / maximum);
            for (double target : factorTargets) roots.emplace_back(target, 0);
            description = "Riemann witness: scaled moment-node polynomial root basins";
        } else {
            for (double q : nodes) {
                const long double sum = 1 + std::sqrt(1 - static_cast<long double>(q));
                const double u = static_cast<double>(static_cast<long double>(q) / (sum * sum));
                if (!(u > 0 && u < 1) || !std::isfinite(u)) {
                    statusText = "Cayley factor construction reached a floating-point limit"; return;
                }
                factorTargets.push_back(u);
                if (gen == 3) factorTargets.push_back(1 / u);
            }
            for (double target : factorTargets) if (!appendCubicRoots(target)) {
                statusText = "Cubic root labelling reached a floating-point limit"; return;
            }
            description = gen == 2
                ? "Riemann witness: cubic pullback of the inner Cayley branch"
                : "Riemann witness: cubic pullback of both reciprocal Cayley branches";
        }
        if (gen <= 3)
            degree = static_cast<int>(factorTargets.size()) * (gen == 1 ? 1 : 3);
        if (roots.size() != static_cast<std::size_t>(degree) ||
            !std::all_of(roots.begin(), roots.end(), [](C z) { return projects::finite(z); })) {
            statusText = "Root labelling is incomplete or non-finite"; return;
        }
        std::sort(roots.begin(), roots.end(), [](C a, C b) {
            return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag());
        });
        validModel = true;
        statusText = std::string(gen >= 4 ? "Five standardized point moments with H2 witness diagnostics" :
            s.order == 1 ? "H1 two-atom" : "H2 anchored Radau three-atom") +
            " finite rounded-input witness; " + (method == 0 ? "Newton" : "Halley") +
            " iteration and residuals are floating-point diagnostics, not an RH proof";
    }

    projects::Result orbit(C initial, int budget, std::vector<C> *trail = nullptr) const {
        projects::Result result; result.endpoint = initial;
        if (trail) trail->clear();
        if (!validModel) { result.status = projects::Status::InvalidModel; return result; }
        budget = std::max(1, std::min(1024, budget));
        C z = initial;
        bool trailOpen = projects::finite(z) && std::abs(z) < 1e12;
        if (trail && trailOpen) trail->push_back(z);
        double structure = 0, weightSum = 0, sampleWeight = 1;
        for (int n = 0; n <= budget; ++n) {
            result.endpoint = z; result.iterations = n;
            const Jet jet = evaluate(z);
            if (!jet.valid) { result.status = projects::Status::NumericalLimit; break; }
            result.residual = jet.residual;
            structure += sampleWeight * (.5 + .5 * std::cos(8 * std::arg(z) +
                2 * std::log1p(std::min(1e100, std::abs(z)))));
            weightSum += sampleWeight; sampleWeight *= .97;
            // This criterion uses the unrelaxed Newton correction even for
            // Halley and complex relaxation, so a tiny step alone is no proof.
            if (jet.exactRoot || (!jet.critical && projects::finite(jet.newton) &&
                jet.residual <= 1e-11 && std::abs(jet.newton) <= 1e-8 * (1 + std::abs(z)))) {
                result.status = projects::Status::Converged; result.period = 1; break;
            }
            if (n == budget) break;
            C correction;
            if (method == 0) {
                if (jet.critical) { result.status = projects::Status::Pole; break; }
                correction = jet.newton;
            } else {
                const C denominator = 2.0 * jet.first * jet.first - jet.value * jet.second;
                if (denominator == C{}) { result.status = projects::Status::Pole; break; }
                correction = (2.0 * jet.value * jet.first) / denominator;
            }
            if (!projects::finite(correction)) {
                result.status = projects::Status::NumericalLimit; break;
            }
            const C next = z - relaxation * correction;
            if (!projects::finite(next)) {
                result.status = projects::Status::NumericalLimit; break;
            }
            z = next;
            if (std::abs(z) >= 1e12) trailOpen = false;
            if (trail && trailOpen) trail->push_back(z);
        }
        result.structure = weightSum > 0 ? structure / weightSum : 0;
        result.phase = rootPhase(result.endpoint);
        return result;
    }
};

inline std::shared_ptr<const Model> model(int generation, const Settings &settings,
                                         C relaxation = {1, 0}, int method = 0) {
    struct Entry {
        Settings settings;
        C relaxation{1, 0};
        int method = 0;
        std::shared_ptr<const Model> value;
    };
    static std::array<Entry, 7> cache;
    static std::mutex guard;
    if (generation < 1 || generation > 7)
        return std::make_shared<const Model>(generation, settings, relaxation, method);
    std::lock_guard<std::mutex> lock(guard);
    Entry &entry = cache[static_cast<std::size_t>(generation - 1)];
    if (entry.value && sameSettings(entry.settings, settings) &&
        entry.relaxation == relaxation && entry.method == method) return entry.value;
    auto value = std::make_shared<const Model>(generation, settings, relaxation, method);
    entry = {settings, relaxation, method, value};
    return value;
}
} // namespace riemann
