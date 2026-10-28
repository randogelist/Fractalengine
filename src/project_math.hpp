#pragma once
#include "project_grh.hpp"
#include "project_pi_bits.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Five project-derived complex iterations.  All orbit classifications below
// are finite floating-point diagnostics, not certified Julia-set membership.
namespace projects {
using C = std::complex<double>;
using Color = std::array<double, 3>;
inline constexpr double pi = 3.141592653589793238462643383279502884;

struct Settings {
    int piTerms = 64, piOffset = 0, nseModes = 6, grhDegree = 24;
    double nseTime = .0005, nseSlope = 1, pcdProbability = 272.0 / 761.0,
           grhRadius = 12;
};
inline bool valid(const Settings &s) {
    return s.piTerms >= 16 && s.piTerms <= 256 && s.piOffset >= 0 &&
           s.piOffset <= 1024 - s.piTerms && s.nseModes >= 4 && s.nseModes <= 7 &&
           s.grhDegree >= 12 && s.grhDegree <= 48 && s.grhDegree % 2 == 0 &&
           std::isfinite(s.nseTime) && s.nseTime >= 0 && s.nseTime <= .01 &&
           std::isfinite(s.nseSlope) && s.nseSlope >= 0 && s.nseSlope <= 3 &&
           std::isfinite(s.pcdProbability) && s.pcdProbability >= .02 &&
           s.pcdProbability <= .9 && std::isfinite(s.grhRadius) &&
           s.grhRadius >= 2 && s.grhRadius <= 24;
}
inline bool finite(C z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
inline double fraction(double x) { return x - std::floor(x); }

enum class Status { Converged, Cycle, Unresolved, Pole, NumericalLimit, InvalidModel };
inline const char *statusName(Status s) {
    switch (s) {
    case Status::Converged: return "numerical root convergence";
    case Status::Cycle: return "numerical periodic recurrence";
    case Status::Unresolved: return "unresolved within iteration budget";
    case Status::Pole: return "pole / projective infinity at final step";
    case Status::NumericalLimit: return "floating-point evaluation limit (not escape)";
    default: return "model unavailable";
    }
}
struct Result {
    Status status = Status::Unresolved;
    C endpoint{};
    int iterations = 0, period = 0;
    double residual = 1, structure = 0, phase = 0;
};

struct PolynomialValue {
    C correction{}, logarithmic{};
    double residual = 1;
    bool valid = false, pole = false, critical = false, exactRoot = false;
};

// Coefficients are in ascending powers.  Reciprocal Horner avoids z^degree
// overflow and retains the rational map's large-z behaviour.  residual is the
// coefficientwise relative backward residual |P| / sum |c_j z^j|.
inline PolynomialValue polynomialValue(const std::vector<C> &c, C z) {
    PolynomialValue r;
    if (c.empty() || !finite(z)) return r;
    int d = static_cast<int>(c.size()) - 1;
    C f{}, df{};
    double bound = 0;
    if (std::abs(z) <= 1) {
        f = c.back(); bound = std::abs(c.back());
        for (int j = d - 1; j >= 0; --j) {
            df = df * z + f;
            f = f * z + c[static_cast<std::size_t>(j)];
            bound = bound * std::abs(z) + std::abs(c[static_cast<std::size_t>(j)]);
        }
        if (!finite(f) || !finite(df) || !std::isfinite(bound) || bound == 0) return r;
        r.residual = std::abs(f) / bound;
        r.exactRoot = f == C{};
        r.pole = r.exactRoot;
        r.critical = df == C{};
        if (!r.critical) r.correction = f / df;
        if (!r.pole) r.logarithmic = df / f;
    } else {
        C w = C(1) / z;
        f = c.front(); bound = std::abs(c.front());
        for (int j = 1; j <= d; ++j) {
            df = df * w + f;
            f = f * w + c[static_cast<std::size_t>(j)];
            bound = bound * std::abs(w) + std::abs(c[static_cast<std::size_t>(j)]);
        }
        C denominator = double(d) * f - w * df;
        if (!finite(f) || !finite(denominator) || !std::isfinite(bound) || bound == 0) return r;
        r.residual = std::abs(f) / bound;
        r.exactRoot = f == C{};
        r.pole = r.exactRoot;
        r.critical = denominator == C{};
        if (!r.critical) r.correction = z * (f / denominator);
        if (!r.pole) r.logarithmic = (double(d) - w * (df / f)) * w;
    }
    r.valid = std::isfinite(r.residual);
    return r;
}

struct SpherePoint { C z{}; bool infinity = false; };
inline double chordal(SpherePoint a, SpherePoint b) {
    if (a.infinity && b.infinity) return 0;
    double sa = std::max({1.0, std::abs(a.z.real()), std::abs(a.z.imag())}),
           sb = std::max({1.0, std::abs(b.z.real()), std::abs(b.z.imag())});
    C ua = a.z / sa, ub = b.z / sb;
    double va = 1 / sa, vb = 1 / sb;
    if (a.infinity) return vb / std::hypot(std::abs(ub), vb);
    if (b.infinity) return va / std::hypot(std::abs(ua), va);
    return std::abs(ua * vb - ub * va) /
           (std::hypot(std::abs(ua), va) * std::hypot(std::abs(ub), vb));
}
inline double pointPhase(C z, int period = 0) {
    return fraction(.5 + std::arg(z) / (2 * pi) +
                    .19 * std::log1p(std::min(1e100, std::abs(z))) + .073 * period);
}

class Model {
    std::vector<C> coefficients;
    Settings settings;
    C a{}, b{1};
    double coordinateScale = 1;

    static bool collatzStep(C z, C &next) {
        // cos(pi*z) grows exponentially off the real axis.  This guard is an
        // evaluation limit, not a mathematical basin or escape criterion.
        if (!finite(z) || std::abs(z.imag()) > 180 || std::abs(z) > 1e100) return false;
        next = z + .25 - (z * .5 + .25) * std::cos(pi * z);
        return finite(next) && std::abs(next) <= 1e100;
    }
    static double structureSample(SpherePoint z) {
        if (z.infinity) return .5;
        return .5 + .5 * std::cos(8 * std::arg(z.z) +
                                 2 * std::log1p(std::min(1e100, std::abs(z.z))));
    }
    double recurrenceDistance(SpherePoint u, SpherePoint v) const {
        if (generation == 1) return chordal(u, v);
        // Infinity is an essential singularity of the entire Collatz map.
        // Chordal proximity there does not mean finite periodic recurrence.
        if (u.infinity || v.infinity) return 1;
        double scale = std::max({1.0, std::abs(u.z.real()), std::abs(u.z.imag()),
                                 std::abs(v.z.real()), std::abs(v.z.imag())});
        return std::abs(u.z / scale - v.z / scale);
    }
    Result newtonOrbit(C initial, int budget, std::vector<C> *trail) const {
        Result r; C z = initial;
        double structure = 0, weightSum = 0, weight = 1;
        bool trailOpen = true;
        if (trail) { trail->clear(); trail->push_back(z); }
        for (int n = 0; n < budget; ++n) {
            r.endpoint = z; r.iterations = n;
            if (!finite(z)) { r.status = Status::NumericalLimit; break; }
            auto v = polynomialValue(coefficients, z);
            if (!v.valid) { r.status = Status::NumericalLimit; break; }
            r.residual = v.residual;
            structure += weight * structureSample({z, false}); weightSum += weight; weight *= .97;
            // Residual alone can misclassify large or ill-conditioned values.
            // Also require a small Newton correction unless P is exactly zero.
            if (v.exactRoot || (!v.critical && finite(v.correction) &&
                v.residual <= 1e-11 && std::abs(v.correction) <= 1e-8 * (1 + std::abs(z)))) {
                r.status = Status::Converged; r.period = 1; break;
            }
            if (v.critical) { r.status = Status::Pole; break; }
            if (!finite(v.correction)) { r.status = Status::NumericalLimit; break; }
            C next = z - v.correction;
            if (!finite(next)) { r.status = Status::NumericalLimit; break; }
            z = next; r.endpoint = z; r.iterations = n + 1;
            if (std::abs(z) >= 1e12) trailOpen = false;
            if (trail && trailOpen) trail->push_back(z);
        }
        r.structure = weightSum > 0 ? structure / weightSum : 0;
        r.phase = finite(r.endpoint) ? pointPhase(r.endpoint) : 0;
        return r;
    }
    Result directOrbit(C initial, int budget, std::vector<C> *trail) const {
        Result r;
        SpherePoint z{initial, false};
        std::array<SpherePoint, 64> history{};
        history[0] = z;
        double structure = 0, weightSum = 0, weight = 1;
        bool trailOpen = true;
        if (trail) { trail->clear(); trail->push_back(initial); }
        for (int n = 0; n < budget; ++n) {
            r.iterations = n;
            if (!z.infinity && !finite(z.z)) { r.status = Status::NumericalLimit; break; }
            structure += weight * structureSample(z); weightSum += weight; weight *= .97;
            SpherePoint next;
            if (generation == 1) {
                if (z.infinity) next = {a, false}; // rational extension: infinity -> a
                else {
                    auto v = polynomialValue(coefficients, z.z / coordinateScale);
                    if (!v.valid) { r.status = Status::NumericalLimit; break; }
                    if (v.pole) next = {{}, true};
                    else {
                        next = {a + b * (v.logarithmic / coordinateScale), false};
                        if (!finite(next.z)) { r.status = Status::NumericalLimit; break; }
                    }
                }
            } else {
                C intermediate;
                if (!collatzStep(z.z, intermediate) || !collatzStep(intermediate, next.z)) {
                    r.status = Status::NumericalLimit; break;
                }
            }
            z = next;
            int count = n + 1;
            history[static_cast<std::size_t>(count) % history.size()] = z;
            r.iterations = count;
            r.endpoint = z.infinity ? C{} : z.z;
            if (z.infinity || (!z.infinity && std::abs(z.z) >= 1e12)) trailOpen = false;
            if (trail && trailOpen) trail->push_back(z.z);
            // Require three returns, and check three neighbouring phases of
            // each candidate period.  This is deliberately a finite diagnostic.
            int maxPeriod = std::min(16, (count - 2) / 3);
            if (count % 4 != 0) continue;
            for (int p = 1; p <= maxPeriod; ++p) {
                double error = recurrenceDistance(z,
                    history[static_cast<std::size_t>(count - p) % history.size()]);
                if (error > 1e-8) continue;
                for (int lag = 0; lag < 3 && error <= 1e-8; ++lag)
                    for (int ret = 1; ret <= 3 && error <= 1e-8; ++ret) {
                        const auto &u = history[static_cast<std::size_t>(count - lag) % history.size()];
                        const auto &v = history[static_cast<std::size_t>(count - lag - ret * p) % history.size()];
                        error = std::max(error, recurrenceDistance(u, v));
                    }
                if (error <= 1e-8) {
                    r.status = Status::Cycle; r.period = p; r.residual = error;
                    // Canonical representative makes hue independent of which
                    // cyclic phase was the last one visited.
                    bool have = false; C representative{};
                    for (int j = 0; j < p; ++j) {
                        auto q = history[static_cast<std::size_t>(count - j) % history.size()];
                        if (!q.infinity && (!have || q.z.real() < representative.real() ||
                            (q.z.real() == representative.real() && q.z.imag() < representative.imag()))) {
                            representative = q.z; have = true;
                        }
                    }
                    r.phase = have ? pointPhase(representative, p) : 0;
                    r.structure = weightSum > 0 ? structure / weightSum : 0;
                    return r;
                }
            }
        }
        if (r.status == Status::Unresolved && z.infinity) r.status = Status::Pole;
        r.endpoint = z.infinity ? C{} : z.z;
        r.phase = finite(r.endpoint) ? pointPhase(r.endpoint) : 0;
        r.structure = weightSum > 0 ? structure / weightSum : 0;
        return r;
    }

public:
    int generation = 1, degree = 0;
    bool validModel = false, piCertified = false;
    std::string description, statusText;

    Model(int gen, const Settings &s, C offset = {}, C feedback = {1, 0})
        : settings(s), a(offset), b(feedback), generation(gen) {
        if (!valid(s) || gen < 1 || gen > 5) { statusText = "Invalid project settings"; return; }
        try {
            if (gen == 1) {
                if (!finite(a) || !finite(b) || b == C{}) {
                    statusText = "GRH pole iteration requires finite a and nonzero finite b"; return;
                }
                auto carrier = project_grh::makeCarrier(s.grhDegree, s.grhRadius);
                if (!carrier.valid) { statusText = carrier.status; return; }
                coefficients = std::move(carrier.coefficients);
                coordinateScale = carrier.radius;
                description = "GRH: finite even completed-L Taylor carrier; rational pole Julia";
                statusText = carrier.status;
            } else if (gen == 2) {
                piCertified = project_pi::certified();
                coefficients.reserve(static_cast<std::size_t>(s.piTerms));
                for (int t = 0; t < s.piTerms; ++t) {
                    int n = 1024 + s.piOffset + t;
                    double c = project_pi::actualPiBit(static_cast<std::size_t>(2 * n)) -
                               project_pi::actualPiBit(static_cast<std::size_t>(2 * n + 1)) -
                               1024.0 * project_pi::drift(n);
                    coefficients.emplace_back(c, 0);
                }
                description = "PI: actual binary-pair discrepancy polynomial Newton basins";
                statusText = "Actual pi bits certified by exact integer interval; drift and orbits floating point";
            } else if (gen == 3) {
                for (int j = 0; j < s.nseModes; ++j) {
                    double energy = std::pow(double(j + 1), -s.nseSlope) *
                                    std::exp(-2 * s.nseTime * std::ldexp(1.0, 2 * j));
                    coefficients.emplace_back(energy, 0);
                }
                description = "NSE: exact periodic-shear energy polynomial Newton basins";
                statusText = "Time parameter is viscosity*time; fractal of complex Newton iteration";
            } else if (gen == 4) {
                description = "COLLATZ: Chamberland entire extension, two-step return Julia";
                statusText = "F squared; integer shortcut agrees exactly; numerical cutoffs are not escape";
                degree = 0; validModel = true; return;
            } else {
                coefficients = {C(-15), C(-3), C(-1), C(0), C(-1),
                                C(s.pcdProbability / (2 - s.pcdProbability))};
                description = "P_c(d,q): q=2 GP-4 critical-dimension quintic Newton basins";
                statusText = "z=2d-1 complexified; analytic GP-4 formula model, not exact lattice threshold";
            }
            while (coefficients.size() > 1 && coefficients.back() == C{}) coefficients.pop_back();
            degree = static_cast<int>(coefficients.size()) - 1;
            validModel = degree >= 2 && std::all_of(coefficients.begin(), coefficients.end(),
                [](C z) { return finite(z); });
            if (!validModel) statusText = "Carrier has insufficient degree or non-finite coefficients";
        } catch (const std::exception &e) {
            statusText = e.what(); validModel = false;
        }
    }
    Result orbit(C initial, int budget, std::vector<C> *trail = nullptr) const {
        if (!validModel) {
            if (trail) trail->clear();
            Result r; r.endpoint = initial; r.status = Status::InvalidModel; return r;
        }
        budget = std::max(1, std::min(4096, budget));
        return generation == 1 || generation == 4 ? directOrbit(initial, budget, trail)
                                                  : newtonOrbit(initial, budget, trail);
    }
};

inline bool sameSettings(const Settings &x, const Settings &y) {
    return x.piTerms == y.piTerms && x.piOffset == y.piOffset &&
           x.nseModes == y.nseModes && x.grhDegree == y.grhDegree &&
           x.nseTime == y.nseTime && x.nseSlope == y.nseSlope &&
           x.pcdProbability == y.pcdProbability && x.grhRadius == y.grhRadius;
}
inline std::shared_ptr<const Model> model(int generation, const Settings &s, C a = {}, C b = {1, 0}) {
    struct Entry { Settings settings; C a{}, b{}; std::shared_ptr<const Model> value; };
    static std::array<Entry, 5> cache;
    static std::mutex guard;
    if (generation < 1 || generation > 5) return std::make_shared<const Model>(generation, s, a, b);
    if (generation != 1) { a = {}; b = {1, 0}; }
    std::lock_guard<std::mutex> lock(guard);
    Entry &entry = cache[static_cast<std::size_t>(generation - 1)];
    if (entry.value && sameSettings(entry.settings, s) && entry.a == a && entry.b == b) return entry.value;
    auto value = std::make_shared<const Model>(generation, s, a, b);
    entry = {s, a, b, value};
    return value;
}

inline Color hsv(double hue, double saturation, double value) {
    double h = fraction(hue) * 6, f = h - std::floor(h);
    double p = value * (1 - saturation), q = value * (1 - saturation * f),
           t = value * (1 - saturation * (1 - f));
    switch (static_cast<int>(h) % 6) {
    case 0: return {value, t, p}; case 1: return {q, value, p};
    case 2: return {p, value, t}; case 3: return {p, q, value};
    case 4: return {t, p, value}; default: return {value, p, q};
    }
}
inline Color colour(const Result &r, int palette, double hue, double bands) {
    if (r.status == Status::InvalidModel) return {.32, .035, .06};
    if (r.status == Status::Pole) return {.97, .93, .76};
    bool captured = r.status == Status::Converged || r.status == Status::Cycle;
    double time = std::log1p(std::max(0, r.iterations));
    double stripes = .5 + .5 * std::cos(bands * (time * 4.7 + r.structure * 3));
    if (palette == 1) {
        double value = captured ? .47 + .49 * stripes : .10 + .55 * stripes;
        return hsv(hue + .64 - .10 * bands * time, .75, value);
    }
    if (palette == 2) {
        double value = captured ? .38 + .59 * stripes : .18 + .61 * stripes;
        return hsv(hue + .68 * r.structure + .16 * r.phase, .73, value);
    }
    if (captured) {
        double edge = .45 + .55 * std::exp(-.016 * std::max(0, r.iterations));
        return hsv(hue + r.phase, .76, edge * (.84 + .16 * stripes));
    }
    if (r.status == Status::NumericalLimit)
        return hsv(hue + .61 + .027 * bands * time, .78, .12 + .62 * stripes);
    return hsv(hue + .65 + .09 * r.structure, .55, .035 + .16 * stripes);
}

inline int selfTest(std::ostream &out) {
    int checks = 0;
    auto check = [&](bool ok, const char *name) {
        ++checks;
        if (!ok) throw std::runtime_error(std::string("project math self-test: ") + name);
    };
    check(valid(Settings{}), "default settings");
    Settings bad; bad.piOffset = 1024;
    check(!valid(bad), "pi prefix range");
    auto q = polynomialValue({C(-15), C(-3), C(-1), C(0), C(-1), C(136.0 / 625)}, C(5));
    check(q.valid && q.residual < 1e-14, "GP-4 real dimension anchor");
    auto n = polynomialValue({C(-1), C(0), C(1)}, C(2));
    check(n.valid && std::abs(n.correction - C(.75)) < 1e-14, "reciprocal Newton ratio");
    check(std::abs(chordal({C(0), false}, {C{}, true}) - 1) < 1e-14,
          "projective infinity distance");
    out << "PROJECT_ITERATION_IDENTITIES=PASS CHECKS=" << checks << '\n';
    return checks;
}
} // namespace projects
