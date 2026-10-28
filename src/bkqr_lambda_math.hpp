#pragma once
#include "rh_master.hpp"
#include "qqrh_math.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// CPU port of lambda_map_only.hip, BKQR lambda-only v0.4.7.
// Preserve the original maps, cell-index seed hashes, shadow-orbit rules,
// integer accumulation, nine metrics, percentile contrast and colour ramp.
// CPU libm and GPU transcendental implementations need not agree bitwise.
// "escape" intentionally retains the legacy finite/numerical stop test.
// "lyap"/"twin" are the same finite, unrenormalized shadow-separation metric.
// No metric is a proof of stability, escape, a phase boundary, or RH.
namespace bkqrlambda {
using C = std::complex<double>;
using Color = std::array<double, 3>;
using Cancel = std::function<bool()>;
inline constexpr double PI = 3.141592653589793238462643383279502884;
inline constexpr int RAW = 0, SELFDUAL = 1, COMPACT = 2;
inline constexpr int JULIA = 0, NEWTON = 1, TRICORN = 2, BURNING = 3, NOVA = 4;
inline constexpr int GREEDY = 0, SUBTRACT = 1, ADD = 2;
inline constexpr int METRICS = 9;

struct Settings {
    C c{-.8, .156};
    int host = JULIA;
    int branch = GREEDY;
    int coupling = RAW;
    int power = 2;
    int quality = 0;
    int precision = 0;
    int iterations = 96;
    int metric = 0;
    // Zero preserves the archived host and arithmetic.  Positive values
    // place the corresponding RH-derived finite polynomial inside BKQR.
    int rhChannel = 0;
    int rhN = 16;
    double rhX = 0.0, rhB = .9, rhH = .28;
    C eta{.45, 0.0};
    int qqrhChannel = 0;
    double qqrhSlope = 11.0;
};

inline const char* rhChannelName(int n) {
    static const char* names[] = {"LEGACY", "MASTER", "SUB_1", "SUB_2"};
    return names[std::clamp(n, 0, 3)];
}

inline const char* qqrhChannelName(int n) {
    static const char* names[] = {"NONE", "RADIAL CURVATURE", "DYADIC RATIO",
        "LAPLACIAN 1", "LAPLACIAN 2", "MOVING BOUNDARY"};
    return names[std::clamp(n, 0, 5)];
}

inline const char* hostName(int n) {
    static const char* names[] = {"JULIA", "NEWTON", "TRICORN", "BURNING", "NOVA"};
    return names[std::clamp(n, 0, 4)];
}
inline const char* branchName(int n) {
    static const char* names[] = {"GREEDY", "SUB", "ADD"};
    return names[std::clamp(n, 0, 2)];
}
inline const char* couplingName(int n) {
    static const char* names[] = {"RAW", "SELFDUAL", "COMPACT"};
    return names[std::clamp(n, 0, 2)];
}
inline const char* metricName(int n) {
    static const char* names[] = {"candidate", "escape", "survival", "plus", "switch",
                                 "lyap", "radius", "twin", "desync"};
    return names[std::clamp(n, 0, 8)];
}
inline const char* precisionName(int n) {
    static const char* names[] = {"BASE", "HIGH", "ULTRA", "EXTREME", "MAX"};
    return names[std::clamp(n, 0, 4)];
}

inline int gridSize(const Settings& s) {
    static constexpr int base[] = {97, 129, 161};
    static constexpr double multiplier[] = {1.0, 1.25, 1.50, 2.00, 2.75};
    int n = static_cast<int>(std::llround(base[std::clamp(s.quality, 0, 2)] *
                                         multiplier[std::clamp(s.precision, 0, 4)]));
    if ((n & 1) == 0) ++n;
    return n;
}
inline int seedCount(const Settings& s) {
    static constexpr int base[] = {8, 12, 16};
    static constexpr double multiplier[] = {1.0, 1.50, 2.00, 3.00, 4.00};
    return std::max(1, static_cast<int>(std::llround(base[std::clamp(s.quality, 0, 2)] *
                                                   multiplier[std::clamp(s.precision, 0, 4)])));
}

// Keep scalar arithmetic in the same order as the archived HIP source.
struct Cx { double x, y; };
inline Cx add(Cx a, Cx b) { return {a.x + b.x, a.y + b.y}; }
inline Cx sub(Cx a, Cx b) { return {a.x - b.x, a.y - b.y}; }
inline Cx mul(Cx a, Cx b) { return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x}; }
inline Cx scl(Cx a, double s) { return {a.x * s, a.y * s}; }
inline double n2(Cx a) { return a.x * a.x + a.y * a.y; }
inline double dot(Cx a, Cx b) { return a.x * b.x + a.y * b.y; }
inline bool finiteScalar(double x) { return x == x && x > -1e300 && x < 1e300; }
inline bool finiteCx(Cx a) { return finiteScalar(a.x) && finiteScalar(a.y); }
inline Cx cdiv(Cx a, Cx b) {
    const double d = n2(b);
    if (d < 1e-30 || !finiteScalar(d)) return {1e300, 1e300};
    return {(a.x * b.x + a.y * b.y) / d, (a.y * b.x - a.x * b.y) / d};
}
inline Cx cpow(Cx z, int power) {
    Cx result{1, 0};
    for (int k = 0; k < power; ++k) result = mul(result, z);
    return result;
}
inline Cx conjc(Cx z) { return {z.x, -z.y}; }

struct Params {
    int host = JULIA, branch = GREEDY, coupling = RAW, power = 2, steps = 96;
    double cre = -.8, cim = .156;
    int nx = 97, ny = 97, seeds = 8;
    double re_min = -2, re_max = 2, im_min = -2, im_max = 2;
    int rhChannel = 0;
    C eta{.45, 0.0};
    std::shared_ptr<const ff_rh::Channel> rh;
    int qqrhChannel = 0;
    std::shared_ptr<const qqrh::Channel> qqrhCarrier;
};

inline Cx hval(Cx lambda, int mode) {
    if (mode == RAW) return lambda;
    const Cx denominator = add({1, 0}, lambda);
    if (n2(denominator) < 1e-20) return {1e300, 1e300};
    if (mode == COMPACT) return cdiv(scl(lambda, 2.0), denominator);
    return cdiv(scl(lambda, 4.0), mul(denominator, denominator));
}
inline int sigma_for(Cx F, Cx eta, int branch) {
    if (branch == SUBTRACT) return -1;
    if (branch == ADD) return 1;
    // Strict inequality: an exact zero dot product selects ADD, as before.
    return dot(F, eta) > 0 ? -1 : 1;
}
inline Cx host_map(Cx z, const Params& p) {
    const Cx c{p.cre, p.cim};
    if (p.qqrhChannel > 0) {
        const C point(z.x, z.y);
        const bool direct = p.qqrhChannel == 1 || p.qqrhChannel == 2;
        const C baseline = direct ? C(p.cre, p.cim) : point * point + C(p.cre, p.cim);
        // Zero eta bypasses every carrier, including poles at the direct
        // centering point: direct hosts reduce to c, the others to z^2+c.
        if (p.eta == C{}) return {baseline.real(), baseline.imag()};
        if (!p.qqrhCarrier || !p.qqrhCarrier->valid()) return {1e300, 1e300};
        const C response = p.qqrhCarrier->response(point);
        const C feedback = p.eta.imag() == 0.0
            ? p.eta.real() * response : p.eta * response;
        const C result = baseline + feedback;
        return {result.real(), result.imag()};
    }
    if (p.rhChannel > 0) {
        // The inner host is the existing RH feedback map, with its own eta.
        // The lambda-plane coordinate never enters this host construction.
        const C point(z.x, z.y);
        const C quadratic = point * point + C(p.cre, p.cim);
        if (p.eta == C{}) return {quadratic.real(), quadratic.imag()};
        if (!p.rh || !p.rh->valid()) return {1e300, 1e300};
        const C response = p.rh->response(point);
        const C feedback = p.eta.imag() == 0.0
            ? p.eta.real() * response : p.eta * response;
        const C result = quadratic + feedback;
        return {result.real(), result.imag()};
    }
    if (p.host == JULIA) return add(cpow(z, p.power), c);
    if (p.host == TRICORN) return add(cpow(conjc(z), p.power), c);
    if (p.host == BURNING)
        return add(cpow({std::fabs(z.x), std::fabs(z.y)}, p.power), c);
    const Cx zp = cpow(z, p.power), zm = cpow(z, p.power - 1);
    const Cx fp = scl(zm, static_cast<double>(p.power));
    if (n2(fp) < 1e-24) return {1e300, 1e300};
    const Cx newton = sub(z, cdiv(sub(zp, {1, 0}), fp));
    return p.host == NOVA ? add(newton, c) : newton;
}
inline Cx step(Cx z, Cx lambda, const Params& p, int& sigma) {
    if (p.rhChannel > 0 || p.qqrhChannel > 0) {
        const Cx H = hval(lambda, p.coupling);
        const bool forced = p.branch == SUBTRACT || p.branch == ADD;
        sigma = forced ? (p.branch == SUBTRACT ? -1 : 1) : 0;
        if (forced) {
            const Cx factor = add({1, 0}, scl(H, static_cast<double>(sigma)));
            // For a forced branch, cancel a zero outer factor first and
            // extend this fixed parameter slice as identity through host poles.
            // In particular forced SUB is exactly identity at lambda=1.
            if (factor.x == 0.0 && factor.y == 0.0) return z;
        }
        const Cx F = host_map(z, p);
        // An undefined greedy host has no branch decision.  Never call it a
        // frozen point just because an unselected subtract factor is zero.
        if (!finiteCx(F)) return F;
        const Cx displacement = sub(F, z);
        if (!forced) sigma = sigma_for(F, displacement, p.branch);
        // Avoid a round-trip F-z+z so lambda=0 is exactly the inner RH host.
        // The branch statistic has still been evaluated above.
        if (lambda.x == 0.0 && lambda.y == 0.0) return F;
        const Cx factor = add({1, 0}, scl(H, static_cast<double>(sigma)));
        if (factor.x == 0.0 && factor.y == 0.0) return z;
        return add(z, mul(factor, displacement));
    }
    // Archived branch is deliberately unchanged.
    const Cx F = host_map(z, p), eta = sub(F, z);
    sigma = sigma_for(F, eta, p.branch);
    const Cx factor = add({1, 0}, scl(hval(lambda, p.coupling), static_cast<double>(sigma)));
    return add(z, mul(factor, eta));
}

struct Accum {
    unsigned long long samples = 0, escaped = 0, steps = 0, plus = 0, switches = 0;
    unsigned long long bounded = 0, radius = 0, desync = 0, lyap_shifted = 0;
};
inline uint32_t mix32(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return x;
}
inline double u01(uint32_t x) { return (static_cast<double>(x & 0x00ffffffu) + .5) / 16777216.0; }

inline Accum scanAccum(const Params& p, int ix, int iy, const Cancel& cancel,
                       bool& aborted) {
    Accum total;
    const double xu = p.nx > 1 ? static_cast<double>(ix) / (p.nx - 1) : 0;
    const double yu = p.ny > 1 ? static_cast<double>(iy) / (p.ny - 1) : 0;
    const Cx lambda{p.re_min + (p.re_max - p.re_min) * xu,
                    p.im_min + (p.im_max - p.im_min) * yu};
    for (int seed = 0; seed < p.seeds; ++seed) {
        if (cancel && cancel()) { aborted = true; return total; }
        const uint32_t hash = mix32(static_cast<uint32_t>(seed) * 0x9e3779b9u +
            static_cast<uint32_t>(ix) * 0x85ebca6bu + static_cast<uint32_t>(iy) * 0xc2b2ae35u);
        const double angle = 2 * PI * u01(hash);
        const double radius = .25 + 1.65 * std::sqrt(u01(mix32(hash ^ 0xa511e9b3u)));
        Cx z{radius * std::cos(angle), radius * std::sin(angle)}, shadow{z.x + 1e-7, z.y};
        int previous = 0, plus = 0, switches = 0, desync = 0, steps = 0, escaped = 0;
        double logarithmicSeparation = 0;
        for (int n = 0; n < p.steps; ++n) {
            if ((n & 31) == 0 && cancel && cancel()) { aborted = true; return total; }
            int a = 0, b = 0;
            const Cx next = step(z, lambda, p, a), shadowNext = step(shadow, lambda, p, b);
            if (a > 0) ++plus;
            if (previous && a != previous) ++switches;
            previous = a;
            if (a != b) desync = 1;
            const double before = std::sqrt(n2(sub(shadow, z))) + 1e-30;
            const double after = std::sqrt(n2(sub(shadowNext, next))) + 1e-30;
            logarithmicSeparation += std::log(after / before);
            z = next; shadow = shadowNext; ++steps;
            const double r2 = n2(z), rt = n2(shadow);
            if (!finiteCx(z) || !finiteCx(shadow) || !finiteScalar(r2) ||
                !finiteScalar(rt) || r2 > 256 || rt > 256) {
                escaped = 1;
                break;
            }
            if (p.rhChannel == 0 && p.qqrhChannel == 0 && (p.host == NEWTON || p.host == NOVA) &&
                std::sqrt(n2(sub(cpow(z, p.power), {1, 0}))) < 1e-8) break;
        }
        ++total.samples;
        total.escaped += static_cast<unsigned long long>(escaped);
        total.steps += static_cast<unsigned long long>(steps);
        total.plus += static_cast<unsigned long long>(plus);
        total.switches += static_cast<unsigned long long>(switches);
        if (!escaped) {
            ++total.bounded;
            total.radius += static_cast<unsigned long long>(
                std::fmin(std::sqrt(n2(z)), 1e6) * 1e6 + .5);
        }
        // fmin/fmax deliberately retain the archived treatment of a NaN
        // shadow-log sum; replacing these with std::clamp would change it.
        const double average = std::fmax(-16.0, std::fmin(16.0,
            logarithmicSeparation / static_cast<double>(steps ? steps : 1)));
        total.lyap_shifted += static_cast<unsigned long long>((average + 16.0) * 1e6 + .5);
        total.desync += static_cast<unsigned long long>(desync);
    }
    return total;
}

struct Cell { std::array<double, METRICS> v{}; };
inline Cell metrics(const Accum& a, int budget) {
    Cell cell;
    const double samples = a.samples ? static_cast<double>(a.samples) : 1;
    const double steps = a.steps ? static_cast<double>(a.steps) : 1;
    const double bounded = a.bounded ? static_cast<double>(a.bounded) : 1;
    cell.v[1] = a.escaped / samples;
    cell.v[2] = a.steps / (samples * budget);
    cell.v[3] = a.plus / steps;
    cell.v[4] = a.switches / steps;
    cell.v[5] = a.lyap_shifted / (samples * 1e6) - 16.0;
    cell.v[6] = a.bounded ? a.radius / (bounded * 1e6) : 0;
    cell.v[7] = cell.v[5];
    cell.v[8] = a.desync / samples;
    return cell;
}

struct Grid {
    int side = 0, nx = 0, ny = 0, seeds = 0, steps = 0;
    bool cancelled = false;
    bool valid = false;
    C center{};
    double span = 4.0;
    Params params{};
    std::vector<Cell> cells;
    std::array<std::pair<double, double>, METRICS> ranges{};
    std::string diagnostic;
};

inline bool validSettings(const Settings& s) {
    const bool legacyFields = s.host >= 0 && s.host <= 4 && s.branch >= 0 && s.branch <= 2 &&
        s.coupling >= 0 && s.coupling <= 2 && s.power >= 2 && s.power <= 8 &&
        s.quality >= 0 && s.quality <= 2 && s.precision >= 0 && s.precision <= 4 &&
        s.iterations >= 32 && s.iterations <= 4096 && s.metric >= 0 && s.metric < METRICS &&
        finiteCx({s.c.real(), s.c.imag()});
    if (!legacyFields || s.rhChannel < 0 || s.rhChannel > 3 ||
        s.qqrhChannel < 0 || s.qqrhChannel > 5 ||
        (s.rhChannel > 0 && s.qqrhChannel > 0)) return false;
    if (s.rhChannel == 0 && s.qqrhChannel == 0) return true;
    if (s.qqrhChannel > 0 && (!finiteScalar(s.qqrhSlope) ||
        s.qqrhSlope < 0.0 || s.qqrhSlope > 64.0)) return false;
    return s.host == JULIA && s.power == 2 &&
        (s.rhN == 8 || s.rhN == 16 || s.rhN == 32) &&
        finiteScalar(s.rhX) && s.rhX >= -1.0 && s.rhX <= 1.0 &&
        finiteScalar(s.rhB) && s.rhB > 0.0 && finiteScalar(s.rhH) &&
        finiteCx({s.eta.real(), s.eta.imag()});
}

inline Grid scanGrid(const Settings& settings, C center, double span,
                     int requestedGrid = 0, Cancel cancelled = {},
                     int seedCap = 0, int stepCap = 0, unsigned requestedThreads = 0) {
    Grid grid;
    grid.center = center; grid.span = span;
    if (!validSettings(settings) || !finiteCx({center.real(), center.imag()}) ||
        !finiteScalar(span) || !(span > 0)) {
        grid.diagnostic = "Invalid legacy lambda-map settings or view";
        return grid;
    }
    const int target = gridSize(settings);
    int side = requestedGrid > 0 ? std::min(target, std::max(17, requestedGrid)) : target;
    if ((side & 1) == 0) ++side;
    grid.side = grid.nx = grid.ny = side;
    grid.seeds = seedCap > 0 ? std::min(seedCount(settings), std::max(1, seedCap)) : seedCount(settings);
    grid.steps = stepCap > 0 ? std::min(settings.iterations, std::max(16, stepCap)) : settings.iterations;
    Params& p = grid.params;
    p.host = settings.host; p.branch = settings.branch; p.coupling = settings.coupling;
    p.power = settings.power; p.steps = grid.steps; p.cre = settings.c.real(); p.cim = settings.c.imag();
    p.rhChannel = settings.rhChannel; p.eta = settings.eta;
    p.qqrhChannel = settings.qqrhChannel;
    if (p.rhChannel > 0) {
        // One immutable channel per scan, shared by every cell and worker.
        p.rh = std::make_shared<const ff_rh::Channel>(p.rhChannel - 1, settings.rhN,
            p.rhChannel == 3 ? -1.0 : settings.rhX, settings.rhB, settings.rhH);
        if (!p.rh->valid()) {
            grid.diagnostic = "Invalid RH polynomial channel";
            return grid;
        }
    }
    if (p.qqrhChannel > 0) {
        p.qqrhCarrier = std::make_shared<const qqrh::Channel>(p.qqrhChannel,
            settings.rhN, settings.rhX, settings.rhB, settings.rhH, settings.qqrhSlope);
        // Structural validation only.  A direct centering pole is deferred
        // until active host evaluation so eta=0 and the forced identity slice
        // remain available without evaluating the inactive rational carrier.
        if (!p.qqrhCarrier->valid()) {
            grid.diagnostic = "Invalid finite QQRH carrier";
            return grid;
        }
    }
    p.nx = p.ny = side; p.seeds = grid.seeds;
    p.re_min = center.real() - .5 * span; p.re_max = center.real() + .5 * span;
    p.im_min = center.imag() - .5 * span; p.im_max = center.imag() + .5 * span;
    if (!finiteScalar(p.re_min) || !finiteScalar(p.re_max) || !finiteScalar(p.im_min) ||
        !finiteScalar(p.im_max) || !(p.re_max > p.re_min) || !(p.im_max > p.im_min)) {
        grid.diagnostic = "View endpoints are not distinguishable finite doubles";
        return grid;
    }
    const std::size_t count = static_cast<std::size_t>(side) * side;
    grid.cells.resize(count);
    std::atomic<std::size_t> next{0};
    std::atomic<bool> aborted{false};
    const Cancel stop = [&] {
        if (aborted.load(std::memory_order_relaxed)) return true;
        if (cancelled && cancelled()) {
            aborted.store(true, std::memory_order_relaxed);
            return true;
        }
        return false;
    };
    auto worker = [&] {
        while (!stop()) {
            const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
            if (index >= count) return;
            const int ix = static_cast<int>(index % side), iy = static_cast<int>(index / side);
            bool interrupted = false;
            const Accum raw = scanAccum(p, ix, iy, stop, interrupted);
            if (interrupted) {
                aborted.store(true, std::memory_order_relaxed);
                return;
            }
            grid.cells[index] = metrics(raw, p.steps);
        }
    };
    // Independent cells preserve exact integer reduction while using native
    // CPU cores.  Every seed and every 32 orbit steps check cancellation.
    const unsigned reported = std::thread::hardware_concurrency();
    const unsigned threads = requestedThreads ? std::clamp(requestedThreads, 1u, 64u)
                                             : std::min(16u, std::max(1u, reported));
    std::vector<std::thread> pool;
    pool.reserve(threads > 0 ? threads - 1 : 0);
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(worker);
    worker();
    for (auto& thread : pool) thread.join();
    if (stop()) {
        grid.cancelled = true;
        grid.diagnostic = "Cancelled before a complete metric grid was available";
        return grid;
    }

    // Candidate is the original eight-channel neighbor gradient, without
    // coordinate-spacing division.  The duplicate lyap/twin terms are kept.
    for (int y = 0; y < side; ++y) {
        if (stop()) { grid.cancelled = true; return grid; }
        for (int x = 0; x < side; ++x) {
            const int xm = std::max(0, x - 1), xp = std::min(side - 1, x + 1);
            const int ym = std::max(0, y - 1), yp = std::min(side - 1, y + 1);
            double value = 0;
            for (int m = 1; m < METRICS; ++m) {
                const double dx = grid.cells[static_cast<std::size_t>(y) * side + xp].v[m] -
                                  grid.cells[static_cast<std::size_t>(y) * side + xm].v[m];
                const double dy = grid.cells[static_cast<std::size_t>(yp) * side + x].v[m] -
                                  grid.cells[static_cast<std::size_t>(ym) * side + x].v[m];
                value += dx * dx + dy * dy;
            }
            grid.cells[static_cast<std::size_t>(y) * side + x].v[0] = std::sqrt(value / 16.0);
        }
    }
    std::vector<double> values(count);
    for (int m = 0; m < METRICS; ++m) {
        if (stop()) { grid.cancelled = true; return grid; }
        for (std::size_t i = 0; i < count; ++i) values[i] = grid.cells[i].v[m];
        std::sort(values.begin(), values.end());
        double lo = values[static_cast<std::size_t>(.02 * (values.size() - 1))];
        double hi = values[static_cast<std::size_t>(.98 * (values.size() - 1))];
        if (!(hi > lo)) hi = lo + 1;
        grid.ranges[m] = {lo, hi};
    }
    grid.valid = true;
    grid.diagnostic = p.qqrhChannel > 0
        ? "Finite QQRH carrier inside original BKQR deformation; direct channels centered at y=h; eta inner, lambda outer; finite diagnostics only"
        : p.rhChannel == 0
            ? "Legacy ensemble diagnostics; finite/numerical stop counts retained, no stability certificate"
            : "RH finite-polynomial host inside original BKQR deformation; eta is inner feedback, lambda is outer coordinate; finite diagnostics only";
    return grid;
}

inline Color viridis(double t) {
    static constexpr Color stops[] = {
        {{.267, .005, .329}}, {{.283, .141, .458}}, {{.254, .265, .530}},
        {{.207, .372, .553}}, {{.164, .471, .558}}, {{.128, .567, .551}},
        {{.135, .659, .518}}, {{.267, .749, .441}}, {{.478, .821, .318}},
        {{.741, .873, .150}}
    };
    t = std::clamp(t, 0.0, 1.0);
    const double position = t * 9;
    const int index = std::clamp(static_cast<int>(std::floor(position)), 0, 8);
    const double fraction = std::clamp(position - index, 0.0, 1.0);
    Color color{};
    for (int channel = 0; channel < 3; ++channel)
        color[channel] = stops[index][channel] +
            (stops[index + 1][channel] - stops[index][channel]) * fraction;
    return color;
}

inline Color metricCell(const Grid& grid, int ix, int iy, int metric) {
    if (!grid.valid || grid.cancelled || grid.cells.empty() || grid.side < 1)
        return {{.01, .012, .018}};
    metric = std::clamp(metric, 0, METRICS - 1);
    ix = std::clamp(ix, 0, grid.nx - 1); iy = std::clamp(iy, 0, grid.ny - 1);
    const auto range = grid.ranges[metric];
    const double value = grid.cells[static_cast<std::size_t>(iy) * grid.nx + ix].v[metric];
    return viridis((value - range.first) / (range.second - range.first));
}

} // namespace bkqrlambda
