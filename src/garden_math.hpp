#pragma once
#include "crown_math.hpp"
#include "restart_math.hpp"
#include "cayley_math.hpp"
#include "lfunction_math.hpp"
#include "grh_cayley_math.hpp"
#include "primitive_pair_math.hpp"
#include "grh_phase_feedback_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace garden {
using C = std::complex<double>;
struct Term {
    int exponent = 0, a = 0, b = 0, sign = 1;
};
struct Data {
    int family = 0, generation = 1, degree = 0;
    std::vector<Term> terms;
    std::vector<crown::Cylinder> cylinders;
    restart::Transfer transfer;
    std::vector<uint64_t> connected, choose;
    // Restart: monic coefficient, return-delay weight, unused, unused.
    // CROWN: odd residue, number of variable bits, unused, unused.
    // GP: exponent, power(a), power(b), sign.
    // Percolation: count, binomial, first derivative, second derivative, divided
    // by 2^M.
    // Cayley and L-function families: one unused record keeps GPU upload valid.
    std::vector<std::array<double, 4>> records;
};
inline C power(C z, int n) {
    C p{1, 0};
    while (n) {
        if (n & 1)
            p *= z;
        z *= z;
        n >>= 1;
    }
    return p;
}
using Partition = std::array<int, 5>; // permanent source + four current-layer vertices
inline Partition canonical(Partition p) {
    std::array<int, 8> seen;
    seen.fill(-1);
    int n = 0;
    for (auto &x : p) {
        if (seen[x] < 0)
            seen[x] = n++;
        x = seen[x];
    }
    return p;
}
inline Partition join(Partition p, int a, int b) {
    int from = p[b], to = p[a];
    for (auto &x : p)
        if (x == from)
            x = to;
    return canonical(p);
}
inline Partition replaceVertex(Partition p, int at, bool open) {
    // The old vertex has no further incident edges. Introduce its successor;
    // an open longitudinal edge retains its component, a closed edge starts a new
    // one.
    if (!open)
        p[at] = 6;
    return canonical(p);
}
inline std::vector<uint64_t> connectivityCounts(int generation) {
    if (generation < 1 || generation > 6)
        throw std::runtime_error("Percolation generation must be 1..6");
    using Counts = std::vector<uint64_t>;
    std::map<Partition, Counts> states;
    states[Partition{0, 0, 1, 2, 3}] = {1};
    int edges = 0;
    auto edge = [&](int a, int b, bool longitudinal) {
        std::map<Partition, Counts> next;
        for (const auto &item : states)
            for (int open = 0; open <= 1; open++) {
                Partition p = longitudinal ? replaceVertex(item.first, a, open != 0)
                                           : (open ? join(item.first, a, b) : item.first);
                auto &v = next[p];
                if (v.empty())
                    v.resize(edges + 2);
                for (size_t k = 0; k < item.second.size(); k++)
                    v[k + open] += item.second[k];
            }
        states = std::move(next);
        edges++;
    };
    auto square = [&] {
        edge(1, 2, false);
        edge(1, 3, false);
        edge(2, 4, false);
        edge(3, 4, false);
    };
    square();
    for (int layer = 1; layer <= generation; layer++) {
        for (int j = 1; j <= 4; j++)
            edge(j, 0, true);
        square();
    }
    Counts total(edges + 1);
    for (const auto &item : states)
        if (item.first[0] == item.first[4])
            for (int k = 0; k <= edges; k++)
                total[k] += item.second[k];
    return total;
}
inline Data build(int family, int generation) {
    Data d;
    d.family = family;
    d.generation = generation;
    if (family == 1) {
        if (generation < 1 || generation > 4)
            throw std::runtime_error("GP generation must be 1..4");
        std::vector<Term> p{{0, 0, 0, 1}};
        for (int n = 0; n < generation; n++) {
            std::vector<Term> next{{0, 0, 0, 1}};
            for (auto t : p) {
                next.push_back({1 + 2 * t.exponent, t.a + 1, t.b, 1});
                next.push_back({4 + 4 * t.exponent, t.a, t.b + 1, 1});
            }
            p = std::move(next);
        }
        for (auto t : p) {
            d.terms.push_back(t);
            t.exponent++;
            t.sign = -1;
            d.terms.push_back(t);
            d.degree = std::max(d.degree, t.exponent);
        }
        for (auto t : d.terms)
            d.records.push_back({double(t.exponent), double(t.a), double(t.b), double(t.sign)});
    } else if (family == 2) {
        d.connected = connectivityCounts(generation);
        d.degree = int(d.connected.size()) - 1;
        int m = d.degree;
        d.choose.resize(m + 1);
        d.choose[0] = 1;
        for (int k = 1; k <= m; k++)
            d.choose[k] = d.choose[k - 1] * uint64_t(m - k + 1) / uint64_t(k);
        std::vector<int64_t> first(m + 1), second(m + 1);
        for (int k = 0; k < m; k++)
            first[k] =
                int64_t(k + 1) * int64_t(d.connected[k + 1]) - int64_t(m - k) * int64_t(d.connected[k]);
        for (int k = 0; k < m - 1; k++)
            second[k] = int64_t(k + 1) * first[k + 1] - int64_t(m - 1 - k) * first[k];
        double norm = std::ldexp(1., m);
        for (int k = 0; k <= m; k++)
            d.records.push_back(
                {d.connected[k] / norm, d.choose[k] / norm, first[k] / norm, second[k] / norm});
    } else if (family == 3) {
        d.cylinders = crown::cylinders(generation);
        d.degree = int(d.cylinders.size());
        for (const auto &w : d.cylinders)
            d.records.push_back({double(w.residue), double(w.bits), 0, 0});
    } else if (family == 4) {
        if (generation < 1 || generation > 4)
            throw std::runtime_error("Restart generation must be 1..4");
        d.transfer = restart::transfer(1 << (generation + 3));
        d.degree = d.transfer.degree;
        for (int k = 0; k <= d.degree; k++)
            d.records.push_back(
                {k == d.degree ? 1. : 0.,
                 k == d.degree ? 0. : double(d.transfer.counts[d.degree - k]) / d.transfer.states, 0, 0});
    } else if (family == 5 || family == 6 || family == 7 || family == 8 || family == 9 || family == 10) {
        if (generation < 1 || generation > 4)
            throw std::runtime_error("Cayley / L-function generation must be 1..4");
        if ((family == 8 || family == 9 || family == 10) && generation != 1)
            throw std::runtime_error(family == 8 ? "Primitive pair generation must be 1" : family == 9 ? "GRH observable generation must be 1" : "GRH phase-feedback generation must be 1");
        d.degree = family == 5 ? 1 << generation : 0;
        d.records.push_back({0, 0, 0, 0});
    } else
        throw std::runtime_error("Unknown garden family");
    return d;
}
inline std::shared_ptr<const Data> data(int family, int generation) {
    static std::mutex mutex;
    static std::map<std::pair<int, int>, std::shared_ptr<const Data>> cache;
    std::lock_guard<std::mutex> lock(mutex);
    auto key = std::make_pair(family, generation);
    auto it = cache.find(key);
    if (it != cache.end())
        return it->second;
    auto d = std::make_shared<Data>(build(family, generation));
    cache[key] = d;
    return d;
}
struct Evaluation {
    C f{}, first{}, second{}, base{1, 0};
    double magnitude = 0;
};
inline Evaluation evaluate(C z, C a, C b, const Data &d) {
    Evaluation e;
    if (d.family == 5) {
        const auto jets = cayley::preimage(z, d.generation);
        e.f = jets.p.f - a * jets.q.f;
        e.first = jets.p.first - a * jets.q.first;
        e.second = jets.p.second - a * jets.q.second;
        // This projective residual remains defined for a=0 at repeated roots.
        // A denominator zero Q=0 has error one, and cannot look like a root.
        e.magnitude = std::abs(jets.p.f) + (1. + std::abs(a)) * std::abs(jets.q.f);
    } else if (d.family == 6) {
        const auto value = lfunction::evaluate(z, d.generation);
        if (!value.valid)
            return e;
        e.f = value.f - a;
        e.first = value.first;
        e.second = value.second;
        e.magnitude = 1.; // L-function termination uses an absolute residual.
    } else if (d.family == 7) {
        const auto value = grh_cayley::evaluate(z, a, d.generation);
        if (!value.valid)
            return e;
        e.f = value.f;
        e.first = value.first;
        e.second = value.second;
        e.magnitude = 1.; // normalized completed-L residual = |L-target|.
    } else if (d.family == 8) {
        const auto value = primitive_pair::evaluate(z, a);
        if (!value.valid)
            return e;
        e.f = value.f;
        e.first = value.first;
        e.second = value.second;
        e.magnitude = 1.; // normalized completed-pair residual.
    } else if (d.family == 1) {
        int degree = d.degree;
        if (b == C(0, 0))
            degree = a == C(0, 0) ? 1 : (1 << d.generation);
        bool reverse = std::norm(z) > 1.;
        C w = reverse ? 1. / z : z;
        std::array<C, 5> ap, bp;
        ap[0] = bp[0] = {1, 0};
        for (int i = 1; i <= d.generation; i++) {
            ap[i] = ap[i - 1] * a;
            bp[i] = bp[i - 1] * b;
        }
        bool tiny = std::norm(z) < 1e-12;
        for (auto t : d.terms) {
            C coefficient = double(t.sign) * ap[t.a] * bp[t.b];
            if (coefficient == C(0, 0))
                continue;
            int k = reverse ? degree - t.exponent : t.exponent;
            C v = coefficient * power(w, k);
            e.f += v;
            e.magnitude += std::abs(v);
            if (tiny) {
                if (t.exponent)
                    e.first += double(t.exponent) * coefficient * power(z, t.exponent - 1);
                if (t.exponent > 1)
                    e.second +=
                        double(t.exponent * (t.exponent - 1)) * coefficient * power(z, t.exponent - 2);
            } else {
                e.first += double(t.exponent) * v;
                e.second += double(t.exponent * (t.exponent - 1)) * v;
            }
        }
        e.base = tiny ? C(1, 0) : z;
    } else if (d.family == 4) {
        bool reverse = std::norm(z) > 1., tiny = std::norm(z) < 1e-12;
        C w = reverse ? 1. / z : z;
        for (int i = d.degree; i >= 0; i--) {
            int k = reverse ? d.degree - i : i;
            C coefficient = d.records[k][0] - a * d.records[k][1];
            e.f = e.f * w + coefficient;
            e.magnitude = e.magnitude * std::abs(w) + std::abs(coefficient);
            if (tiny) {
                if (k)
                    e.first += double(k) * coefficient * power(z, k - 1);
                if (k > 1)
                    e.second += double(k * (k - 1)) * coefficient * power(z, k - 2);
            } else {
                e.first = e.first * w + double(k) * coefficient;
                e.second = e.second * w + double(k * (k - 1)) * coefficient;
            }
        }
        e.base = tiny ? C(1, 0) : z;
    } else {
        int m = d.degree;
        bool reverse = std::norm(z) > std::norm(1. - z);
        e.base = reverse ? z : 1. - z;
        C w = reverse ? (1. - z) / z : z / (1. - z);
        for (int i = m; i >= 0; i--) {
            int k = reverse ? m - i : i;
            C value = d.records[k][0] - a * d.records[k][1];
            e.f = e.f * w + value;
            e.magnitude = e.magnitude * std::abs(w) + std::abs(value);
        }
        for (int i = m - 1; i >= 0; i--) {
            int k = reverse ? m - 1 - i : i;
            e.first = e.first * w + d.records[k][2];
        }
        for (int i = m - 2; i >= 0; i--) {
            int k = reverse ? m - 2 - i : i;
            e.second = e.second * w + d.records[k][3];
        }
    }
    return e;
}
using Step = crown::Correction;
inline Step correction(C z, C a, C b, const Data &d, int method) {
    if (d.family == 3)
        return crown::correction(z, crown::roots(d.cylinders, a), b, method);
    if (d.family == 4 && a == C(0, 0)) {
        Step s;
        s.error = std::abs(z) / (1 + std::abs(z));
        s.delta = b * z * (method ? 2. / (d.degree + 1) : 1. / d.degree);
        s.valid = std::isfinite(s.delta.real()) && std::isfinite(s.delta.imag());
        return s;
    }
    if (d.family == 10) {
        auto f = grh_phase_feedback::correction(z, a, b, method);
        Step s;
        s.error = f.residual;
        s.delta = f.delta;
        s.label = f.phase;
        s.valid = f.valid;
        return s;
    }
    auto e = evaluate(z, a, b, d);
    Step s;
    if (!(e.magnitude > 0) || !std::isfinite(e.magnitude))
        return s;
    s.error = std::abs(e.f) / e.magnitude;
    // Exact multiple roots have f=f'=f''=0. Accept them before derivative
    // normalization, while leaving the established families unchanged.
    if ((d.family == 5 || d.family == 6 || d.family == 7 || d.family == 8 || d.family == 10) && e.f == C(0, 0)) {
        s.valid = true;
        return s;
    }
    double norm = std::max({std::abs(e.f), std::abs(e.first), std::abs(e.second)});
    if (!(norm > 0) || !std::isfinite(norm))
        return s;
    e.f /= norm;
    e.first /= norm;
    e.second /= norm;
    C denominator = method == 1 ? 2. * e.first * e.first - e.f * e.second : e.first;
    if (std::norm(denominator) < 1e-26) {
        s.valid = s.error < 3e-6;
        return s;
    }
    s.delta = e.base * (method == 1 ? 2. * e.f * e.first : e.f) / denominator;
    if (d.family == 2 || d.family == 4 || d.family == 5 || d.family == 6 || d.family == 7 || d.family == 8)
        s.delta *= b;
    s.valid = std::isfinite(s.delta.real()) && std::isfinite(s.delta.imag());
    return s;
}
struct Orbit {
    C z;
    int status = 1, iterations = 0;
    double error = 1, label = 0, path = 0, phase = 0;
};
inline Orbit orbit(C z, C a, C b, const Data &d, int method, int iterations, std::vector<C> *trail = nullptr,
                   const std::vector<C> *cachedRoots = nullptr) {
    Orbit o;
    o.z = z;
    std::vector<C> localRoots;
    if (d.family == 3 && !cachedRoots) {
        localRoots = crown::roots(d.cylinders, a);
        cachedRoots = &localRoots;
    }
    for (int i = 0; i < iterations; i++) {
        if (trail && trail->size() < 192)
            trail->push_back(z);
        if (!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z) > 1e7) {
            o.status = 3;
            break;
        }
        if ((d.family == 6 || d.family == 8) && (z.real() < -8. || z.real() > 12. || std::abs(z.imag()) > 80.)) {
            o.status = 3;
            break;
        }
        if (d.family == 7 || d.family == 10) {
            C mapped;
            if (!grh_cayley::sheetPoint(z, d.family == 10 ? grh_phase_feedback::sheet : d.generation, mapped)) {
                o.status = 3;
                break;
            }
        }
        auto step =
            d.family == 3 ? crown::correction(z, *cachedRoots, b, method) : correction(z, a, b, d, method);
        o.label = step.label;
        o.error = step.error;
        o.iterations = i;
        if (step.error < 3e-6 && step.valid) {
            if (d.family == 4 && a == C(0, 0))
                z = {0, 0};
            o.status = 0;
            break;
        }
        if (!step.valid) {
            o.status = 2;
            break;
        }
        if (d.family == 10) {
            o.path += std::abs(step.delta);
            o.phase += step.label;
        }
        z -= step.delta;
        o.z = z;
        o.iterations = i + 1;
    }
    o.z = z;
    return o;
}
inline std::array<double, 3> colour(const Orbit &o, int palette, double hue, double bands, int family) {
    auto rainbow = [](double t) {
        constexpr double tau = 6.283185307179586;
        return std::array<double, 3>{.52 + .46 * std::cos(tau * t), .52 + .46 * std::cos(tau * (t + .33)),
                                     .52 + .46 * std::cos(tau * (t + .67))};
    };
    if (o.status == 2)
        return {.19, .095, .065};
    if (o.status == 3)
        return {.09, .10, .13};
    if (o.status != 0) {
        auto c = rainbow(hue + .72 + .015 * std::log1p(o.error) + (family == 10 ? .025*std::sin(.11*o.phase) : 0));
        for (auto &x : c)
            x *= family == 10 ? (.07 + .10*(1-std::exp(-.3*o.path))) : .12;
        return c;
    }
    if (family == 10) {
        double basin = .08 + .53*(.5 + .5*std::tanh(3.0*o.z.imag())) + .06*std::sin(3.0*std::arg(o.z));
        double ribbon = .5 + .5*std::cos(bands*(2.4*o.path + .75*o.iterations + .20*o.phase));
        double t=hue+basin;
        double light=.30 + .64*(1-std::exp(-.16*std::max(0,64-o.iterations)))*(.72+.28*ribbon);
        if (palette==1) {
            t=hue+.03+bands*.12*std::log1p(o.iterations)+.025*std::sin(.08*o.phase);
            light=.76*(.82+.18*ribbon);
        } else if (palette==2) {
            t=hue+.10-bands*.065*std::log10(std::max(o.error,3e-7))+.035*std::sin(.12*o.phase);
            light=.70*(.80+.20*ribbon);
        }
        auto c=rainbow(t);
        for(auto &x:c) x*=light;
        return c;
    }
    C root = o.z - C(family == 2 ? .5 : 0, 0);
    double angle = std::arg(root) / 6.283185307179586;
    double t = hue + (family == 3 ? o.label : angle + .06 * std::log1p(std::abs(root)));
    if ((family == 6 || family == 7 || family == 8) && palette == 0)
        t = hue + std::fmod(root.imag() * .071 + root.real() * .035, 1.);
    double light = .37 + .44 * std::exp(-o.iterations * .035);
    if (palette == 0)
        light *= .84 + .16 * std::cos(std::log1p(o.iterations) * bands * 3.8);
    if (palette == 1) {
        t = hue + .04 + bands * std::log1p(o.iterations) * .13;
        light = .78;
    }
    if (palette == 2) {
        t = hue + .08 - bands * std::log10(std::max(o.error, family >= 3 ? 3e-7 : 1e-16)) * .065;
        light = .72;
    }
    auto c = rainbow(t);
    for (auto &x : c)
        x *= light;
    return c;
}
} // namespace garden
