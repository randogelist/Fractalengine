#pragma once
#include "amp_core.hpp"
#include <iostream>

namespace garden_tests {
using Z = std::complex<long double>;
struct Jet {
    Z f{}, d{}, dd{};
};
inline Jet add(Jet a, Jet b) { return {a.f + b.f, a.d + b.d, a.dd + b.dd}; }
inline Jet mul(Jet a, Jet b) {
    return {a.f * b.f, a.d * b.f + a.f * b.d, a.dd * b.f + 2.L * a.d * b.d + a.f * b.dd};
}
inline Jet scale(Jet a, Z b) { return {a.f * b, a.d * b, a.dd * b}; }
inline Jet pow(Jet a, int n) {
    Jet p{{1, 0}};
    while (n) {
        if (n & 1)
            p = mul(p, a);
        a = mul(a, a);
        n >>= 1;
    }
    return p;
}
inline Jet recursion(Jet z, int n, Z a, Z b) {
    if (n == 0)
        return {{1, 0}};
    return add(Jet{{1, 0}}, add(scale(mul(z, recursion(pow(z, 2), n - 1, a, b)), a),
                                scale(mul(pow(z, 4), recursion(pow(z, 4), n - 1, a, b)), b)));
}
inline std::vector<uint64_t> bruteConnectivity(int length) {
    std::vector<std::pair<int, int>> edges;
    for (int x = 0; x <= length; x++) {
        for (auto pair : {std::pair<int, int>{0, 1}, {0, 2}, {1, 3}, {2, 3}})
            edges.push_back({4 * x + pair.first, 4 * x + pair.second});
        if (x)
            for (int j = 0; j < 4; j++)
                edges.push_back({4 * (x - 1) + j, 4 * x + j});
    }
    int m = int(edges.size());
    std::vector<uint64_t> counts(m + 1);
    for (uint32_t mask = 0; mask < (1u << m); mask++) {
        std::array<int, 12> parent;
        for (int i = 0; i < 4 * (length + 1); i++)
            parent[i] = i;
        auto find = [&](int n) {
            while (parent[n] != n)
                n = parent[n];
            return n;
        };
        int k = 0;
        for (int i = 0; i < m; i++)
            if (mask & (1u << i)) {
                k++;
                auto e = edges[i];
                parent[find(e.first)] = find(e.second);
            }
        if (find(0) == find(4 * length + 3))
            counts[k]++;
    }
    return counts;
}
inline uint64_t modularPower(uint64_t a, int n, uint64_t prime) {
    uint64_t result = 1;
    while (n) {
        if (n & 1)
            result = result * a % prime;
        a = a * a % prime;
        n >>= 1;
    }
    return result;
}
inline uint64_t transferDeterminant(const restart::Transfer &t, uint64_t z, uint64_t q, uint64_t prime) {
    // Independent dense determinant of M*(zI-T(q)), over a prime field.
    int m = t.states;
    std::vector<std::vector<uint64_t>> a(m, std::vector<uint64_t>(m));
    for (int i = 0; i < m; i++)
        a[i][i] = uint64_t(m) * z % prime;
    for (int j = 0; j < m; j++) {
        uint64_t n = 2 * j + 1, next = 3 * n + 1;
        while (!(next & 1))
            next >>= 1;
        if (n == 1 || next > uint64_t(2 * m - 1)) {
            for (int i = 0; i < m; i++)
                a[i][j] = (a[i][j] + prime - q) % prime;
        } else {
            int i = int((next - 1) / 2);
            a[i][j] = (a[i][j] + prime - m) % prime;
        }
    }
    uint64_t determinant = 1;
    for (int i = 0; i < m; i++) {
        int pivot = i;
        while (pivot < m && a[pivot][i] == 0)
            pivot++;
        if (pivot == m)
            return 0;
        if (pivot != i) {
            std::swap(a[pivot], a[i]);
            determinant = (prime - determinant) % prime;
        }
        determinant = determinant * a[i][i] % prime;
        uint64_t inverse = modularPower(a[i][i], int(prime - 2), prime);
        for (int k = i + 1; k < m; k++) {
            uint64_t factor = a[k][i] * inverse % prime;
            for (int j = i + 1; j < m; j++)
                a[k][j] = (a[k][j] + prime - factor * a[i][j] % prime) % prime;
            a[k][i] = 0;
        }
    }
    return determinant;
}
inline int run(std::ostream &out) {
    int checks = 0;
    auto check = [&](bool yes, const char *name) {
        checks++;
        if (!yes)
            throw std::runtime_error(std::string("Garden test: ") + name);
    };
    for (int n = 1; n <= 4; n++) {
        auto data = garden::data(1, n);
        int expected[] = {5, 21, 85, 341};
        check(data->degree == expected[n - 1], "GP polynomial degree");
        for (int i = 0; i < 48; i++) {
            Z z{1.43L * std::sin(.73L * i) + .07L, 1.12L * std::cos(.37L * i)}, a{.82L, .13L}, b{.70L, -.12L};
            Jet variable{z, {1, 0}, {0, 0}};
            Jet g = mul(add(Jet{{1, 0}}, scale(variable, {-1, 0})), recursion(variable, n, a, b));
            for (int method = 0; method <= 1; method++) {
                Z expected = method ? 2.L * g.f * g.d / (2.L * g.d * g.d - g.f * g.dd) : g.f / g.d;
                auto step = garden::correction({double(z.real()), double(z.imag())},
                                               {double(a.real()), double(a.imag())},
                                               {double(b.real()), double(b.imag())}, *data, method);
                Z actual{step.delta.real(), step.delta.imag()};
                check(step.valid && std::abs(actual - expected) < 2e-8L * (1 + std::abs(expected)),
                      "Sparse GP derivative vs independent recursive jets");
            }
        }
        for (amp::C z : {amp::C(0, 0), amp::C(.2, .1), amp::C(3, 2)}) {
            auto step = garden::correction(z, {0, 0}, {0, 0}, *data, 0);
            check(step.valid && std::abs(step.delta - (z - amp::C(1, 0))) < 1e-10,
                  "Zero GP weights reduce to 1-z");
        }
    }
    for (int n = 1; n <= 2; n++) {
        auto expected = bruteConnectivity(n);
        auto actual = garden::connectivityCounts(n);
        check(expected == actual, "Transfer counts equal exhaustive bond enumeration");
    }
    for (int n = 1; n <= 6; n++) {
        auto data = garden::data(2, n);
        int m = data->degree;
        check(m == 4 + 8 * n, "Cubic-strip bond count");
        check(data->connected[0] == 0 && data->connected[m] == 1, "Connectivity endpoints");
        for (int k = 0; k <= m; k++)
            check(data->connected[k] <= data->choose[k], "Counts bounded by all k-edge subsets");
        if (n > 4)
            continue;
        for (int i = 0; i < 24; i++) {
            Z z{.5L + .73L * std::sin(i * .83L), .47L * std::cos(i * .61L)}, target{.43L, .07L};
            Jet variable{z, {1, 0}, {0, 0}}, other = add(Jet{{1, 0}}, scale(variable, {-1, 0})), f{};
            for (int k = 0; k <= m; k++)
                f = add(f, scale(mul(pow(variable, k), pow(other, m - k)),
                                 {static_cast<long double>(data->connected[k]), 0}));
            f = add(f, Jet{-target});
            for (int method = 0; method <= 1; method++) {
                Z expected = method ? 2.L * f.f * f.d / (2.L * f.d * f.d - f.f * f.dd) : f.f / f.d;
                auto step =
                    garden::correction({double(z.real()), double(z.imag())},
                                       {double(target.real()), double(target.imag())}, {1, 0}, *data, method);
                Z actual{step.delta.real(), step.delta.imag()};
                check(step.valid && std::abs(actual - expected) < 2e-7L * (1 + std::abs(expected)),
                      "Bernstein derivatives vs independent polynomial jets");
            }
        }
    }
    for (int depth = 1, count = 3; depth <= 4; depth++, count *= 3) {
        auto d = garden::data(3, depth);
        check(d->degree == count, "CROWN cylinder count is 3^depth");
        auto prefix = depth > 1 ? crown::cylinders(depth - 1) : std::vector<crown::Cylinder>{};
        for (size_t i = 0; i < d->cylinders.size(); i++) {
            const auto &w = d->cylinders[i];
            for (uint64_t lift = 0; lift < 4; lift++) {
                uint64_t n = w.residue + (lift << (w.bits + 1));
                bool agrees = true;
                for (int j = 0; j < depth; j++) {
                    n = 3 * n + 1;
                    int a = 0;
                    while (!(n & 1)) {
                        n >>= 1;
                        a++;
                    }
                    agrees = agrees && a == w.word[j];
                }
                check(agrees, "CROWN residue and its lifts have the exact Syracuse valuation word");
            }
            if (depth > 1) {
                auto p = prefix[i / 3];
                check((w.residue & ((1u << (p.bits + 1)) - 1)) == p.residue,
                      "CROWN child cylinder refines its parent");
            }
        }
        // Independent exhaustive classification of every odd residue at the largest modulus.
        bool partition = true;
        for (uint32_t n = 1; n < (1u << (3 * depth + 1)); n += 2) {
            uint64_t orbit = n;
            bool allowed = true;
            int wordIndex = 0;
            for (int j = 0; j < depth; j++) {
                orbit = 3 * orbit + 1;
                int a = 0;
                while (!(orbit & 1)) {
                    orbit >>= 1;
                    a++;
                }
                if (a > 3) {
                    allowed = false;
                    break;
                }
                wordIndex = wordIndex * 3 + a - 1;
            }
            int hits = 0, index = -1;
            for (int j = 0; j < count; j++) {
                auto w = d->cylinders[j];
                if ((n & ((1u << (w.bits + 1)) - 1)) == w.residue) {
                    hits++;
                    index = j;
                }
            }
            partition = partition && hits == int(allowed) && (!allowed || index == wordIndex);
        }
        check(partition, "CROWN cylinders equal exhaustive odd-residue classification without overlaps");
        for (amp::C q : {amp::C(.45, .35), amp::C(.30, .55), amp::C(-.48, .34), amp::C(0, 0)}) {
            auto roots = crown::roots(d->cylinders, q);
            for (size_t j = 0; j < roots.size(); j++) {
                Z expected{}, power{1, 0}, qq{q.real(), q.imag()};
                auto w = d->cylinders[j];
                for (uint32_t bit = 1; bit <= w.bits; bit++) {
                    expected += static_cast<long double>(2 * int((w.residue >> bit) & 1u) - 1) * power;
                    power *= qq;
                }
                check(std::abs(expected - Z(roots[j].real(), roots[j].imag())) < 1e-13L,
                      "CROWN Horner embedding equals independent digit sum");
            }
            for (int i = 0; i < 24; i++) {
                Z z{1.73L * std::sin(.83L * i) + .09L, 1.41L * std::cos(.47L * i)}, relaxation{.86L, .16L};
                Jet f{{1, 0}};
                for (auto root : roots)
                    f = mul(f, Jet{z - Z(root.real(), root.imag()), {1, 0}, {0, 0}});
                for (int method = 0; method < 2; method++) {
                    Z expected =
                        relaxation * (method ? 2.L * f.f * f.d / (2.L * f.d * f.d - f.f * f.dd) : f.f / f.d);
                    auto step =
                        crown::correction({double(z.real()), double(z.imag())}, roots, {.86, .16}, method);
                    check(step.valid && std::abs(Z(step.delta.real(), step.delta.imag()) - expected) <
                                            2e-9L * (1 + std::abs(expected)),
                          "CROWN logarithmic step equals independent product derivatives, including repeated "
                          "roots");
                }
            }
            auto exact = crown::correction(roots[0], roots, {1, 0}, 0);
            check(exact.valid && exact.error == 0 && exact.delta == amp::C(0, 0),
                  "Exact cylinder root remains fixed");
        }
    }
    const std::vector<std::vector<int>> expectedDelays = {
        {0, 4, 3, 2, 1, 1, 2, 2, 1},
        {0, 7, 5, 4, 3, 3, 4, 3, 2, 1},
        {0, 12, 7, 6, 6, 6, 8, 6, 5, 3, 2, 2, 1},
        {0, 23, 12, 11, 10, 6, 8, 8, 9, 6, 6, 7, 5, 4, 5, 4, 3, 1}};
    for (int generation = 1; generation <= 4; generation++) {
        auto d = garden::data(4, generation);
        auto &t = d->transfer;
        check(t.counts == expectedDelays[generation - 1], "Exact restart delay histogram");
        int total = 0;
        for (int h : t.counts)
            total += h;
        check(total == t.states && t.counts.back() > 0,
              "Restart delays partition all states and attain degree");
        for (int j = 0; j < t.states; j++) {
            int at = j, steps = 0;
            while (at >= 0 && steps <= t.states) {
                at = t.next[at];
                steps++;
            }
            check(at == -1 && steps == t.delay[j],
                  "Every finite transfer column exits within the nilpotency bound");
        }
        for (uint64_t prime : {1000000007ull, 998244353ull})
            for (auto parameters : {std::pair<int, int>{0, 1}, {1, 1}, {2, 3}, {3, 7}, {5, 0}}) {
                uint64_t z = parameters.first, q = parameters.second, sum = 0;
                for (int l = 1; l <= t.degree; l++)
                    sum = (sum + uint64_t(t.counts[l]) * modularPower(z, t.degree - l, prime)) % prime;
                uint64_t polynomial =
                    (uint64_t(t.states) * modularPower(z, t.degree, prime) + prime - q * sum % prime) % prime;
                uint64_t expected = modularPower(t.states, t.states - 1, prime) *
                                    modularPower(z, t.states - t.degree, prime) % prime * polynomial % prime;
                check(transferDeterminant(t, z, q, prime) == expected,
                      "Full transfer determinant equals renewal spectral factor modulo a prime");
            }
        check(std::abs(garden::evaluate({1, 0}, {1, 0}, {1, 0}, *d).f) < 1e-15,
              "Markov feedback q=1 has eigenvalue 1");
        for (int i = 0; i < 48; i++) {
            Z z{i == 0 ? 0.L : 1.43L * std::sin(.73L * i), i == 0 ? 0.L : 1.21L * std::cos(.39L * i)},
                q{.85L, .35L}, relaxation{.86L, .18L};
            Jet variable{z, {1, 0}, {0, 0}}, f = pow(variable, t.degree);
            for (int l = 1; l <= t.degree; l++)
                f = add(f, scale(pow(variable, t.degree - l),
                                 -q * (static_cast<long double>(t.counts[l]) / t.states)));
            for (int method = 0; method < 2; method++) {
                auto step = garden::correction({double(z.real()), double(z.imag())}, {.85, .35}, {.86, .18},
                                               *d, method);
                Z expected =
                    relaxation * (method ? 2.L * f.f * f.d / (2.L * f.d * f.d - f.f * f.dd) : f.f / f.d);
                check(step.valid && std::abs(Z(step.delta.real(), step.delta.imag()) - expected) <
                                        2e-9L * (1 + std::abs(expected)),
                      "Restart scaled Horner corrections equal independent polynomial jets");
            }
        }
        for (int method = 0; method < 2; method++) {
            auto step = garden::correction({2, .3}, {0, 0}, {.86, .18}, *d, method);
            amp::C expected =
                amp::C(.86, .18) * amp::C(2, .3) * (method ? 2. / (t.degree + 1) : 1. / t.degree);
            check(step.valid && std::abs(step.delta - expected) < 1e-14,
                  "Zero feedback uses the exact monomial Newton/Halley step");
            auto o = garden::orbit({2, .3}, {0, 0}, {1, 0}, *d, method, 512);
            check(o.status == 0 && o.z == amp::C(0, 0),
                  "Zero-feedback nilpotent spectrum converges to its sole root");
        }
    }
    out << "PASS: " << checks << " garden checks, including exact transfer determinant comparisons\n";
    return checks;
}
} // namespace garden_tests
