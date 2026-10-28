#pragma once
#include <algorithm>
#include <cmath>
#include <complex>

namespace cayley {
using C = std::complex<double>;
struct Jet {
    C f{}, first{}, second{};
};
inline Jet add(const Jet &a, const Jet &b) {
    return {a.f + b.f, a.first + b.first, a.second + b.second};
}
inline Jet multiply(const Jet &a, const Jet &b) {
    return {a.f * b.f, a.first * b.f + a.f * b.first,
            a.second * b.f + 2. * a.first * b.first + a.f * b.second};
}
inline Jet scale(const Jet &a, double s) {
    return {a.f * s, a.first * s, a.second * s};
}
struct Pair {
    Jet p, q;
};
inline Pair preimage(C z, int generation) {
    Pair v{{z, {1, 0}, {}}, {{1, 0}, {}, {}}};
    for (int k = 0; k < generation; k++) {
        // A common, locally constant scale for all six jet components preserves
        // the derivatives of the homogeneous polynomial, up to one common
        // factor. Do not differentiate this numerical normalization.
        const Jet sum = add(v.p, v.q);
        const Jet p = scale(multiply(v.p, v.q), 4.);
        const Jet q = multiply(sum, sum);
        double norm = std::max({std::abs(p.f), std::abs(p.first), std::abs(p.second),
                                std::abs(q.f), std::abs(q.first), std::abs(q.second)});
        if (!(norm > 0) || !std::isfinite(norm))
            return {p, q};
        v = {scale(p, 1. / norm), scale(q, 1. / norm)};
    }
    return v;
}
} // namespace cayley
