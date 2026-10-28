#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

namespace crown {
using C = std::complex<double>;
struct Cylinder {
    std::array<int, 4> word{};
    uint32_t residue = 0, bits = 0;
};
inline std::vector<Cylinder> cylinders(int depth) {
    if (depth < 1 || depth > 4)
        throw std::runtime_error("CROWN depth must be 1..4");
    std::vector<Cylinder> result;
    std::array<int, 4> word{};
    std::function<void(int, uint32_t, uint32_t, uint32_t)> visit = [&](int n, uint32_t sum, uint32_t offset,
                                                                       uint32_t multiplier) {
        if (n == depth) {
            // S^N(n) = (3^N*n+B)/2^A must be odd. Invert 3^N modulo 2^(A+1).
            uint32_t inverse = 1;
            for (int j = 0; j < 5; j++)
                inverse *= 2u - multiplier * inverse;
            uint32_t mask = (1u << (sum + 1)) - 1;
            result.push_back({word, (((1u << sum) - offset) * inverse) & mask, sum});
            return;
        }
        for (int a = 1; a <= 3; a++) {
            word[n] = a;
            visit(n + 1, sum + a, 3 * offset + (1u << sum), 3 * multiplier);
        }
    };
    visit(0, 0, 0, 1);
    return result;
}
inline std::vector<C> roots(const std::vector<Cylinder> &words, C q) {
    std::vector<C> result;
    result.reserve(words.size());
    for (const auto &w : words) {
        C root{};
        // Discard the shared odd bit; embed the remaining digits least-significant
        // first.
        for (int j = int(w.bits); j >= 1; j--)
            root = root * q + double((w.residue >> j) & 1u) * 2. - 1.;
        result.push_back(root);
    }
    return result;
}
struct Correction {
    C delta{};
    double error = 1, label = 0;
    bool valid = false;
};
inline Correction correction(C z, const std::vector<C> &roots, C relaxation, int method) {
    Correction s;
    C first{}, second{};
    double nearest = 1e300;
    for (size_t i = 0; i < roots.size(); i++) {
        C diff = z - roots[i];
        double distance = std::abs(diff) / (1 + std::abs(roots[i]));
        if (!std::isfinite(distance))
            return s;
        if (distance < nearest) {
            nearest = distance;
            s.label = (double(i) + .5) / roots.size();
        }
    }
    s.error = nearest;
    if (nearest < 3e-6) {
        s.valid = true;
        return s;
    }
    for (C root : roots) {
        C inverse = 1. / (z - root);
        first += inverse;
        second += inverse * inverse;
    }
    // Logarithmic derivatives of the product, retaining any repeated roots.
    C denominator = method ? first * first + second : first;
    if (std::norm(denominator) < 1e-26)
        return s;
    s.delta = relaxation * (method ? 2. * first : C(1, 0)) / denominator;
    s.valid = std::isfinite(s.delta.real()) && std::isfinite(s.delta.imag());
    return s;
}
} // namespace crown
