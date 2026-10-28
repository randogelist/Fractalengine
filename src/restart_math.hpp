#pragma once
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace restart {
struct Transfer {
    int states = 0, cutoff = 0, degree = 0;
    std::vector<int> next, delay, counts;
};
inline uint32_t syracuse(uint32_t n) {
    uint32_t value = 3 * n + 1;
    while (!(value & 1u))
        value >>= 1;
    return value;
}
inline Transfer transfer(int states) {
    if (states < 8 || states > 128 || (states & (states - 1)))
        throw std::runtime_error("Restart state count must be 8, 16, 32, 64 or 128");
    Transfer t;
    t.states = states;
    t.cutoff = 2 * states - 1;
    t.next.resize(states, -1);
    t.delay.resize(states);
    for (int j = 0; j < states; j++) {
        uint32_t n = 2 * j + 1, next = syracuse(n);
        // n=1 and a jump outside the finite set are restart columns.
        if (n != 1 && next <= uint32_t(t.cutoff))
            t.next[j] = int((next - 1) / 2);
    }
    for (int j = 0; j < states; j++) {
        int at = j, delay = 1;
        while (t.next[at] >= 0) {
            at = t.next[at];
            if (++delay > states)
                throw std::runtime_error("Finite transfer has a cycle; renewal reduction is invalid");
        }
        t.delay[j] = delay;
        t.degree = std::max(t.degree, delay);
    }
    t.counts.resize(t.degree + 1);
    for (int delay : t.delay)
        t.counts[delay]++;
    return t;
}
} // namespace restart
