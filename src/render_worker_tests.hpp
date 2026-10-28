#pragma once
#include "render_workers.hpp"
#include <array>
#include <atomic>
#include <ostream>
#include <stdexcept>

namespace render_worker_tests {
inline bool run(std::ostream& out) {
    renderworkers::Pool pool;
    bool good = true;
    int cases = 0;
    for (unsigned requested : {1u, 4u, 2u, 6u, 1u, 3u}) {
        std::array<std::atomic<int>, 257> visits{};
        for (auto& visit : visits) visit.store(0);
        std::atomic<unsigned> next{0};
        pool.run(requested, [&] {
            for (unsigned index; (index = next.fetch_add(1)) < visits.size();)
                visits[index].fetch_add(1);
        });
        for (auto& visit : visits) good = good && visit.load() == 1;
        ++cases;
    }
    // A failing task must still join its helpers and permit a later frame.
    bool caught = false;
    try { pool.run(3, [] { throw std::runtime_error("worker fixture"); }); }
    catch (const std::runtime_error&) { caught = true; }
    good = good && caught;
    ++cases;
    std::atomic<int> completed{0};
    pool.run(2, [&] { completed.fetch_add(1); });
    good = good && completed.load() >= 1 && completed.load() <= 2;
    ++cases;
    out << "RENDER_WORKER_CASES=" << cases << '\n'
        << "RENDER_WORKERS=" << (good ? "PASS" : "FAIL") << '\n';
    return good;
}
} // namespace render_worker_tests
