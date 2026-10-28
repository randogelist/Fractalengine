#include "self_tests.hpp"
#include "fresh_classic_tests.hpp"
#include "lab_tests.hpp"
#include "focus_tests.hpp"
#include "cpu_optimisation_tests.hpp"
#include "render_worker_tests.hpp"
#include "stable_curvature_tests.hpp"
#include "rh_branch_extension_tests.hpp"
#include "new_structure_tests.hpp"
#include "frk_extra_tests.hpp"
#include "qqrh_hilbert_tests.hpp"
#include <iostream>

int main() {
    bool pass = true;
    pass = rhfresh::runSelfTests(std::cout) && pass;
    pass = fresh_classic_tests::run(std::cout) && pass;
    pass = lab_tests::run(std::cout) && pass;
    pass = focus_tests::run(std::cout) && pass;
    pass = cpu_optimisation_tests::run(std::cout) && pass;
    pass = render_worker_tests::run(std::cout) && pass;
    pass = stable_curvature_tests::run(std::cout) && pass;
    pass = rh_branch_extension_tests::run(std::cout) && pass;
    pass = new_structure_tests::run(std::cout) && pass;
    pass = frk_extra_tests::run(std::cout) && pass;
    pass = qqrh_hilbert_tests::run(std::cout) && pass;
    std::cout << "PORTABLE_CPU_SUITES=11 RESULT=" << (pass ? "PASS" : "FAIL") << '\n';
    return pass ? 0 : 1;
}
