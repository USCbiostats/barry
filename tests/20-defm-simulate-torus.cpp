#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

// Regression test for the support-vs-array index confusion in
// Model::sample(const Array_Type &, params) and DEFM::simulate().
//
// The model describes a "torus" walk over ten mutually-exclusive states:
// the single active state y_k deterministically moves to y_{k+1}, wrapping
// around from y9 back to y0. With a strong reward on each transition and a
// strong penalty on extra ones, the walk is deterministic, so every
// simulated row must have exactly one active state that advances by one.
//
// Before the fix, the walk broke the first time a mid-simulation
// conditioning state was revisited (row 13 in the README example): the
// sampler resolved the reused support through arrays2support a second time
// and drew from the wrong (all-zeros) support, collapsing the process.

BARRY_TEST_CASE("DEFM simulate torus keeps moving", "[DEFM simulate torus]") {

    using namespace defm;

    const size_t n  = 20u;
    const size_t ny = 10u;
    const size_t nx = 1u;

    std::vector< int >    id(n, 1);
    std::vector< int >    y(n * ny, 0);
    std::vector< double > x(n * nx, 0.0);

    // Column-major storage; the walk starts with y0 active.
    y[0] = 1;

    DEFM model(&id[0], &y[0], &x[0], n, ny, nx, 1u, true, true);

    // Ten transition terms: y_i -> y_{i+1}, wrapping at the end.
    for (size_t i = 0u; i < (ny - 1u); ++i)
    {
        std::string f =
            "{y"  + std::to_string(i) + ", 0y" + std::to_string(i + 1u) +
            "} > {0y" + std::to_string(i) + ", y"  + std::to_string(i + 1u) + "}";
        counter_formula(model.get_counters(), f, 1u, ny);
    }
    counter_formula(model.get_counters(), "{0y0, y9} > {y0, 0y9}", 1u, ny);

    counter_ones(model.get_counters());

    model.init(false);
    model.set_seed(33);

    // Strong transition reward, strong ones penalty => deterministic torus
    // (independent of the RNG, so the assertions are portable).
    std::vector< double > par(ny, 200.0);
    par.push_back(-20.0);

    std::vector< int > out(n * ny, 0);
    out[0] = 1; // baseline row (simulate() only writes rows 2..n)

    model.simulate(par, &out[0]);

    // Each row must have exactly one active state, advancing by one.
    bool ok = true;
    for (size_t r = 0u; r < n; ++r)
    {
        size_t active     = ny;      // sentinel: "none"
        size_t n_active   = 0u;
        for (size_t c = 0u; c < ny; ++c)
            if (out[r * ny + c] == 1)
            {
                active = c;
                ++n_active;
            }

        size_t expected = r % ny;    // y0, y1, ..., y9, y0, ...
        if ((n_active != 1u) || (active != expected))
        {
            ok = false;
            std::cout << "Row " << (r + 1u) << ": expected single active state y"
                      << expected << ", got " << n_active << " active (first at "
                      << (active == ny ? -1 : static_cast<int>(active)) << ")\n";
        }
    }

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(ok);
    #else
    if (!ok)
        std::cout << "FAIL: torus walk did not advance as expected." << std::endl;
    else
        std::cout << "OK: torus walk advanced through all rows." << std::endl;
    #endif

}
