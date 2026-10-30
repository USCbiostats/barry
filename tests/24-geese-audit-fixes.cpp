#include "tests.hpp"
#include "../include/barry/models/geese.hpp"

// Regression tests for the GEESE correctness audit
// (plans/geese-correctness-audit.md). Finding numbers refer to that document.

inline void geese_audit_init(geese::Geese & g)
{
    geese::counter_gains(g.get_counters(), {0});
    geese::counter_loss(g.get_counters(), {0});
    g.init(0);
}

BARRY_TEST_CASE("GEESE audit fixes", "[geese-audit]")
{

    using namespace geese;

    // Copies are not allowed (findings 4 and 11)
    static_assert(
        !std::is_copy_constructible< Geese >::value,
        "Geese should not be copy constructible."
    );

    std::vector< size_t > id     = {0, 1, 2, 3, 4};
    std::vector< int >    parent = {-1, 0, 0, 1, 1};
    std::vector< bool >   dupl(5, true);
    std::vector< double > par    = {-2, -1, 0.3};

    using Ann = std::vector< std::vector< size_t > >;
    auto fresh = [&](Ann a) {
        Geese f(a, id, parent, dupl);
        geese_audit_init(f);
        return f.likelihood(par);
    };

    // Findings 1, 2, and 10: updates match fresh models -----------------------
    Ann ann = Ann(5, {9});
    Geese g(ann, id, parent, dupl);
    geese_audit_init(g);

    REQUIRE(g.likelihood(par) == Approx(1.0));
    REQUIRE(g.likelihood(par, true) == Approx(0.0).margin(1e-12));

    // only_annotated prediction on a tree without annotations
    REQUIRE_NOTHROW(g.predict(par, nullptr, false, true));

    size_t narrays = g.get_model()->get_arrays2support()->size();

    // Value computed before the fixes (plans/geese-audit/results.txt)
    Ann ann_tip3 = ann;
    ann_tip3[3] = {1};
    REQUIRE(fresh(ann_tip3) == Approx(0.4071904237454682).epsilon(1e-12));

    // Missing -> observed
    g.update_annotations(2, {1});
    ann[2] = {1};
    REQUIRE(g.likelihood(par) == Approx(fresh(ann)).epsilon(1e-12));

    // Annotating a previously excluded branch
    g.update_annotations(4, {0});
    ann[4] = {0};
    REQUIRE(g.likelihood(par) == Approx(fresh(ann)).epsilon(1e-12));
    REQUIRE(g.likelihood(par) == Approx(g.likelihood(par, false, false)).epsilon(1e-12));

    // Observed -> missing (node 1 drops out of the reduced sequence)
    g.update_annotations(4, {9});
    ann[4] = {9};
    REQUIRE(g.likelihood(par) == Approx(fresh(ann)).epsilon(1e-12));

    g.update_annotations(2, {9});
    REQUIRE(g.likelihood(par, true) == Approx(0.0).margin(1e-12));

    // Leave-one-out matches freshly masked models, and does not grow the model
    Ann ann_obs = {{9}, {9}, {1}, {0}, {1}};
    Geese h(ann_obs, id, parent, dupl);
    geese_audit_init(h);

    size_t narrays_h = h.get_model()->get_arrays2support()->size();
    auto loo = h.predict(par, nullptr, true);
    REQUIRE(h.get_model()->get_arrays2support()->size() == narrays_h);

    // Values computed before the fixes (plans/geese-audit/results.txt)
    REQUIRE(h.likelihood(par) == Approx(0.07809990947097682).epsilon(1e-12));
    REQUIRE(loo[2][0] == Approx(0.5272873804128089).epsilon(1e-12));
    REQUIRE(loo[3][0] == Approx(0.6846375092305400).epsilon(1e-12));
    REQUIRE(loo[4][0] == Approx(0.3501804795632941).epsilon(1e-12));

    for (size_t tip : {2u, 3u, 4u})
    {
        Ann masked = ann_obs;
        masked[tip] = {9};
        Geese f(masked, id, parent, dupl);
        geese_audit_init(f);
        REQUIRE(loo[tip][0] == Approx(f.predict(par)[tip][0]).epsilon(1e-12));
    }

    g.update_annotations(3, {1});
    REQUIRE(g.get_model()->get_arrays2support()->size() == narrays);

    // Finding 9: internal annotations are ignored by all likelihoods
    Ann ann_internal = {{1}, {0}, {1}, {0}, {1}};
    Geese gi(ann_internal, id, parent, dupl);
    geese_audit_init(gi);
    double l_obs = h.likelihood(par);
    REQUIRE(gi.likelihood(par) == Approx(l_obs).epsilon(1e-12));
    REQUIRE(gi.likelihood_exhaust(par) == Approx(l_obs).epsilon(1e-12));

    h.update_annotations(1, {1});
    REQUIRE(h.likelihood(par) == Approx(l_obs).epsilon(1e-12));
    REQUIRE(h.likelihood_exhaust(par) == Approx(l_obs).epsilon(1e-12));
    h.update_annotations(1, {9});

    REQUIRE_THROWS(h.update_annotations(3, {2}));

    // Finding 5: large root logits ---------------------------------------------
    for (double sgn : {1.0, -1.0})
    {
        double l_large = h.likelihood({-2, -1, sgn * 1000});
        double l_limit = h.likelihood({-2, -1, sgn * 40});
        REQUIRE(std::isfinite(l_large));
        REQUIRE(l_large == Approx(l_limit).epsilon(1e-12));
        REQUIRE(
            h.likelihood_exhaust({-2, -1, sgn * 1000}) ==
            Approx(l_large).epsilon(1e-12)
        );
    }

    REQUIRE(h.simulate({-2, -1, 1000})[0][0] == 1u);

    // Finding 3: simulation after likelihood at other parameters ---------------
    std::vector< double > par_a = {-8, -8, 0.3};
    std::vector< double > par_b = {8, 8, 0.3};

    Geese s1(ann_obs, id, parent, dupl);
    geese_audit_init(s1);
    Geese s2(ann_obs, id, parent, dupl);
    geese_audit_init(s2);

    for (size_t i = 0u; i < 10u; ++i)
        s1.simulate(par_a);

    s1.likelihood(par_b);
    s1.likelihood_exhaust(par_b);

    s1.set_seed(123);
    s2.set_seed(123);

    std::vector< std::vector< std::vector< size_t > > > sims1, sims2;
    for (size_t i = 0u; i < 200u; ++i)
    {
        sims1.push_back(s1.simulate(par_a));
        sims2.push_back(s2.simulate(par_a));
    }

    REQUIRE(sims1 == sims2);

    // Finding 7: distinct parent states get distinct support keys --------------
    PhyloArray a1(17, 2), a2(17, 2);
    std::vector< bool > x(17, false);
    std::vector< bool > y = x;
    x[16] = true;
    y[16] = true;
    y[0]  = true;
    a1.set_data(new NodeData({1, 1}, x), true);
    a2.set_data(new NodeData({1, 1}, y), true);
    REQUIRE(keygen_full(a1, nullptr) != keygen_full(a2, nullptr));

    // Finding 8: exhaustive likelihood with sparse and permuted ids ------------
    Ann ann_sparse = {{9}, {1}, {0}};
    std::vector< size_t > id_sparse     = {10, 20, 30};
    std::vector< int >    parent_sparse = {-1, 10, 10};
    std::vector< bool >   dupl3(3, true);

    Geese e1(ann_sparse, id_sparse, parent_sparse, dupl3);
    geese_audit_init(e1);
    double l_sparse = e1.likelihood(par);
    REQUIRE(e1.likelihood_exhaust(par) == Approx(l_sparse).epsilon(1e-12));

    Ann ann_perm = {{0}, {9}, {1}};
    std::vector< size_t > id_perm     = {30, 10, 20};
    std::vector< int >    parent_perm = {10, -1, 10};

    Geese e2(ann_perm, id_perm, parent_perm, dupl3);
    geese_audit_init(e2);
    REQUIRE(e2.likelihood(par) == Approx(l_sparse).epsilon(1e-12));
    REQUIRE(e2.likelihood_exhaust(par) == Approx(l_sparse).epsilon(1e-12));

    // Finding 6: no underflow on large trees -----------------------------------
    // With all parameters at 0, each of the n observed tips is a fair coin.
    size_t n = 1100u;
    Ann ann_big(2 * n - 1, {9});
    std::vector< size_t > id_big(2 * n - 1);
    std::vector< int >    parent_big(2 * n - 1, -1);
    std::vector< bool >   dupl_big(2 * n - 1, true);
    for (size_t i = 0u; i < id_big.size(); ++i)
    {
        id_big[i] = i;
        if (i > 0u)
            parent_big[i] = static_cast< int >((i - 1) / 2);
        if (i >= (n - 1))
            ann_big[i] = {1};
    }

    Geese big(ann_big, id_big, parent_big, dupl_big);
    geese_audit_init(big);

    std::vector< double > par0 = {0, 0, 0};
    REQUIRE(
        big.likelihood(par0, true) ==
        Approx(-static_cast< double >(n) * std::log(2.0)).epsilon(1e-10)
    );
    REQUIRE(big.likelihood(par0, true, false) == Approx(big.likelihood(par0, true)));

    // Posteriors are still probabilities
    auto pred_big = big.predict(par0);
    bool all_probs = true;
    for (const auto & p : pred_big)
        all_probs = all_probs && (p[0] >= 0.0) && (p[0] <= 1.0);

    REQUIRE(all_probs);

}
