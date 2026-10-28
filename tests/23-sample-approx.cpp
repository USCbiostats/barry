#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

// Approximate (MCMC) sampling with Model::sample_approx() (tie/no-tie
// Metropolis-Hastings). The chains are compared against the exact expected
// statistics, E[s(Y); theta] = s(y_obs) - grad, where grad is the exact
// gradient of the log-likelihood of a single observed array.

#ifdef BARRY_VALGRIND
    #define SAMPLE_APPROX_N    500u
    #define SAMPLE_APPROX_TOL  0.35
#else
    #define SAMPLE_APPROX_N    20000u
    #define SAMPLE_APPROX_TOL  0.1
#endif

// Column means of a row-major (n x k) matrix
inline std::vector< double > sample_approx_colmeans(
    const std::vector< double > & x,
    size_t k
) {

    std::vector< double > ans(k, 0.0);
    size_t n = x.size() / k;
    for (size_t i = 0u; i < n; ++i)
        for (size_t j = 0u; j < k; ++j)
            ans[j] += x[i * k + j] / static_cast< double >(n);

    return ans;

}

// Largest absolute difference between the chain means and the exact
// expectation, relative to max(1, |E[s]|)
template< typename Model_Type >
inline double sample_approx_maxdiff(
    Model_Type & model,
    const std::vector< double > & par,
    unsigned int seed
) {

    auto expected = model.likelihood_total_grad(par);
    auto & target = model.get_stats_target()->at(0u);
    for (size_t k = 0u; k < expected.size(); ++k)
        expected[k] = target[k] - expected[k];

    model.set_seed(seed);
    auto state = model.sample_approx_init(model.get_pset(0u)->at(0u));

    // Burn-in, then the sample
    std::vector< double > stats;
    model.sample_approx(state, par, 1000u, 1u, nullptr);
    model.sample_approx(state, par, SAMPLE_APPROX_N, 10u, &stats);

    auto means = sample_approx_colmeans(stats, par.size());

    double ans = 0.0;
    for (size_t k = 0u; k < means.size(); ++k)
    {
        double d = std::fabs(means[k] - expected[k]) /
            std::max(1.0, std::fabs(expected[k]));

        if (d > ans)
            ans = d;
    }

    return ans;

}

BARRY_TEST_CASE("Approximate sampling (TNT)", "[sample approx]") {

    using namespace barry::counters::network;

    // ERGM vs exact -----------------------------------------------------------
    NetModel<Network> model;
    model.store_psets();
    counter_edges<>(model.get_counters());
    counter_mutual<>(model.get_counters());
    counter_ttriads<>(model.get_counters());
    rules_zerodiag<>(model.get_rules());

    Network net(4, 4);
    net.set_data(new NetworkData({0, 0, 1, 0}), true);
    net(0, 1) = 1;
    net(1, 0) = 1;
    net(2, 3) = 1;
    model.add_array(net);

    // Sparse and dense regimes (the latter exercises the proposal ratio)
    std::vector< double > par_sparse = {-2.0, 1.0, 0.2};
    std::vector< double > par_dense  = {2.0, -0.5, 0.1};

    double err_sparse = sample_approx_maxdiff(model, par_sparse, 123u);
    double err_dense  = sample_approx_maxdiff(model, par_dense, 321u);

    // Same with a dense array
    NetModel<NetworkDense> model_dense;
    model_dense.store_psets();
    counter_edges<NetworkDense>(model_dense.get_counters());
    counter_mutual<NetworkDense>(model_dense.get_counters());
    counter_ttriads<NetworkDense>(model_dense.get_counters());
    rules_zerodiag<NetworkDense>(model_dense.get_rules());

    NetworkDense net_dense(4, 4);
    net_dense.set_data(new NetworkData({0, 0, 1, 0}), true);
    net_dense(0, 1) = 1;
    net_dense(1, 0) = 1;
    net_dense(2, 3) = 1;
    model_dense.add_array(net_dense);

    double err_dense_arr = sample_approx_maxdiff(model_dense, par_sparse, 55u);

    // Nonlinear transformation (edges x mutual interaction) -------------------
    NetModel<Network> model_t;
    model_t.store_psets();
    counter_edges<>(model_t.get_counters());
    counter_mutual<>(model_t.get_counters());
    rules_zerodiag<>(model_t.get_rules());
    model_t.add_array(net);

    std::function<std::vector<double>(double *, size_t)> tfun = []
    (double * dat, size_t) {
        return std::vector< double >({dat[0u], dat[1u], dat[0u] * dat[1u] / 10.0});
    };
    model_t.set_transform_model(tfun, {"edges", "mutual", "edges x mutual"});

    std::vector< double > par_t = {-1.0, 0.5, 0.8};
    double err_t = sample_approx_maxdiff(model_t, par_t, 99u);

    printf_barry(
        "sample_approx: max rel. error sparse: %.3f, dense: %.3f, "
        "dense array: %.3f, transformed: %.3f\n",
        err_sparse, err_dense, err_dense_arr, err_t
    );

    // Persistence and incremental statistics ----------------------------------
    model.set_seed(10u);
    auto state = model.sample_approx_init(net);
    model.sample_approx(state, par_sparse, 300u);
    model.sample_approx(state, par_dense, 300u);
    model.sample_approx(state, {0.0, 0.0, 0.0}, 300u);

    barry::StatsCounter<Network, NetCounterData> recount(&state.array);
    recount.set_counters(model.get_counters());
    auto stats_recount = recount.count_all();

    bool ties_ok = true;
    size_t nties = 0u;
    for (size_t k = 0u; k < state.free.size() / 2u; ++k)
    {
        size_t i = state.free[k * 2u];
        size_t j = state.free[k * 2u + 1u];
        bool on = static_cast< double >(
            static_cast< const Network & >(state.array)(i, j)
        ) != 0.0;

        if (on != (state.pos[k] >= 0))
            ties_ok = false;

        if (on)
            ++nties;
    }

    bool diag_ok = true;
    for (size_t i = 0u; i < 4u; ++i)
        if (static_cast< double >(static_cast< const Network & >(state.array)(i, i)) != 0.0)
            diag_ok = false;

    // Two calls == one call (same seed)
    std::vector< double > stream_a, stream_b;
    model.set_seed(77u);
    auto state_a = model.sample_approx_init(net);
    model.sample_approx(state_a, par_sparse, 50u, 3u, &stream_a);
    model.sample_approx(state_a, par_sparse, 50u, 3u, &stream_a);

    model.set_seed(77u);
    auto state_b = model.sample_approx_init(net);
    model.sample_approx(state_b, par_sparse, 100u, 3u, &stream_b);

    // Single draw
    std::vector< double > stream_single;
    Network draw = model.sample_approx(net, par_sparse, 100u, &stream_single);

    // Locked cells keep their starting values, as in the exact support ------
    NetModel<Network> model_lock;
    model_lock.store_psets();
    counter_edges<>(model_lock.get_counters());
    rules_zerodiag<>(model_lock.get_rules());

    Network net_lock(3, 3);
    net_lock.set_data(new NetworkData({0, 0, 0}), true);
    net_lock(0, 0) = 1;
    net_lock(0, 1) = 1;
    model_lock.add_array(net_lock);
    model_lock.set_seed(8u);

    bool lock_exact_ok  = true;
    bool lock_approx_ok = true;
    auto state_lock = model_lock.sample_approx_init(net_lock);
    for (size_t s = 0u; s < 50u; ++s)
    {
        Network x = model_lock.sample(0u, {-1.0});
        if (static_cast< double >(static_cast< const Network & >(x)(0, 0)) == 0.0)
            lock_exact_ok = false;

        model_lock.sample_approx(state_lock, {-1.0}, 10u);
        if (static_cast< double >(
            static_cast< const Network & >(state_lock.array)(0, 0)
        ) == 0.0)
            lock_approx_ok = false;
    }

    // Dynamic rules: at most 3 edges -------------------------------------------
    NetModel<Network> model_dyn;
    counter_edges<>(model_dyn.get_counters());
    counter_mutual<>(model_dyn.get_counters());
    rules_zerodiag<>(model_dyn.get_rules());

    std::vector< double > * curstats =
        model_dyn.get_support_fun()->get_current_stats();
    model_dyn.add_rule_dyn(
        [curstats](const Network &, size_t, size_t, bool &) {
            return curstats->at(0u) <= 3.0;
        },
        false
    );

    model_dyn.set_seed(5u);
    auto state_dyn = model_dyn.sample_approx_init(net);
    std::vector< double > stream_dyn;
    model_dyn.sample_approx(state_dyn, {3.0, 0.0}, 2000u, 1u, &stream_dyn);

    double max_edges = 0.0;
    for (size_t s = 0u; s < stream_dyn.size(); s += 2u)
        max_edges = std::max(max_edges, stream_dyn[s]);

    Network net_bad(4, 4);
    net_bad.set_data(new NetworkData({0, 0, 1, 0}), true);
    for (size_t i = 0u; i < 4u; ++i)
        for (size_t j = 0u; j < 4u; ++j)
            if (i != j)
                net_bad(i, j) = 1;

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(err_sparse < SAMPLE_APPROX_TOL);
    REQUIRE(err_dense < SAMPLE_APPROX_TOL);
    REQUIRE(err_dense_arr < SAMPLE_APPROX_TOL);
    REQUIRE(err_t < SAMPLE_APPROX_TOL);

    REQUIRE(stats_recount == state.stats);
    REQUIRE(ties_ok);
    REQUIRE(nties == state.ties.size());
    REQUIRE(diag_ok);
    REQUIRE(state.nsteps == 900u);
    REQUIRE(state.naccepted > 0u);

    REQUIRE(stream_a.size() == 100u * 3u);
    REQUIRE(stream_a == stream_b);
    REQUIRE(state_a.stats == state_b.stats);

    REQUIRE(stream_single.size() == 100u * 3u);
    REQUIRE(draw.nrow() == 4u);

    REQUIRE(lock_exact_ok);
    REQUIRE(lock_approx_ok);
    REQUIRE(state_lock.free.size() / 2u == 6u);

    REQUIRE(max_edges <= 3.0);
    REQUIRE(max_edges == 3.0);
    REQUIRE_THROWS_AS(model_dyn.sample_approx_init(net_bad), std::logic_error);

    REQUIRE_THROWS_AS(
        model.sample_approx(state, {1.0, 2.0}, 1u),
        std::length_error
    );
    #endif

}

BARRY_TEST_CASE("Approximate sampling with DEFM", "[sample approx DEFM]") {

    using namespace defm;

    // Independent cells: Pr(y = 1) = logistic(par), so the column means of
    // the simulated rows should be close to logistic(par).
    const size_t nind = 20u;
    const size_t per  = 50u;
    const size_t n    = nind * per;
    const size_t ny   = 2u;
    const size_t nx   = 1u;

    std::vector< int >    id(n);
    std::vector< int >    y(n * ny, 0);
    std::vector< double > x(n * nx, 0.0);
    for (size_t r = 0u; r < n; ++r)
        id[r] = static_cast< int >(r / per) + 1;

    DEFM model(&id[0], &y[0], &x[0], n, ny, nx, 1u, true, true);
    counter_ones(model.get_counters());
    model.init(false);
    model.set_seed(2024);

    std::vector< double > par = {-1.0};
    std::vector< int > out(n * ny, 0);
    model.simulate_approx(par, &out[0], 200u);

    double mean = 0.0;
    size_t cnt  = 0u;
    for (size_t r = 0u; r < n; ++r)
    {
        // The first row of each individual is not simulated
        if ((r % per) == 0u)
            continue;

        for (size_t c = 0u; c < ny; ++c)
            mean += static_cast< double >(out[r * ny + c]);

        cnt += ny;
    }
    mean /= static_cast< double >(cnt);

    double expected = 1.0 / (1.0 + std::exp(1.0));
    printf_barry(
        "DEFM simulate_approx: mean %.4f, expected %.4f\n", mean, expected
    );

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::fabs(mean - expected) < 0.04);
    #endif

}
