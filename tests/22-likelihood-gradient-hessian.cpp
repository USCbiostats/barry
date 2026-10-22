#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

// Largest relative difference between two vectors
inline double grad_test_maxreldiff(
    const std::vector< double > & a,
    const std::vector< double > & b
) {

    if (a.size() != b.size())
        throw std::length_error("-a- and -b- should be of the same length.");

    double ans = 0.0;
    for (size_t i = 0u; i < a.size(); ++i)
    {
        double d = std::fabs(a[i] - b[i]) / std::max(1.0, std::fabs(b[i]));
        if (d > ans)
            ans = d;
    }

    return ans;

}

// Central finite differences of the total log-likelihood
template< typename Model_Type >
inline std::vector< double > grad_test_fd_gradient(
    Model_Type & model,
    std::vector< double > par,
    double h = 1e-5
) {

    std::vector< double > ans(par.size());
    for (size_t j = 0u; j < par.size(); ++j)
    {

        double p0 = par[j];
        par[j] = p0 + h;
        double l_up = model.likelihood_total(par, true, 1u);
        par[j] = p0 - h;
        double l_down = model.likelihood_total(par, true, 1u);
        par[j] = p0;

        ans[j] = (l_up - l_down) / (2.0 * h);

    }

    return ans;

}

// Central finite differences of the (analytic) gradient, row-major
template< typename Model_Type >
inline std::vector< double > grad_test_fd_hessian(
    Model_Type & model,
    std::vector< double > par,
    double h = 1e-5
) {

    size_t K = par.size();
    std::vector< double > ans(K * K);
    for (size_t j1 = 0u; j1 < K; ++j1)
    {

        double p0 = par[j1];
        par[j1] = p0 + h;
        auto g_up = model.likelihood_total_grad(par);
        par[j1] = p0 - h;
        auto g_down = model.likelihood_total_grad(par);
        par[j1] = p0;

        for (size_t j0 = 0u; j0 < K; ++j0)
            ans[j0 * K + j1] = (g_up[j0] - g_down[j0]) / (2.0 * h);

    }

    return ans;

}

// Solves A x = b (A row-major, K x K) with Gaussian elimination
inline std::vector< double > grad_test_solve(
    std::vector< double > A,
    std::vector< double > b
) {

    size_t K = b.size();
    for (size_t c = 0u; c < K; ++c)
    {

        size_t piv = c;
        for (size_t r = c + 1u; r < K; ++r)
            if (std::fabs(A[r * K + c]) > std::fabs(A[piv * K + c]))
                piv = r;

        for (size_t j = 0u; j < K; ++j)
            std::swap(A[c * K + j], A[piv * K + j]);
        std::swap(b[c], b[piv]);

        for (size_t r = c + 1u; r < K; ++r)
        {
            double f = A[r * K + c] / A[c * K + c];
            for (size_t j = c; j < K; ++j)
                A[r * K + j] -= f * A[c * K + j];
            b[r] -= f * b[c];
        }

    }

    std::vector< double > x(K);
    for (size_t c = K; c-- > 0u;)
    {
        double s = b[c];
        for (size_t j = c + 1u; j < K; ++j)
            s -= A[c * K + j] * x[j];
        x[c] = s / A[c * K + c];
    }

    return x;

}

// Whether a symmetric matrix (row-major) is positive definite (Cholesky)
inline bool grad_test_is_posdef(const std::vector< double > & A, size_t K)
{

    std::vector< double > L(K * K, 0.0);
    for (size_t i = 0u; i < K; ++i)
        for (size_t j = 0u; j <= i; ++j)
        {

            double s = A[i * K + j];
            for (size_t k = 0u; k < j; ++k)
                s -= L[i * K + k] * L[j * K + k];

            if (i == j)
            {
                if (s <= 0.0)
                    return false;
                L[i * K + i] = std::sqrt(s);
            }
            else
                L[i * K + j] = s / L[j * K + j];

        }

    return true;

}

BARRY_TEST_CASE("Gradient and Hessian of the likelihood", "[likelihood derivatives]") {

    using namespace defm;

    // Simulated panel: 60 individuals, 4 visits, 3 outcomes, 2 covariates
    // (one binary, one continuous so most arrays get their own support).
    std::mt19937 rengine(1231);
    std::uniform_real_distribution< double > runif(0.0, 1.0);

    size_t n_id = 60u, n_visits = 4u, n_y = 3u, n_x = 2u;
    size_t n = n_id * n_visits;

    std::vector< int > id(n), Y(n * n_y);
    std::vector< double > X(n * n_x);
    for (size_t i = 0u; i < n_id; ++i)
    {

        double x0 = runif(rengine) < .4 ? 1.0 : 0.0;
        double x1 = runif(rengine) * 2.0 - 1.0;

        for (size_t v = 0u; v < n_visits; ++v)
        {

            size_t r = i * n_visits + v;
            id[r] = static_cast< int >(i);
            X[r]     = x0;
            X[r + n] = x1;
            for (size_t j = 0u; j < n_y; ++j)
                Y[r + j * n] = runif(rengine) < (.3 + .1 * j + .2 * x0) ? 1 : 0;

        }

    }

    // Model with covariate interactions, a motif, and a transition term
    DEFM model(&id[0u], &Y[0u], &X[0u], n, n_y, n_x, 1);
    model.set_names({"A", "B", "C"}, {"X1", "X2"});

    counter_logit_intercept(model.get_counters(), n_y);
    counter_logit_intercept(model.get_counters(), n_y, {}, 0, &model.get_X_names());
    counter_logit_intercept(model.get_counters(), n_y, {}, 1, &model.get_X_names());
    counter_formula(
        model.get_counters(), "{y0, y1}", 1, n_y,
        &model.get_X_names(), &model.get_Y_names()
    );
    counter_formula(
        model.get_counters(), "{y0_0} > {y2} x X2", 1, n_y,
        &model.get_X_names(), &model.get_Y_names()
    );

    model.init();

    auto & m = model.get_model();
    size_t K = m.nterms();

    std::vector< double > par(K);
    for (size_t j = 0u; j < K; ++j)
        par[j] = runif(rengine) - .5;

    // Log-likelihood from the fused function equals likelihood_total
    std::vector< double > grad;
    double ll_fused = m.likelihood_total_and_grad(par, grad);
    double ll_total = m.likelihood_total(par, true, 1u);

    auto grad_fd = grad_test_fd_gradient(m, par);
    auto hess    = m.likelihood_total_hess(par);
    auto hess_fd = grad_test_fd_hessian(m, par);

    // Multi-threaded versions
    std::vector< double > grad_2;
    double ll_fused_2 = m.likelihood_total_and_grad(par, grad_2, 2u);
    auto hess_2 = m.likelihood_total_hess(par, 2u);

    double err_grad = grad_test_maxreldiff(grad, grad_fd);
    double err_hess = grad_test_maxreldiff(hess, hess_fd);
    printf_barry("DEFM: K = %i, supports = %i\n", static_cast<int>(K), static_cast<int>(m.get_stats_support_sizes()->size()));
    printf_barry("DEFM: max rel. error gradient: %.2e, Hessian: %.2e\n", err_grad, err_hess);

    bool hess_symmetric = true;
    for (size_t j0 = 0u; j0 < K; ++j0)
        for (size_t j1 = 0u; j1 < K; ++j1)
            if (hess[j0 * K + j1] != hess[j1 * K + j0])
                hess_symmetric = false;

    std::vector< double > neg_hess(hess);
    for (auto & h : neg_hess)
        h = -h;

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(m.get_stats_support_sizes()->size() >= 128u);
    REQUIRE(ll_fused == Approx(ll_total).epsilon(1e-10));
    REQUIRE(err_grad < 1e-5);
    REQUIRE(err_hess < 1e-5);
    REQUIRE(hess_symmetric);
    REQUIRE(grad_test_is_posdef(neg_hess, K));
    REQUIRE(ll_fused_2 == Approx(ll_fused).epsilon(1e-12));
    REQUIRE(grad_test_maxreldiff(grad_2, grad) < 1e-12);
    REQUIRE(grad_test_maxreldiff(hess_2, hess) < 1e-12);
    #endif

    // Order-0 model with logit intercepts only: independent logistic
    // regressions, so we can compare with the closed-form score and Hessian.
    std::vector< int > id0(n);
    for (size_t r = 0u; r < n; ++r)
        id0[r] = static_cast< int >(r);

    DEFM model_logit(&id0[0u], &Y[0u], &X[0u], n, n_y, n_x, 0);
    counter_logit_intercept(model_logit.get_counters(), n_y);
    counter_logit_intercept(model_logit.get_counters(), n_y, {}, 0);
    counter_logit_intercept(model_logit.get_counters(), n_y, {}, 1);
    model_logit.init();

    auto & m0 = model_logit.get_model();
    size_t K0 = m0.nterms();

    // Parameter layout: [b_y (n_y), g_{y,x0} (n_y), g_{y,x1} (n_y)]
    auto logit_score_hess = [&](
        const std::vector< double > & p,
        std::vector< double > & score,
        std::vector< double > & H
    ) -> void {

        score.assign(K0, 0.0);
        H.assign(K0 * K0, 0.0);
        for (size_t r = 0u; r < n; ++r)
        {

            double design[3] = {1.0, X[r], X[r + n]};
            for (size_t j = 0u; j < n_y; ++j)
            {

                double eta = 0.0;
                for (size_t c = 0u; c < 3u; ++c)
                    eta += p[c * n_y + j] * design[c];

                double prob = 1.0 / (1.0 + std::exp(-eta));
                double resid = Y[r + j * n] - prob;
                for (size_t c0 = 0u; c0 < 3u; ++c0)
                {

                    score[c0 * n_y + j] += design[c0] * resid;
                    for (size_t c1 = 0u; c1 < 3u; ++c1)
                        H[(c0 * n_y + j) * K0 + (c1 * n_y + j)] -=
                            prob * (1.0 - prob) * design[c0] * design[c1];

                }

            }

        }

    };

    std::vector< double > par0(K0, 0.1), score_cf, hess_cf;
    logit_score_hess(par0, score_cf, hess_cf);

    auto grad0 = m0.likelihood_total_grad(par0);
    auto hess0 = m0.likelihood_total_hess(par0);

    double err_logit_grad = grad_test_maxreldiff(grad0, score_cf);
    double err_logit_hess = grad_test_maxreldiff(hess0, hess_cf);
    printf_barry("Logit: max rel. error gradient: %.2e, Hessian: %.2e\n", err_logit_grad, err_logit_hess);

    // Newton-Raphson from zero: the objective is concave, so it should
    // converge to a point where the gradient is (numerically) zero.
    std::vector< double > theta(K0, 0.0);
    double ll_start = m0.likelihood_total(theta, true, 1u);
    for (size_t iter = 0u; iter < 25u; ++iter)
    {

        auto g = m0.likelihood_total_grad(theta);
        auto H = m0.likelihood_total_hess(theta);
        for (auto & h : H)
            h = -h;

        auto step = grad_test_solve(H, g);
        for (size_t j = 0u; j < K0; ++j)
            theta[j] += step[j];

    }

    auto grad_mle = m0.likelihood_total_grad(theta);
    double grad_mle_max = 0.0;
    for (auto g : grad_mle)
        grad_mle_max = std::max(grad_mle_max, std::fabs(g));

    double ll_mle = m0.likelihood_total(theta, true, 1u);
    printf_barry("Logit: Newton |grad|_max at MLE: %.2e\n", grad_mle_max);

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(err_logit_grad < 1e-10);
    REQUIRE(err_logit_hess < 1e-10);
    REQUIRE(grad_mle_max < 1e-8);
    REQUIRE(ll_mle > ll_start);
    #endif

    // ERGM (small networks), before and after transforming the model
    using namespace barry::counters::network;

    NetModel<NetworkDense> model_net;
    counter_edges<NetworkDense>(model_net.get_counters());
    counter_mutual<NetworkDense>(model_net.get_counters());
    counter_ttriads<NetworkDense>(model_net.get_counters());
    counter_nodematch<NetworkDense>(model_net.get_counters(), 0u);
    rules_zerodiag<NetworkDense>(model_net.get_rules());

    for (size_t net_i = 0u; net_i < 3u; ++net_i)
    {

        NetworkDense net(4, 4);
        net.set_data(
            new NetworkData({0.0, 0.0, 1.0, static_cast< double >(net_i % 2)}),
            true
        );

        for (size_t i = 0u; i < 4u; ++i)
            for (size_t j = 0u; j < 4u; ++j)
                if ((i != j) && (runif(rengine) < .4))
                    net(i, j) = 1;

        model_net.add_array(net);

    }

    std::vector< double > par_net = {-1.0, .5, .2, .3};
    auto err_net_grad = grad_test_maxreldiff(
        model_net.likelihood_total_grad(par_net),
        grad_test_fd_gradient(model_net, par_net)
    );
    auto err_net_hess = grad_test_maxreldiff(
        model_net.likelihood_total_hess(par_net),
        grad_test_fd_hessian(model_net, par_net)
    );

    // Transformation: drop mutual, rescale edges
    std::function<std::vector<double>(double *, size_t)> tfun = []
    (double * dat, size_t) {
        return std::vector< double >({dat[0u] / 2.0, dat[2u], dat[3u]});
    };

    model_net.set_transform_model(tfun, {"edges/2", "ttriads", "nodematch"});

    std::vector< double > par_net_t = {-2.0, .2, .3};
    auto err_net_t_grad = grad_test_maxreldiff(
        model_net.likelihood_total_grad(par_net_t),
        grad_test_fd_gradient(model_net, par_net_t)
    );
    auto err_net_t_hess = grad_test_maxreldiff(
        model_net.likelihood_total_hess(par_net_t),
        grad_test_fd_hessian(model_net, par_net_t)
    );

    printf_barry(
        "ERGM: max rel. error gradient: %.2e, Hessian: %.2e (transformed: %.2e, %.2e)\n",
        err_net_grad, err_net_hess, err_net_t_grad, err_net_t_hess
    );

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(err_net_grad < 1e-5);
    REQUIRE(err_net_hess < 1e-5);
    REQUIRE(err_net_t_grad < 1e-5);
    REQUIRE(err_net_t_hess < 1e-5);
    REQUIRE_THROWS_AS(model_net.likelihood_total_grad(par_net), std::length_error);
    #endif

}
