// Benchmark: analytic gradient/Hessian vs finite differences in a DEFM with
// many covariate interactions (5 outcomes, Markov order 1, 21 binary
// covariates, K = 122 terms). Data are simulated.
#include "../../include/barry/barry.hpp"
#include "../../include/barry/models/defm.hpp"
#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> clocktime;

#define TIME_START(a) clocktime a##_start = std::chrono::steady_clock::now();
#define TIME_END(a) std::chrono::duration< double, std::milli > a##_diff = \
    std::chrono::steady_clock::now() - a##_start;
#define ELAPSED(a) a##_diff.count()

int main(int argc, char * argv[]) {

    size_t ncores = (argc > 1) ? std::strtoul(argv[1], nullptr, 10) : 1u;

    using namespace defm;

    std::mt19937 rengine(9123);
    std::uniform_real_distribution< double > runif(0.0, 1.0);

    // Categorical baseline covariates (constant within individual), expanded
    // into dummies (reference level dropped)
    std::vector< std::pair< std::string, size_t > > cats = {
        {"sex", 2}, {"site", 3}, {"age", 4}, {"race", 7}, {"branch", 5},
        {"rank", 6}
    };

    std::vector< std::string > x_names;
    for (auto & c : cats)
        for (size_t l = 2u; l <= c.second; ++l)
            x_names.push_back(c.first + std::to_string(l));

    size_t n_id = 2000u, n_y = 5u, n_x = x_names.size();

    std::vector< int > id;
    std::vector< std::vector< double > > x_rows;
    for (size_t i = 0u; i < n_id; ++i)
    {

        // Level probabilities decreasing with the level
        std::vector< double > xi;
        for (auto & c : cats)
        {

            double total = c.second * (c.second + 1) / 2.0, u = runif(rengine);
            size_t lvl = 1u;
            double acc = c.second / total;
            while ((u > acc) && (lvl < c.second))
                acc += (c.second - lvl++) / total;

            for (size_t l = 2u; l <= c.second; ++l)
                xi.push_back(lvl == l ? 1.0 : 0.0);

        }

        size_t n_visits = 2u + static_cast< size_t >(runif(rengine) * 6.0);
        for (size_t v = 0u; v < n_visits; ++v)
        {
            id.push_back(static_cast< int >(i));
            x_rows.push_back(xi);
        }

    }

    size_t n = id.size();
    std::vector< double > X(n * n_x);
    std::vector< int > Y(n * n_y);
    for (size_t r = 0u; r < n; ++r)
    {

        for (size_t c = 0u; c < n_x; ++c)
            X[r + c * n] = x_rows[r][c];

        for (size_t j = 0u; j < n_y; ++j)
            Y[r + j * n] = runif(rengine) < .3 ? 1 : 0;

    }

    // Model: intercepts, intercepts x every covariate, symptom-treatment
    // motifs, and motifs x site.
    TIME_START(build)
    DEFM model(&id[0u], &Y[0u], &X[0u], n, n_y, n_x, 1);
    model.set_names({"y0", "y1", "y2", "y3", "y4"}, x_names);

    auto * counters = model.get_counters();
    counter_logit_intercept(counters, n_y);
    for (size_t c = 0u; c < n_x; ++c)
        counter_logit_intercept(
            counters, n_y, {}, static_cast< int >(c), &model.get_X_names(),
            &model.get_Y_names()
        );

    for (size_t k = 1u; k < n_y; ++k)
    {

        std::string motif = "{y0, y" + std::to_string(k) + "}";
        counter_formula(
            counters, motif, 1, n_y, &model.get_X_names(), &model.get_Y_names()
        );
        for (std::string s : {"site2", "site3"})
            counter_formula(
                counters, motif + " x " + s, 1, n_y, &model.get_X_names(),
                &model.get_Y_names()
            );

    }

    model.init();
    TIME_END(build)

    auto & m = model.get_model();
    size_t K = m.nterms();

    std::vector< double > par(K);
    for (auto & p : par)
        p = (runif(rengine) - .5) / 5.0;

    size_t reps = 20u;

    // Log-likelihood alone (perturbing the parameters so the cache is not
    // reused between calls)
    TIME_START(loglik)
    for (size_t r = 0u; r < reps; ++r)
    {
        par[0u] += 1e-4;
        m.likelihood_total(par, true, ncores);
    }
    TIME_END(loglik)

    std::vector< double > grad;
    TIME_START(grad)
    for (size_t r = 0u; r < reps; ++r)
    {
        par[0u] += 1e-4;
        m.likelihood_total_and_grad(par, grad, ncores);
    }
    TIME_END(grad)

    TIME_START(hess)
    auto H = m.likelihood_total_hess(par, ncores);
    TIME_END(hess)

    double t_ll   = ELAPSED(loglik) / reps;
    double t_grad = ELAPSED(grad) / reps;

    printf("Rows: %i, arrays: %i, unique supports: %i, K = %i, ncores = %i\n",
        static_cast<int>(n), static_cast<int>(m.size()),
        static_cast<int>(m.get_stats_support_sizes()->size()), static_cast<int>(K),
        static_cast<int>(ncores));
    printf("Model build (init)                : %10.1f ms\n", ELAPSED(build));
    printf("Log-likelihood                    : %10.2f ms\n", t_ll);
    printf("Log-likelihood + analytic gradient: %10.2f ms\n", t_grad);
    printf("Gradient by central differences   : %10.2f ms (2K evaluations)\n",
        2.0 * K * t_ll);
    printf("Analytic Hessian                  : %10.2f ms\n", ELAPSED(hess));
    printf("Hessian by finite differences     : %10.2f ms (4K^2 evaluations)\n",
        4.0 * K * K * t_ll);

    return 0;

}
