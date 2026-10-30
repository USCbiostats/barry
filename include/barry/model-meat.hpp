#ifndef BARRY_MODEL_MEAT_HPP 
#define BARRY_MODEL_MEAT_HPP 1

/**
 * @defgroup stat-models Statistical Models
 * @brief Statistical models available in `barry`.
 */

inline double update_normalizing_constant(
    const std::vector<double> & params,
    const double * support,
    size_t k,
    size_t n
)
{
    double res = 0.0;

    std::vector< double > resv(n, 0.0);

    for (size_t j = 0u; j < (k - 1u); ++j)
    {

        const double p = params[j];
        
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd 
        #elif defined(__GNUC__) && !defined(__clang__)
            #pragma GCC ivdep
        #endif
        for (size_t i = 0u; i < n; ++i)
            resv[i] += (*(support + i * k + 1u + j)) * p;

    }

    // Accumulate resv to a double res        
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd reduction(+:res)
    #elif defined(__GNUC__) && !defined(__clang__)
        #pragma GCC ivdep
    #endif
    for (size_t i = 0u; i < n; ++i)
    {
        res += std::exp(resv[i] BARRY_SAFE_EXP) * (*(support + i * k));
    }




    #ifdef BARRY_DEBUG
    if (std::isnan(res))
        throw std::overflow_error(
            std::string("NaN in update_normalizing_constant. ") +
            std::string("res = ") + std::to_string(res) +
            std::string(", k = ") + std::to_string(k) +
            std::string(", n = ") + std::to_string(n)
            );
    if (std::isinf(res))
        throw std::overflow_error(
            std::string("Inf in update_normalizing_constant. ") +
            std::string("res = ") + std::to_string(res) +
            std::string(", k = ") + std::to_string(k) +
            std::string(", n = ") + std::to_string(n)
            );

    #endif

    return res;
    
}

inline double likelihood_(
        const double * stats_target,
        const std::vector< double > & params,
        const double normalizing_constant,
        size_t n_params,
        bool log_ = false
) {
    
    if (n_params != params.size())
        throw std::length_error("-stats_target- and -params- should have the same length.");
        
    double numerator = 0.0;
    
    // Computing the numerator
    #ifdef __INTEL_LLVM_COMPILER
    #pragma code_align 32
    #endif
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd reduction(+:numerator)
    #endif
    for (size_t j = 0u; j < n_params; ++j)
        numerator += *(stats_target + j) * params[j];

    if (!log_)
        numerator = std::exp(numerator BARRY_SAFE_EXP);
    else
        return numerator BARRY_SAFE_EXP - std::log(normalizing_constant);

    double ans = numerator/normalizing_constant;

    #ifdef BARRY_DEBUG
    if (std::isnan(ans))
        throw std::overflow_error(
            std::string("NaN in likelihood_. ") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );
    if (std::isinf(ans))
        throw std::overflow_error(
            std::string("Inf in likelihood_. ") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );

    if (ans > 1.0)
        throw std::overflow_error(
            std::string("Likelihood > 1.0") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );
    #endif

    return ans;

}

/**
 * @brief Probabilities and expected statistics of a single support
 * @details Computes, in a numerically stable way (shifting the exponent by its
 * maximum), the probability of each row of the support and the expected value
 * of the sufficient statistics:
 * \f[
 * p_i(\theta) = \mathrm{Pr}_{\mathcal{Y}}(Y = y_i; \theta) = \frac{w_i\exp{\left(\theta^{\mathbf{t}}s(y_i)\right)}}{\eta},\quad
 * \mu = \sum_i p_i s(y_i)
 * \f]
 * @param params Vector of model parameters (length `k - 1`).
 * @param support Pointer to the first element of the support. Each row has
 * `k` elements: the weight followed by the `k - 1` statistics.
 * @param k Number of elements per row (number of terms + 1).
 * @param n Number of rows in the support.
 * @param probs Output vector, resized to `n` and filled with \f$p_i\f$.
 * @param mu Output pointer to `k - 1` elements, filled with \f$\mu\f$.
 * @return The log of the normalizing constant, \f$\log \eta\f$, or `-Inf` if the
 * support is empty.
 * @see `doc/likelihood-derivatives.md` for the derivation.
 */
inline double update_expected_stats(
    const std::vector< double > & params,
    const double * support,
    size_t k,
    size_t n,
    std::vector< double > & probs,
    double * mu
)
{

    const size_t nterms = k - 1u;
    std::fill(mu, mu + nterms, 0.0);
    probs.resize(n);

    // Linear predictor of each row, and its maximum (among rows with weight)
    double eta_max = -std::numeric_limits< double >::infinity();
    for (size_t i = 0u; i < n; ++i)
    {

        const double * row = support + i * k + 1u;
        double eta = 0.0;
        for (size_t j = 0u; j < nterms; ++j)
            eta += row[j] * params[j];

        probs[i] = eta;
        if ((*(support + i * k) > 0.0) && (eta > eta_max))
            eta_max = eta;

    }

    if (std::isinf(eta_max) && (eta_max < 0.0))
        return -std::numeric_limits< double >::infinity();

    // Shifted (unnormalized) probabilities
    double z = 0.0;
    for (size_t i = 0u; i < n; ++i)
    {
        probs[i] = *(support + i * k) * std::exp(probs[i] - eta_max);
        z += probs[i];
    }

    for (size_t i = 0u; i < n; ++i)
    {

        probs[i] /= z;

        const double * row = support + i * k + 1u;
        for (size_t j = 0u; j < nterms; ++j)
            mu[j] += probs[i] * row[j];

    }

    return eta_max + std::log(z);

}

/**
 * @brief Adds the (weighted) covariance of the statistics of a support
 * @details Adds \f$-m \sum_i p_i (s(y_i) - \mu)(s(y_i) - \mu)^{\mathbf{t}}\f$ to the
 * lower triangle of `H`. The inputs `probs` and `mu` are the outputs of
 * `update_expected_stats()` for the same support and parameters.
 * @param support,k,n See `update_expected_stats()`.
 * @param probs Probabilities of each row of the support.
 * @param mu Expected statistics (`k - 1` elements).
 * @param weight Multiplier \f$m\f$ (number of arrays sharing the support).
 * @param H Pointer to a row-major `(k - 1) x (k - 1)` matrix.
 * @param dev Buffer of length `k - 1` used to store \f$s(y_i) - \mu\f$.
 * @see `doc/likelihood-derivatives.md` for the derivation.
 */
inline void add_covariance(
    const double * support,
    size_t k,
    size_t n,
    const std::vector< double > & probs,
    const double * mu,
    double weight,
    double * H,
    double * dev
)
{

    const size_t nterms = k - 1u;
    for (size_t i = 0u; i < n; ++i)
    {

        if (probs[i] == 0.0)
            continue;

        const double * row = support + i * k + 1u;
        for (size_t j = 0u; j < nterms; ++j)
            dev[j] = row[j] - mu[j];

        const double w = weight * probs[i];
        for (size_t j0 = 0u; j0 < nterms; ++j0)
        {

            const double wd = w * dev[j0];
            double * H_row = H + j0 * nterms;
            for (size_t j1 = 0u; j1 <= j0; ++j1)
                H_row[j1] -= wd * dev[j1];

        }

    }

    return;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type, Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::update_normalizing_constants(
    const std::vector< double > & params,
    size_t ncores,
    int i
) {

    const size_t n = stats_support_sizes.size();

    // Barrier to make sure paralelization makes sense
    if ((ncores > 1u) && (n < 128u))
        ncores = 1u;

    
    if (i >= 0)
        ncores = 1u;
    
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for firstprivate(params) num_threads(ncores) \
        shared(n, normalizing_constants, first_calc_done, \
            stats_support, stats_support_sizes, stats_support_sizes_acc, i) \
        default(none)
    #endif
    for (size_t s = 0u; s < n; ++s)
    {

        if ((i > -1) && (i != static_cast<int>(s)))
            continue;

        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[s];

        first_calc_done[s] = true;
        normalizing_constants[s] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[s] * k
                ], k, n
        );

    }

    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::update_likelihoods(
    const std::vector< double > & params,
    size_t ncores
) {
    
    update_normalizing_constants(params, ncores);

    size_t n_params = params.size();

    if (stats_likelihood.size() != stats_target.size())
        stats_likelihood.resize(stats_target.size());

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for simd num_threads(ncores) \
        shared(n_params, stats_target, normalizing_constants, arrays2support, \
            params) \
        default(none)
    #endif
    for (size_t s = 0u; s < stats_target.size(); ++s)
    {
        stats_likelihood[s] = 0.0;
        for (size_t j = 0u; j < n_params; ++j)
            stats_likelihood[s] += stats_target[s][j] * params[j];

        stats_likelihood[s] =
            std::exp(stats_likelihood[s] BARRY_SAFE_EXP)/
            normalizing_constants[arrays2support[s]];
    }
    
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::update_pset_probs(
    const std::vector< double > & params,
    size_t ncores,
    int i
) {

    update_normalizing_constants(params, ncores, i);
    
    if (i > -1)
        params_last[i] = params;

    size_t n_params = params.size();
    pset_probs.resize(
        pset_locations.back() + 
        pset_sizes.back()
        );

    // No need to paralelize if there is only one core
    if (i >= 0)
       ncores = 1u; 

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for num_threads(ncores) collapse(1) \
        shared(n_params, pset_stats, pset_probs, normalizing_constants, pset_sizes, \
            params, i) \
        default(none)
    #endif
    for (size_t s = 0u; s < pset_sizes.size(); ++s)
    {

        if ((i >= 0) && (i != static_cast<int>(s)))
            continue;

        // When does the pset starts
        size_t pset_start = pset_locations[s];

        // Looping over observations of the pset
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd 
        #endif
        for (size_t a = 0u; a < pset_sizes[s]; ++a)
        {

            // Start location in the array
            size_t start_loc = pset_start * n_params + a * n_params;
            
            pset_probs[pset_start + a] = 0.0;

            // Looping over the parameters
            for (size_t j = 0u; j < n_params; ++j)
                pset_probs[pset_start + a] +=
                    pset_stats[start_loc + j] * params[j];

            // Now turning into a probability
            pset_probs[pset_start + a] =
                std::exp(pset_probs[pset_start + a] BARRY_SAFE_EXP)/
                normalizing_constants[s];
        }

        #ifdef BARRY_DEBUG
        // Making sure the probabilities add to one
        double totprob = 0.0;
        for (size_t i_ = 0u; i_ < pset_sizes[s]; ++i)
            totprob =+ pset_probs[pset_start + i_];

        if (std::abs(totprob - 1) > 1e-6)
            throw std::runtime_error(
                std::string("Probabilities do not add to one! ") +
                std::string("totprob = ") + std::to_string(totprob)
            );

        #endif
    }

    // Recording which parameters each updated block corresponds to
    pset_params_last.resize(pset_sizes.size());
    for (size_t s = 0u; s < pset_sizes.size(); ++s)
        if ((i < 0) || (i == static_cast<int>(s)))
            pset_params_last[s] = params;

    return;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model() :
    stats_support(0u),
    stats_support_sizes(0u),
    stats_support_sizes_acc(0u),
    stats_support_n_arrays(0u),
    stats_target(0u),
    stats_likelihood(0u),
    arrays2support(0u),
    keys2support(0u),
    pset_arrays(0u), pset_stats(0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(nullptr),
    transform_model_term_names(0u)
{  

    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);
    
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model(
    size_t size_
    ) :
    stats_support(0u),
    stats_support_sizes(0u),
    stats_support_sizes_acc(0u),
    stats_support_n_arrays(0u),
    stats_target(0u),
    stats_likelihood(0u),
    arrays2support(0u), keys2support(0u), 
    pset_arrays(0u), pset_stats(0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(nullptr),
    transform_model_term_names(0u)
{
    
    stats_target.reserve(size_);
    arrays2support.reserve(size_);

    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);
        
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model(
    const Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & Model_
    ) : 
    stats_support(Model_.stats_support),
    stats_support_sizes(Model_.stats_support_sizes),
    stats_support_sizes_acc(Model_.stats_support_sizes_acc),
    stats_support_n_arrays(Model_.stats_support_n_arrays),
    stats_target(Model_.stats_target),
    stats_likelihood(Model_.stats_likelihood),
    arrays2support(Model_.arrays2support),
    keys2support(Model_.keys2support),
    pset_arrays(Model_.pset_arrays),
    pset_stats(Model_.pset_stats),
    pset_probs(Model_.pset_probs),
    pset_sizes(Model_.pset_sizes),
    pset_locations(Model_.pset_locations),
    counters(new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters))),
    rules(new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules))),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn))),
    support_fun(),
    counter_fun(),
    params_last(Model_.params_last),
    normalizing_constants(Model_.normalizing_constants),
    first_calc_done(Model_.first_calc_done),
    delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(Model_.transform_model_fun),
    transform_model_term_names(Model_.transform_model_term_names)
    {

    pset_params_last = Model_.pset_params_last;
    with_pset        = Model_.with_pset;

    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);

    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);

    return;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & 
    Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::operator=(
    const Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & Model_
) {
    
    // Clearing
    if (this != &Model_) {

        if (delete_counters)
            delete counters;

        if (delete_rules)
            delete rules;
        
        if (delete_rules_dyn)
            delete rules_dyn;
        
        stats_support              = Model_.stats_support;
        stats_support_sizes        = Model_.stats_support_sizes;
        stats_support_sizes_acc    = Model_.stats_support_sizes_acc;
        stats_support_n_arrays     = Model_.stats_support_n_arrays;
        stats_target               = Model_.stats_target;
        stats_likelihood           = Model_.stats_likelihood;
        arrays2support             = Model_.arrays2support;
        keys2support               = Model_.keys2support;
        pset_arrays                = Model_.pset_arrays;
        pset_stats                 = Model_.pset_stats;
        pset_probs                 = Model_.pset_probs;
        pset_sizes                 = Model_.pset_sizes;
        pset_locations             = Model_.pset_locations;
        counters                   = new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters));
        rules                      = new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules));
        rules_dyn                  = new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn));
        delete_counters            = true;
        delete_rules               = true;
        delete_rules_dyn           = true;
        params_last                = Model_.params_last;
        normalizing_constants      = Model_.normalizing_constants;
        first_calc_done            = Model_.first_calc_done;
        pset_params_last           = Model_.pset_params_last;
        with_pset                  = Model_.with_pset;
        transform_model_fun        = Model_.transform_model_fun;
        transform_model_term_names = Model_.transform_model_term_names;

        // Counters are shared
        support_fun.set_counters(counters);
        counter_fun.set_counters(counters);
        
        // Rules are shared
        support_fun.set_rules(rules);
        support_fun.set_rules_dyn(rules_dyn);
        
    }
        
    return *this;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: store_psets() noexcept {
    with_pset = true;
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: gen_key(
    const Array_Type & Array_
) {
    return this->counters->gen_hash(Array_);   
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_counter(
        Counter<Array_Type, Data_Counter_Type> & counter
) {
    
    counters->add_counter(counter, Data_Counter_Type());
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_counter(
    Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_,
    Data_Counter_Type                              data_
) {
    
    counters->add_counter(
        count_fun_,
        init_fun_,
        data_
    );
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_counters(
    Counters<Array_Type,Data_Counter_Type> * counters_
) {

    if (delete_counters) {
        delete counters;
        delete_counters = false;
    }
    
    this->counters = counters_;
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_hasher(
    Hasher_fun_type<Array_Type,Data_Counter_Type> fun_
) {

    counters->add_hash(fun_);

}

////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule(
    Rule<Array_Type, Data_Rule_Type> & rules
) {
    
    rules->add_rule(rules, Data_Rule_Type());
    return;
}


template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_rules(
    Rules<Array_Type,Data_Rule_Type> * rules_
) {

    if (delete_rules)
        delete rules;

    this->rules = rules_;
    this->delete_rules = false;

    support_fun.set_rules(rules);
    return;

}

////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule_dyn(
    Rule<Array_Type, Data_Rule_Dyn_Type> & rules_
) {
    
    rules_dyn->add_rule(rules_, Data_Rule_Dyn_Type());
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule_dyn(
    Rule_fun_type<Array_Type,Data_Rule_Dyn_Type> rule_fun_,
    Data_Rule_Dyn_Type                           data_
) {
    
    rules_dyn->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_rules_dyn(
    Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {

    if (delete_rules_dyn)
        delete rules_dyn;

    this->rules_dyn = rules_;
    this->delete_rules_dyn = false;
    support_fun.set_rules_dyn(rules_dyn);
    return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::add_array(
    const Array_Type & Array_,
    bool force_new
) {
    
    // Array counts (target statistics)
    counter_fun.reset_array(&Array_);
    
    if (transform_model_fun)
    {
        
        auto tmpcounts = counter_fun.count_all();
        stats_target.emplace_back(
            transform_model_fun(&tmpcounts[0u], tmpcounts.size())
            );

    } else
        stats_target.push_back(counter_fun.count_all());
    
    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = counters->gen_hash(Array_);
    MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
    if (force_new | (locator == keys2support.end()))
    {

        // Current size of the support stats
        size_t stats_support_size = stats_support.size();
        
        // Adding to the map
        keys2support[key] = stats_support_sizes.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support_sizes.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);

            // Current size of the powerset
            size_t pset_stats_size = pset_stats.size();
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &pset_stats
                );
                
            }
            catch (const std::exception& e)
            {
                
                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                printf_barry("Here is the array that generated the error.\n");
                Array_.print();
                throw std::logic_error("");
                
            }

            // Recording the number of elements
            pset_locations.push_back(
                pset_locations.size() == 0u ?
                    0u :
                    pset_locations.back() + pset_sizes.back()
                );

            pset_sizes.push_back(
                (pset_stats.size() - pset_stats_size) / (counter_fun.size())
                );
                
                

            
        }
        else
        {
            try
            {

                support_fun.calc();
                
            }
            catch (const std::exception& e)
            {

                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                printf_barry("Here is the array that generated the error.\n");
                Array_.print();
                throw std::logic_error("");

            }
        }
        
        if (transform_model_fun)
        {
            auto tmpsupport = support_fun.get_counts();
            size_t k = counter_fun.size();
            size_t n = tmpsupport.size() / (k + 1);

            std::vector< double > s_new(0u);            
            s_new.reserve(tmpsupport.size());

            for (size_t i = 0u; i < n; ++i)
            {

                // Appending size
                s_new.push_back(tmpsupport[i * (k + 1u)]);

                // Applying transformation and adding to the new set
                auto res = transform_model_fun(&tmpsupport[i * (k + 1u) + 1u], k);
                std::copy(res.begin(), res.end(), std::back_inserter(s_new));

            }

            for (auto & s : s_new)
                stats_support.push_back(s);
            

        } else {
            for (const auto & s: support_fun.get_counts())
                stats_support.push_back(s);
        }
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);

        // Incrementing the size of the support set
        if (stats_support_sizes.size() == 0u)
        {
            stats_support_sizes_acc.push_back(0u);    
        } else {
            stats_support_sizes_acc.push_back(
                stats_support_sizes.back() + 
                stats_support_sizes_acc.back()
            );
        }


        stats_support_sizes.push_back(
            
            (stats_support.size() - stats_support_size)/
                (counter_fun.size() + 1u)

            );
        
        return arrays2support.size() - 1u;
        
    }
    
    // Increasing the number of arrays in that stat
    ++stats_support_n_arrays[locator->second];
    
    // Adding the corresponding map
    arrays2support.push_back(locator->second);
    
    return arrays2support.size() - 1u;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t idx = arrays2support[i];

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[idx] == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[idx] || !vec_equal_approx(params, params_last[idx])))
    {
        
        first_calc_done[idx] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[idx];

        normalizing_constants[idx] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[idx] * k
                ], k, n
        );
        
        params_last[idx] = params;
        
    }
    
    return likelihood_(
        &stats_target[i],
        params,
        normalizing_constants[idx],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const Array_Type & Array_,
    int i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Key of the support set to use
    int loc;

    if (i < 0)
    {

        std::vector< double > key = counters->gen_hash(Array_);
        MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
        if (locator == keys2support.end()) 
            throw std::range_error(
                "This type of array has not been included in the model."
                );

        loc = locator->second;

    }
    else
    {

        if (static_cast<size_t>(i) >= arrays2support.size())
            throw std::range_error(
                "This type of array has not been included in the model."
                );

        loc = arrays2support[i];

    }

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Counting stats_target
    StatsCounter< Array_Type, Data_Counter_Type> tmpstats(&Array_);

    tmpstats.set_counters(this->counters);
    
    std::vector< double > target_ = tmpstats.count_all();

    if (transform_model_fun)
        target_ = transform_model_fun(&target_[0u], target_.size());

    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc])) )
    {
        
        first_calc_done[loc] = true;

        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];
        
        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
                ], k, n
        );
        
        params_last[loc] = params;
        
    }

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    return likelihood_(
        &target_[0u],
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const std::vector<double> & target_,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t loc = arrays2support[i];

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
    {

        // Concatenating the elements of target_ into aa single string
        std::string target_str = "";
        for (size_t i = 0u; i < target_.size(); ++i)
            target_str += std::to_string(target_[i]) + " ";

        throw std::range_error(
            "The array is not in the support set. The array's statistics are: " +
            target_str +
            std::string(".")
            );
    }
        

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
    {
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc])) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];

        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
                ], k, n
        );
        
        params_last[loc] = params;
        
    }
    
    return likelihood_(
        &target_[0u],
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const double * target_,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t loc = arrays2support[i];

    // Checking if passes the rules
    if (support_fun.get_rules_dyn()->size() > 0u)
    {

        std::vector< double > tmp_target;
        tmp_target.reserve(nterms());
        for (size_t t = 0u; t < nterms(); ++t)
            tmp_target.push_back(*(target_ + t));

        if (!support_fun.eval_rules_dyn(tmp_target, 0u, 0u))
        {
            // Concatenating the elements of target_ into aa single string
            std::string target_str = "";
            for (size_t i = 0u; i < nterms(); ++i)
                target_str += std::to_string((*target_ + i)) + " ";

            throw std::range_error(
                "The array is not in the support set. The array's statistics are: " + target_str + std::string(".")
                );
        }

    }

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
    {
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) )) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];

        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
            ], k, n
        );
        
        params_last[loc] = params;
        
    }
    
    return likelihood_(
        target_,
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood_total(
    const std::vector<double> & params,
    bool as_log,
    BARRY_NCORES_ARG(),
    bool no_update_normconst
) {
    
    size_t params_last_size = params_last.size();

    if (!no_update_normconst)
    {
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp parallel for num_threads(ncores) \
            shared(normalizing_constants, params_last, first_calc_done, \
                stats_support, stats_support_sizes, stats_support_sizes_acc) \
            firstprivate(params)
        #endif
        for (size_t i = 0u; i < params_last_size; ++i)
        {

            if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) )
            {

                size_t k = params.size() + 1u;
                size_t n = stats_support_sizes[i];
                
                first_calc_done[i] = true;
                normalizing_constants[i] = update_normalizing_constant(
                    params, &stats_support[
                        stats_support_sizes_acc[i] * k
                    ], k, n
                );
                
                params_last[i] = params;
                
            }

        }
    }
    
    double res = 0.0;
    if (as_log)
    {

        for (size_t i = 0; i < stats_target.size(); ++i) 
            res += vec_inner_prod(
                &stats_target[i][0u],
                &params[0u],
                params.size()
                ) BARRY_SAFE_EXP;
        
        #if defined(__OPENMP) || defined(_OPENMP) 
        #pragma omp simd reduction(-:res)
        #endif
        for (size_t i = 0u; i < params_last_size; ++i)
            res -= (std::log(normalizing_constants[i]) * this->stats_support_n_arrays[i]);

    } else {
        
        res = 1.0;
        size_t stats_target_size = stats_target.size();
        #if defined(__OPENMP) || defined(_OPENMP) 
        #pragma omp simd reduction(*:res)
        #endif
        for (size_t i = 0; i < stats_target_size; ++i)
            res *= std::exp(
                vec_inner_prod(
                    &stats_target[i][0u],
                    &params[0u],
                    params.size()
                ) BARRY_SAFE_EXP) / 
                normalizing_constants[arrays2support[i]];
        
    }
    
    return res;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood_total_and_grad(
    const std::vector<double> & params,
    std::vector<double> & grad,
    size_t ncores
) const {

    const size_t K = nterms();
    if (params.size() != K)
        throw std::length_error(
            "The length of -params- (" + std::to_string(params.size()) +
            ") doesn't match the number of terms in the model (" +
            std::to_string(K) + ")."
            );

    // Observed part of the gradient: sum of the target statistics
    grad.assign(K, 0.0);
    for (const auto & s : stats_target)
        for (size_t j = 0u; j < K; ++j)
            grad[j] += s[j];

    double res = 0.0;
    for (size_t j = 0u; j < K; ++j)
        res += grad[j] * params[j];

    // Expected part: one pass over the unique supports
    const size_t n_supports = stats_support_sizes.size();
    if ((ncores > 1u) && (n_supports < 128u))
        ncores = 1u;

    double log_z_total = 0.0;

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel num_threads(ncores) reduction(+:log_z_total)
    #endif
    {

        std::vector< double > probs;
        std::vector< double > mu(K);
        std::vector< double > mu_total(K, 0.0);

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp for schedule(static)
        #endif
        for (size_t u = 0u; u < n_supports; ++u)
        {

            const double m_u = static_cast< double >(stats_support_n_arrays[u]);
            const double log_z = update_expected_stats(
                params,
                stats_support.data() + stats_support_sizes_acc[u] * (K + 1u),
                K + 1u,
                stats_support_sizes[u],
                probs,
                mu.data()
            );

            log_z_total += m_u * log_z;

            // Empty support: the log-likelihood is -Inf
            if (std::isinf(log_z))
                continue;

            for (size_t j = 0u; j < K; ++j)
                mu_total[j] += m_u * mu[j];

        }

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp critical
        #endif
        for (size_t j = 0u; j < K; ++j)
            grad[j] -= mu_total[j];

    }

    return res - log_z_total;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< double > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood_total_grad(
    const std::vector<double> & params,
    size_t ncores
) const {

    std::vector< double > grad;
    likelihood_total_and_grad(params, grad, ncores);
    return grad;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< double > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood_total_hess(
    const std::vector<double> & params,
    size_t ncores
) const {

    const size_t K = nterms();
    if (params.size() != K)
        throw std::length_error(
            "The length of -params- (" + std::to_string(params.size()) +
            ") doesn't match the number of terms in the model (" +
            std::to_string(K) + ")."
            );

    std::vector< double > H(K * K, 0.0);

    const size_t n_supports = stats_support_sizes.size();
    if ((ncores > 1u) && (n_supports < 128u))
        ncores = 1u;

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel num_threads(ncores)
    #endif
    {

        std::vector< double > probs;
        std::vector< double > mu(K);
        std::vector< double > dev(K);
        std::vector< double > H_local(K * K, 0.0);

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp for schedule(static)
        #endif
        for (size_t u = 0u; u < n_supports; ++u)
        {

            const double * support =
                stats_support.data() + stats_support_sizes_acc[u] * (K + 1u);

            const double log_z = update_expected_stats(
                params, support, K + 1u, stats_support_sizes[u], probs,
                mu.data()
            );

            if (std::isinf(log_z))
                continue;

            add_covariance(
                support, K + 1u, stats_support_sizes[u], probs, mu.data(),
                static_cast< double >(stats_support_n_arrays[u]),
                H_local.data(), dev.data()
            );

        }

        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp critical
        #endif
        for (size_t j = 0u; j < (K * K); ++j)
            H[j] += H_local[j];

    }

    // Only the lower triangle was filled
    for (size_t j0 = 0u; j0 < K; ++j0)
        for (size_t j1 = 0u; j1 < j0; ++j1)
            H[j1 * K + j0] = H[j0 * K + j1];

    return H;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline const std::vector< double > &
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_normalizing_constants() const {
    
    return normalizing_constants;
    
}

template<
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline const std::vector< double > &
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_likelihoods() const {
    
    return stats_likelihood;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const std::vector< Array_Type > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset(
    const size_t & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");


    return &pset_arrays[arrays2support[i]];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const double *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats(
    const size_t & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    return &pset_stats[pset_locations[arrays2support[i]] * counter_fun.size()];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: print_stats(size_t i) const
{
    
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t k       = nterms();
    size_t nunique = stats_support_sizes.size();

    for (size_t l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5i ", static_cast<int>(l));

        printf_barry("counts: %.0f motif: ", stats_support[
            stats_support_sizes_acc[l] * (k + 1u) 
            // l * (k + 1u)
            ]);
        
        for (size_t j = 0u; j < k; ++j)
        {
            printf_barry(
                "%.2f, ",
                stats_support[
                    stats_support_sizes_acc[l] * (k + 1u) + j + 1u
                    ]);
        }

        printf_barry("\n");

    }
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::print() const
{

    // Relevant information:
    // - Number of arrays involved
    // - Size of the support
    // - Terms involved

    int min_v = std::numeric_limits<int>::max();
    int max_v = 0;

    for (const auto & stat : this->stats_support_sizes)
    {

        if (static_cast<int>(stat) > max_v)
            max_v = static_cast<int>(stat);
        
        if (static_cast<int>(stat) < min_v)
            min_v = static_cast<int>(stat);

    }  

    // The vectors in the support reflec the size of nterms x entries
    // max_v /= static_cast<int>(nterms() + 1);
    // min_v /= static_cast<int>(nterms() + 1);

    if (this->size() > 0u)
    {
        printf_barry(
            "Num. of Arrays       : %i\n",
            static_cast<int>(this->size())
        );
        printf_barry(
            "Support size         : %i\n",
            static_cast<int>(this->size_unique())
        );
        printf_barry("Support size range   : [%i, %i]\n", min_v, max_v);
    }
    else 
    {
        printf_barry("Num. of Arrays       : 0\n");
        printf_barry("Support size         : -\n");
        printf_barry("Support size range   : -\n");
    }
    

    if (with_pset)
    {
        printf_barry("Arrays in powerset   : %i\n",
            static_cast<int>(std::accumulate(pset_sizes.begin(), pset_sizes.end(), 0u))
        );
    }


    printf_barry("Transform. Fun.      : %s\n", transform_model_fun ? "yes": "no");
    printf_barry("Model terms (%i)    :\n", static_cast<int>(this->nterms()));
    for (auto & cn : this->colnames())
    {
        printf_barry(" - %s\n", cn.c_str());
    }

    if (this->nrules() > 0u)
    {
        printf_barry(
            "Model rules (%i)     :\n",
            static_cast<int>(this->nrules())
        );
    
        for (auto & rn : rules->get_names())
        {
            printf_barry(" - %s\n", rn.c_str());
        }
    }

    if (this->nrules_dyn() > 0u)
    {
        printf_barry(
            "Model rules dyn (%i):\n",
            static_cast<int>(this->nrules_dyn())
        );
    
        for (auto & rn : rules_dyn->get_names())
        {
            printf_barry(" - %s\n", rn.c_str());
        }
    }

    return;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: size() const noexcept
{
    // INITIALIZED()
    return this->stats_target.size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: size_unique() const noexcept
{

    // INITIALIZED()
    return this->stats_support_sizes.size();

} 

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nterms() const noexcept
{
 
    if (transform_model_fun)
        return transform_model_term_names.size();
    else
        return this->counters->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nrules() const noexcept
{
 
    return this->rules->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nrules_dyn() const noexcept
{
 
    return this->rules_dyn->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: support_size() const noexcept
{

    // INITIALIZED()
    return stats_support_sizes_acc.back();
    // size_t tot = 0u;
    // for (auto& a : stats_support)
    //     tot += a.size();

    // return tot;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::string > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: colnames() const
{
    
    if (transform_model_fun)
        return transform_model_term_names;
    else
        return counters->get_names();

}
    
template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Array_Type
Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::sample(
    const size_t & i,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // Getting the index
    size_t a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    // Making sure the powerset probabilities are current
    if (!pset_probs_current(a, params))
        update_pset_probs(params, 1u, static_cast<int>(a));

    // Sampling an array (bounded, in case roundoff leaves cumprob < r)
    size_t j = 0u;
    const double * probs = &pset_probs[pset_locations[a]];
    while ((j < pset_sizes[a]) && (cumprob < r))
        cumprob += *(probs + j++);

    if (j > 0u)
        j--;
    
    #ifdef BARRY_DEBUG
    return this->pset_arrays.at(a).at(j);   
    #else
    return this->pset_arrays[a][j];   
    #endif

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Array_Type Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: sample(
    const Array_Type & Array_,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    size_t i;

    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = counters->gen_hash(Array_);
    MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
    if (locator == keys2support.end())
    {
        size_t stats_support_size = stats_support.size();

        // Adding to the map
        keys2support[key] = stats_support_sizes.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support_sizes.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Current size of the powerset
            size_t pset_stats_size = pset_stats.size();

            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);
            // pset_stats.resize(pset_stats.size() + 1u);
            // pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &pset_stats
                );
                
            }
            catch (const std::exception& e)
            {
                
                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                throw std::logic_error("");
                
            }

            // Recording the number of elements
            pset_locations.push_back(
                pset_locations.size() == 0u ?
                    0u :
                    pset_locations.back() + pset_sizes.back()
                );

            pset_sizes.push_back(
                (pset_stats.size() - pset_stats_size) / (counter_fun.size())
                );

            // Increasing the space to store probabilities
            pset_probs.resize(pset_probs.size() + pset_sizes.back());

            
        }
        else
        {
            support_fun.calc();
        }
        
        if (transform_model_fun)
        {
            auto tmpsupport = support_fun.get_counts();
            size_t k = counter_fun.size();
            size_t n = tmpsupport.size() / (k + 1);

            std::vector< double > s_new(0u);            
            s_new.reserve(tmpsupport.size());

            for (size_t i = 0u; i < n; ++i)
            {

                // Appending size
                s_new.push_back(tmpsupport[i * (k + 1u)]);

                // Applying transformation and adding to the new set
                auto res = transform_model_fun(&tmpsupport[i * (k + 1u) + 1u], k);
                std::copy(res.begin(), res.end(), std::back_inserter(s_new));

            }

            for (auto & s : s_new)
                stats_support.push_back(s);
            // stats_support.push_back(s_new);

        } else {
            for (auto & s : support_fun.get_counts())
                stats_support.push_back(s);

            // stats_support.push_back(support_fun.get_counts());
        }
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);

        // Incrementing the size of the support set
        if (stats_support_sizes.size() == 0u)
        {
            stats_support_sizes_acc.push_back(0u);    
        } else {
            stats_support_sizes_acc.push_back(
                stats_support_sizes.back() + 
                stats_support_sizes_acc.back()
            );
        }


        stats_support_sizes.push_back(
            
            (stats_support.size() - stats_support_size)/
                (counter_fun.size() + 1u)

            );

        
        i = arrays2support.size() - 1u;
    } else
        // Retrieving the corresponding position in the support
        i = locator->second;

    // Getting the support index. In the new-array branch above, `i` is an
    // array index and must be resolved through `arrays2support`. In the
    // reuse branch, `keys2support` already stores a SUPPORT index, so
    // `locator->second` (assigned to `i`) is the support index directly and
    // must NOT be mapped through `arrays2support` a second time.
    size_t a = (locator == keys2support.end()) ?
        arrays2support[i] : i;

    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    // Making sure the powerset probabilities for this support are current
    // before sampling from them (mirrors the sample(size_t, params) path).
    // `pset_params_last` is only set by update_pset_probs(), so this also
    // covers an empty `pset_probs` and a prior likelihood() call that updated
    // only the normalizing constants.
    if (!pset_probs_current(a, params))
        update_pset_probs(params, 1u, static_cast<int>(a));

    // Sampling an array (bounded, in case roundoff leaves cumprob < r)
    size_t j = 0u;
    double * probs = &pset_probs[ pset_locations[a] ];
    while ((j < pset_sizes[a]) && (cumprob < r))
        cumprob += *(probs + j++);

    if (j > 0u)
        j--;


    #ifdef BARRY_DEBUG
    return this->pset_arrays.at(a).at(j);   
    #else
    return this->pset_arrays[a][j];   
    #endif

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: conditional_prob(
    const Array_Type & Array_,
    const std::vector< double > & params,
    size_t i,
    size_t j
) {

    // Generating a copy of the array so we can update
    Array_Type A(Array_, true);

    // Making sure we add it first
    A.insert_cell(i, j, A.default_val(), true, false);

    // Computing the change stats_target
    std::vector< double > tmp_counts;
    tmp_counts.reserve(counters->size());
    for (size_t ii = 0u; ii < counters->size(); ++ii)
        tmp_counts.push_back(counters->operator[](ii).count(A, i, j));

    // If there is a transformation function, it needs to be
    // applied before dealing with the likelihood.
    if (transform_model_fun)
        tmp_counts = transform_model_fun(&tmp_counts[0u], tmp_counts.size());

    return 1.0/
        (1.0 + std::exp(-vec_inner_prod<double>(
            &params[0u], &tmp_counts[0u], params.size()
            )));

    
}

#define MODEL_TEMPLATE_ARGS() <Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>

#define MODEL_TYPE() Model MODEL_TEMPLATE_ARGS()

#define MODEL_TEMPLATE(a,b) \
    template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type> \
    inline a MODEL_TYPE()::b

MODEL_TEMPLATE(double, runif_)() {

    if (rengine == nullptr)
    {
        rengine        = new std::mt19937(std::random_device{}());
        delete_rengine = true;
    }

    std::uniform_real_distribution<> urand(0, 1);
    return urand(*rengine);

}

MODEL_TEMPLATE(std::vector< double >, transform_stats_)(
    const std::vector< double > & stats
) {

    if (!transform_model_fun)
        return stats;

    std::vector< double > tmp(stats);
    return transform_model_fun(&tmp[0u], tmp.size());

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline typename MODEL_TYPE()::SampleApproxState MODEL_TYPE()::sample_approx_init(
    const Array_Type & Array_
) {

    SampleApproxState state;
    state.array = Array_;

    // Free cells, as in the support. Locked cells keep their values from
    // Array_ (the support conditions on them too, e.g., DEFM's Markov rows),
    // so non-zero locked cells are not an error.
    rules->get_seq(state.array, &state.free, nullptr);

    // Current ties (non-zero free cells)
    size_t nfree = state.free.size() / 2u;
    state.pos.resize(nfree, -1);
    for (size_t k = 0u; k < nfree; ++k)
    {

        const Array_Type & A = state.array;
        if (static_cast< double >(A(state.free[k * 2u], state.free[k * 2u + 1u], false)) != 0.0)
        {
            state.pos[k] = static_cast< int >(state.ties.size());
            state.ties.push_back(k);
        }

    }

    // Initial (raw) sufficient statistics
    counter_fun.reset_array(&state.array);
    state.stats = counter_fun.count_all();

    if (!support_fun.eval_rules_dyn(state.stats, 0u, 0u))
        throw std::logic_error(
            "The starting array does not satisfy the dynamic rules."
        );

    return state;

}

MODEL_TEMPLATE(void, sample_approx_tnt_step_)(
    SampleApproxState & state,
    const std::vector< double > & params
) {

    size_t D = state.free.size() / 2u;
    size_t E = state.ties.size();

    state.nsteps++;
    if (D == 0u)
        return;

    // Picking a cell: with probability 1/2 (if there are ties) a tie, otherwise
    // any free cell. Indices are drawn from uniforms so results do not depend
    // on the standard library's integer distributions.
    size_t k;
    if ((E > 0u) && (runif_() < 0.5))
        k = state.ties[std::min(static_cast< size_t >(runif_() * E), E - 1u)];
    else
        k = std::min(static_cast< size_t >(runif_() * D), D - 1u);

    size_t i = state.free[k * 2u];
    size_t j = state.free[k * 2u + 1u];
    bool is_tie = state.pos[k] >= 0;

    // Change statistics are computed with the cell set to one, as when
    // computing the support.
    Array_Type & A = state.array;
    if (!is_tie)
        A.insert_cell(i, j, A.default_val(), false, false);

    size_t nterms_raw = counters->size();
    std::vector< double > stats_new(state.stats);
    for (size_t n = 0u; n < nterms_raw; ++n)
    {

        double delta = counters->operator[](n).count(A, i, j);
        stats_new[n] += is_tie ? -delta : delta;

    }

    // Change in the linear predictor
    double logratio;
    if (transform_model_fun)
    {

        std::vector< double > t_old = transform_stats_(state.stats);
        std::vector< double > t_new = transform_stats_(stats_new);
        logratio = 0.0;
        for (size_t n = 0u; n < params.size(); ++n)
            logratio += params[n] * (t_new[n] - t_old[n]);

    } else {

        logratio = 0.0;
        for (size_t n = 0u; n < params.size(); ++n)
            logratio += params[n] * (stats_new[n] - state.stats[n]);

    }

    // Proposal ratio q(y' -> y) / q(y -> y'). The probability of proposing
    // cell c from a state with e ties is 1/2 * 1{c is a tie} / e + 1/2 / D
    // if e > 0, and 1 / D otherwise.
    auto q = [D](bool tie, size_t e) -> double {
        if (e == 0u)
            return 1.0 / static_cast< double >(D);
        return (tie ? 0.5 / static_cast< double >(e) : 0.0) +
            0.5 / static_cast< double >(D);
    };

    size_t E_new = is_tie ? E - 1u : E + 1u;
    logratio += std::log(q(!is_tie, E_new)) - std::log(q(is_tie, E));

    // Accepting or rejecting
    bool accept = std::log(runif_()) < logratio;
    if (accept && (support_fun.get_rules_dyn()->size() > 0u))
        accept = support_fun.eval_rules_dyn(stats_new, i, j);

    if (!accept)
    {

        if (!is_tie)
            A.rm_cell(i, j, false, false);

        return;

    }

    state.naccepted++;
    state.stats = std::move(stats_new);

    if (is_tie)
    {

        // Swap-with-last removal from the list of ties
        size_t p    = static_cast< size_t >(state.pos[k]);
        size_t last = state.ties.back();
        state.ties[p]    = last;
        state.pos[last]  = static_cast< int >(p);
        state.ties.pop_back();
        state.pos[k] = -1;

        A.rm_cell(i, j, false, false);

    } else {

        state.pos[k] = static_cast< int >(state.ties.size());
        state.ties.push_back(k);

    }

    return;

}

MODEL_TEMPLATE(void, sample_approx)(
    SampleApproxState & state,
    const std::vector< double > & params,
    size_t n,
    size_t thin,
    std::vector< double > * stats_out,
    SampleApproxMethod method
) {

    if (params.size() != nterms())
        throw std::length_error(
            "The length of params (" + std::to_string(params.size()) +
            ") does not match the number of terms (" +
            std::to_string(nterms()) + ")."
        );

    if (thin == 0u)
        throw std::invalid_argument("thin must be greater than zero.");

    if (stats_out != nullptr)
        stats_out->reserve(stats_out->size() + n * nterms());

    for (size_t s = 0u; s < n; ++s)
    {

        for (size_t t = 0u; t < thin; ++t)
        {
            switch (method)
            {
            case SampleApproxMethod::tnt:
                sample_approx_tnt_step_(state, params);
                break;
            default:
                throw std::invalid_argument("Unknown sampling method.");
            }
        }

        if (stats_out != nullptr)
        {
            std::vector< double > tmp = transform_stats_(state.stats);
            stats_out->insert(stats_out->end(), tmp.begin(), tmp.end());
        }

    }

    return;

}

MODEL_TEMPLATE(Array_Type, sample_approx)(
    const Array_Type & Array_,
    const std::vector< double > & params,
    size_t nsteps,
    std::vector< double > * stats_out,
    SampleApproxMethod method
) {

    SampleApproxState state = sample_approx_init(Array_);
    sample_approx(state, params, nsteps, 1u, stats_out, method);

    return state.array;

}

#undef MODEL_TEMPLATE
#undef MODEL_TYPE
#undef MODEL_TEMPLATE_ARGS

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const std::mt19937 * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_rengine() const {
    return this->rengine;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Counters<Array_Type,Data_Counter_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_counters() {
    return this->counters;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules() {
    return this->rules;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Dyn_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules_dyn() {
    return this->rules_dyn;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_support_fun() {
    return &this->support_fun;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::vector< double > > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_stats_target()
{
    return &stats_target;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support()
{
    return &stats_support;
}

// Implementation of get_stats_support_sizes()
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support_sizes()
{
    return &stats_support_sizes;
}

// Implementation of get_stats_support_sizes_acc()
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support_sizes_acc()
{
    return &stats_support_sizes_acc;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_arrays2support()
{
    return &arrays2support;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::vector< Array_Type > > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_arrays() {
    return &pset_arrays;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector<double> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats() {
    return &pset_stats;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector<double> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_probs() {
    return &pset_probs;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< size_t > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_pset_sizes()
{
    return &pset_sizes;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< size_t > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_pset_locations()
{
    return &pset_locations;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::set_transform_model(
    std::function<std::vector<double>(double *,size_t)> fun,
    std::vector< std::string > names
    )
{

    if (transform_model_fun)
        throw std::logic_error("A transformation function for the model has already been established.");
    
    transform_model_fun = fun;
    transform_model_term_names = names;

    size_t k = counters->size(); 

    auto stats_support_old = stats_support;

    // Applying over the support
    for (size_t nsupport = 0u; nsupport < stats_support_sizes.size(); ++nsupport)
    {

        // How many observations in the support
        size_t n = stats_support_sizes[nsupport];

        // Iterating through each observation in the nsupport'th 
        for (size_t i = 0; i < n; ++i)
        {

            // Applying transformation and adding to the new set
            auto res = transform_model_fun(
                &stats_support_old[
                    stats_support_sizes_acc[nsupport] * (k + 1u) +
                    i * (k + 1u) + 1u
                    ],
                k
                );

            if (res.size() != transform_model_term_names.size())
                throw std::length_error(
                    std::string("The transform vector from -transform_model_fun- ") +
                    std::string(" does not match the size of ") + 
                    std::string("-transform_model_term_names-.")
                    );

            // Resizing stats_support if the transform stats do not match the
            // previous size
            if ((nsupport == 0u) && (i == 0u) && (res.size() != k))
                stats_support.resize(
                    (res.size() + 1) * (
                        stats_support_sizes_acc.back() +
                        stats_support_sizes.back()
                        )
                );

            // Weigth
            stats_support[
                stats_support_sizes_acc[nsupport] * (res.size() + 1u) +
                (res.size() + 1u) * i
                ] = stats_support_old[
                    stats_support_sizes_acc[nsupport] * (k + 1u) +
                    i * (k + 1u)
                ];

            // Copying the rest of the elements
            for (size_t j = 0u; j < res.size(); ++j)
                stats_support[
                    stats_support_sizes_acc[nsupport] * (res.size() + 1u) +
                    (res.size() + 1u) * i + j + 1u
                    ] = res[j];

        }

    }

    // Applying over the target statistics
    for (auto & s : stats_target)
        s = transform_model_fun(&s[0u], k);

    // Checking if there is support included
    if (with_pset)
    {

        // Applying it to the support
        for (auto s = 0u; s < pset_arrays.size(); ++s)
        {
            std::vector< double > new_stats;
            size_t pset_stats_loc = pset_locations[s] * k;

            for (auto a = 0u; a < pset_arrays[s].size(); ++a)
            {
                // Computing the transformed version of the data
                auto tmpstats = transform_model_fun(
                    &pset_stats[pset_stats_loc + a * k], k
                    );

                // Storing the new values
                for (auto p = 0u; p < k; ++p)
                    new_stats.push_back(tmpstats[p]);
            }

            // Updating the dataset
            for (size_t stat = 0u; stat < new_stats.size(); ++stat)
                pset_stats[pset_stats_loc + stat] = new_stats[stat];

        }

    }

    // And, resizing the last set of parameters
    for (auto & p : params_last)
        p.resize(transform_model_term_names.size());

    // The pset stats changed, so stored probabilities are no longer valid
    pset_params_last.clear();

    return;

}

#undef MODEL_TEMPLATE
#undef MODEL_TEMPLATE_ARGS
#undef MODEL_TYPE

#endif
