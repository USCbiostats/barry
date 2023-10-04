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
    #pragma code_align 32
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

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type, Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::update_normalizing_constants(
    const std::vector< double > & params,
    size_t ncores
) {

    const size_t n = stats_support_sizes.size();

    // Barrier to make sure paralelization makes sense
    if ((ncores > 1u) && (n < 128u))
        ncores = 1u;
    
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for firstprivate(params) num_threads(ncores) \
        shared(n, normalizing_constants, first_calc_done, \
            stats_support, stats_support_sizes, stats_support_sizes_acc) \
        default(none)
    #endif
    for (size_t i = 0u; i < n; ++i)
    {

        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[i];

        first_calc_done[i] = true;
        normalizing_constants[i] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[i] * k
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
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model() :
    stats_support(0u),
    stats_support_sizes(0u),
    stats_support_sizes_acc(0u),
    stats_support_n_arrays(0u),
    stats_target(0u), arrays2support(0u),
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
    stats_target(0u), arrays2support(0u), keys2support(0u), 
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
    arrays2support(Model_.arrays2support),
    keys2support(Model_.keys2support),
    pset_arrays(Model_.pset_arrays),
    pset_stats(Model_.pset_stats),
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
        arrays2support             = Model_.arrays2support;
        keys2support               = Model_.keys2support;
        pset_arrays                = Model_.pset_arrays;
        pset_stats                 = Model_.pset_stats;
        counters                   = new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters));
        rules                      = new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules));
        rules_dyn                  = new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn));
        delete_counters            = true;
        delete_rules               = true;
        delete_rules_dyn           = true;
        params_last                = Model_.params_last;
        normalizing_constants      = Model_.normalizing_constants;
        first_calc_done            = Model_.first_calc_done;
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
            pset_stats.resize(pset_stats.size() + 1u);
            pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &(pset_stats[pset_stats.size() - 1u])
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

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > &
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_normalizing_constants() {
    
    return normalizing_constants;
    
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
inline const std::vector< double > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats(
    const size_t & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    return &pset_stats[arrays2support[i]];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: print_stats(size_t i) const
{
    
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // const auto & S = stats_support[arrays2support[i]];
    size_t array_id = arrays2support[i];

    size_t k       = nterms();
    size_t nunique = stats_support_sizes.size();

    for (size_t l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5li ", l);

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
    max_v /= static_cast<int>(nterms() + 1);
    min_v /= static_cast<int>(nterms() + 1);

    printf_barry("Num. of Arrays       : %li\n", this->size());
    printf_barry("Support size         : %li\n", this->size_unique());
    printf_barry("Support size range   : [%i, %i]\n", min_v, max_v);
    printf_barry("Transform. Fun.      : %s\n", transform_model_fun ? "yes": "no");
    printf_barry("Model terms (%li)    :\n", this->nterms());
    for (auto & cn : this->colnames())
    {
        printf_barry(" - %s\n", cn.c_str());
    }

    if (this->nrules() > 0u)
    {
        printf_barry("Model rules (%li)     :\n", this->nrules());
    
        for (auto & rn : rules->get_names())
        {
            printf_barry(" - %s\n", rn.c_str());
        }
    }

    if (this->nrules_dyn() > 0u)
    {
        printf_barry("Model rules dyn (%li):\n", this->nrules_dyn());
    
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

    size_t k = params.size();

    // Sampling an array
    size_t j = 0u;
    std::vector< double > & probs = pset_probs[a];
    if ((probs.size() > 0u) && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += probs[j++];

        if (j > 0u)
            j--;

    } else { 
       
        probs.resize(pset_arrays[a].size());
        std::vector< double > temp_stats(params.size());
        const std::vector< double > & stats = pset_stats[a];

        int i_matches = -1;
        for (size_t array = 0u; array < probs.size(); ++array)
        {

            // Filling out the parameters
            for (auto p = 0u; p < params.size(); ++p)
                temp_stats[p] = stats[array * k + p];

            probs[array] = this->likelihood(params, temp_stats, i, false);
            cumprob += probs[array];

            if (i_matches == -1 && cumprob >= r)
                i_matches = array;
        }

        #ifdef BARRY_DEBUG
        if (i_matches < 0)
            throw std::logic_error(
                std::string(
                    "Something went wrong when sampling from a different set of.") +
                std::string("parameters. Please report this bug: ") +
                std::string(" cumprob: ") + std::to_string(cumprob) +
                std::string(" r: ") + std::to_string(r)
                );
        #endif

        j = i_matches;
        
    }
    
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
            
            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);
            pset_stats.resize(pset_stats.size() + 1u);
            pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &(pset_stats[pset_stats.size() - 1u])
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

    // Getting the index
    size_t a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    size_t k = params.size();

    // Sampling an array
    size_t j = 0u;
    std::vector< double > & probs = pset_probs[a];
    if ((probs.size() > 0u) && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += probs[j++];

        if (j > 0u)
            j--;

    } else { 
       
        probs.resize(pset_arrays[a].size());
        std::vector< double > temp_stats(params.size());
        const std::vector< double > & stats = pset_stats[a];

        int i_matches = -1;
        for (size_t array = 0u; array < probs.size(); ++array)
        {

            // Filling out the parameters
            for (auto p = 0u; p < params.size(); ++p)
                temp_stats[p] = stats[array * k + p];

            probs[array] = this->likelihood(params, temp_stats, i, false);
            cumprob += probs[array];

            if (i_matches == -1 && cumprob >= r)
                i_matches = array;
        }

        #ifdef BARRY_DEBUG
        if (i_matches < 0)
            throw std::logic_error(
                std::string(
                    "Something went wrong when sampling from a different set of.") +
                std::string("parameters. Please report this bug: ") +
                std::string(" cumprob: ") + std::to_string(cumprob) +
                std::string(" r: ") + std::to_string(r)
                );
        #endif

        j = i_matches;
        
    }
    

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
inline std::vector< std::vector<double> > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats() {
    return &pset_stats;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::vector<double> > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_probs() {
    return &pset_probs;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
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

            for (auto a = 0u; a < pset_arrays[s].size(); ++a)
            {
                // Computing the transformed version of the data
                auto tmpstats = transform_model_fun(
                    &pset_stats[s][a * k], k
                    );

                // Storing the new values
                for (auto p = 0u; p < k; ++p)
                    new_stats.push_back(tmpstats[p]);
            }

            // Updating the dataset
            std::swap(pset_stats[s], new_stats);

        }

    }

    // And, resizing the last set of parameters
    for (auto & p : params_last)
        p.resize(transform_model_term_names.size());

    return;

}

#undef MODEL_TEMPLATE
#undef MODEL_TEMPLATE_ARGS
#undef MODEL_TYPE

#endif