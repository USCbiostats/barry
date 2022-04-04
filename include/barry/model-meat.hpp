#ifndef BARRY_MODEL_MEAT_HPP 
#define BARRY_MODEL_MEAT_HPP 1

/**
 * @defgroup stat-models Statistical Models
 * @brief Statistical models available in `barry`.
 */

inline double update_normalizing_constant(
    const double * params,
    const double * support,
    size_t k,
    size_t n
)
{
    
    double res = 0.0;
    
    #ifdef __OPENMP
    #pragma omp simd reduction(+:res) 
    #else
    #pragma GCC ivdep
    #endif
    for (unsigned int i = 0u; i < n; ++i)
    {

        double tmp = 0.0;
        const double * support_n = support + i * k + 1u;
        
        for (unsigned int j = 0u; j < (k - 1u); ++j)
            tmp += (*(support_n + j)) * (*(params + j));
        
        res += std::exp(tmp BARRY_SAFE_EXP) * (*(support + i * k));

    }
    
    // This will only evaluate if the option BARRY_CHECK_FINITE
    // is defined
    BARRY_ISFINITE(res)

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
    #ifdef __OPENMP
    #pragma omp simd reduction(+:numerator)
    #else
    #pragma GCC ivdep
    #endif
    for (unsigned int j = 0u; j < params.size(); ++j)
        numerator += *(stats_target + j) * params[j];

    if (!log_)
        numerator = exp(numerator BARRY_SAFE_EXP);
    else
        return numerator BARRY_SAFE_EXP - log(normalizing_constant);

    double ans = numerator/normalizing_constant;

    if (ans > 1.0)
        printf_barry("ooo\n");

    return ans;
    
}

#define MODEL_TYPE() Model<Array_Type,Data_Counter_Type,Data_Rule_Type,\
    Data_Rule_Dyn_Type>

#define MODEL_TEMPLATE_ARGS() <typename Array_Type, typename Data_Counter_Type,\
    typename Data_Rule_Type, typename Data_Rule_Dyn_Type>

#define MODEL_TEMPLATE(a,b) \
    template MODEL_TEMPLATE_ARGS() inline a MODEL_TYPE()::b


MODEL_TEMPLATE(,Model)() :
    stats_support(0u),
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

    // Checking with the hasher function: Is this present?
    keygen = keygen_default<Array_Type>;
    
    return;
    
}

MODEL_TEMPLATE(,Model)(uint size_) :
    stats_support(0u),
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
    
    // Checking with the hasher function: Is this present?
    keygen = keygen_default<Array_Type>;
    
    return;
    
}

MODEL_TEMPLATE(,Model)(
    const MODEL_TYPE() & Model_) : 
    stats_support(Model_.stats_support),
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

    keygen = Model_.keygen;
    
    return;
    
}

MODEL_TEMPLATE(MODEL_TYPE() &, operator)=(
    const MODEL_TYPE() & Model_
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

MODEL_TEMPLATE(void, store_psets)() noexcept {
    // if (with_pset)
    //   throw std::logic_error("Powerset storage alreay activated.");
    with_pset = true;
    return;
}

MODEL_TEMPLATE(void, set_keygen)(
    std::function<std::vector<double>(const Array_Type &)> keygen_
) {
    keygen = keygen_;
    return;
}

MODEL_TEMPLATE(std::vector< double >, gen_key)(
    const Array_Type & Array_
) {
    return this->keygen(Array_);   
}

MODEL_TEMPLATE(void, add_counter)(
        Counter<Array_Type, Data_Counter_Type> & counter
) {
    
    counters->add_counter(counter, Data_Counter_Type());
    return;
}

MODEL_TEMPLATE(void, add_counter)(
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

MODEL_TEMPLATE(void, set_counters)(
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

////////////////////////////////////////////////////////////////////////////////

MODEL_TEMPLATE(void, add_rule)(
    Rule<Array_Type, Data_Rule_Type> & rules
) {
    
    rules->add_rule(rules, Data_Rule_Type());
    return;
}


MODEL_TEMPLATE(void, set_rules)(
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

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule<Array_Type, Data_Rule_Dyn_Type> & rules_
) {
    
    rules_dyn->add_rule(rules_, Data_Rule_Dyn_Type());
    return;
}

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule_fun_type<Array_Type,Data_Rule_Dyn_Type> rule_fun_,
    Data_Rule_Dyn_Type                           data_
) {
    
    rules_dyn->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

MODEL_TEMPLATE(void, set_rules_dyn)(
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

MODEL_TEMPLATE(uint, add_array)(
    const Array_Type & Array_,
    bool force_new
) {
    
    // Array counts (target statistics)
    counter_fun.reset_array(&Array_);
    
    if (transform_model_fun)
    {
        auto tmpcounts = counter_fun.count_all();
        stats_target.push_back(transform_model_fun(&tmpcounts[0u], tmpcounts.size()));
    } else
        stats_target.push_back(counter_fun.count_all());
    
    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = keygen(Array_);
    MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
    if (force_new | (locator == keys2support.end()))
    {
        
        // Adding to the map
        keys2support[key] = stats_support.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support.size()); // Map of the array id to the support
        
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
                    // "with error: " << e.what();
                throw std::logic_error("");
                
            }
            
        }
        else
        {
            
            // try
            // {

                support_fun.calc();

            // }
            // catch (const std::exception& e)
            // {
                
            //     printf_barry("A problem ocurred while trying to add the array. ");
            //     printf_barry("with error: %s", e.what());
            //     throw std::logic_error("");
                
            // }
            
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

            stats_support.push_back(s_new);

        } else 
            stats_support.push_back(support_fun.get_counts());
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);
        
        return arrays2support.size() - 1u;
        
    }
    
    // Increasing the number of arrays in that stat
    ++stats_support_n_arrays[locator->second];
    
    // Adding the corresponding map
    arrays2support.push_back(locator->second);
    
    return arrays2support.size() - 1u;

}

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    unsigned int idx = arrays2support[i];

    // Checking if this actually has a change of happening
    if (this->stats_support[idx].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[idx] || !vec_equal_approx(params, params_last[idx]) )
    {
        
        first_calc_done[idx] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[idx].size() / k;

        normalizing_constants[idx] = update_normalizing_constant(
            &params[0u], &stats_support[idx][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const Array_Type & Array_,
    int i,
    bool as_log
) {
    
    // Key of the support set to use
    int loc;

    if (i < 0)
    {

        std::vector< double > key = keygen(Array_);
        MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
        if (locator == keys2support.end()) 
            throw std::range_error("This type of array has not been included in the model.");

        loc = locator->second;

    }
    else
    {

        if (static_cast<uint>(i) >= arrays2support.size())
            throw std::range_error("This type of array has not been included in the model.");

        loc = arrays2support[i];

    }

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Counting stats_target
    StatsCounter< Array_Type, Data_Counter_Type> tmpstats(&Array_);

    tmpstats.set_counters(this->counters);
    
    std::vector< double > target_ = tmpstats.count_all();

    if (transform_model_fun)
        target_ = transform_model_fun(&target_[0u], target_.size());

    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) )
    {
        
        first_calc_done[loc] = true;

        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;
        
        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const std::vector<double> & target_,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    uint loc = arrays2support[i];

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;

        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const double * target_,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    uint loc = arrays2support[i];

    // Checking if passes the rules
    if (support_fun.get_rules_dyn()->size() > 0u)
    {

        std::vector< double > tmp_target(nterms(), 0.0);
        for (size_t t = 0u; t < nterms(); ++t)
            tmp_target[t] = *(target_ + t);

        if (!support_fun.eval_rules_dyn(tmp_target, 0u, 0u))
            return as_log ? -std::numeric_limits<double>::infinity() : 0.0;

    }

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;

        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood_total)(
    const std::vector<double> & params,
    bool as_log
) {
    
    size_t params_last_size = params_last.size();

    for (uint i = 0u; i < params_last_size; ++i)
    {

        if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) )
        {

            size_t k = params.size() + 1u;
            size_t n = stats_support[i].size() / k;
            
            first_calc_done[i] = true;
            normalizing_constants[i] = update_normalizing_constant(
                &params[0u], &stats_support[i][0u], k, n
            );
            
            params_last[i] = params;
            
        }

    }
    
    double res = 0.0;
    if (as_log)
    {

        for (uint i = 0; i < stats_target.size(); ++i) 
            res += vec_inner_prod(
                &stats_target[i][0u],
                &params[0u],
                params.size()
                ) BARRY_SAFE_EXP;
        
        #ifdef __OPENM 
        #pragma omp simd reduction(-:res)
        #else
        #pragma GCC ivdep
        #endif
        for (unsigned int i = 0u; i < params_last_size; ++i)
            res -= (std::log(normalizing_constants[i]) * this->stats_support_n_arrays[i]);

    } else {
        
        res = 1.0;
        size_t stats_target_size = stats_target.size();
        #ifdef __OPENM 
        #pragma omp simd reduction(*:res)
        #else
        #pragma GCC ivdep
        #endif
        for (unsigned int i = 0; i < stats_target_size; ++i)
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

MODEL_TEMPLATE(double, get_norm_const)(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    const auto id = arrays2support[i];
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[id] || !vec_equal_approx(params, params_last[id]) )
    {
        
        first_calc_done[id] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[id].size() / k;

        normalizing_constants[id] = update_normalizing_constant(
            &params[0u], &stats_support[id][0u], k, n
        );
        
        params_last[id] = params;
        
    }
    
    return as_log ? 
        std::log(normalizing_constants[id]) :
        normalizing_constants[id]
        ;
    
}

MODEL_TEMPLATE(const std::vector< Array_Type > *, get_pset)(
    const uint & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");


    return &pset_arrays[arrays2support[i]];

}

MODEL_TEMPLATE(const std::vector< double > *, get_pset_stats)(
    const uint & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    return &pset_stats[arrays2support[i]];

}

MODEL_TEMPLATE(void, print_stats)(uint i) const
{
    
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    const auto & S = stats_support[arrays2support[i]];

    size_t k       = nterms();
    size_t nunique = S.size() / (k + 1u);

    for (uint l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5i ", l);

        printf_barry("counts: %.0f motif: ", S[l * (k + 1u)]);
        
        for (unsigned int j = 0u; j < k; ++j)
            printf_barry("%.2f, ", S[l * (k + 1) + j + 1]);

        printf_barry("\n");

    }
    
    return;
    
}

MODEL_TEMPLATE(void, print)() const
{

    // Relevant information:
    // - Number of arrays involved
    // - Size of the support
    // - Terms involved

    uint min_v = std::numeric_limits<uint>::infinity();
    uint max_v = 0u;

    for (const auto & stat : this->stats_support)
    {
        if (stat.size() > max_v)
            max_v = stat.size();
        
        if (stat.size() < min_v)
            min_v = stat.size();

    }  

    printf_barry("Num. of Arrays     : %i\n", this->size());
    printf_barry("Support size       : %i\n", this->size_unique());
    printf_barry("Support size range : [%i, %i]\n", min_v, max_v);
    printf_barry("Transform. Fun.    : %s\n", transform_model_fun ? "yes": "no");
    printf_barry("Model terms (%i)   :\n", this->nterms());
    
    for (auto & cn : this->colnames())
        printf_barry(" - %s\n", cn.c_str());

    return;

}

MODEL_TEMPLATE(uint, size)() const noexcept
{
    // INITIALIZED()
    return this->stats_target.size();

}

MODEL_TEMPLATE(uint, size_unique)() const noexcept
{

    // INITIALIZED()
    return this->stats_support.size();

} 

MODEL_TEMPLATE(uint, nterms)() const noexcept
{
 
    if (transform_model_fun)
        return transform_model_term_names.size();
    else
        return this->counters->size();

}

MODEL_TEMPLATE(uint, support_size)() const noexcept
{

    // INITIALIZED()
    uint tot = 0u;
    for (auto& a : stats_support)
        tot += a.size();

    return tot;

}

MODEL_TEMPLATE(std::vector< std::string >, colnames)() const
{
    
    if (transform_model_fun)
        return transform_model_term_names;
    else
        return counters->get_names();

}
    
MODEL_TEMPLATE(Array_Type, sample)(
    const unsigned int & i,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // Getting the index
    unsigned int a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    size_t k = params.size();

    // Sampling an array
    unsigned int j = 0u;
    std::vector< double > & probs = pset_probs[a];
    if ((probs.size() > 0u) && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += probs[j++];

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

        j = i_matches;
        
    }
    

    return this->pset_arrays[a][j];   

}

MODEL_TEMPLATE(double, conditional_prob)(
    const Array_Type & Array_,
    const std::vector< double > & params,
    unsigned int i,
    unsigned int j
) {

    // Generating a copy of the array so we can update
    Array_Type A(Array_, true);

    // Making sure we add it first
    A.insert_cell(i, j, A.default_val(), true, false);

    // Computing the change stats_target
    std::vector< double > tmp_counts(counters->size());
    for (unsigned int ii = 0u; ii < tmp_counts.size(); ++ii)
        tmp_counts[ii] = counters->operator[](ii).count(A, i, j);

    // If there is a transformation function, it needs to be
    // applied before dealing with the likelihood.
    if (transform_model_fun)
        tmp_counts = transform_model_fun(&tmp_counts[0u], tmp_counts.size());

    return 1.0/
        (1.0 + std::exp(-vec_inner_prod<double>(
            &params[0u], &tmp_counts[0u], params.size()
            )));

    
}

MODEL_TEMPLATE(const std::mt19937 *, get_rengine)() const {
    return this->rengine;
}

template MODEL_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Counter_Type> * MODEL_TYPE()::get_counters() {
    return this->counters;
}

template MODEL_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Type> * MODEL_TYPE()::get_rules() {
    return this->rules;
}

template MODEL_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Dyn_Type> * MODEL_TYPE()::get_rules_dyn() {
    return this->rules_dyn;
}

template MODEL_TEMPLATE_ARGS()
inline Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> *
MODEL_TYPE()::get_support_fun() {
    return &this->support_fun;
}

MODEL_TEMPLATE(std::vector< std::vector< double > > *, get_stats_target)()
{
    return &stats_target;
}

MODEL_TEMPLATE(std::vector< std::vector< double > > *, get_stats_support)()
{
    return &stats_support;
}

MODEL_TEMPLATE(std::vector< unsigned int > *, get_arrays2support)()
{
    return &arrays2support;
}

MODEL_TEMPLATE(std::vector< std::vector< Array_Type > > *, get_pset_arrays)() {
    return &pset_arrays;
}

MODEL_TEMPLATE(std::vector< std::vector<double> > *, get_pset_stats)() {
    return &pset_stats;
}

MODEL_TEMPLATE(std::vector< std::vector<double> > *, get_pset_probs)() {
    return &pset_probs;
}

MODEL_TEMPLATE(void, set_transform_model)(
    std::function<std::vector<double>(double *,unsigned int)> fun,
    std::vector< std::string > names
    )
{

    if (transform_model_fun)
        throw std::logic_error("A transformation function for the model has already been established.");
    
    transform_model_fun = fun;
    transform_model_term_names = names;

    size_t k = counters->size(); 

    // Applying over the support
    for (auto & s : stats_support)
    {

        // Making room for the new support
        std::vector< double > s_new(0u);
        s_new.reserve(s.size());

        size_t n = s.size() / (k + 1u);

        // Iterating through the unique sets
        for (size_t i = 0; i < n; ++i)
        {

            // Appending size
            s_new.push_back(s[i * (k + 1u)]);

            // Applying transformation and adding to the new set
            auto res = transform_model_fun(&s[i * (k + 1u) + 1u], k);

            if (res.size() != transform_model_term_names.size())
                throw std::length_error("The transform vector from -transform_model_fun- does not match the size of -transform_model_term_names-.");

            std::copy(res.begin(), res.end(), std::back_inserter(s_new));

        }

        // Exchanging with the original
        std::swap(s, s_new);

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
            std::vector< double > new_stats(0u);

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
