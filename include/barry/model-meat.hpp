// #include "model-bones.hpp"

#ifndef BARRY_MODEL_MEAT_HPP 
#define BARRY_MODEL_MEAT_HPP 1

/**
 * @defgroup stat-models Statistical Models
 * @brief Statistical models available in `barry`.
 */

inline double update_normalizing_constant(
    const std::vector< double > & params,
    const std::vector< double > & support
)
{
    
    double res = 0.0;
    size_t k   = params.size();
    size_t n   = support.size() / (k + 1u);

    double tmp;
    for (unsigned int i = 0u; i < n; ++i)
    {

        tmp = 0.0;
        
        for (unsigned int j = 0u; j < params.size(); ++j)
            tmp += support[i * (k + 1u) + j + 1u] * params[j];
        
        res += exp(tmp BARRY_SAFE_EXP) * support[i * (k + 1u)];

    }
    
    // This will only evaluate if the option BARRY_CHECK_FINITE
    // is defined
    BARRY_ISFINITE(res)

    return res;
    
}

inline double likelihood_(
        const std::vector< double > & target_stats,
        const std::vector< double > & params,
        const double normalizing_constant,
        bool log_ = false
) {
    
    if (target_stats.size() != params.size())
        throw std::length_error("-target_stats- and -params- should have the same length.");
        
    double numerator = 0.0;
    
    // Computing the numerator
    for (unsigned int j = 0u; j < target_stats.size(); ++j)
        numerator += target_stats[j] * params[j];

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
    stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
    target_stats(0u), arrays2support(0u), keys2support(0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true)
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
    stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
    target_stats(0u), arrays2support(0u), keys2support(0u), 
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true)
{
    
    target_stats.reserve(size_);
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
    stats(Model_.stats),
    n_arrays_per_stats(Model_.n_arrays_per_stats),
    pset_arrays(Model_.pset_arrays),
    pset_stats(Model_.pset_stats),
    target_stats(Model_.target_stats),
    arrays2support(Model_.arrays2support),
    keys2support(Model_.keys2support),
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
    delete_rules_dyn(true)
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
        
        stats                 = Model_.stats;
        n_arrays_per_stats    = Model_.n_arrays_per_stats;
        pset_arrays           = Model_.pset_arrays;
        pset_stats            = Model_.pset_stats;
        target_stats          = Model_.target_stats;
        arrays2support        = Model_.arrays2support;
        keys2support          = Model_.keys2support;
        counters              = new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters));
        rules                 = new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules));
        rules_dyn             = new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn));
        delete_counters       = true;
        delete_rules          = true;
        delete_rules_dyn      = true;
        params_last           = Model_.params_last;
        normalizing_constants = Model_.normalizing_constants;
        first_calc_done       = Model_.first_calc_done;

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
    
    counters->add_counter(counter);
    return;
}

MODEL_TEMPLATE(void, add_counter)(
        Counter<Array_Type, Data_Counter_Type> * counter
) {
    
    counters->add_counter(counter);
    return;
    
}

MODEL_TEMPLATE(void, add_counter)(
    Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_,
    Data_Counter_Type *                            data_,
    bool                                           delete_data_
) {
    
    counters->add_counter(
        count_fun_,
        init_fun_,
        data_,
        delete_data_
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
    
    rules->add_rule(rules);
    return;
}

MODEL_TEMPLATE(void, add_rule)(
    Rule<Array_Type, Data_Rule_Type> * rule
) {
    
    rules->add_rule(rule);
    return;
    
}

MODEL_TEMPLATE(void, add_rule)(
    Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
    Data_Rule_Type *                         data_,
    bool                                     delete_data_
) {
    
    rules->add_rule(
        rule_fun_,
        data_,
        delete_data_
    );
    
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
    
    rules_dyn->add_rule(rules_);
    return;
}

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule<Array_Type, Data_Rule_Dyn_Type> * rules_
) {
    
    rules_dyn->add_rule(rules_);
    return;
    
}

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule_fun_type<Array_Type,Data_Rule_Dyn_Type> rule_fun_,
    Data_Rule_Dyn_Type *                         data_,
    bool                                     delete_data_
) {
    
    rules_dyn->add_rule(
        rule_fun_,
        data_,
        delete_data_
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
    target_stats.push_back(counter_fun.count_all());
    
    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = keygen(Array_);
    MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
    if (force_new | (locator == keys2support.end()))
    {
        
        // Adding to the map
        keys2support[key] = stats.size();
        n_arrays_per_stats.push_back(1u);       // How many elements now
        arrays2support.push_back(stats.size()); // Map of the array id to the support
        
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
            
            try
            {

                support_fun.calc();

            }
            catch (const std::exception& e)
            {
                
                printf_barry("A problem ocurred while trying to add the array. ");
                printf_barry("with error: %s", e.what());
                throw std::logic_error("");
                
            }
            
        }
        
        stats.push_back(support_fun.get_counts());
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(target_stats[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);
        
        return arrays2support.size() - 1u;
        
    }
    
    // Increasing the number of arrays in that stat
    ++n_arrays_per_stats[locator->second];
    
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
    if (this->stats[idx].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[idx] || !vec_equal_approx(params, params_last[idx]) )
    {
        
        first_calc_done[idx] = true;
        
        normalizing_constants[idx] = update_normalizing_constant(
            params, stats[idx]
        );
        
        params_last[idx] = params;
        
    }
    
    return likelihood_(
        target_stats[i],
        params,
        normalizing_constants[idx],
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
    if (this->stats[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Counting stats
    StatsCounter< Array_Type, Data_Counter_Type> tmpstats(&Array_);

    tmpstats.set_counters(this->counters);
    
    std::vector< double > target_ = tmpstats.count_all();

    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) )
    {
        
        first_calc_done[loc] = true;
        
        normalizing_constants[loc] = update_normalizing_constant(
            params, stats[loc]
        );
        
        params_last[loc] = params;
        
    }

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    return likelihood_(
        target_,
        params,
        normalizing_constants[loc],
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
    if (this->stats[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) ) {
        
        first_calc_done[loc] = true;
        
        normalizing_constants[loc] = update_normalizing_constant(
            params, stats[loc]
        );
        
        params_last[loc] = params;
        
    }
    
    return likelihood_(
        target_,
        params,
        normalizing_constants[loc],
        as_log
    );
    
}

MODEL_TEMPLATE(double, likelihood_total)(
    const std::vector<double> & params,
    bool as_log
) {
    
    for (uint i = 0u; i < params_last.size(); ++i)
    {

        if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) )
        {
            
            first_calc_done[i] = true;
            normalizing_constants[i] = update_normalizing_constant(
                params, stats[i]
            );
            
            params_last[i] = params;
            
        }

    }
    
    double res = 0.0;
    if (as_log)
    {

        for (uint i = 0; i < target_stats.size(); ++i) 
            res += vec_inner_prod(target_stats[i], params) BARRY_SAFE_EXP;
        
        for (uint i = 0u; i < params_last.size(); ++i)
            res -= (std::log(normalizing_constants[i]) * this->n_arrays_per_stats[i]);

    } else {
        
        res = 1.0;
        for (uint i = 0; i < target_stats.size(); ++i)
            res *= std::exp(vec_inner_prod(target_stats[i], params) BARRY_SAFE_EXP) / 
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
        
        normalizing_constants[id] = update_normalizing_constant(
                params, stats[id]
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

MODEL_TEMPLATE(const std::vector< std::vector< double > > *, get_pset_stats)(
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

    const auto & S = stats[arrays2support[i]];

    size_t k       = params_last.size();
    size_t nunique = S.size() / (k + 1u);

    for (uint l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5i ", l);

        printf_barry("counts: %i motif: ", S[l * (k + 1)]);
        
        for (unsigned int j = 0u; j < k; ++j)
            printf_barry("%.2f, ", S[l * (k + 1) + k + 1]);

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

    for (const auto & stat : this->stats)
    {
        if (stat.size() > max_v)
            max_v = stat.size();
        
        if (stat.size() < min_v)
            min_v = stat.size();

    }  

    printf_barry("Num. of Arrays     : %i\n", this->size());
    printf_barry("Support size       : %i\n", this->size_unique());
    printf_barry("Support size range : [%i, %i]\n", min_v, max_v);
    printf_barry("Model terms (%i)   :\n", this->nterms());
        
    for (auto & cn : this->colnames())
        printf_barry(" - %s\n", cn.c_str());

    return;

}

MODEL_TEMPLATE(uint, size)() const noexcept
{
    // INITIALIZED()
    return this->target_stats.size();

}

MODEL_TEMPLATE(uint, size_unique)() const noexcept
{

    // INITIALIZED()
    return this->stats.size();

} 

MODEL_TEMPLATE(uint, nterms)() const noexcept
{

    return this->counters->size();

}

MODEL_TEMPLATE(uint, support_size)() const noexcept
{

    // INITIALIZED()
    uint tot = 0u;
    for (auto& a : stats)
        tot += a.size();

    return tot;

}

MODEL_TEMPLATE(std::vector< std::string >, colnames)() const
{
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
    // if (pset_probs.at(a).size() == 0u) {
    //   pset_probs.at(a).resize(pset_arrays.at(a).size(), 0u);
    // }
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    // Updating until reach above
    unsigned int j = 0u;
    while (cumprob < r) {

        cumprob += this->likelihood(params, this->pset_stats[a][j], i, false);
        ++j;
    }

    return this->pset_arrays[a][j-1u];   

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

    // Computing the change stats
    std::vector< double > tmp_counts(counters->size());
    for (unsigned int ii = 0u; ii < tmp_counts.size(); ++ii)
        tmp_counts[ii] = counters->operator[](ii).count(A, i, j);

    return 1.0/
        (1.0 + std::exp(-vec_inner_prod<double>(params, tmp_counts)));

    
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
MODEL_TYPE()::get_support() {
    return &this->support_fun;
}

#undef MODEL_TEMPLATE
#undef MODEL_TEMPLATE_ARGS
#undef MODEL_TYPE

#endif
