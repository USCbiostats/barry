#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

#include "formula.hpp"

/**
 * @ingroup counting 
 * @details Details on the available counters for `DEFMworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{

///@}

/**
 * @brief Factory to create a hasher function for DEFMArray
 * @param covar_index If >= 0, then the hasher will include the
 * covariate at that index as part of the hash.
 * @return A hasher function for DEFMArray
 */
inline barry::Hasher_fun_type<DEFMArray, DEFMCounterData> 
defm_hasher_factory(int covar_index = -1) {
    
    // With no covariate index, we skip adding that
    // layer to the hasher
    if (covar_index >= 0)
    {

        return [covar_index](
            const DEFMArray & array,
            DEFMCounterData * d
        ) -> std::vector< double > {

            std::vector< double > res;

            // Adding the column feature
            for (size_t i = 0u; i < array.nrow(); ++i)
            {
                res.push_back(array.D()(i, covar_index));
            }

            // Adding the fixed dims
            for (size_t i = 0u; i < (array.nrow() - 1); ++i)
            {
                for (size_t j = 0u; j < array.ncol(); ++j)
                {
                    res.push_back(array(i, j));
                }
            }

            return res;
        };

    } else {
        
        return [](
            const DEFMArray & array,
            DEFMCounterData * d
        ) -> std::vector< double > {

            std::vector< double > res;

            // Adding the fixed dims
            for (size_t i = 0u; i < (array.nrow() - 1); ++i)
            {
                for (size_t j = 0u; j < array.ncol(); ++j)
                {
                    res.push_back(array(i, j));
                }
            }

            return res;
        };

    }

}

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_COUNTER(a) \
inline double (a) (const DEFMArray & Array, size_t i, size_t j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a) \
barry::Counter_fun_type<DEFMArray, DEFMCounterData> a = \
    [](const DEFMArray & Array, size_t i, size_t j, DEFMCounterData & data) -> double

///@}

/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_RULE(a) \
inline bool (a) (const DEFMArray & Array, size_t i, size_t j, bool & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_RULE_LAMBDA(a) \
barry::Rule_fun_type<DEFMArray, DEFMRuleData> a = \
[](const DEFMArray & Array, size_t i, size_t j, DEFMRuleData & data) -> bool
///@}

/**Lambda function for definition of a network counter function*/
#define DEFM_RULEDYN_LAMBDA(a) \
barry::Rule_fun_type<DEFMArray, DEFMRuleDynData> a = \
[](const DEFMArray & Array, size_t i, size_t j, DEFMRuleDynData & data) -> bool
///@}

/**
  * @weakgroup  counters-network DEFMArray counters
  * @brief Counters for network models
  * @param counters A pointer to a `DEFMCounters` object (`Counters`<`DEFMArray`, `DEFMCounterData`>).
  */
///@{
// -----------------------------------------------------------------------------
/**
 * @brief Prevalence of ones
 * 
 * This is similar to the `edges` term in ERGMs.
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_ones(
    DEFMCounters * counters,
    int covar_index   = -1,
    const std::vector< std::string > * x_names = nullptr
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {   

        auto hasher = defm_hasher_factory(covar_index);

        DEFM_COUNTER_LAMBDA(counter_tmp)
        {

            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return Array.D()(i, data.idx(0u));

        };

        std::string vname;
        if (x_names != nullptr)
            vname = x_names->operator[](covar_index);
        else
            vname = std::string("attr")+ std::to_string(covar_index);

        counters->add_counter(
            counter_tmp, nullptr, hasher,
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}, true), 
            "Num. of ones x " + vname, 
            "Overall number of ones" + (
                covar_index >= 0 ? (" weighted by " + vname) : std::string("")
            )
        );

    }
    else
    {

        DEFM_COUNTER_LAMBDA(count_ones)
        {
            
            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return 1.0;
        };

        DEFMCounterData dat;
        dat.is_motif = true;

        counters->add_counter(
            count_ones, nullptr, nullptr,
            dat, // DEFMCounterData(),
            "Num. of ones", 
            "Overall number of ones"
        );
    }

    return;

}


/**
 * Calculates the logit intercept for the DEFM model.
 *
 * @param counters A pointer to the DEFMCounters object.
 * @param n_y The number of response variables.
 * @param which A vector of indices indicating which response variables to use. If empty, all response variables are used.
 * @param covar_index The index of the covariate to use as the intercept. 
 * @param x_names A pointer to a vector of strings containing the names of the covariates.
 * @param y_names A pointer to a vector of strings containing the names of the response variables.
 */
inline void counter_logit_intercept(
    DEFMCounters * counters,
    size_t n_y,
    std::vector< size_t > which = {},
    int covar_index = -1,
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
) {


    if (which.size() == 0u)
    {
        which.resize(n_y, 0u);
        std::iota(which.begin(), which.end(), 0u);
    } else {
        for (auto w : which)
            if (w >= n_y)
                throw std::logic_error("Values in `which` are out of range.");
    }

    // Case when no interaction happens, whatsoever.
    std::string vname;
    if (covar_index < 0)
    {

        DEFM_COUNTER_LAMBDA(tmp_counter)
        {
            if (i != (Array.nrow() - 1))
                return 0.0;

            if (j != data.idx(0u))
                return 0.0;

            return 1.0;
        };

        for (auto i : which)
        {

            if (y_names != nullptr)
                vname = y_names->operator[](i);
            else
                vname = std::to_string(i);

            counters->add_counter(
                tmp_counter, nullptr, nullptr,
                DEFMCounterData({i}, {}, {}, false), 
                "Logit intercept " + vname, 
                "Equal to one if the outcome " + vname + " is one. Equivalent to the logistic regression intercept."
            );

        }

    } else {

        DEFM_COUNTER_LAMBDA(tmp_counter)
        {
            if (i != Array.nrow() - 1)
                return 0.0;

            if (j != data.idx(0u))
                return 0.0;

            return Array.D()(i, data.idx(1u));
        };

        auto hasher = defm_hasher_factory(covar_index);
        bool hasher_added = false;

        std::string yname, vname;
        for (auto i : which)
        {

            if (y_names != nullptr)
            {
                yname = y_names->operator[](i);
            }
            else
            {
                yname = std::to_string(i);
            }

            if (x_names != nullptr)
                vname = x_names->operator[](covar_index);
            else
                vname = std::string("attr")+ std::to_string(covar_index);

            if (hasher_added)
                counters->add_counter(
                    tmp_counter, nullptr, nullptr,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}, false), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );
            else {

                hasher_added = true;

                counters->add_counter(
                    tmp_counter, nullptr, hasher,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}, false), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );

            }

        }

    }
    

}

/**
 * @brief Generalized counter for DEFM models
 * 
 * @param counters Pointer ot a vector of counters
 * @param signs A vector of signs for the motif variables. If empty, all
 * are assumed to be positive. If the size does not match the size of
 * -coords-, an error is thrown.
 * @param coords A vector of coordinates for the motif variables.
 * Each coordinate must be between 0 and (m_order + 1) * n_y - 1.
 * The coordinates are specified in column-major order.
 * @param m_order The Markov order of the data.
 * @param n_y The number of response variables.
 * @param covar_index If >= than 0, then the interaction with the covariate
 * at that index will be included.
 * @param x_names A pointer to a vector of strings containing the names of
 * the covariates.
 * @param y_names A pointer to a vector of strings containing the names of
 * the response variables.
 * 
 * @details
 * This function adds a counter to the DEFM model that can compute 
 * either a motif-based count (combinations) or transitions. These
 * can also be interacted with a covariate.
 *
 * If either `x_names` or `y_names` is nullptr, then generic names will be used.
 * If both are provided, then the names will be used in the counter description.
 * If only one is provided, then only that name will be used.
 */
inline void counter_generic(
    DEFMCounters * counters,
    std::vector< size_t > coords,
    std::vector< bool > signs,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
)
{

    // A vector to store the type of dat
    if (signs.size() == 0u)
        signs.resize(coords.size(), true);
    else if (signs.size() != coords.size())
        throw std::length_error("Size of -coords- and -signs- must match.");

    if (covar_index >= 0)
        coords.push_back(static_cast<size_t>(covar_index));
    else
        coords.push_back(1000u);

    DEFM_COUNTER_LAMBDA(count_init)
    {

        auto indices = data.indices;
        auto sgn     = data.logical;
        int covaridx = indices[indices.size() - 1u];

        // Notice that the indices vector contains:
        // - 1st, the indices of the motif. That's why we set the lenght
        //   using -1.
        // - the last is, the covariate index
        for (size_t k = 0u; k < (indices.size() - 1u); ++k)
        {
            if (indices[k] >= (Array.ncol()* Array.nrow()))
                throw std::range_error("The motif includes entries out of range.");
        }

        // Counting
        const auto & array = Array.get_data();
        for (size_t k = 0u; k < (indices.size() - 1); ++k)
        {
            auto cellv = array[indices[k]];
            if (sgn[k] && (cellv != 1))
                return 0.0;
        }
            
        // If nothing happens, then is one or the covaridx
        return (covaridx < 1000) ? Array.D()(Array.nrow() - 1u, covaridx) : 1.0;
        
    };

    DEFM_COUNTER_LAMBDA(count_ones)
    {
        
        auto dat = data.indices;
        auto sgn = data.logical;
        int covaridx = dat[dat.size() - 1u];

        // Checking if the observation is in the stat. We 
        const auto & array = Array.get_data();
        size_t loc = i + j * Array.nrow();
        size_t n_cells = dat.size() - 1u;

        // Only one currently needs to be a zero for it
        // to change
        size_t n_now = 0;
        bool baseline_value = false;
        bool i_in_array = false;
        for (size_t e = 0u; e < n_cells; ++e)
        {

            // Is the current cell in the list?
            if (dat[e] == loc)
            {
                i_in_array = true;
                baseline_value = sgn[e];
            }

            if ((sgn[e] && (array[dat[e]] == 1)) || (!sgn[e] && (array[dat[e]] == 0)))
                n_now++;
            
        }

        // If i in array still false, then no change
        if (!i_in_array)
            return 0.0;
        
        size_t n_prev = n_now;
        if (baseline_value)
            n_prev--;
        else
            n_prev++;

        // Computing stats
        if (covaridx < 1000)
        {
            
            double val = Array.D()(Array.nrow() - 1u, covaridx);
            double value_now  = n_now == n_cells ?  val : 0.0;
            double value_prev = n_prev == n_cells ? val : 0.0;

            return value_now - value_prev;

        } 
        else
        {

            double value_now  = n_now == n_cells ? 1.0 : 0.0;
            double value_prev = n_prev == n_cells ? 1.0 : 0.0;

            return value_now - value_prev;

        }

    };

    // Creating name of the structure
    std::string name, vname;
    if (coords.size() == 1u)
        name = "";
    else
        name = "Motif ";

    // Creating an empty motif filled with zeros
    barry::BArrayDense<int> motif(m_order + 1u, n_y, 0);

    // Filling the matrix in, negative values are 0s and 1s are... 1s.
    // Zero are values not used.
    size_t n_cells = coords.size() - 1u;
    for (size_t d = 0u; d < n_cells; ++d)
    {
        size_t c = std::floor(coords[d] / (m_order + 1u));
        size_t r = coords[d] - c * (m_order + 1u);
        motif(r, c) = signs[d] ? 1 : -1;
        
    }

    #ifdef BARRY_WITH_LATEX
        name += "$";
    #endif

    // Generate name by showing each time point separately
    // Loop through all time points (0 to m_order)
    for (size_t i = 0u; i <= m_order; ++i)
    {
        // Check if this time point has any variables
        bool has_vars = false;
        for (size_t j = 0u; j < n_y; ++j)
        {
            if (motif(i, j) != 0)
            {
                has_vars = true;
                break;
            }
        }

        if (!has_vars)
            continue;

        // Open bracket for this time point
        #ifdef BARRY_WITH_LATEX
            name += "(";
        #else
            name += "{";
        #endif

        // Add variables for this time point
        bool row_start = true;
        for (size_t j = 0u; j < n_y; ++j)
        {
            // Is it included?
            if (motif(i, j) == 0)
                continue;

            // Is not the first?
            if (row_start)
                row_start = false;
            else
                name += ", ";

            if (y_names != nullptr)
                name += y_names->operator[](j);
            else
                name += (std::string("y") + std::to_string(j));

            #ifdef BARRY_WITH_LATEX
                name += (motif(i, j) < 0 ? "^-" : "^+");
            #else
                name += (motif(i, j) < 0 ? "-" : "+");
            #endif
        }

        // Close bracket for this time point
        #ifdef BARRY_WITH_LATEX
            name += ")";
        #else
            name += "}";
        #endif

        // Add arrow if not the last time point
        if (i < m_order)
        {
            name += ">";
        }
    }

    #ifdef BARRY_WITH_LATEX
    name += "$";
    #endif

    if (covar_index >= 0)
    {

        auto hasher = defm_hasher_factory(covar_index);

        if (x_names != nullptr)
            vname = x_names->operator[](covar_index);
        else
            vname = std::string("attr")+ std::to_string(covar_index);

        counters->add_counter(
            count_ones, count_init, hasher,
            DEFMCounterData(coords, {}, signs, coords.size() > 1u ? true : false), 
            name + " x " + vname, 
            "Motif weighted by " + vname
        );

    } else {

        counters->add_counter(
            count_ones, count_init, nullptr,
            DEFMCounterData(coords, {}, signs, coords.size() > 1u ? true : false), 
            name, 
            "Motif"
        );

    }
    

    return;

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_formula(
    DEFMCounters * counters,
    std::string formula,
    size_t m_order,
    size_t n_y,
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
) {

    std::vector< size_t > coords;
    std::vector< bool > signs;
    std::string covar_name = "";

    defm_motif_parser(
        formula, coords, signs, m_order, n_y, covar_name
    );

    int covar_index = -1;
    if (covar_name != "")
    {

        if (x_names != nullptr)
        {
            for (size_t i = 0u; i < x_names->size(); ++i)
                if (x_names->operator[](i) == covar_name)
                {
                    covar_index = static_cast<int>(i);
                    break;
                }
        }

        if (covar_index < 0)
            throw std::logic_error(
                std::string("The covariate name '") +
                covar_name +
                std::string("' was not found in the list of covariates.")
                );

    }

    // Checking the number of coords, could be single intercept
    if (coords.size() == 1u)
    {

        // Getting the column
        size_t coord = static_cast< size_t >(
            std::floor(
            static_cast<double>(coords[0u]) / static_cast<double>(m_order + 1)
            ));

        counter_logit_intercept(
            counters, n_y, {coord},
            covar_index,
            x_names,
            y_names
        );

    }
    else 
    {

        counter_generic(
            counters, coords, signs, m_order, n_y, covar_index,
            x_names, y_names
        );

    }


}

/**
 * @name Returns true if the cell is free
 * @param rules A pointer to a `DEFMRules` object (`Rules`<`DEFMArray`, `bool`>).
 */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
inline void rules_markov_fixed(
    DEFMRules * rules,
    size_t markov_order
    ) {
    
    DEFM_RULE_LAMBDA(no_self_tie) {
        return i >= data.idx(0u);
    };
    
    rules->add_rule(
        no_self_tie,
        DEFMRuleData({},{markov_order}),
        std::string("Markov model of order ") + std::to_string(markov_order),
        std::string("Blocks the first m-order cells of the array.")
        );
    
    return;
}

/**
 * @brief Blocks switching a one to zero.
 * 
 * @param rules 
 * @param term_indices Ids of the variables that will follow this rule.
 */
inline void rules_dont_become_zero(
    DEFMSupport * support,
    std::vector<size_t> term_indices
    ) {
    
    DEFM_RULE_LAMBDA(rule) {

        if (!data.init)
        {
            std::vector< size_t > tmp(Array.ncol(), 0u);

            for (auto v : data.indices)
            {
                if (v >= Array.ncol())
                    throw std::range_error("The specified id for `dont_become_zero` is out of range.");

                tmp[v] = 1u;
            }

            data.indices.resize(Array.ncol());
            for (size_t v = 0u; v < tmp.size(); ++v)
                data.indices[v] = tmp[v];

            data.init = true;
        }

        // If not considered, then continue
        if (data.indices[j] == 0u)
            return true;

        // The data outside of the markov chain is checked by other rule
        if (i != (Array.nrow() - 1))
            return true;

        // This is now one, is the next different zero? If so,
        // we can include it (1->1)
        return (Array(i - 1, j) != 1) || (Array(i, j) != 1);

    };
    
    support->get_rules()->add_rule(
        rule,
        DEFMRuleData({}, {term_indices}),
        std::string("Ones can't become zero"),
        std::string("Blocks cells that have became equal to one.")
        );
    
    return;
}

/**
 * @brief Overall functional gains
 * @param support Support of a model.
 * @param term_index Position of the focal statistic.
 * @param lb Lower bound
 * @param ub Upper bound
 * @details 
 * @return (void) adds a rule limiting the support of the model.
 */
inline void rule_constrain_support(
    DEFMSupport * support,
    size_t term_index,
    double lb,
    double ub
)
{
  
    DEFM_RULEDYN_LAMBDA(tmp_rule)
    {

        if (data() < data.lb)
            return false;
        else if (data() > data.ub)
            return false;
        else
            return true;
      
    };

    
    support->get_rules_dyn()->add_rule(
        tmp_rule,
        DEFMRuleDynData(
            support->get_current_stats(),
            term_index, lb, ub
            ),
        support->get_counters()->get_names()[term_index] +
            "' within [" + std::to_string(lb) + ", " +
            std::to_string(ub) + std::string("]"),
        std::string("When the support is enumerated, only states where the statistic '") +
            support->get_counters()->get_names()[term_index] +
            std::to_string(term_index) + "' falls within [" + std::to_string(lb) + ", " +
            std::to_string(ub) + "] are included."
    );
    
    return;
  
}


///@}

///@}

#endif
