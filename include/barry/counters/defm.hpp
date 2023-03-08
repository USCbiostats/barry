#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

#include "defm-formula.hpp"

/**
 * @ingroup counting 
 * @details Details on the available counters for `DEFMworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{

/**
 * @brief Data class for DEFM arrays.
 * 
 * This holds information pointing to the data array, including information
 * regarding the number of observations, the time slices of the observation,
 * and the number of covariates in the data.
 * 
 */

class DEFMData;

typedef BArrayDense<int, DEFMData> DEFMArray;

class DEFMData {
public:
    
    DEFMArray * array; // Pointer to the owner of this data
    const double * covariates; ///< Vector of covariates (complete vector)
    size_t obs_start;    ///< Index of the observation in the data.
    size_t X_ncol; ///< Number of columns in the array of covariates.
    size_t X_nrow; ///< Number of rows in the array of covariates.
    std::vector< size_t > covar_sort; /// Value where the sorting of the covariates is stored.
    std::vector< size_t > covar_used; /// Vector indicating which covariates are included in the model
    
    DEFMData() {};
    
    /**
     * @brief Constructor
     * @param covariates_ Pointer to the attribute data.
     * @param obs_start_ Location of the current observation in the covariates
     *  vector
     * @param X_ncol_ Number of columns (covariates.)
     */
    DEFMData(
        DEFMArray * array_,
        const double * covariates_,
        size_t obs_start_,
        size_t X_ncol_,
        size_t X_nrow_
    ) : array(array_), covariates(covariates_), obs_start(obs_start_),
    X_ncol(X_ncol_), X_nrow(X_nrow_) {}; 

    /**
     * @brief Access to the row (i) colum (j) data
     * 
     * @param i 
     * @param j 
     * @return double 
     */
    double operator()(size_t i, size_t j) const;
    double at(size_t i, size_t j) const;
    size_t ncol() const;
    size_t nrow() const;
    void print() const;
    
    ~DEFMData() {};

};

/**
  * @brief Data class used to store arbitrary uint or double vectors */
class DEFMCounterData {
public:
    
    std::vector< size_t > indices;
    std::vector< double > numbers;
    std::vector< bool >   logical;
    
    DEFMCounterData() : indices(0u), numbers(0u) {};
    DEFMCounterData(
        const std::vector< size_t > indices_,
        const std::vector< double > numbers_,
        const std::vector< bool > logical_
    ): indices(indices_), numbers(numbers_), 
        logical(logical_) {};

    size_t idx(size_t i) const {return indices[i];};
    double num(size_t i) const {return numbers[i];};
    bool is_true(size_t i) const {return logical[i];};
    
    ~DEFMCounterData() {};
    
};

class DEFMRuleData {
public:

    std::vector< double > numbers;
    std::vector< size_t > indices;
    std::vector< bool >   logical;

    bool init = false;

    double num(size_t i) const {return numbers[i];};
    size_t idx(size_t i) const {return indices[i];};
    bool is_true(size_t i) const {return logical[i];};

    DEFMRuleData() {};

    DEFMRuleData(
        std::vector< double > numbers_,
        std::vector< size_t > indices_,
        std::vector< bool > logical_
    ) : numbers(numbers_), indices(indices_), logical(logical_) {};

    DEFMRuleData(
        std::vector< double > numbers_,
        std::vector< size_t > indices_
    ) : numbers(numbers_), indices(indices_), logical(numbers_.size()) {};

};

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{

class DEFMRuleDynData : public DEFMRuleData {
public:
    const std::vector< double > * counts;
    
    DEFMRuleDynData(
        const std::vector< double > * counts_,
        std::vector< double > numbers_ = {},
        std::vector< size_t > indices_ = {},
        std::vector< bool > logical_ = {}
        ) : DEFMRuleData(numbers_, indices_, logical_), counts(counts_) {};
    
    ~DEFMRuleDynData() {};
    
};

/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef Counter<DEFMArray, DEFMCounterData > DEFMCounter;
typedef Counters<DEFMArray, DEFMCounterData> DEFMCounters;
typedef Support<DEFMArray, DEFMCounterData, DEFMRuleData,DEFMRuleDynData> DEFMSupport;
typedef StatsCounter<DEFMArray, DEFMCounterData> DEFMStatsCounter;
typedef Model<DEFMArray, DEFMCounterData,DEFMRuleData,DEFMRuleDynData> DEFMModel;


typedef Rule<DEFMArray, DEFMRuleData> DEFMRule;
typedef Rules<DEFMArray, DEFMRuleData> DEFMRules;
typedef Rule<DEFMArray, DEFMRuleDynData> DEFMRuleDyn;
typedef Rules<DEFMArray, DEFMRuleDynData> DEFMRulesDyn;



///@}

inline double DEFMData::operator()(size_t i, size_t j) const
{
    return *(covariates + (obs_start + j * X_nrow + i));
}

inline size_t DEFMData::ncol() const {
    return X_ncol;
}

inline size_t DEFMData::nrow() const {
    return X_nrow;
}

inline void DEFMData::print() const {

    for (size_t i = 0u; i < array->nrow(); ++i)
    {

        printf_barry("row %li (%li): ", i, obs_start + i);
        for (size_t j = 0u; j < X_ncol; ++j)
            printf_barry("% 5.2f, ", operator()(i, j));
        printf_barry("\n");
        
    }

}

#define MAKE_DEFM_HASHER(hasher,a,cov) Hasher_fun_type<DEFMArray,DEFMCounterData> hasher = [cov](const DEFMArray & array, DEFMCounterData * d) { \
            std::vector< double > res; \
            /* Adding the column feature */ \
            for (size_t i = 0u; i < array.nrow(); ++i) \
                res.push_back(array.D()(i, cov)); \
            /* Adding the fixed dims */ \
            for (size_t i = 0u; i < (array.nrow() - 1); ++i) \
                for (size_t j = 0u; j < array.ncol(); ++j) \
                    res.push_back(array(i, j)); \
            return res;\
        };
    

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_COUNTER(a) \
inline double (a) (const DEFMArray & Array, uint i, uint j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a) \
Counter_fun_type<DEFMArray, DEFMCounterData> a = \
    [](const DEFMArray & Array, uint i, uint j, DEFMCounterData & data) -> double

///@}

/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_RULE(a) \
inline bool (a) (const DEFMArray & Array, uint i, uint j, bool & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_RULE_LAMBDA(a) \
Rule_fun_type<DEFMArray, DEFMRuleData> a = \
[](const DEFMArray & Array, uint i, uint j, DEFMRuleData & data) -> bool
///@}

/**Lambda function for definition of a network counter function*/
#define DEFM_RULEDYN_LAMBDA(a) \
Rule_fun_type<DEFMArray, DEFMRuleDynData> a = \
[](const DEFMArray & Array, uint i, uint j, DEFMRuleDynData & data) -> bool
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
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_ones(
    DEFMCounters * counters,
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {   

        MAKE_DEFM_HASHER(hasher, array, covar_index)

        DEFM_COUNTER_LAMBDA(counter_tmp)
        {

            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return Array.D()(i, data.idx(0u));

        };


        if (vname == "")
        {
            if (x_names != nullptr)
                vname = x_names->operator[](covar_index);
            else
                vname = std::string("attr")+ std::to_string(covar_index);
        }

        counters->add_counter(
            counter_tmp, nullptr, hasher,
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}), 
            "Num. of ones x " + vname, 
            "Overall number of ones"
        );



    } else {

        DEFM_COUNTER_LAMBDA(count_ones)
        {
            
            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return 1.0;
        };

        counters->add_counter(
            count_ones, nullptr, nullptr,
            DEFMCounterData(),
            "Num. of ones", 
            "Overall number of ones"
        );
    }

    return;

}

inline void counter_logit_intercept(
    DEFMCounters * counters,
    size_t n_y,
    std::vector< size_t > which = {},
    int covar_index = -1,
    std::string vname = "",
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
                DEFMCounterData({i}, {}, {}), 
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

        MAKE_DEFM_HASHER(hasher, array, covar_index)
        bool hasher_added = false;

        std::string yname;
        for (auto i : which)
        {

            if (y_names != nullptr)
                yname = y_names->operator[](i);
            else
                yname = std::to_string(i);

            if (vname == "")
            {
                if (x_names != nullptr)
                    vname = x_names->operator[](covar_index);
                else
                    vname = std::string("attr")+ std::to_string(covar_index);
            }

            if (hasher_added)
                counters->add_counter(
                    tmp_counter, nullptr, nullptr,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );
            else {

                hasher_added = true;

                counters->add_counter(
                    tmp_counter, nullptr, hasher,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );

            }

        }

    }
    

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_transition(
    DEFMCounters * counters,
    std::vector< size_t > coords,
    std::vector< bool > signs,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    std::string vname = "",
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

        for (size_t i = 0u; i < (indices.size() - 1u); ++i)
        {
            if (
                std::floor(indices[i] / Array.nrow()) >= 
                static_cast<int>(Array.ncol())
                )
                throw std::range_error("The motif includes entries out of range.");
        }
            
        return 0.0;
        
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

            if ((sgn[e] & (array[dat[e]] == 1)) | (!sgn[e] & (array[dat[e]] == 0)))
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
    std::string name;
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

    // Checking if any prior to the event
    bool any_before_event = false;
    
    for (size_t i = 0u; i < m_order; ++i)
    {
        for (size_t j = 0u; j < n_y; ++j)
        {
            if (motif(i,j) != 0)
            {
                any_before_event = true;
                break;
            }

        }
    }
    
    #ifdef BARRY_WITH_LATEX
        name += "$";
    #endif

    if (any_before_event)
        #ifdef BARRY_WITH_LATEX
            name += "(";
        #else
            name += "{";
        #endif

    #ifdef BARRY_WITH_LATEX
        #define UNI_SUB(a) \
            (\
                ((a) == 0) ? "_0" : (\
                ((a) == 1) ? "_1" : (\
                ((a) == 2) ? "_2" : (\
                ((a) == 3) ? "_3" : (\
                ((a) == 4) ? "_4" : (\
                ((a) == 5) ? "_5" : (\
                ((a) == 6) ? "_6" : (\
                ((a) == 7) ? "_7" : (\
                ((a) == 8) ? "_8" : \
                "_9"))))))))\
            )
    #else
        #define UNI_SUB(a) \
            (\
                ((a) == 0) ? "\u2080" : (\
                ((a) == 1) ? "\u2081" : (\
                ((a) == 2) ? "\u2082" : (\
                ((a) == 3) ? "\u2083" : (\
                ((a) == 4) ? "\u2084" : (\
                ((a) == 5) ? "\u2085" : (\
                ((a) == 6) ? "\u2086" : (\
                ((a) == 7) ? "\u2087" : (\
                ((a) == 8) ? "\u2088" : \
                "\u2089"))))))))\
            )
    #endif

    // If order is greater than zero, the starting point of the transtion
    for (size_t i = 0u; i < m_order; ++i)
    {

        bool row_start = true;
        for (size_t j = 0u; j < n_y; ++j)
        {

            // Is it included?
            if (motif(i,j) == 0)
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
                name += (motif(i,j) < 0 ? "^-" : "^+");
            #else
                name += (motif(i,j) < 0 ? "\u207B" : "\u207A");
            #endif

        }

    }

    // If it has starting point, then need to close.
    if (any_before_event & (m_order > 0u))
        #ifdef BARRY_WITH_LATEX
            name += ") -> (";
        #else
            name += "} \u21E8 {";
        #endif
    else
        #ifdef BARRY_WITH_LATEX
            name += "(";
        #else
            name += "{";
        #endif

    // Looking onto the transtions
    bool row_start = true;
    for (size_t j = 0u; j < n_y; ++j)
    {

        if (motif(m_order, j) == 0)
            continue;

        if (row_start)
            row_start = false;
        else
            name += ", ";

        if (y_names != nullptr)
            name += y_names->operator[](j);
        else
            name += (std::string("y") + std::to_string(j));

        #ifdef BARRY_WITH_LATEX
        name += (motif(m_order, j) < 0 ? "^-" : "^+" );
        #else
        name += (motif(m_order, j) < 0 ? "\u207B" : "\u207A" );
        #endif


    }

    #undef UNI_SUB

    #ifdef BARRY_WITH_LATEX
    name += ")$";
    #else
    name += "}";
    #endif

    if (covar_index >= 0)
    {

        MAKE_DEFM_HASHER(hasher, array, covar_index)

        if (vname == "")
        {
            if (x_names != nullptr)
                vname = x_names->operator[](covar_index);
            else
                vname = std::string("attr")+ std::to_string(covar_index);
        }

        counters->add_counter(
            count_ones, count_init, hasher,
            DEFMCounterData(coords, {}, signs), 
            name + " x " + vname, 
            "Motif weighted by single attribute"
        );

    } else {

        counters->add_counter(
            count_ones, count_init, nullptr,
            DEFMCounterData(coords, {}, signs), 
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
inline void counter_transition_formula(
    DEFMCounters * counters,
    std::string formula,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
) {

    std::vector< size_t > coords;
    std::vector< bool > signs;

    defm_motif_parser(
        formula, coords, signs, m_order, n_y
    );

    counter_transition(
        counters, coords, signs, m_order, n_y, covar_index, vname,
        x_names, y_names
    );

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_fixed_effect(
    DEFMCounters * counters,
    int covar_index,
    double k,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr
)
{

    DEFM_COUNTER_LAMBDA(count_init)
    {
        return std::pow(Array.D()((size_t) i, data.idx(0u)), data.num(0u));
    };

    DEFM_COUNTER_LAMBDA(count_tmp)
    {
        return 0.0;
    };

    MAKE_DEFM_HASHER(hasher, array, covar_index)

    if (x_names != nullptr)
        vname = x_names->operator[](covar_index);
    else
        vname = std::string("attr")+ std::to_string(covar_index);

    counters->add_counter(
        count_tmp, count_init, hasher,
        DEFMCounterData({static_cast<size_t>(covar_index)}, {k}, {}), 
        "Fixed effect feature (" + vname + ")^" + std::to_string(k)
    );

    return;

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
        std::string("Blocks the first morder cells of the array.")
        );
    
    return;
}

/**
 * @brief Blocks switching a one to zero.
 * 
 * @param rules 
 * @param ids Ids of the variables that will follow this rule.
 */
inline void rules_dont_become_zero(
    DEFMSupport * support,
    std::vector<size_t> ids
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
        return (Array(i - 1, j) != 1) | (Array(i, j) != 1);

    };
    
    support->get_rules()->add_rule(
        rule,
        DEFMRuleData({}, {ids}),
        std::string("Ones can't become zero"),
        std::string("Blocks cells that have became equal to one.")
        );
    
    return;
}

///@}

///@}

#endif
