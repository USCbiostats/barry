#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

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
    size_t X_ncol; ///< Number of covariates included in the model.
    size_t X_nrow; ///< Number of covariates included in the model.
    
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

    size_t idx(size_t i) {return indices[i];};
    double num(size_t i) {return numbers[i];};
    bool is_true(size_t i) {return logical[i];};
    
    ~DEFMCounterData() {};
    
};

class DEFMRuleData {
private: 
    std::vector< double > numbers;
    std::vector< size_t > indices;

public:

    double num(size_t i) {return numbers[i];};
    size_t idx(size_t i) {return indices[i];};

    DEFMRuleData() {};

    DEFMRuleData(
        std::vector< double > numbers_,
        std::vector< size_t > indices_
    ) : numbers(numbers_), indices(indices_) {};

};

/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef Counter<DEFMArray, DEFMCounterData > DEFMCounter;
typedef Counters<DEFMArray, DEFMCounterData> DEFMCounters;
typedef Support<DEFMArray, DEFMCounterData, DEFMRuleData> DEFMSupport;
typedef StatsCounter<DEFMArray, DEFMCounterData> DEFMStatsCounter;
typedef Model<DEFMArray, DEFMCounterData,DEFMRuleData,DEFMRuleData> DEFMModel;
typedef Rule<DEFMArray, DEFMRuleData> DEFMRule;
typedef Rules<DEFMArray, DEFMRuleData> DEFMRules;
///@}

inline double DEFMData::operator()(size_t i, size_t j) const
{
    return *(covariates + (obs_start + j * X_nrow + i));
}

inline size_t DEFMData::ncol() const {
    return X_ncol;
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

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_COUNTER(a) \
inline double (a) (const DEFMArray & Array, uint i, uint j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a) \
Counter_fun_type<DEFMArray, DEFMCounterData> a = \
    [](const DEFMArray & Array, uint i, uint j, DEFMCounterData & data)

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
[](const DEFMArray & Array, uint i, uint j, DEFMRuleData & data)
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
    int covar_index = -1
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {      

        DEFM_COUNTER_LAMBDA(counter_tmp)
        {
            return Array.D()(Array.nrow() - 1u, data.idx(0u));

        };

        counters->add_counter(
            counter_tmp, nullptr,
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}), 
            "# of ones x attr" + std::to_string(covar_index), 
            "Overall number of ones"
        );

    } else {

        DEFM_COUNTER_LAMBDA(count_ones) {return 1.0;};

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData(),
            "# of ones", 
            "Overall number of ones"
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
inline void counter_transition(
    DEFMCounters * counters,
    std::vector< size_t > coords,
    std::vector< bool > signs,
    size_t m_order,
    size_t n_y,
    int covar_index = -1
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
            int c = std::floor(indices[i] / Array.nrow());
            int r = indices[i] - c * Array.nrow();

            if (c >= Array.ncol())
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
    std::string name = "Motif {";


    barry::BArrayDense<int> motif(m_order + 1u, n_y);
    for (size_t i = 0u; i < m_order; ++i)
        for (size_t j = 0u; j < n_y; ++j)
            motif(i, j) = 0;

    size_t n_cells = coords.size() - 1u;
    for (size_t d = 0u; d < n_cells; ++d)
    {
        size_t c = std::floor(coords[d] / (m_order + 1u));
        size_t r = coords[d] - c * (m_order + 1u);
        motif(r, c) = signs[d] ? 1 : -1;
        
    }

    // From
    for (size_t i = 0u; i < m_order; ++i)
    {
        for (size_t j = 0u; j < n_y; ++j)
        {

            // Is it included?
            if (motif(i,j) == 0)
                continue;

            // Is not the first?
            if ((i != 0u) | (j == 0u))
                name += ", ";

            name += (
                (motif(i,j) < 0 ? "-y" : "+y") +
                std::to_string(i) + std::to_string(j) + ")"
                );
        }
    }

    if (m_order > 0u)
        name += "} > ";

    for (size_t j = 0u; j < n_y; ++j)
    {

        if (motif(m_order, j) == 0)
            continue;

        if (j != 0u)
            name += ", ";

        name += (
            (motif(m_order, j) < 0 ? "-y" : "+y") +
            std::to_string(m_order) + std::to_string(j) + ")"
            );

    }

    if (covar_index >= 0)
    {

        counters->add_counter(
            count_ones, count_init,
            DEFMCounterData(coords, {}, signs), 
            name + " with attr " + std::to_string(covar_index), 
            "Motif weighted by single attribute"
        );

    } else {

        counters->add_counter(
            count_ones, count_init,
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
inline void counter_fixed_effect(
    DEFMCounters * counters,
    int covar_index,
    double k
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

    counters->add_counter(
        count_tmp, count_init,
        DEFMCounterData({static_cast<size_t>(covar_index)}, {k}, {}), 
        "Fixed effect feature " + std::to_string(covar_index) + "^" + std::to_string(k)
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
        DEFMRuleData({},{markov_order})
        );
    
    return;
}

///@}

///@}

#endif
