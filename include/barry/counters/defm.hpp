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
            return Array.D()(static_cast<size_t>(i), data.idx(0u));

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
    int covar_index = -1
)
{

    // A vector to store the type of dat
    if (signs.size() == 0u)
        signs.resize(coords.size(), true);
    else if (signs.size() != coords.size())
        throw std::length_error("Size of -coords- and -signs- must match.");

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {      

        coords.push_back(static_cast<size_t>(covar_index));

        DEFM_COUNTER_LAMBDA(count_ones)
        {
            
            auto dat = data.indices;
            auto sgn = data.logical;

            // Checking if the observation is in the stat. We 
            const auto & array = Array.get_data();
            size_t loc = i + j * Array.nrow();
            size_t n_cells = dat.size() - 1u;


            // Only one currently needs to be a zero for it
            // to change
            size_t n_present = 0;
            bool baseline_value = 0;
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
                    n_present++;
                
            }

            // If i in array still false, then no change
            if (!i_in_array)
                return 0.0;

            // If the difference is greater than one, then nothing
            // happens
            if (std::fabs(n_present - n_cells) > 1)
                return 0.0;

            double val = Array.D()(static_cast<size_t>(i), static_cast<size_t>(dat[n_cells]));
            if (n_present == n_cells) // We now match (regardless)
                return val;

            // We know we added one now, so we have two cases:
            // false -> Now disagreen, so removed a counter
            //   n_present > n_cells: Was above alreadu => 0.0;
            //   n_present < n_cells: Used to match => -val;
            // true -> Now agree so adding a counter
            //   n_present > n_cells: Used to match => -val;
            //   n_present < n_cells: Was below already => 0.0;


            if (!baseline_value)
                return (n_present > n_cells) ? 0.0 : -val;
            else
                return (n_present > n_cells) ? -val: 0.0;

        };

        // Creating name of the structure
        std::string name = "Motif";
        for (size_t d = 0u; d < (coords.size() - 1u); ++d)
            name += (" "+ std::to_string(coords[d]));

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData(coords, {}, signs), 
            name + " with attr " + std::to_string(covar_index), 
            "Motif weighted by single attribute"
        );

    } else {

        DEFM_COUNTER_LAMBDA(count_ones)
        {

            auto dat = data.indices;
            auto sgn = data.logical;

            // Checking if the observation is in the stat. We 
            const auto & array = Array.get_data();
            size_t loc = i + j * Array.nrow();
            size_t n_cells = dat.size();

            // Only one currently needs to be a zero for it
            // to change
            size_t n_present = 0;
            bool baseline_value = 0;
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
                    n_present++;
                
            }

            // If i in array still false, then no change
            if (!i_in_array)
                return 0.0;

            // If the difference is greater than one, then nothing
            // happens
            if (std::fabs(n_present - n_cells) > 1)
                return 0.0;

            if (n_present == n_cells) // We now match (regardless)
                return 1.0;

            // We know we added one now, so we have two cases:
            // false -> Now disagreen, so removed a counter
            //   n_present > n_cells: Was above alreadu => 0.0;
            //   n_present < n_cells: Used to match => -val;
            // true -> Now agree so adding a counter
            //   n_present > n_cells: Used to match => -val;
            //   n_present < n_cells: Was below already => 0.0;


            if (!baseline_value)
                return (n_present > n_cells) ? 0.0 : -1.0;
            else
                return (n_present > n_cells) ? -1.0: 0.0;   

        };

        // Creating name of the structure
        std::string name = "Motif";
        for (size_t d = 0u; d < coords.size(); ++d)
            name += (" "+ std::to_string(coords[d]));

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData({coords}, {}, signs), 
            name, 
            "Structural motif"
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
