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
class DEFMData {
public:
    
    std::vector< double > * covariates; ///< Vector of covariates (complete vector)
    size_t obs_start;    ///< Index of the observation in the data.
    size_t obs_n_times;  ///< Number of records of the observation in the model.
    size_t n_covariates; ///< Number of covariates included in the model.
    
    DEFMData() {};
    
    /**
     * @brief Constructor
     * @param covariates_ Pointer to the attribute data.
     * @param obs_start_ Location of the current observation in the covariates
     *  vector
     * @param obs_n_times_ Number of observations in the model
     * @param n_covariates_ Number of columns (covariates.)
     */
    DEFMData(
        std::vector< double > * covariates_,
        size_t obs_start_,
        size_t obs_n_times_,
        size_t n_states_,
        size_t n_covariates_
    ) : obs_start(obs_start_), obs_n_times(obs_n_times_),
    n_states(n_states_), n_covariates(n_covariates_) {}; 

    /**
     * @brief Access to the row (i) colum (j) data
     * 
     * @param i 
     * @param j 
     * @return double 
     */
    double operator()(size_t i, size_t j);
    double at(size_t i, size_t j);
    
    ~DEFMData() {};

};

inline double DEFMData::operator()(size_t i, size_t j)
{
    return covariates->operator[](obs_start + i * n_covariates + j);
}

inline double DEFMData::at(size_t i, size_t j)
{
    return covariates->at(obs_start + i * n_covariates + j);
}

/**
  * @brief Data class used to store arbitrary uint or double vectors */
class DEFMCounterData {
public:
    
    std::vector< size_t > indices;
    std::vector< double > numbers;
    size_t markov_order; ///< Order of the markov process
    
    DEFMCounterData() : indices(0u), numbers(0u) {};
    DEFMCounterData(
        const std::vector< uint > indices_,
        const std::vector< double > numbers_,
        size_t markov_order_
    ): indices(indices_), numbers(numbers_), markov_order(markov_order_) {};

    size_t idx(size_t i) {return indices[i];};
    double num(size_t i) {return numbers[i];};
    
    ~DEFMCounterData() {};
    
};

/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef BArrayDense<int, DEFMData> DEFMArray;

template <typename Tarray = DEFMArray>
using DEFMCounter =  Counter<Tarray, DEFMCounterData >;

template <typename Tarray = DEFMArray>
using DEFMCounters =  Counters<Tarray, DEFMCounterData>;

template <typename Tarray = DEFMArray>
using DEFMSupport =  Support<Tarray, DEFMCounterData >;

template <typename Tarray = DEFMArray>
using DEFMStatsCounter =  StatsCounter<Tarray, DEFMCounterData>;

template <typename Tarray>
using DEFMModel =  Model<Tarray, DEFMCounterData>;

template <typename Tarray = DEFMArray>
using DEFMRule =  Rule<Tarray, bool>;

template <typename Tarray = DEFMArray>
using DEFMRules =  Rules<Tarray, bool>;
///@}

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_COUNTER(a) \
template<typename Tarray = DEFMArray>\
inline double (a) (const Tarray & Array, uint i, uint j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a,b) \
Counter_fun_type<Tarray, DEFMCounterData> a = \
    [b](const Tarray & Array, uint i, uint j, DEFMCounterData & data)

#define NETWORKDENSE_COUNTER_LAMBDA(a) \
Counter_fun_type<DEFMworkDense, DEFMCounterData> a = \
    [](const DEFMworkDense & Array, uint i, uint j, DEFMCounterData & data)
///@}


/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_RULE(a) \
template<typename Tarray = DEFMArray>\
inline bool (a) (const Tarray & Array, uint i, uint j, bool & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_RULE_LAMBDA(a) \
Rule_fun_type<Tarray, bool> a = \
[](const Tarray & Array, uint i, uint j, bool & data)
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
 * @tparam Tarray 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
template<typename Tarray = DEFMArray>
inline void counter_ones(
    DEFMCounters<Tarray> * counters,
    int covar_index = -1
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {      

        DEFM_COUNTER_LAMBDA(count_ones, vid = covar_index)
            {
                return Array.D()(i, vid);

            };

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData(), 
            "# of ones x attr" + std::to_str(covar_index), 
            "Overall number of ones"
        );

    } else {

        DEFM_COUNTER_LAMBDA(count_ones,) {return 1.0;};

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
 * @tparam Tarray 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
template<typename Tarray = DEFMArray>
inline void counter_transition(
    DEFMCounters<Tarray> * counters,
    std::vector< int > coords,
    int covar_index = -1
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {      

        coords.push_back(covar_index);

        DEFM_COUNTER_LAMBDA(count_ones, dat = coords)
        {
            
            // Checking if the observation is in the stat. We 
            const auto & array = Array.get_data();
            int loc = i + j * Array.nrow();

            // Only one currently needs to be a zero for it
            // to change
            int n_present = 0;
            bool i_in_array = false;
            for (size_t e = 0u; e < (dat.size() - 1); ++e)
                if (array[dat[e]] != 1)
                {
                    n_present++;
                    
                    if (dat[e] == loc)
                        i_in_array = true;
                }

            // If i in array still false, then no change
            if (!i_in_array)
                return 0.0;

            // Measuring the number of elements
            int nele = static_cast<int>(dat.size()) - 1;

            // We now match, so adding a new one
            if (nele == n_present)
                return Array.D()(i, dat[nele]);

            return 0.0;            

        };

        // Creating name of the structure
        std::string name = "Motif";
        for (size_t d = 0u; d < (coords.size() - 1u); ++d)
            name += (" "+ std::to_str(coords[d]));

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData(), 
            name + " with attr " + std::to_str(covar_index), 
            "Motif weigher by single attribute"
        );

    } else {

        DEFM_COUNTER_LAMBDA(count_ones, dat = coords)
        {
            
            // Checking if the observation is in the stat. We 
            const auto & array = Array.get_data();
            int loc = i + j * Array.nrow();

            // Only one currently needs to be a zero for it
            // to change
            int n_present = 0;
            bool i_in_array = false;
            for (size_t e = 0u; e < dat.size(); ++e)
                if (array[dat[e]] != 1)
                {
                    n_present++;
                    
                    if (dat[e] == loc)
                        i_in_array = true;
                }

            // If i in array still false, then no change
            if (!i_in_array)
                return 0.0;

            // Measuring the number of elements
            int nele = static_cast<int>(dat.size()) - 1;

            // We now match, so adding a new one
            if (nele == n_present)
                return 1.0;

            return 0.0;            

        };

        // Creating name of the structure
        std::string name = "Motif";
        for (size_t d = 0u; d < (coords.size() - 1u); ++d)
            name += (" "+ std::to_str(coords[d]));

        counters->add_counter(
            count_ones, nullptr,
            DEFMCounterData(), 
            name, 
            "Structural motif"
        );
    }

    return;

}

/**
 * @name Rules for network models
 * @param rules A pointer to a `DEFMRules` object (`Rules`<`DEFMArray`, `bool`>).
 */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
template<typename Tarray = DEFMArray>
inline void rules_zerodiag(DEFMRules<Tarray> * rules) {
    
    DEFM_RULE_LAMBDA(no_self_tie) {
        return i != j;
    };
    
    rules->add_rule(no_self_tie);
    
    return;
}

///@}

///@}

#endif
