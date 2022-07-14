#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

/**
 * @brief Parses a motif formula
 * 
 * @details This function will take the formula and generate the corresponding
 * input for defm::counter_transition(). Formulas can be specified in the
 * following ways:
 * 
 * - Intercept effect: {...} No transition, only including the current state.
 * - Transition effect: {...} > {...} Includes current and previous states.
 * 
 * The general notation is `[0]y[column id]_[row id]`. A preceeding zero
 * means that the value of the cell is considered to be zero. The column
 * id goes between 0 and the number of columns in the array - 1 (so it
 * is indexed from 0,) and the row id goes from 0 to m_order.
 * 
 * ## Intercept effects
 * 
 * Intercept effects only involve a single set of curly brackets. Using the
 * 'greater-than' symbol (i.e., '<') is only for transition effects. When
 * specifying intercept effects, users can skip the `row_id`, e.g.,
 * `y0_0` is equivalent to `y0`. If the passed `row id` is different from
 * the Markov order, i.e., `row_id != m_order`, then the function returns
 * with an error. 
 * 
 * Examples:
 * 
 * - `"{y0, 0y1}"` is equivalent to set a motif with the first element equal
 * to one and the second to zero. 
 * 
 * ## Transition effects
 * 
 * Transition effects can be specified using two sets of curly brackets and
 * an greater-than symbol, i.e., `{...} > {...}`. The first set of brackets,
 * which we call LHS, can only hold `row id` that are less than `m_order`.
 * 
 * 
 * 
 * @param formula 
 * @param locations 
 * @param signs 
 * @param m_order 
 * @param y_ncol 
 */
inline void defm_motif_parser(
    std::string formula,
    std::vector< size_t > & locations,
    std::vector< bool > & signs,
    size_t m_order,
    size_t y_ncol
)
{
    // Resetting the results
    locations.clear();
    signs.clear();

    std::regex pattern_intercept(
        "\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\s*\\}"
        );
    std::regex pattern_transition(
        std::string("\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\}\\s*(>)\\s*") +
        std::string("\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\s*\\}")
        );

    auto empty = std::sregex_iterator();

    // This column-major vector indicates true if the variable has already been
    // selected
    std::vector< bool > selected((m_order + 1) * y_ncol, false);

    std::smatch match;
    std::regex_match(formula, match, pattern_transition);
    if (!match.empty())
    {

        if (m_order == 0)
            throw std::logic_error("Transition effects are only valid when the data is a markov process.");

        // Will indicate where the arrow is located at
        size_t arrow_position = match.position(4u);

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {

            // Baseline position
            size_t current_location = i->position(0u);

            // First value true/false
            bool is_positive;
            if (i->operator[](1u).str() == "")
                is_positive = true;
            else if (i->operator[](1u).str() == "0")
                is_positive = false;
            else
                throw std::logic_error("The number preceding y should be either none or zero.");

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            std::string tmp_str = i->operator[](4u).str();
            if (m_order > 1)
            {
                // If missing, we replace with the location 
                if (tmp_str == "")
                {

                    if (current_location > arrow_position)
                        y_row = m_order;
                    else
                        throw std::logic_error("LHS of transition must specify time when m_order > 1");

                } else
                    y_row = std::stoul(tmp_str);

                if (y_row > m_order)
                    throw std::logic_error("The proposed row is out of range.");


            } else {

                // If missing, we replace with the location 
                if (tmp_str != "")
                    y_row = std::stoul(tmp_str);
                else
                    y_row = (current_location < arrow_position ? 0u: 1u);

            }

            if (selected[y_col * (m_order + 1) + y_row])
                throw std::logic_error(
                    "The term " + i->str() + " shows more than once in the formula.");

            // Only the end of the chain can be located at position after the
            // arrow
            if ((current_location > arrow_position) && (y_row != m_order))
                throw std::logic_error(
                    "Only the row " + std::to_string(m_order) +
                    " can be specified at the RHS of the motif."
                    );

            selected[y_col * (m_order + 1) + y_row] = true;

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    std::regex_match(formula, match, pattern_intercept);
    if (!match.empty()){

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {
            
            // First value true/false
            bool is_positive;
            if (i->operator[](1u).str() == "")
                is_positive = true;
            else if (i->operator[](1u).str() == "0")
                is_positive = false;
            else
                throw std::logic_error("The number preceding y should be either none or zero.");

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            if (i->operator[](4u).str() == "") // Assume is the last
                y_row = m_order;
            else {

                y_row = std::stoul(i->operator[](4u).str());

                if (y_row != m_order)
                    throw std::logic_error(
                        std::string("Intercept motifs cannot feature past events. ") +
                        std::string("Only transition motifs can: {...} > {...}.")
                        );

            }

            if (selected[y_col * (m_order + 1) + y_row])
                throw std::logic_error(
                    "The term " + i->str() + " shows more than once in the formula.");

            selected[y_col * (m_order + 1) + y_row] = true;

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    throw std::logic_error(
        "The motif specified in the formula: " + formula +
        " has the wrong syntax."
        );
    
}

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
    int covar_index = -1,
    std::string vname = ""
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {      

        DEFM_COUNTER_LAMBDA(counter_tmp)
        {

            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return Array.D()(i, data.idx(0u));

        };

        counters->add_counter(
            counter_tmp, nullptr,
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}), 
            "# of ones x " + ((vname != "")? vname : ("attr" + std::to_string(covar_index))), 
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
            count_ones, nullptr,
            DEFMCounterData(),
            "# of ones", 
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
    std::string vname = ""
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

            counters->add_counter(
                tmp_counter, nullptr,
                DEFMCounterData({i}, {}, {}), 
                "Logit intercept " + std::to_string(i), 
                "Equal to one if the outcome " + std::to_string(i) + " is one. Equivalent to the logistic regression intercept."
            );

        }

    } else {

        DEFM_COUNTER_LAMBDA(tmp_counter)
        {
            if (i != Array.nrow())
                return 0.0;

            if (j != data.idx(0u))
                return 0.0;

            return Array.D()(i, j);
        };

        for (auto i : which)
        {

            counters->add_counter(
                tmp_counter, nullptr,
                DEFMCounterData({i}, {}, {}), 
                "Logit intercept " + std::to_string(i) + " x " + ((vname != "")? vname : ("attr" + std::to_string(covar_index))), 
                "Equal to one if the outcome " + std::to_string(i) + " is one. Equivalent to the logistic regression intercept."
            );

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
    std::string vname = ""
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
    

    if (any_before_event)
        name += "{";

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

            name += (motif(i,j) < 0 ? "y\u207B" : "y\u207A");

            if (m_order > 1)
                name += UNI_SUB(i);
                
            name += UNI_SUB(j);
        }
    }

    #define UNI0S \u2080
    #define UNI1S \u2081
    #define UNI2S \u2082
    #define UNI3S \u2083
    #define UNI4S \u2084
    #define UNI5S \u2085
    #define UNI6S \u2086
    #define UNI7S \u2087
    #define UNI8S \u2088
    #define UNI9S \u2089

    // If it has starting point, then need to close.
    if (any_before_event & (m_order > 0u))
        name += "} \u21E8 {";
    else
        name += "{";

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

        name += (motif(m_order, j) < 0 ? "y\u207B" : "y\u207A" );

        if (m_order > 1)
            name += UNI_SUB(m_order);

        name += UNI_SUB(j);

    }

    #undef UNI_SUB

    name += "}";

    if (covar_index >= 0)
    {

        counters->add_counter(
            count_ones, count_init,
            DEFMCounterData(coords, {}, signs), 
            name + " x " + ((vname != "")? vname : ("attr" + std::to_string(covar_index))), 
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
inline void counter_transition_formula(
    DEFMCounters * counters,
    std::string formula,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    std::string vname = ""
) {

    std::vector< size_t > coords;
    std::vector< bool > signs;

    defm_motif_parser(
        formula, coords, signs, m_order, n_y
    );

    counter_transition(
        counters, coords, signs, m_order, n_y, covar_index, vname
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
    std::string vname = ""
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
        "Fixed effect feature (" +
        ((vname != "")? vname : ("attr" + std::to_string(covar_index))) + ")^" + std::to_string(k)
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
