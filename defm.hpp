
#ifndef DEFM_HPP
#define DEFM_HPP 1

// #include "../barry.hpp"

#include <iterator>
#include <regex>

namespace defm {

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/defm/defm-types.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef DEFM_TYPES_HPP
#define DEFM_TYPES_HPP
class DEFMData;

typedef barry::BArrayDense<int, DEFMData> DEFMArray;

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
  * @brief Data class used to store arbitrary size_t or double vectors */
class DEFMCounterData {
public:

    std::vector< size_t > indices;
    std::vector< double > numbers;
    std::vector< bool >   logical;
    bool is_motif; ///< If false, then is a logit intercept.
    
    DEFMCounterData() : indices(0u), numbers(0u), logical(0u), is_motif(true) {};
    DEFMCounterData(
        const std::vector< size_t > indices_,
        const std::vector< double > numbers_,
        const std::vector< bool > logical_,
        bool is_motif_ = true
    ): indices(indices_), numbers(numbers_), 
        logical(logical_), is_motif(is_motif_) {};

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
typedef barry::Counter<DEFMArray, DEFMCounterData > DEFMCounter;
typedef barry::Counters<DEFMArray, DEFMCounterData> DEFMCounters;
typedef barry::Support<DEFMArray, DEFMCounterData, DEFMRuleData,DEFMRuleDynData> DEFMSupport;
typedef barry::StatsCounter<DEFMArray, DEFMCounterData> DEFMStatsCounter;
typedef barry::Model<DEFMArray, DEFMCounterData,DEFMRuleData,DEFMRuleDynData> DEFMModel;


typedef barry::Rule<DEFMArray, DEFMRuleData> DEFMRule;
typedef barry::Rules<DEFMArray, DEFMRuleData> DEFMRules;
typedef barry::Rule<DEFMArray, DEFMRuleDynData> DEFMRuleDyn;
typedef barry::Rules<DEFMArray, DEFMRuleDynData> DEFMRulesDyn;
///@}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/defm/defm-types.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/defm/defm-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef DEFM_BONES_HPP
#define DEFM_BONES_HPP 1

class DEFM : public defmcounters::DEFMModel {
private:

    // std::shared_ptr< std::mt19937 > rengine = nullptr;
    // std::shared_ptr< defmcounters::DEFMModel > model = nullptr;

    /**
     * @brief Model data
     */
    ///@{
    const int * Y = nullptr;    ///< Outcome variable
    const int * ID = nullptr;   ///< Individual ids
    const double * X = nullptr; ///< Covariates
    
    size_t N;             ///< Number of agents/individuals
    size_t ID_length;     ///< Length of the vector IDs
    size_t Y_ncol;        ///< Number of columns in the response
    size_t Y_length;      ///< Length of the vector Y
    size_t X_ncol;        ///< Number of columns in the features
    size_t X_length;      ///< Length of the vector X
    size_t M_order;       ///< Markov order of the model

    std::vector< std::string > Y_names;
    std::vector< std::string > X_names;
    std::vector< size_t > start_end;
    std::vector< size_t > model_ord;
    ///@}

public:

    DEFM(
        const int * id,
        const int * y,
        const double * x,
        size_t id_length,
        size_t y_ncol,
        size_t x_ncol,
        size_t m_order
    );

    // ~DEFM() {
    //     defmcounters::DEFMModel::~Model();
    // };

    defmcounters::DEFMModel & get_model() {
        return *this;
    };

    void init();

    double likelihood(std::vector< double > & par, bool as_log = false);
    void simulate(std::vector< double > par, int * y_out);

    size_t get_n_y() const;
    size_t get_n_obs() const;
    size_t get_n_covars() const;
    size_t get_m_order() const;
    size_t get_n_rows() const;

    const int * get_Y() const;
    const int * get_ID() const;
    const double * get_X() const;

    barry::FreqTable<int> motif_census(
        std::vector< size_t > idx
    );

    std::vector< double > logodds(
        const std::vector< double > & par,
        size_t i,
        size_t j
    );

    void set_names(
        std::vector< std::string > Y_names_,
        std::vector< std::string > X_names_
    );

    const std::vector< std::string > & get_Y_names() const;
    const std::vector< std::string > & get_X_names() const;

    void print() const;

    std::vector< bool > is_motif();

};

#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/defm/defm-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/defm/defm-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef DEFM_MEAT_HPP
#define DEFM_MEAT_HPP 1

inline std::vector< double > keygen_defm(
    const defmcounters::DEFMArray & Array_,
    defmcounters::DEFMCounterData * data
    ) {
    
    size_t nrow = Array_.nrow();
    size_t ncol = Array_.ncol();

    std::vector< double > res(
        2u +                // Rows + cols
        ncol * (nrow - 1u) // Markov cells
        );

    res[0u] = static_cast<double>(nrow);
    res[1u] = static_cast<double>(ncol);

    size_t iter = 2u;
    // Adding the cells
    for (size_t i = 0u; i < (nrow - 1); ++i)
        for (size_t j = 0u; j < ncol; ++j)
            res[iter++] = Array_(i, j);

    return res;

}

#define DEFM_RANGES(a) \
    size_t start_i = start_end[a * 2u];\
    size_t end_i   = start_end[a * 2u + 1u];\
    size_t nobs_i  = end_i - start_i + 1u;

#define DEFM_LOOP_ARRAYS(a) \
    for (size_t a = 0u; a < (nobs_i - M_order); ++a)

inline void DEFM::simulate(
    std::vector< double > par,
    int * y_out
) {

    size_t model_num = 0u; 
    size_t n_entry = M_order * Y_ncol;
    auto idx = this->get_arrays2support();
    defmcounters::DEFMArray last_array;
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(proc_n)
        {

            // In the first process, we take the data as is
            if (proc_n == 0u)
            {
                last_array = this->sample(idx->at(model_num++), par);
                for (size_t y = 0u; y < Y_ncol; ++y)
                    *(y_out + n_entry++) = last_array(M_order, y, false);

                // last_array.print("i: %li, proc_n: %li\n", i, proc_n);

            }
            else
            // Otherwise, we need to continue using the previous data!
            {
                // Removing the previous row
                defmcounters::DEFMArray tmp_array(M_order + 1u, Y_ncol);
                for (size_t t_i = 1u; t_i < (M_order + 1u); ++t_i)
                    for (size_t t_j = 0u; t_j < Y_ncol; ++t_j)
                        tmp_array(t_i - 1u, t_j) = last_array(t_i, t_j);

                // Setting the data
                tmp_array.set_data(
                    new defmcounters::DEFMData(&tmp_array, X, (start_i + proc_n), X_ncol, ID_length),
                    true // Delete the data
                );

                // Baseline
                // tmp_array.print("baseline i: %li, proc_n: %li\n", i, proc_n);
                // tmp_array.D().print();

                model_num++;
                last_array = this->sample(tmp_array, par);
                for (size_t y = 0u; y < Y_ncol; ++y)
                    *(y_out + n_entry++) = last_array(M_order, y, false);

                // last_array.print("generated i: %li, proc_n: %li\n", i, proc_n);

            }


            
        }

        n_entry += M_order * Y_ncol;

    }

}

inline DEFM::DEFM(
    const int * id,
    const int * y,
    const double * x,
    size_t id_length,
    size_t y_ncol,
    size_t x_ncol,
    size_t m_order
) {

    // Pointers
    ID = id;
    Y  = y;
    X  = x;

    // Overall dimmensions
    ID_length = id_length;

    Y_ncol    = y_ncol;
    Y_length  = y_ncol * id_length;

    X_ncol    = x_ncol;
    X_length  = x_ncol * id_length;

    M_order   = m_order;

    // Creating the model and engine
    this->rengine = new std::mt19937();
    this->delete_rengine = true;

    this->store_psets();
    auto kgen = keygen_defm;
    this->add_hasher(kgen);

    // Iterating for adding observations
    start_end.reserve(id_length);
    start_end.push_back(0);

    // Identifying the start and end of each observation
    N = 0u;
    for (size_t row = 1u; row < id_length; ++row)
    {

        // Still in the individual
        if (*(id + row) != *(id + row - 1u))
        {

            // End of the previous observation
            start_end.push_back(row - 1u);

            // In the case that the start and end do not fit
            // within the markov process order, then it should fail
            size_t n_rows_i = (row - 1u) - start_end[N++ * 2u] + 1;
            if (n_rows_i < (M_order + 1u))
                throw std::length_error(
                    "Obs. id: " + std::to_string(*(id + row - 1u)) + " (row " +
                    std::to_string(row) + ") has fewer rows (" +
                    std::to_string(n_rows_i) + ") than those needed (" +
                    std::to_string(M_order + 1) + ") for the Markov Model."
                );

            // Beginning of the current
            start_end.push_back(row);

        }
        
    }

    start_end.push_back(id_length - 1u);

    N++;

    // Creating the names
    for (auto i = 0u; i < Y_ncol; ++i)
        Y_names.push_back(std::string("y") + std::to_string(i));

    for (auto i = 0u; i < X_ncol; ++i)
        X_names.push_back(std::string("X") + std::to_string(i));

    return;    

}


inline void DEFM::init() 
{

    // Adding the rule
    defmcounters::rules_markov_fixed(this->get_rules(), M_order);

    // Creating the arrays
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        size_t start_i = start_end[i * 2u];
        size_t end_i   = start_end[i * 2u + 1u];
        size_t nobs_i  = end_i - start_i + 1u;

        // Creating the observations.
        // Number of processes : (N rows) - (Process size)
        for (size_t n_proc = 0u; n_proc < (nobs_i - M_order); ++n_proc)
        {

            // Creating the array for process n_proc and setting the data
            defmcounters::DEFMArray array(M_order + 1u, Y_ncol);
            array.set_data(
                new defmcounters::DEFMData(&array, X, (start_i + n_proc), X_ncol, ID_length),
                true // Delete the data
            );

            // Filling-out the array
            for (size_t k = 0u; k < Y_ncol; ++k)
                for (size_t o = 0u; o < (M_order + 1u); ++o)
                    array(o, k) = *(Y + k * ID_length + start_i + n_proc + o);

            // Adding to the model
            model_ord.push_back( this->add_array(array, true) );

        }

    }

}

inline size_t DEFM::get_n_y() const
{
    return Y_ncol;
}

inline size_t DEFM::get_n_obs() const
{
    return N;
}

inline size_t DEFM::get_n_covars() const
{
    return X_ncol;
}

inline size_t DEFM::get_m_order() const
{
    return M_order;
}

inline size_t DEFM::get_n_rows() const
{
    return ID_length;
}

inline const int * DEFM::get_Y() const
{
    return Y;
}

inline const int * DEFM::get_ID() const
{
    return ID;
}

inline const double * DEFM::get_X() const
{
    return X;
}


inline barry::FreqTable<int> DEFM::motif_census(
        std::vector< size_t > idx
) {

    // Checking all sizes
    for (const auto & i : idx)
        if (i >= Y_ncol)
            throw std::range_error("The -idx- for motif accounting is out of range.");

    barry::FreqTable<int> ans;
    std::vector<int> array(idx.size() * (M_order + 1));

    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(proc_n)
        {

            // Generating an integer array between the parts
            size_t nele = 0u;

            for (size_t o = 0u; o < (M_order + 1u); ++o)
                for (auto & k : idx)
                    array[nele++] = *(Y + k * ID_length + start_i + proc_n + o);

            ans.add(array, nullptr);

        }

    }

    return ans;

}

inline std::vector< double > DEFM::logodds(
    const std::vector< double > & par,
    size_t i_,
    size_t j_
) {
    

    std::vector< double > res(ID_length, std::nan(""));

    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(n_proc)
        {

            // Creating the array for process n_proc and setting the data
            defmcounters::DEFMArray array(M_order + 1u, Y_ncol);
            array.set_data(
                new defmcounters::DEFMData(&array, X, (start_i + n_proc), X_ncol, ID_length),
                true // Delete the data
            );

            // Filling-out the array
            for (size_t k = 0u; k < Y_ncol; ++k)
                for (size_t o = 0u; o < (M_order + 1u); ++o)
                    array(o, k) = *(Y + k * ID_length + start_i + n_proc + o);

            double p_1 = this->conditional_prob(array, par, i_, j_);
            res[M_order + start_i + n_proc] = std::log(p_1/(1.0 - p_1));

        }

    }

    return res;


}

inline void DEFM::set_names(
    std::vector< std::string > Y_names_,
    std::vector< std::string > X_names_
) {

    // Checking the length
    if (Y_names_.size() != Y_ncol)
        throw std::length_error("The length of Y_names_ doesn't match the number of dependent variables.");

    if (X_names_.size() != X_ncol)
        throw std::length_error("The length of X_names_ doesn't match the number of dependent variables.");

    Y_names = Y_names_;
    X_names = X_names_;

}

inline const std::vector<std::string > & DEFM::get_Y_names() const {
    return Y_names;
}

inline const std::vector<std::string > & DEFM::get_X_names() const {
    return X_names;
}

inline void DEFM::print() const
{
    defmcounters::DEFMModel::print();
    printf_barry("Model Y variables (%i):\n", static_cast<int>(get_n_y()));
    int ny = 0;
    for (const auto & y : get_Y_names())
    {

        printf_barry(" % 2i) %s\n", ny++, y.c_str());

    }
}

inline std::vector< bool > DEFM::is_motif()
{
    std::vector< bool > res(0u);
    auto * counterss = defmcounters::DEFMModel::get_counters();
    for (size_t i = 0u; i < counters->size(); ++i)
        res.push_back(counterss->operator[](i).data.is_motif);

    return res;
}

#undef DEFM_RANGES
#undef DEFM_LOOP_ARRAYS

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/defm/defm-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/defm/counters.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models//defm/formula.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_DEFM_MOTIF_FORMULA_HPP
#define BARRY_DEFM_MOTIF_FORMULA_HPP
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
#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models//defm/formula.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/**
 * @ingroup counting 
 * @details Details on the available counters for `DEFMworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{





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
inline double (a) (const DEFMArray & Array, size_t i, size_t j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a) \
Counter_fun_type<DEFMArray, DEFMCounterData> a = \
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
Rule_fun_type<DEFMArray, DEFMRuleData> a = \
[](const DEFMArray & Array, size_t i, size_t j, DEFMRuleData & data) -> bool
///@}

/**Lambda function for definition of a network counter function*/
#define DEFM_RULEDYN_LAMBDA(a) \
Rule_fun_type<DEFMArray, DEFMRuleDynData> a = \
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
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}, true), 
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
            DEFMCounterData(coords, {}, signs, coords.size() > 1u ? true : false), 
            name + " x " + vname, 
            "Motif weighted by single attribute"
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
        return (Array(i - 1, j) != 1) || (Array(i, j) != 1);

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
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/defm/counters.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



}

#endif
