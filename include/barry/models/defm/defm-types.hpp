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
    bool column_major;
    
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
        size_t X_nrow_,
        bool column_major_
    ) : array(array_), covariates(covariates_), obs_start(obs_start_),
    X_ncol(X_ncol_), X_nrow(X_nrow_), column_major(column_major_) {}; 

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


inline double DEFMData::operator()(size_t i, size_t j) const
{

    if (column_major)
        return *(covariates + (obs_start + i + X_nrow * j));
    else
        return *(covariates + ((obs_start + i) * X_ncol + j));

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

        printf_barry(
            "row %i (%i): ",
            static_cast<int>(i),
            static_cast<int>(obs_start + i)
        );
        for (size_t j = 0u; j < X_ncol; ++j)
        {
            printf_barry("% 5.2f, ", operator()(i, j));
        }
        printf_barry("\n");
        
    }

}

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{

class DEFMRuleDynData {
public:
    const std::vector< double > * counts;
    size_t pos;
    size_t lb;
    size_t ub;
    
    DEFMRuleDynData(
        const std::vector< double > * counts_,
        size_t pos_,
        size_t lb_,
        size_t ub_
        ) : counts(counts_), pos(pos_), lb(lb_), ub(ub_) {};
    
    ~DEFMRuleDynData() {};

    const double operator()() const
    {
        return (*counts)[pos];
    }
    
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