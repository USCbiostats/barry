// #include <vector>
// #include <unordered_map>
#include "barray-bones.hpp"
#include "support-bones.hpp"
#include "statscounter-bones.hpp"
#include "rules-bones.hpp"

#ifndef BARRY_MODEL_BONES_HPP 
#define BARRY_MODEL_BONES_HPP 1


/**
 * @brief Array Hasher class (used for computing support)
 * 
 */
template<typename Array_Type>
inline std::vector< double > keygen_default(const Array_Type & Array_) {
    return {static_cast<double>(Array_.nrow()), static_cast<double>(Array_.ncol())};
}


/**
  * @ingroup stat-models
  * @brief General framework for discrete exponential models.
  * This class allows generating discrete exponential models in the form of a linear
  * exponential model:
  * \f[
  * \frac{
  *    \exp{\left(\theta^{\mbox{t}}c(A)\right)}
  *  }{
  *    \sum_{A'\in \mathcal{A}}\exp{\left(\theta^{\mbox{t}}c(A')\right)}
  *  }
  * \f]
  * 
  * This implementation aims to reduce the number of times that the support
  * needs to be computed. Models included here use more than a single array, and
  * thus allow the function to recycle support sets as needed. For example,
  * if we are looking at directed graphs all of the same size and without
  * vertex level features, i.e. a model that only counts edges, triangles, etc.
  * then the support needs to be fully computed only once.
  * 
  * @tparam Array_Type Class of `BArray` object.
  * @tparam Data_Counter_Type Any type.
  * @tparam Data_Rule_Type Any type.
  */
template <
    typename Array_Type         = BArray<>,
    typename Data_Counter_Type  = bool,
    typename Data_Rule_Type     = bool,
    typename Data_Rule_Dyn_Type = bool
    >
class Model {

private:
    /**
      * @name Random number generation
      * @brief Random number generation
      */
    ///@{
    std::mt19937 * rengine = nullptr;
    bool delete_rengine    = false;

    std::vector< Counts_type >         stats;              ///< Sufficient statistics of the model (support)
    std::vector< uint >                n_arrays_per_stats; ///< Number of arrays included per support.

    /**
      * @name Container space for the powerset (and its sufficient stats)
      * @details This is useful in the case of using simulations or evaluating
      * functions that need to account for the full set of states.
      */
    ///@{
    bool with_pset = false;
    std::vector< std::vector< Array_Type > >          pset_arrays; ///< Arrays of the support(s)
    std::vector< std::vector< std::vector<double> > > pset_stats;  ///< Statistics of the support(s)
    std::vector< std::vector<double> >                pset_probs;  ///< Probabilities of the support(s)
    ///@}

    /**
      * @name Information about the arrays used in the model 
      * @details `target_stats` holds the observed sufficient statistics for each
      * array in the dataset. `array_frequency` contains the frequency with which
      * each of the target stats (arrays) shows in the support. `array2support` 
      * maps array indices (0, 1, ...) to the corresponding support.
      */
    ///@{
    std::vector< std::vector< double > > target_stats;
    std::vector< uint >                 array_frequency;
    std::vector< uint >                 arrays2support;
    ///@}

    /**
      * @brief Map of types of arrays to support sets
      * @details This is of the same length as the vector `stats`.
      */
    MapVec_type< double, uint > keys2support;
    
    /**
      * @name Functions to compute statistics
      * @details Arguments are recycled to save memory and computation.
      */
    ///@{
    Counters<Array_Type,Data_Counter_Type> *                                counters;
    Rules<Array_Type,Data_Rule_Type> *                                      rules;
    Rules<Array_Type,Data_Rule_Dyn_Type> *                                  rules_dyn;
    Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> support_fun;
    StatsCounter<Array_Type,Data_Counter_Type>                              counter_fun;
    ///@}
    
    /**@brief Vector of the previously used parameters */
    std::vector< std::vector<double> > params_last;
    std::vector< double > normalizing_constants;
    std::vector< bool > first_calc_done;

    /**@brief Function to extract features of the array to be hash
    */
    std::function<std::vector<double>(const Array_Type &)> keygen = nullptr;  

    bool delete_counters  = false;
    bool delete_rules     = false;
    bool delete_rules_dyn = false;

public:
    
    void set_rengine(std::mt19937 * rengine_, bool delete_ = false) {

        if (delete_rengine)
            delete rengine;

        rengine        = rengine_;
        delete_rengine = delete_;
        
    };

    void set_seed(unsigned int s) {

        if (rengine == nullptr) {
            rengine = new std::mt19937;
            delete_rengine = true;
        }

        rengine->seed(s);

    };
    ///@}
        
    Model();
    Model(uint size_);
    Model(const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_);
    Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & operator=(
        const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_
    );
    ~Model() {
        if (delete_counters)
            delete counters;

        if (delete_rules)
            delete rules;

        if (delete_rules_dyn)
            delete rules_dyn;

        if (delete_rengine)
            delete rengine;
    };
    
    void store_psets() noexcept;
    void set_keygen(std::function<std::vector<double>(const Array_Type &)> keygen_);
    
    /**
      * @name Wrappers for the `Counters` member. 
      * @details These will add counters to the model, which are shared by the
      * support and the actual counter function.
      */
    ///@{
    void add_counter(Counter<Array_Type, Data_Counter_Type> & counter);
    void add_counter(Counter<Array_Type, Data_Counter_Type> * counter);
    void add_counter(
        Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_    = nullptr,
        Data_Counter_Type *                            data_        = nullptr,
        bool                                           delete_data_ = false
    );
    void set_counters(Counters<Array_Type,Data_Counter_Type> * counters_);
    ///@}
    
    /**
      * @name Wrappers for the `Rules` member. 
      * @details These will add rules to the model, which are shared by the
      * support and the actual counter function.
      */
    ///@{
    void add_rule(Rule<Array_Type, Data_Rule_Type> & rule);
    void add_rule(Rule<Array_Type, Data_Rule_Type> * rule);
    void add_rule(
        Rule_fun_type<Array_Type, Data_Rule_Type> count_fun_,
        Data_Rule_Type *                          data_        = nullptr,
        bool                                      delete_data_ = false
    );
    
    void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);

    void add_rule_dyn(Rule<Array_Type, Data_Rule_Dyn_Type> & rule);
    void add_rule_dyn(Rule<Array_Type, Data_Rule_Dyn_Type> * rule);
    void add_rule_dyn(
        Rule_fun_type<Array_Type, Data_Rule_Dyn_Type> count_fun_,
        Data_Rule_Dyn_Type *                          data_        = nullptr,
        bool                                      delete_data_ = false
    );
    
    void set_rules_dyn(Rules<Array_Type,Data_Rule_Dyn_Type> * rules_);
    ///@}
    

    /**
      * @brief Adds an array to the support of not already included.
      * @param Array_ array to be added
      * @param force_new If `false`, it will use `keygen` to obtain a double vector
      * and create a hash of it. If the hash has been computed earlier, the support
      * is recycled.
      * 
      * @return The number of the array.
      */
    uint add_array(const Array_Type & Array_, bool force_new = false);
    
    
    /**
      * @name Likelihood functions.
      * @details Calculation of likelihood functions is done reusing normalizing
      * constants. Before recalculating the normalizing constant, the function 
      * checks whether `params` matches the last set vector of parameters used
      * to compute it.
      * 
      * 
      * @param params Vector of parameters
      * @param as_log When `true`, the function returns the log-likelihood.
      */
    ///@{
    double likelihood(
        const std::vector<double> & params,
        const uint & i,
        bool as_log = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const Array_Type & Array_,
        int i = -1,
        bool as_log = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const std::vector<double> & target_,
        const uint & i,
        bool as_log = false
    );

    double likelihood_total(
        const std::vector<double> & params,
        bool as_log = false
    );
    ///@}

    /**
      * @name Extract elements by index 
      * @param i Index relative to the array in the model.
      * @param params A new vector of model parameters to compute the normalizing
      * constant.
      * @param as_log When `true` returns the logged version of the normalizing
      * constant.
      */
    ///@{
    double get_norm_const(
        const std::vector< double > & params,
        const uint & i,
        bool as_log = false
    );

    const std::vector< Array_Type > * get_pset(
        const uint & i
    );

    const std::vector< std::vector< double > > * get_stats(
        const uint & i
    );
    ///@}
    
    void print_stats(uint i) const;
    
    Array_Type sample(const Array_Type & Array_, const std::vector<double> & params = {});
    Array_Type sample(const uint & i, const std::vector<double> & params);
    
    /**
      * @name Size of the model
      * 
      * @brief Number of different supports included in the model
      * 
      * This will return the size of `stats`.
      * 
      * @return `size()` returns the number of arrays in the model.
      * @return `size_unique()` returns the number of unique arrays (according to
      * the hasher) in the model.
      * @return `nterms()` returns the number of terms in the model.
      */
    ///@{
    unsigned int size() const noexcept;
    unsigned int size_unique() const noexcept;
    unsigned int nterms() const noexcept;
    unsigned int support_size() const noexcept;
    std::vector< std::string > colnames() const;
    ///@}

    const std::mt19937 * get_rengine() const;

    Counters<Array_Type,Data_Counter_Type> * get_counters();
    Rules<Array_Type,Data_Rule_Type>       * get_rules();
    Rules<Array_Type,Data_Rule_Dyn_Type>   * get_rules_dyn();
    Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> * get_support();

};


#endif
