// #include <vector>
// #include <unordered_map>
#include "barray-bones.hpp"
#include "support-bones.hpp"
#include "statscounter-bones.hpp"
#include "rules-bones.hpp"

#ifndef BARRY_MODEL_BONES_HPP 
#define BARRY_MODEL_BONES_HPP 1

inline double update_normalizing_constant(
  const std::vector< double > & params,
  const Counts_type & support
) {
  
  double res = 0.0;
  double tmp;
  for (unsigned int n = 0u; n < support.size(); ++n) {
    
    tmp = 0.0;
    for (unsigned int j = 0u; j < params.size(); ++j)
      tmp += support[n].first[j] * params[j];
    
    res += exp(tmp) * support[n].second;
  }
  
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
  
  double numerator=0.0;
  
  // Computing the numerator
  for (unsigned int j = 0u; j < target_stats.size(); ++j)
    numerator += target_stats[j] * params[j];
  if (!log_)
    numerator = exp(numerator);
  
  if (log_)
    return numerator - log(normalizing_constant);
  
  return numerator/normalizing_constant;
  
}

/**@brief Array Hasher class (used for computing support)
 * 
 */
template<typename Array_Type>
inline std::vector< double > keygen_default(const Array_Type & Array_) {
  return {(double) Array_.N, (double) Array_.M};
}


/**
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
  typename Array_Type = BArray<>,
  typename Data_Counter_Type = bool,
  typename Data_Rule_Type = bool>
class Model {

public:
  
  /**@name Random number generation*/
  ///@{
  
  ///@}
  
  /**@brief */
  std::vector< Counts_type >         stats;
  std::vector< uint >                n_arrays_per_stats;
  
  /**
   * @name Container space for the powerset (and its sufficient stats)
   * @details This is useful in the case of using simulations or evaluating
   * functions that need to account for the full set of states.
   */
  ///@{
  bool with_pset = false;
  std::vector< std::vector< Array_Type > >          pset_arrays;
  std::vector< std::vector< std::vector<double> > > pset_stats;
  std::vector< std::vector<double> >                pset_probs;
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
  
  /**@name Functions to compute statistics
   * @details Arguments are recycled to save memory and computation.
   */
  ///@[{
  Counters<Array_Type,Data_Counter_Type>          counters;
  Rules<Array_Type,Data_Rule_Type>                     rules;
  Support<Array_Type,Data_Counter_Type,Data_Rule_Type> support_fun;
  StatsCounter<Array_Type,Data_Counter_Type>           counter_fun;
  ///@}
  
  /**@brief Vector of the previously used parameters */
  std::vector< std::vector<double> > params_last;
  std::vector< double > normalizing_constants;
  std::vector< bool > first_calc_done;

  /**@brief Function to extract features of the array to be hash
  */
  std::function<std::vector<double>(const Array_Type &)> keygen = nullptr;  

  Model();
  Model(uint size_);
  Model(const Model<Array_Type, Data_Counter_Type, Data_Rule_Type> & Model_);
  Model<Array_Type, Data_Counter_Type, Data_Rule_Type> & operator=(
    const Model<Array_Type, Data_Counter_Type, Data_Rule_Type> & Model_
    );
  ~Model() {};
  
  void store_psets();
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
    bool as_log = false
    );
  
  double likelihood_total(
    const std::vector<double> & params,
    bool as_log = false
  );
  ///@}
  
  void print_stats(uint i) const;
  
  Array_Type sample(const Array_Type & Array_, const std::vector<double> & params = {});
  Array_Type sample(const uint & i, const std::vector<double> & params = {});
  
  /**
   * @brief Number of different supports included in the model
   * 
   * This will return the size of `stats`.
   * 
   * @return unsigned int 
   */
  ///@{
  unsigned int size() const;
  unsigned int size_unique() const;
  ///@}
};


#endif
