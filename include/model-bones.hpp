// #include <vector>
// #include <unordered_map>
#include "barray-bones.hpp"
#include "support.hpp"
#include "statscounter.hpp"

#ifndef MODEL_BONES_HPP 
#define MODEL_BONES_HPP 1

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
struct Array_Hasher {
  
  virtual std::vector<double> operator()(Array_Type const& dat) const noexcept {
    return {(double) dat.N, (double) dat.M};
  }
  
};

template<>
struct Array_Hasher<BArray<>> {
  std::vector<double> operator()(BArray<> const&dat) const noexcept {
    return {(double) dat.N, (double) dat.M};
  }
};

template<typename Array_Type>
inline std::vector< double > keygen_default(const Array_Type & Array_) {
  return {(double) Array_.N, (double) Array_.M};
}


/**@brief General framework for discrete exponential models.
 * 
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
 */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Model {

public:
  
  /**@brief */
  std::vector< Counts_type >         stats;
  std::vector< uint >                n_arrays_per_stats;
  
  /**@name Container space for the powerset (and its sufficient stats)
   * @details This is useful in the case of using simulations or evaluating
   * functions that need to account for the full set of states.
   */
  ///@{
  bool with_pset = false;
  std::vector< std::vector< Array_Type > >          pset_arrays;
  std::vector< std::vector< std::vector<double> > > pset_stats;
  ///@}
  
  /**@name Information about the arrays used in the model 
   * @details `target_stats` holds the observed sufficient statistics for each
   * array in the dataset. `array_frequency` contains the frequency with which
   * each of the target stats (arrays) shows in the support. `array2support` 
   * maps array indices (0, 1, ...) to the corresponding support.
   */
  ///@{
  std::vector< std::vector< double >> target_stats;
  std::vector< uint >                 array_frequency;
  std::vector< uint >                 arrays2support;
  ///@}
  
  /**@brief Map of types of arrays to support sets
   * @details This is of the same length as the vector `stats`.
   */
  MapVec_type< double, uint > keys2support;
  
  /**@name Functions to compute statistics
   * @details Arguments are recycled to save memory and computation.
   */
  ///@[{
  CounterVector<Array_Type,Data_Type> counters;
  Support<Array_Type,Data_Type>       support_fun;
  StatsCounter<Array_Type,Data_Type>  counter_fun;
  ///@}
  
  /**@brief Vector of the previously used parameters */
  std::vector< std::vector<double> > params_last;
  std::vector< double > normalizing_constants;
  bool first_calc_done = false;

  /**@brief Function to extract features of the array to be hash
  */
  std::function<std::vector<double>(const Array_Type &)> keygen = nullptr;  

  Model();
  Model(uint size_);
  ~Model() {};
  
  void store_psets();
  void set_keygen(std::function<std::vector<double>(const Array_Type &)> keygen_);
  
  /**@name Wrappers for the `CounterVector` member. 
   * @details These will add counters to the model, which are shared by the
   * support and the actual counter function.
   */
  ///@{
  void add_counter(Counter<Array_Type, Data_Type> & counter);
  void add_counter(Counter<Array_Type, Data_Type> * counter);
  void add_counter(
      Counter_fun_type<Array_Type,Data_Type> count_fun_,
      Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
      Data_Type *                            data_        = nullptr,
      bool                                   delete_data_ = false
  );
  ///@}
  
  /**@brief Adds an array to the support of not already included.
   * 
   * @param Array_ 
   */
  uint add_array(const Array_Type & Array_, bool force_new = false);
  
  
  /**@name Likelihood functions.
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
  
};

template <typename Array_Type, typename Data_Type>
inline Model<Array_Type,Data_Type>::Model() :
  stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
  target_stats(0u), arrays2support(0u), keys2support(0u), counters() 
{
  
  // Counters are shared
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline Model<Array_Type,Data_Type>::Model(uint size_) :
  stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
  target_stats(0u), arrays2support(0u), keys2support(0u), counters() 
{
  
  target_stats.reserve(size_);
  arrays2support.reserve(size_);

  // Counters are shared
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Model<Array_Type,Data_Type>::store_psets() {
  if (with_pset)
    throw std::logic_error("Powerset storage alreay activated.");
  with_pset = true;
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Model<Array_Type,Data_Type>::set_keygen(
    std::function<std::vector<double>(const Array_Type &)> keygen_
) {
  keygen = keygen_;
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Model<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> & counter
) {
  
  counters.add_counter(counter);
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Model<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * counter
) {
  
  counters.add_counter(counter);
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Model<Array_Type,Data_Type>::add_counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Data_Type *                            data_,
    bool                                   delete_data_
) {
  
  counters.add_counter(
      count_fun_,
      init_fun_,
      data_,
      delete_data_
  );
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline uint Model<Array_Type,Data_Type>::add_array(
  const Array_Type & Array_,
  bool force_new
) {
  
  // Array counts (target statistics)
  counter_fun.reset_array(&Array_);
  target_stats.push_back(counter_fun.count_all());
  
  // Checking with the hasher function: Is this present?
  if (keygen == nullptr)
    keygen = keygen_default<Array_Type>;
  
  // If the data hasn't been analyzed earlier, then we need to compute
  // the support
  std::vector< double > key = keygen(Array_);
  MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
  if (force_new | (locator == keys2support.end())) {
    
    // Adding to the map
    keys2support[key] = stats.size();
    n_arrays_per_stats.push_back(1u);       // How many elements now
    arrays2support.push_back(stats.size()); // Map of the array id to the support
    
    // Computing support using the counters included in the model
    support_fun.reset_array(&Array_);
    
    /** When computing with the powerset, we need to grow the corresponding
      * vectors on the fly */
    if (with_pset) {
      
      pset_arrays.resize(pset_arrays.size() + 1u);
      pset_stats.resize(pset_stats.size() + 1u);
      support_fun.calc(
        0u,
        true,
        &(pset_arrays[pset_arrays.size() - 1u]),
        &(pset_stats[pset_stats.size() - 1u])
      );
      
    } else
      support_fun.calc();
    
    stats.push_back(support_fun.get_counts());
    
    // Making room for the previous parameters. This will be used to check if
    // the normalizing constant has been updated or not.
    params_last.push_back(target_stats[0u]);
    normalizing_constants.push_back(0.0);
    
    return keys2support.size();
    
  }
  
  // Increasing the number of arrays in that stat
  ++n_arrays_per_stats[locator->second];
  
  // Adding the corresponding map
  arrays2support.push_back(locator->second);
  
  return locator->second;

}

template <typename Array_Type, typename Data_Type>
inline double Model<Array_Type,Data_Type>::likelihood(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
  
  // Checking if the index exists
  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");
  
  // Checking if we have updated the normalizing constant or not
  if (!first_calc_done || !vec_equal_approx(params, params_last[arrays2support[i]]) ) {
    
    first_calc_done = true;
    
    normalizing_constants[arrays2support[i]] = update_normalizing_constant(
      params, stats[arrays2support[i]]
    );
    
    params_last[arrays2support[i]] = params;
    
  }
  
  return likelihood_(
    target_stats[i],
    params,
    normalizing_constants[arrays2support[i]],
    as_log
  );
  
}

template <typename Array_Type, typename Data_Type>
inline double Model<Array_Type,Data_Type>::likelihood(
    const std::vector<double> & params,
    const Array_Type & Array_,
    bool as_log
) {
  
  std::vector< double > key = keygen(Array_);
  MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
  if (locator == keys2support.end()) 
    throw std::range_error("This type of array has not been included in the model.");
  
  return likelihood(
    params,
    locator->second,
    as_log
  );
  
}

template <typename Array_Type, typename Data_Type>
inline double Model<Array_Type,Data_Type>::likelihood_total(
    const std::vector<double> & params,
    bool as_log
) {
  
  for (uint i = 0u; i < params_last.size(); ++i) {
    if (!first_calc_done || !vec_equal_approx(params, params_last[i]) ) {
      
      first_calc_done = true;
      normalizing_constants[i] = update_normalizing_constant(
        params, stats[i]
      );
      
      params_last[i] = params;
      
    }
  }
  
  double res = 0.0;
  if (as_log) {
    for (uint i = 0; i < target_stats.size(); ++i) 
      res += vec_inner_prod(target_stats[i], params);
    
    for (uint i = 0u; i < params_last.size(); ++i) {
      res -= (std::log(normalizing_constants[i]) * this->n_arrays_per_stats[i]);
    }
  } else {
    
    res = 1.0;
    for (uint i = 0; i < target_stats.size(); ++i)
      res *= std::exp(vec_inner_prod(target_stats[i], params)) / 
        normalizing_constants[arrays2support[i]];
    
  }
  
  return res;
  
}

#endif
