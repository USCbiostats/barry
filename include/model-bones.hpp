// #include <vector>
// #include <unordered_map>
#include "barray-bones.hpp"
#include "support.hpp"
#include "statscounter.hpp"

#ifndef MODEL_BONES_HPP 
#define MODEL_BONES_HPP 1

inline double likelihood(
    const std::vector< double > & target_stats,
    const std::vector< double > & params,
    const Counts_type & support,
    bool log_ = false
) {
  
  if (target_stats.size() != params.size())
    throw std::length_error("-target_stats- and -params- should have the same length.");
  
  double numerator=0.0, denominator=0.0;
  
  // Computing the numerator
  for (unsigned int j = 0u; j < target_stats.size(); ++j)
    numerator += target_stats[j] * params[j];
  if (!log_)
    numerator = exp(numerator);
  
  // Computing denominator
  for (unsigned int n = 0u; n < support.size(); ++n) {
    
    double tmp = 0.0;
    
    for (unsigned int j = 0u; j < params.size(); ++j)
      tmp += support[n].first[j] * params[j];
    
    tmp = exp(tmp);
    denominator += tmp * support[n].second;
  }
  
  if (log_)
    return numerator - log(denominator);
  
  return numerator/denominator;
  
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
  
  /**@name Functions to compute statistics
   * @details Arguments are recycled to save memory and computation.
   */
  ///@[{
  Support<Array_Type,Data_Type>       support_fun;
  StatsCounter<Array_Type,Data_Type>  counter_fun;
  CounterVector<Array_Type,Data_Type> counters;
  ///@}
  
  /**@brief */
  std::vector< Counts_type >          stats;
  std::vector< uint >                 n_arrays_per_stats;
  std::vector< std::vector< double >> target_stats;
  std::vector< uint >                 arrays2support;
  
  /**@brief Vector of the previously used parameters */
  std::vector< double > params_last;
  
  /**@brief Map of types of arrays to support sets
   * @details This is of the same length as the vector `stats`.
   */
  MapVec_type< double, uint > keys2support;

  /**@brief Function to extract features of the array to be hash
  */
  std::function<std::vector<double>(const Array_Type &)> keygen = nullptr;  

  Model() {
    
    // Counters are shared
    support_fun.set_counters(&counters);
    counter_fun.set_counters(&counters);
    
    return;
  };
  ~Model() {};
  
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
  uint add_array(const Array_Type & Array_);
  double likelihood(const Array_Type & Array_, bool add = false);
  
};

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
  const Array_Type & Array_
) {
  
  // Checking with the hasher function: Is this present?
  if (keygen == nullptr)
    keygen = keygen_default;
  
  // If the data hasn't been analyzed earlier, then we need to compute
  // the support
  std::vector< double > key = keygen(Array_);
  MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
  if (locator == keys2support.end()) {
    
    // Adding to the map
    keys2support[key] = keys2support.size();
    
    // Array counts
    counter_fun.reset_array(&Array_);
    target_stats.push_back(counter_fun.count_all());
    
    // Computing support using the counters included in the model
    support_fun.reset(&Array_);
    support_fun.calc();
    stats.push_back(support_fun.get_counts());
    
    return keys2support.size();
    
  } 
  
  arrays2support.push_back(locator->second);
  
  return locator->second;

}

#endif
