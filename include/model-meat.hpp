#include "model-bones.hpp"

#ifndef MODEL_MEAT_HPP 
#define MODEL_MEAT_HPP 1

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
    
    return keys2support.size() - 1u;
    
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
