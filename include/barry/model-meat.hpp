#include "model-bones.hpp"

#ifndef BARRY_MODEL_MEAT_HPP 
#define BARRY_MODEL_MEAT_HPP 1

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::Model() :
  stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
  target_stats(0u), arrays2support(0u), keys2support(0u), counters(), rules()
{  

  // Counters are shared
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  // Rules are shared
  support_fun.set_rules(&rules);

  // Checking with the hasher function: Is this present?
  keygen = keygen_default<Array_Type>;
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::Model(uint size_) :
  stats(0u), n_arrays_per_stats(0u), pset_arrays(0u), pset_stats(0u),
  target_stats(0u), arrays2support(0u), keys2support(0u), counters(), rules()
{
  
  target_stats.reserve(size_);
  arrays2support.reserve(size_);

  // Counters are shared
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  // Rules are shared
  support_fun.set_rules(&rules);
  
  // Checking with the hasher function: Is this present?
  keygen = keygen_default<Array_Type>;
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::Model(
  const Model<Array_Type,Data_Counter_Type,Data_Rule_Type> & Model_) : 
  stats(Model_.stats),
  n_arrays_per_stats(Model_.n_arrays_per_stats),
  pset_arrays(Model_.pset_arrays),
  pset_stats(Model_.pset_stats),
  target_stats(Model_.target_stats),
  arrays2support(Model_.arrays2support),
  keys2support(Model_.keys2support),
  counters(Model_.counters),
  rules(Model_.rules),
  params_last(Model_.params_last),
  normalizing_constants(Model_.normalizing_constants),
  first_calc_done(Model_.first_calc_done)
  {
  
  // Counters are shared
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  // Rules are shared
  support_fun.set_rules(&rules);

  keygen = Model_.keygen;
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type> & Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::operator=(
  const Model<Array_Type,Data_Counter_Type,Data_Rule_Type> & Model_
) {
  
  // Clearing
  if (this != &Model_) {
    
    stats                 = Model_.stats;
    n_arrays_per_stats    = Model_.n_arrays_per_stats;
    pset_arrays           = Model_.pset_arrays;
    pset_stats            = Model_.pset_stats;
    target_stats          = Model_.target_stats;
    arrays2support        = Model_.arrays2support;
    keys2support          = Model_.keys2support;
    counters              = Model_.counters;
    rules                 = Model_.rule;
    params_last           = Model_.params_last;
    normalizing_constants = Model_.normalizing_constants;
    first_calc_done       = Model_.first_calc_done;

    // Counters are shared
    support_fun.set_counters(&counters);
    counter_fun.set_counters(&counters);
    
    // Rules are shared
    support_fun.set_rules(&rules);
    
  }
    
  return *this;
  
}


template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::store_psets() {
  if (with_pset)
    throw std::logic_error("Powerset storage alreay activated.");
  with_pset = true;
  return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::set_keygen(
    std::function<std::vector<double>(const Array_Type &)> keygen_
) {
  keygen = keygen_;
  return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_counter(
    Counter<Array_Type, Data_Counter_Type> & counter
) {
  
  counters.add_counter(counter);
  return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_counter(
    Counter<Array_Type, Data_Counter_Type> * counter
) {
  
  counters.add_counter(counter);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_counter(
    Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_,
    Data_Counter_Type *                            data_,
    bool                                           delete_data_
) {
  
  counters.add_counter(
      count_fun_,
      init_fun_,
      data_,
      delete_data_
  );
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type>::set_counters(
    Counters<Array_Type,Data_Counter_Type> * counters_
) {
  
  this->counters = *counters_;
  support_fun.set_counters(&counters);
  counter_fun.set_counters(&counters);
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> & rule
) {
  
  rules.add_rule(rule);
  return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> * rule
) {
  
  rules.add_rule(rule);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_rule(
    Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
    Data_Rule_Type *                         data_,
    bool                                     delete_data_
) {
  
  rules.add_rule(
    rule_fun_,
    data_,
    delete_data_
  );
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type>::set_rules(
    Rules<Array_Type,Data_Rule_Type> * rules_
) {

  this->rules = *rules_;
  support_fun.set_rules(&rules);
  return;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline uint Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::add_array(
  const Array_Type & Array_,
  bool force_new
) {
  
  // Array counts (target statistics)
  counter_fun.reset_array(&Array_);
  target_stats.push_back(counter_fun.count_all());
  
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
      
      // Making space for storing the support
      pset_arrays.resize(pset_arrays.size() + 1u);
      pset_stats.resize(pset_stats.size() + 1u);
      pset_probs.resize(pset_probs.size() + 1u);
      
      try {
        
        support_fun.calc(
          &(pset_arrays[pset_arrays.size() - 1u]),
          &(pset_stats[pset_stats.size() - 1u])
        );
        
      } catch (const std::exception& e) {
        
        std::cout << "A problem ocurred while trying to add the array (and recording the powerset)." <<
          "with error: " << e.what();
        throw std::logic_error("");
        
      }
      
    } else {
      
      try {
        support_fun.calc();
      } catch (const std::exception& e) {
        
        std::cout << "A problem ocurred while trying to add the array." <<
          "with error: " << e.what();
        throw std::logic_error("");
        
      }
      
    }
    
    stats.push_back(support_fun.get_counts());
    
    // Making room for the previous parameters. This will be used to check if
    // the normalizing constant has been updated or not.
    params_last.push_back(target_stats[0u]);
    normalizing_constants.push_back(0.0);
    first_calc_done.push_back(false);
    
    return arrays2support.size() - 1u;
    
  }
  
  // Increasing the number of arrays in that stat
  ++n_arrays_per_stats[locator->second];
  
  // Adding the corresponding map
  arrays2support.push_back(locator->second);
  
  return arrays2support.size() - 1u;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline double Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::likelihood(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
  
  // Checking if the index exists
  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");
  
  // Checking if we have updated the normalizing constant or not
  if (!first_calc_done[arrays2support[i]] || !vec_equal_approx(params, params_last[arrays2support[i]]) ) {
    
    first_calc_done[arrays2support[i]] = true;
    
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

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline double Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::likelihood(
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

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline double Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::likelihood(
    const std::vector<double> & params,
    const std::vector<double> & target_,
    const uint & i,
    bool as_log
) {
  
  // Checking if the index exists
  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");
  
  // Checking if we have updated the normalizing constant or not
  if (!first_calc_done[arrays2support[i]] || !vec_equal_approx(params, params_last[arrays2support[i]]) ) {
    
    first_calc_done[arrays2support[i]] = true;
    
    normalizing_constants[arrays2support[i]] = update_normalizing_constant(
      params, stats[arrays2support[i]]
    );
    
    params_last[arrays2support[i]] = params;
    
  }
  
  return likelihood_(
    target_,
    params,
    normalizing_constants[arrays2support[i]],
    as_log
  );
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline double Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::likelihood_total(
    const std::vector<double> & params,
    bool as_log
) {
  
  for (uint i = 0u; i < params_last.size(); ++i) {
    if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) ) {
      
      first_calc_done[i] = true;
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

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline double Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::get_norm_const(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
  
  // Checking if the index exists
  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");
  
  // Checking if we have updated the normalizing constant or not
  if (!first_calc_done[arrays2support[i]] || !vec_equal_approx(params, params_last[arrays2support[i]]) ) {
    
    first_calc_done[arrays2support[i]] = true;
    
    normalizing_constants[arrays2support[i]] = update_normalizing_constant(
      params, stats[arrays2support[i]]
    );
    
    params_last[arrays2support[i]] = params;
    
  }
  
  return as_log ? 
    std::log(normalizing_constants[arrays2support[i]]) :
    normalizing_constants[arrays2support[i]]
    ;
  
}


template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline const std::vector< Array_Type > * Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::get_pset(
    const uint & i
) {

  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");


  return &pset_arrays[arrays2support[i]];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline const std::vector< std::vector<double> > * Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::get_stats(
    const uint & i
) {

  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");

  return &pset_stats[arrays2support[i]];

}

template<typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::print_stats(uint i) const {
  
  if (i >= arrays2support.size())
    throw std::range_error("The requested support is out of range");

  for (uint l = 0u; l < stats[arrays2support[i]].size(); ++l) {
    printf("% 5i ", l);
    std::cout << "counts " << stats[arrays2support[i]][l].second << " motif: ";
    for (unsigned int k = 0u; k < stats[arrays2support[i]][l].first.size(); ++k) {
      std::cout << stats[arrays2support[i]][l].first[k] << ", ";
    }
    std::cout << std::endl;
  }
  
  return;
  
}

template<typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline uint Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::size() const {
  return this->target_stats.size();
}

template<typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline uint Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::size_unique() const {
  return this->stats.size();
} 
  
template<typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Array_Type Model<Array_Type,Data_Counter_Type,Data_Rule_Type>::sample(
    const unsigned int & i,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // Getting the index
    unsigned int a = arrays2support[i];
    // if (pset_probs.at(a).size() == 0u) {
    //   pset_probs.at(a).resize(pset_arrays.at(a).size(), 0u);
    // }
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    // Updating until reach above
    unsigned int j = 0u;
    while (cumprob < r) {

        cumprob += this->likelihood(params, this->pset_stats[a][j], i, false);
        ++j;
    }

    return this->pset_arrays[a][j-1u];   

}

#endif
