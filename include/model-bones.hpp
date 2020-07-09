// #include <vector>
// #include <unordered_map>
#include "barray-bones.hpp"
#include "support.hpp"

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

template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Model {

public:
  
  std::vector< Support<Array_Type,Data_Type> > supports;
  std::unordered_map< std::vector< double > , uint, vecHasher<double> > keys;

  ~Model() {};
};

#endif
