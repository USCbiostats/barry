// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"

#ifndef SUPPORT_HPP 
#define SUPPORT_HPP 1

/**@brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 */ 
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Support {
public:
  
  const Array_Type * Array;
  Data_Type * data = nullptr;
  Array_Type EmptyArray;
  StatsDB support;
  std::vector< Counter<Array_Type, Data_Type> > counters;
  std::vector< double > current_stats;
  
  uint N, M;
  bool initialized = false;
  
  Support(const Array_Type * Array_) : Array(Array_), EmptyArray(Array_->N, Array_->M),
    N(Array_->N), M(Array_->M) {
    
    if (Array_->data != nullptr) 
      EmptyArray.data = Array_->data;
    
    return;
    
  };
  Support(uint N_, uint M_) : EmptyArray(N_, M_) ,N(N_), M(M_){};
  ~Support() {};
  
  /**@brief Resets the support calculator
   * 
   * If needed, the counters of a support object can be reused.
   * 
   * @param Array_ New array over which the support will be computed.
   */
  void reset();
  void reset(const Array_Type * Array_);
  
  void add_counter(Counter<Array_Type, Data_Type> & f_);  
  
  /**@brief Computes the entire support
   * 
   * Not to be used by the user. Sets the starting point in the array
   * (column-major).
   * 
   * @param diag When `true`, includes the diagonal (i=j) in the counts.
   * 
   * @param array_bank If specified, the counter will add to the vector each 
   * possible state of the array, as it counts.
   * 
   * @param stats_bank If specified, the counter will add to the vector each
   * possible set of statistics, as it counts.
   * 
   */
  void calc(
      uint pos = 0u,
      const bool & diag = false,
      std::vector< Array_Type > * array_bank = nullptr,
      std::vector< std::vector< double > > * stats_bank = nullptr
    );
  
  
};

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::reset() {
  
  support.clear();
  initialized = false;
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::reset(const Array_Type * Array_) {
  
  support.clear();
  initialized = false;
  Array = Array_;
  N = Array_->N;
  M = Array_->M;
  EmptyArray.resize(N, M);
  if (Array_->data != nullptr) 
    EmptyArray.data = Array_->data;
  
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::calc(
    uint pos, const bool & diag,
    std::vector< Array_Type > * array_bank,
    std::vector< std::vector< double > > * stats_bank
  ) {
  
  // Getting the location 
  uint i = (int) pos % (int) N;
  uint j = floor((int) pos / (int) N);
  
  // No self ties, go to the next step and return.
  if (!diag && (i == j))
    return calc(pos + 1u, diag, array_bank, stats_bank);
  
  // If reached the end, also return
  if ((i >= N) || (j >= M))
    return;
  
  // Initializing
  if (!initialized) {
    
    // Initializing
    initialized = true;
    EmptyArray.clear();
    current_stats.resize(counters.size());
    
    for (uint n = 0u; n < counters.size(); ++n) 
      current_stats[n] = counters[n].init(&EmptyArray, i, j);
    
    // Adding to the overall count
    support.add(current_stats);
    
    if (array_bank != nullptr)
      array_bank->push_back(EmptyArray);
    
    if (stats_bank != nullptr)
      stats_bank->push_back(current_stats);
    
  }
  
  // We will pass it to the next step, if the iteration makes sense.
  calc(pos + 1u, diag, array_bank, stats_bank);
  
  // Once we have returned, everything will be back as it used to be, so we
  // treat the data as if nothing has changed.
  
  // Toggle the cell (we will toggle it back after calling the counter)
  EmptyArray.insert_cell(i, j, true, false, false);
  
  // Counting
  std::vector< double > change_stats(counters.size());
  for (uint n = 0u; n < counters.size(); ++n) {
    change_stats[n] = counters[n].count(&EmptyArray, i, j);
    current_stats[n] += change_stats[n];
  }
  
  // Adding to the overall count
  support.add(current_stats);
  
  // Need to save?
  if (array_bank != nullptr)
    array_bank->push_back(EmptyArray);
  
  if (stats_bank != nullptr)
    stats_bank->push_back(current_stats);
  
  // Again, we only pass it to the next level iff the next level is not
  // passed the last step.
  calc(pos + 1, diag, array_bank, stats_bank);
  
  // We need to restore the state of the cell
  EmptyArray.rm_cell(i, j, false, false);
  for (uint n = 0u; n < counters.size(); ++n) 
    current_stats[n] -= change_stats[n];
  
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> & f_
  ) {
  counters.push_back(f_);
  return;
}


#endif
