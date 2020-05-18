#include <vector>
#include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"

#ifndef LBARRAY_BONES_HPP 
#define LBARRAY_BONES_HPP 1

/***
 * ! Users can a list of functions that can be used with this.
 * ! The baseline set of arguments is a pointer to a binary array.
 */ 
template <typename Cell_Type>
class StatsCounter {
public:
  
  // Should receive an array
  const BArray<Cell_Type> * Array;
  BArray<Cell_Type> EmptyArray;
  std::vector< double > current_stats;
  std::vector< double > change_stats;
   
  // We will save the data here
  StatsDB * stats;
  std::vector< Counter<Cell_Type> > counters;
  
  // Functions
  StatsCounter(const BArray<Cell_Type> * data, StatsDB * stats_) :
    counters(0u), Array(data), EmptyArray(data->N, data->M), stats(stats_) {
    // Copying the information
    EmptyArray.meta = Array->meta;
    return;
  }
  
  void add_counter(Counter<Cell_Type> f_);
  
  /***
   * ! This function recurses through the entries of Array and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_init(uint i, uint j);
  void count_current(uint i, uint j);
  void count_all();
  
};

template <typename Cell_Type>
inline void StatsCounter<Cell_Type>::add_counter(Counter<Cell_Type> f_){
  counters.push_back(f_);
  return;
}

template <typename Cell_Type>
inline void StatsCounter<Cell_Type>::count_init(uint i, uint j) {
  // Iterating through the functions, and updating the set of
  // statistics.
  current_stats.resize(counters.size(), 0.0);
  change_stats.resize(counters.size(), 0.0);
  for (uint n = 0u; n < counters.size(); ++n) 
    current_stats.at(n) = counters.at(n).init(&EmptyArray, i, j);
  
  return;
}

template <typename Cell_Type>
inline void StatsCounter<Cell_Type>::count_current(uint i, uint j) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  for (uint n = 0u; n < counters.size(); ++n) {
    change_stats.at(n)   = counters.at(n).count(&EmptyArray, i, j);
    current_stats.at(n) += change_stats.at(n);
  }

  return;
  
}

template <typename Cell_Type>
inline void StatsCounter<Cell_Type>::count_all() {
  
  // Initializing the counter on the empty array
  count_init(0u, 0u);
  
  // Setting it to zero.
  EmptyArray.clear();
  
  // Start iterating through the data
  uint N = Array->N;
  uint M = Array->M;
  for (uint row = 0; row < N; ++row) {
    
    // Any element?
    if (A_ROW(row).size() == 0u)
      continue;
    
    // If there's one, then update the statistic, by iterating
    for (auto col = A_ROW(row).begin(); col != A_ROW(row).end(); ++col) {
      
      // Adding a cell
      EmptyArray.insert_cell(row, col->first);
      
      // Computing the change statistics
      count_current(row, col->first);
     
    } 
    
  }
  
  // Adding to the sufficient statistics
  stats->add(current_stats);
  return;
  
}

/***
 * ! Compute the support of a model by iterating through all
 * ! possible combinations, as we would do in a powerset.
 */

template <typename Cell_Type>
class Support {
public:
  
  const BArray<Cell_Type> * Array;
  BArray<Cell_Type> EmptyArray;
  StatsDB support;
  std::vector< Counter<Cell_Type> > counters;
  std::vector< double > current_stats;
  
  uint N, M;
  bool initialized = false;
  
  Support(const BArray<Cell_Type> * Array_) : Array(Array_), EmptyArray(Array_->N, Array_->M),
    N(Array_->N), M(Array_->M) {
    
    EmptyArray.meta = Array->meta;
    return;
    
  };
  Support(uint N_, uint M_) : EmptyArray(N_, M_) ,N(N_), M(M_) {};
  ~Support() {};

  void add_counter(Counter<Cell_Type> f_);  
  
  void calc(uint pos = 0u) {
    

    // Getting the location 
    uint i = (int) pos % (int) N;
    uint j = floor((int) pos / (int) N);
    
    // No self ties, go to the next step and return.
    if (i == j)
      return calc(pos + 1u);
    
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
        current_stats.at(n) = counters.at(n).init(&EmptyArray, i, j);
      
      // Adding to the overall count
      support.add(current_stats);
      
    }

    // We will pass it to the next step, if the iteration makes sense.
    calc(pos + 1u);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    
    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(i, j, true, false, false);
    
    // Counting
    std::vector< double > change_stats(counters.size());
    for (uint n = 0u; n < counters.size(); ++n) {
      change_stats.at(n) = counters.at(n).count(&EmptyArray, i, j);
      current_stats.at(n) += change_stats.at(n);
    }
    
    // Adding to the overall count
    support.add(current_stats);
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc(pos + 1);

    // We need to restore the state of the cell
    EmptyArray.rm_cell(i, j, false, false);
    for (uint n = 0u; n < counters.size(); ++n) 
      current_stats.at(n) -= change_stats.at(n);
     
    
    return;
    
  }
  
};

template <typename Cell_Type>
inline void Support<Cell_Type>::add_counter(Counter<Cell_Type> f_) {
  counters.push_back(f_);
  return;
}

/***
 * Some example functions
 */

#endif
