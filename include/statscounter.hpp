// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"

#ifndef STATSCOUNTER_HPP 
#define STATSCOUNTER_HPP 1

/**
 * @brief Count stats for a single Array.
 * 
 * Users can a list of functions that can be used with this. The baseline set of
 * arguments is a pointer to a binary array and a dataset to add the counts to.
 */ 
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class StatsCounter {
public:
  
  // Should receive an array
  const Array_Type * Array;
  Array_Type EmptyArray;
  std::vector< double > current_stats;
  std::vector< double > change_stats;
   
  // We will save the data here
  std::vector< Counter<Array_Type,Data_Type> > counters;
  
  /**
   * @brief Creator of a `StatsCounter`
   * 
   * @param data A const pointer to a `BArray`.
   * @param Stats_ A pointer to a dataset of stats `StatsDB`.
   */
  StatsCounter(const Array_Type * data) :
    counters(0u), Array(data), EmptyArray(data->N, data->M) {
    
    if (data->data != nullptr)
      EmptyArray.data = data->data;
    
    return;
  }
  
  void add_counter(Counter<Array_Type,Data_Type> f_);
  
  /***
   * ! This function recurses through the entries of Array and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_init(uint i, uint j);
  void count_current(uint i, uint j);
  std::vector< double > count_all();
  
};

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type,Data_Type> f_
  ) {
  
  counters.push_back(f_);
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_init(
    uint i,
    uint j
  ) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  current_stats.resize(counters.size(), 0.0);
  change_stats.resize(counters.size(), 0.0);
  for (uint n = 0u; n < counters.size(); ++n) 
    current_stats[n] = counters[n].init(&EmptyArray, i, j);
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_current(
    uint i,
    uint j
  ) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  for (uint n = 0u; n < counters.size(); ++n) {
    change_stats[n]   = counters[n].count(&EmptyArray, i, j);
    current_stats[n] += change_stats[n];
  }

  return;
  
}

template <typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type, Data_Type>::count_all() {
  
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
  return current_stats;
  
}

#endif
