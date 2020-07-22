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
  CounterVector<Array_Type,Data_Type> * counters;
  bool                                  counter_deleted  = false;
  
  /**
   * @brief Creator of a `StatsCounter`
   * 
   * @param data A const pointer to a `BArray`.
   * @param Stats_ A pointer to a dataset of stats `StatsDB`.
   */
  StatsCounter(const Array_Type * Array_) :
    Array(Array_), EmptyArray(*Array_),
    counters(new CounterVector<Array_Type,Data_Type>()) {
    
    // We are removing the entries without freeing the memory. This should
    // make the insertion faster.
    EmptyArray.clear(false);
    
    return;
  }
  
  /**@brief Can be created without setting the array.
   * 
   */
  StatsCounter() :
    Array(nullptr), EmptyArray(0u,0u),
    counters(new CounterVector<Array_Type,Data_Type>()) {};
  ~StatsCounter();
  
  /**@brief Changes the reference array for the counting.
   * 
   * @param Array_ A pointer to an array of class `Array_Type`.
   */
  void reset_array(const Array_Type * Array_);
  
  void add_counter(Counter<Array_Type,Data_Type> * f_);
  void add_counter(Counter<Array_Type,Data_Type> f_);
  void set_counters(CounterVector<Array_Type,Data_Type> * counters_);
  
  /***
   * ! This function recurses through the entries of `Array` and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_init(uint i, uint j);
  void count_current(uint i, uint j);
  std::vector< double > count_all();
  
};


template <typename Array_Type, typename Data_Type>
inline StatsCounter<Array_Type,Data_Type>::~StatsCounter() {
  if (!counter_deleted)
    delete counters;
  return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::reset_array(
  const Array_Type * Array_
) {
  
  Array = Array_;
  EmptyArray = *Array_;
  EmptyArray.delete_data = false;
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type,Data_Type> * f_
  ) {
  
  counters->add_counter(f_);
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type,Data_Type> f_
) {
  
  counters->add_counter(f_);
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::set_counters(
    CounterVector<Array_Type,Data_Type> * counters_
) {
  
  // Cleaning up before replacing the memory
  if (!counter_deleted)
    delete counters;
  counter_deleted = true;
  counters = counters_;
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_init(
    uint i,
    uint j
  ) {
  
  // Do we have any counter?
  if (counters->size() == 0u)
    throw std::logic_error("No counters added: Cannot count without knowning what to count!");
  
  // Iterating through the functions, and updating the set of
  // statistics.
  current_stats.resize(counters->size(), 0.0);
  change_stats.resize(counters->size(), 0.0);
  for (uint n = 0u; n < counters->size(); ++n) 
    current_stats[n] = counters->operator[](n)->init(&EmptyArray, i, j);
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_current(
    uint i,
    uint j
  ) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  for (uint n = 0u; n < counters->size(); ++n) {
    change_stats[n]   = counters->operator[](n)->count(&EmptyArray, i, j);
    current_stats[n] += change_stats[n];
  }

  return;
  
}

template <typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type, Data_Type>::count_all() {
  
  // Initializing the counter on the empty array
  count_init(0u, 0u);
  
  // Setting it to zero.
  EmptyArray.clear(false);
  
  // Start iterating through the data
  for (uint row = 0; row < Array->N; ++row) {
    
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
