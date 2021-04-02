#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"

#ifndef BARRY_STATSCOUNTER_BONES_HPP 
#define BARRY_STATSCOUNTER_BONES_HPP 1

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
  // std::vector< double > change_stats;
   
  // We will save the data here
  Counters<Array_Type,Data_Type> * counters;
  bool                                  counter_deleted  = false;
  
  /**
   * @brief Creator of a `StatsCounter`
   * 
   * @param Array_ A const pointer to a `BArray`.
   */
  StatsCounter(const Array_Type * Array_) :
    Array(Array_), EmptyArray(*Array_),
    counters(new Counters<Array_Type,Data_Type>()) {
    
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
    counters(new Counters<Array_Type,Data_Type>()) {};
  ~StatsCounter();
  
  /**@brief Changes the reference array for the counting.
   * 
   * @param Array_ A pointer to an array of class `Array_Type`.
   */
  void reset_array(const Array_Type * Array_);
  
  void add_counter(Counter<Array_Type,Data_Type> * f_);
  void add_counter(Counter<Array_Type,Data_Type> f_);
  void set_counters(Counters<Array_Type,Data_Type> * counters_);
  
  /***
   * ! This function recurses through the entries of `Array` and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_init(uint i, uint j);
  void count_current(uint i, uint j);
  std::vector< double > count_all();
  
};



#endif
