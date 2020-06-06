#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_COUNTERS_H
#define BARRAY_COUNTERS_H

#define COUNTER_FUNCTION_DEFAULT(a) template <typename Array_Type = BArray<>, typename Data_Type = bool> \
  inline double (a) (Array_Type * Array, uint i, uint j, Data_Type * data)\
    
#define COUNTER_FUNCTION(a) template <typename Array_Type, typename Data_Type> \
  inline double (a) (Array_Type * Array, uint i, uint j, Data_Type * data)\
    

/**
 * @brief A counter function based on change statistics.
 * 
 * This class is used by `CountStats` and `StatsCounter` as a way to count
 * statistics using change statistics.
 */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Counter {
public:
  
  Data_Type * data = nullptr;
  Counter_fun_type<Array_Type,Data_Type> count_fun;
  Counter_fun_type<Array_Type,Data_Type> init_fun;

  /***
   * ! Initializers
   */
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  
  /**
   * @brief Creator passing only a counter function
   * @param count_fun The main counter function.
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_
    ) :
    count_fun(count_fun_), init_fun(nullptr) {};
  /**
   * @brief Creator passing a counter and an initializer
   * 
   * @param count_fun The main counter function.
   * @param init_fun The initializer function can also be used to check if the
   *  `BArray` has the required metadata (e.g. is it a directed graph?).
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_
    ): count_fun(count_fun_), init_fun(init_fun_) {};
  
  ~Counter() {};
  
  /***
   * ! Main functions.
   */
  double count(Array_Type * Array, uint i, uint j);
  double init(Array_Type * Array, uint i, uint j);
};

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::count(
    Array_Type * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j, data);
}

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::init(
    Array_Type * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j, data);
}

#endif
