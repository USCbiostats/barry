#include "counters-bones.hpp"

#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

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

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type> * CounterVector<Array_Type,Data_Type>::operator[](uint idx) {
  return data[idx];
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::add_counter(
  Counter<Array_Type, Data_Type> & counter
) {
  
  to_be_deleted.push_back(data.size());
  data.push_back(new Counter<Array_Type, Data_Type>(counter));
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * counter
) {
  
  data.push_back(counter);
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::add_counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Data_Type *                            data_,
    bool                                   delete_data_
) {
 
  /* We still need to delete the counter since we are using the 'new' operator.
   * Yet, the actual data may not need to be deleted.
   */
  to_be_deleted.push_back(data.size());
  
  data.push_back(new Counter<Array_Type,Data_Type>(
    count_fun_,
    init_fun_,
    data_,
    delete_data_
  ));
 
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::clear() {
  
  for (auto iter = to_be_deleted.begin(); iter != to_be_deleted.end(); ++iter)
    delete data[*iter];
  
  to_be_deleted.clear();
  
  return;
  
}

#endif 