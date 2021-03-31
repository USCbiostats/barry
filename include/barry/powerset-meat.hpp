#include "powerset-bones.hpp"

#ifndef BARRY_POWERSET_MEAT_HPP
#define BARRY_POWERSET_MEAT_HPP 1

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::PowerSet(
  const Array_Type & array
) : EmptyArray(array), data(0u),
    rules(new Rules<Array_Type,Data_Rule_Type>()), N(array.N), M(array.M) {

}

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::~PowerSet() {
  if (!this->rules_deleted)
    delete rules;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::init_support() {
  
  // Computing the locations
  coordinates_free.clear();
  coordinates_locked.clear();
  rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);
  
  // Computing initial statistics
  if (EmptyArray.nnozero() > 0u) {
    for (uint i = 0u; i < coordinates_free.size(); ++i) 
      EmptyArray.rm_cell(coordinates_free[i].first, coordinates_free[i].second, false, true);
  }

  // EmptyArray.clear(true);
  // EmptyArray.reserve();
  
  // Resizing support
  data.reserve(pow(2.0, coordinates_free.size())); 

  // Adding the empty array to the set
  data.push_back(EmptyArray);
  
  return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc_backend(uint pos) {
  
  // Did we reached the end??
  if (pos >= coordinates_free.size())
    return;
      
  // We will pass it to the next step, if the iteration makes sense.
  calc_backend(pos + 1u);
    
  // Toggle the cell (we will toggle it back after calling the counter)
  EmptyArray.insert_cell(
    coordinates_free[pos].first,
    coordinates_free[pos].second,
    EmptyArray.Cell_default , false, false
    );

  data.push_back(EmptyArray);
  
  // Again, we only pass it to the next level iff the next level is not
  // passed the last step.
  calc_backend(pos + 1u);
  
  // We need to restore the state of the cell
  EmptyArray.rm_cell(
    coordinates_free[pos].first,
    coordinates_free[pos].second,
    false, false
    );  
  
  return;
  
}


/***
 * Function to generate the powerset of the 
 */
template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc() {

  // Generating sequence
  this->init_support();

  // Recursive function to count
  calc_backend(0u);

  return;
  
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::reset(
    uint N_,
    uint M_
  ) {
  
  data.empty();
  N = N_, M = M_;
  
  return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> & rule
) {
  
  rules->add_rule(rule);
  return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> * rule
) {
  
  rules->add_rule(rule);
  return;
  
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
    Data_Rule_Type *                            data_,
    bool                                   delete_data_
) {
  
  rules->add_rule(
    rule_fun_,
    data_,
    delete_data_
  );
  
  return;
  
}

#endif