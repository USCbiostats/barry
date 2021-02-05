// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "rules-bones.hpp"

#ifndef LBARRAY_POWERSET_HPP 
#define LBARRAY_POWERSET_HPP 1

/* A cell is a fundamental type that holds information about a cell
 * for now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
template <typename Array_Type = BArray<>, typename Data_Rule_Type = bool> 
class PowerSet {
  
public:
  std::vector< Array_Type > data;
  Rules<Array_Type,Data_Rule_Type> rules;
  uint N, M;
  uint counter;
  
  /***
   * ! Creator and destructor functions
   */
  PowerSet() : data(0u), rules(), M(0u), N(0u) {};
  PowerSet(uint N_, uint M_) : data(0u), rules(), M(M_), N(N_), counter(0u) {};
  ~PowerSet() {};
  
  void calc(uint i = 0u, Array_Type * a0 = nullptr);
  void reset(uint N_, uint M_);
  
  /**@name Wrappers for the `Rules` member. 
   * @details These will add rules to the model, which are shared by the
   * support and the actual counter function.
   */
  ///@{
  void add_rule(Rule<Array_Type, Data_Rule_Type> & rule);
  void add_rule(Rule<Array_Type, Data_Rule_Type> * rule);
  void add_rule(
      Rule_fun_type<Array_Type,Data_Rule_Type> count_fun_,
      Data_Rule_Type * data_ = nullptr,
      bool delete_data_ = false
  );
  ///@}
  
  /**
   * Getter functions
   */
  const std::vector< Array_Type > * get_data_ptr() const {return &data;};
  std::vector< Array_Type > get_data() const {return data;};
  
};

/***
 * Function to generate the powerset of the 
 */
template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::calc(
    uint i,
    Array_Type * a0
  ) {
  
  // If it is size 0, then need to recalculate the size, and initialize
  // the first datum
  if (a0 == nullptr) {
    
    this->data.reserve(pow(2.0, this->N * this->M));
    this->data.push_back(Array_Type(this->N, this->M));
    a0 = &this->data[0u];
    
  }
  
  // Here is the deal
  this->data.push_back(*a0); // Making it one
  this->data[++this->counter].insert_cell(i, true, false, false);
  
  // If the sum is even, we increase i, otherwise j
  if ((i + 1) < (this->N * this->M)) {
    calc(i + 1u, &this->data[this->counter]);
    calc(i + 1u, a0);
  }
  
  return;
  
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::reset(
    uint N_,
    uint M_
  ) {
  
  data.empty();
  N = N_, M = M_, counter = 0u;
  
  return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> & rule
) {
  
  rules.add_rule(rule);
  return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> * rule
) {
  
  rules.add_rule(rule);
  return;
  
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
    Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
    Data_Rule_Type *                            data_,
    bool                                   delete_data_
) {
  
  rules.add_rule(
    rule_fun_,
    data_,
    delete_data_
  );
  
  return;
  
}

#endif
