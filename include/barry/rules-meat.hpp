// #include <vector>
// #include <stdexcept>

#include "rules-bones.hpp"

#ifndef BARRY_RULES_MEAT_HPP
#define BARRY_RULES_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline Rules<Array_Type,Data_Type>::Rules(
  const Rules<Array_Type,Data_Type> & rules_
) {

  // Checking which need to be deleted
  std::vector< bool > tobedeleted(rules_.size(), false);
  for (auto i = rules_.to_be_deleted.begin(); i != rules_.to_be_deleted.end(); ++i)
    tobedeleted[*i] = true;

  // Copy all rules, if a rule is tagged as 
  // to be deleted, then copy the value
  for (auto i = 0u; i != rules_.size(); ++i) {
    if (tobedeleted[i]) {
      this->add_rule(*rules_.data[i]);
    } else {
      this->add_rule(rules_.data[i]);
    } 
  }

  return;

};

template <typename Array_Type, typename Data_Type>
Rules<Array_Type,Data_Type> Rules<Array_Type,Data_Type>::operator=(
  const Rules<Array_Type,Data_Type> & rules_
) {

  if (this != &rules_) {

    // Checking which need to be deleted
    std::vector< bool > tobedeleted(rules_.size(), false);
    for (auto i = rules_.to_be_deleted.begin(); i != rules_.to_be_deleted.end(); ++i)
      tobedeleted[*i] = true;

    // Copy all rules, if a rule is tagged as 
    // to be deleted, then copy the value
    for (auto i = 0u; i != rules_.size(); ++i) {
      if (tobedeleted[i]) {
        this->add_rule(*rules_.data[i]);
      } else {
        this->add_rule(rules_.data[i]);
      } 
    }

  }

  return *this;

};

template<typename Array_Type, typename Data_Type>
inline bool Rule<Array_Type,Data_Type>::locked(const Array_Type * a, uint i, uint j) {
  return fun(a, i, j, dat);
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
    Rule<Array_Type, Data_Type> & rule
) {
  
  to_be_deleted.push_back(data.size());
  data.push_back(new Rule<Array_Type, Data_Type>(rule));
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
    Rule<Array_Type, Data_Type> * rule
) {
  
  data.push_back(rule);
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
    Rule_fun_type<Array_Type,Data_Type> rule_,
    Data_Type *                         data_,
    bool                                delete_data_
) {
  
  /* We still need to delete the counter since we are using the 'new' operator.
   * Yet, the actual data may not need to be deleted.
   */
  to_be_deleted.push_back(data.size());
  
  data.push_back(new Rule<Array_Type,Data_Type>(
      rule_,
      data_,
      delete_data_
  ));
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline bool Rules<Array_Type,Data_Type>::locked(
  const Array_Type * a, uint i, uint j
) {
  
  if (data.size()==0u)
    return false;
  
  for (auto iter = data.begin(); iter != data.end(); ++iter)
    if ((*iter)->locked(a, i, j))
      return true;
  
  return false;
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::clear() {
  
  for (auto iter = to_be_deleted.begin(); iter != to_be_deleted.end(); ++iter)
    delete data[*iter];
  
  to_be_deleted.clear();
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::get_seq(
  const Array_Type * a,
  std::vector< std::pair<uint,uint> > * free,
  std::vector< std::pair<uint,uint> > * locked
) {

  
  uint N = a->nrow();
  uint K = a->ncol();
  
  // Reserving some space
  free->empty();
  free->reserve(N*K);
  
  for (uint i = 0u; i < N; ++i) {
    for (uint j = 0u; j < K; ++j) {
      if (!this->locked(a, i, j))
        free->push_back({i, j});
      else if (locked != nullptr)
        locked->push_back({i,j});
    }
  }
  
  free->shrink_to_fit();

  return;

}

#endif
