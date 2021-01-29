// #include <vector>
// #include <stdexcept>

#include "typedefs.hpp"

#ifndef BARRAY_RULES_BONES_HPP
#define BARRAY_RULES_BONES_HPP 1


template <typename Array_Type, typename Data_Type>
inline bool rule_fun_default(const Array_Type * A, uint i, uint j, Data_Type *) {
  return false;
}

/**@brief Sequence of cells indices.
 * 
 */
template<typename Array_Type, typename Data_Type>
class Rule {
  
private:
  Rule_fun_type<Array_Type,Data_Type> fun;
  Data_Type * dat = nullptr;
  bool delete_dat = false;
  
public:
  Rule() : fun(rule_fun_default<Array_Type,Data_Type>) {};
  Rule(
    Rule_fun_type<Array_Type,Data_Type> fun_,
    Data_Type * dat_ = nullptr,
    bool delete_dat_ = false
    ) : fun(fun_), dat(dat_), delete_dat(delete_dat_) {};
  
  ~Rule() {
    if (delete_dat)
      delete dat;
    return;
  }
  
  bool locked(const Array_Type * a, uint i, uint j);
  
};

template<typename Array_Type, typename Data_Type>
inline bool Rule<Array_Type,Data_Type>::locked(const Array_Type * a, uint i, uint j) {
  return fun(a, i, j, dat);
}

template<typename Array_Type, typename Data_Type>
class Rules {

private:
  std::vector< Rule<Array_Type,Data_Type> * > data = {};
  std::vector< uint > to_be_deleted                 = {};
  
public:
  
  Rules() {};

  ~Rules() {
    this->clear();
    return;
  }
  
  // Functions to add rules
  void add_rule(Rule<Array_Type, Data_Type> & rule);
  void add_rule(Rule<Array_Type, Data_Type> * rule);
  void add_rule(
      Rule_fun_type<Array_Type,Data_Type> rule_,
      Data_Type *                         data_        = nullptr,
      bool                                delete_data_ = false
  );
  
  bool is_locked(const Array_Type * a, uint i, uint j);
  
  void clear();
  
  /* Computes the sequence of free and locked cells in 
   * 
   */
  void get_seq(
    const Array_Type * a,
    std::vector< std::pair<uint,uint> > * free,
    std::vector< std::pair<uint,uint> > * locked
  );
  
};

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
inline bool Rules<Array_Type,Data_Type>::is_locked(
  const Array_Type * a, uint i, uint j
) {
  
  if (data.size()==0u)
    return false;
  
  for (auto iter = data.begin(); iter != data.end(); ++iter)
    if (iter->is_locked(a, i, j))
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

  uint N = a->nrows();
  uint K = a->ncols();
  for (uint i = 0u; i < N; ++i) {
    for (uint j = 0u; j < K; ++j) {
      if (!this->is_locked(a, i, j))
        free->push_back({i, j});
      else if (locked != nullptr)
        locked->push_back({i,j});
    }
  }

  return;

}

#endif
