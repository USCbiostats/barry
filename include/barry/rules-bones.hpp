// #include <vector>
// #include <stdexcept>

#include "typedefs.hpp"

#ifndef BARRY_RULES_BONES_HPP
#define BARRY_RULES_BONES_HPP 1

template <typename Array_Type, typename Data_Type>
bool rule_fun_default(const Array_Type * array, uint i, uint j, Data_Type * dat) {
  return false;
}

/**@brief Sequence of cells indices.
 * 
 */
template<typename Array_Type = BArray<>, typename Data_Type = bool>
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
class Rules {

private:
  std::vector< Rule<Array_Type,Data_Type> * > data = {};
  std::vector< uint > to_be_deleted                = {};
  
public:
  Rules() {};

  Rules(const Rules<Array_Type,Data_Type> & rules_);
  Rules<Array_Type,Data_Type> operator=(const Rules<Array_Type,Data_Type> & rules_);

  ~Rules() {
    this->clear();
    return;
  }

  uint size() const {
    return data.size();
  };
  
  // Functions to add rules
  void add_rule(Rule<Array_Type, Data_Type> & rule);
  void add_rule(Rule<Array_Type, Data_Type> * rule);
  void add_rule(
      Rule_fun_type<Array_Type,Data_Type> rule_,
      Data_Type *                         data_        = nullptr,
      bool                                delete_data_ = false
  );
  
  bool locked(const Array_Type * a, uint i, uint j);
  
  void clear();
  
  /* Computes the sequence of free and locked cells in 
   * 
   */
  void get_seq(
    const Array_Type * a,
    std::vector< std::pair<uint,uint> > * free,
    std::vector< std::pair<uint,uint> > * locked = nullptr
  );
  
};


#endif
