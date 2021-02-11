// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "rules-bones.hpp"

#ifndef BARRY_POWERSET_BONES_HPP 
#define BARRY_POWERSET_BONES_HPP 1

/* A cell is a fundamental type that holds information about a cell
 * for now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
template <typename Array_Type = BArray<>, typename Data_Rule_Type = bool> 
class PowerSet {
  
private:
  void calc_backend(uint pos = 0u);  

public:
  Array_Type                         EmptyArray;
  std::vector< Array_Type >          data;
  Rules<Array_Type,Data_Rule_Type> * rules;

  uint N, M;
  bool rules_deleted   = false;

  // Tempvars
  std::vector< std::pair<uint,uint> >  coordinates_free;
  std::vector< std::pair<uint,uint> >  coordinates_locked;
  
  /***
   * ! Creator and destructor functions
   */
  PowerSet() : 
  EmptyArray(), data(0u), rules(new Rules<Array_Type,Data_Rule_Type>()), N(0u), M(0u) {};
  PowerSet(uint N_, uint M_) : EmptyArray(N_, M_), data(0u), rules(new Rules<Array_Type,Data_Rule_Type>()), N(N_), M(M_) {};
  ~PowerSet();
  
  void init_support();
  void calc();
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
  ///@{
  const std::vector< Array_Type > * get_data_ptr() const {return &data;};
  std::vector< Array_Type > get_data() const {return data;};
  typename std::vector< Array_Type >::iterator begin() {return data.begin();};
  typename std::vector< Array_Type >::iterator end() {return data.end();};
  ///@}
  
};

#endif
