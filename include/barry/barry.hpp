// #include <iostream>
// #include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <map>
#include <algorithm>
#include <utility>

#ifndef BARRY_HPP
#define BARRY_HPP 

/**
 * @brief barry: Your go-to motif accountant
 */
namespace barry {
  
  //! Tree class and TreeIterator class
  #include "typedefs.hpp"

  #include "cell-bones.hpp"
  #include "cell-meat.hpp"

  #include "barray-bones.hpp"
  #include "barraycell-bones.hpp"
  #include "barray-meat.hpp"
  #include "barraycell-meat.hpp"
  #include "barray-meat-operators.hpp"
  
  #include "counters-bones.hpp"
  #include "counters-meat.hpp"

  #include "statscounter-bones.hpp"
  #include "statscounter-meat.hpp"

  #include "support-bones.hpp"
  #include "support-meat.hpp"

  #include "block-bones.hpp"
  #include "block-meat.hpp"

  #include "powerset-bones.hpp"
  #include "powerset-meat.hpp"

  #include "model-bones.hpp"
  #include "model-meat.hpp"
  
  #include "rules-bones.hpp"
  #include "rules-meat.hpp"
  
  namespace counters {
    namespace network {
      #include "counters/network.hpp"
    }
    namespace phylo {
      #include "counters/phylo.hpp"
    }
  }
  
}

namespace netcounters = barry::counters::network;
namespace phylocounters = barry::counters::phylo;

#define COUNTER_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  inline double (a) (const Array_Type & Array, uint i, uint j, Data_Type * data)\

#define COUNTER_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  Counter_fun_type<Array_Type, Data_Type> a = \
  [](const Array_Type & Array, uint i, uint j, Data_Type * data)

#define RULE_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  inline bool (a) (const Array_Type & Array, uint i, uint j, Data_Type * data)\

#define RULE_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  Rule_fun_type<Array_Type, Data_Type> a = \
  [](const Array_Type & Array, uint i, uint j, Data_Type * data)

#endif