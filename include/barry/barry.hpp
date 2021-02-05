// #include <iostream>
// #include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <map>

#ifndef BARRY_HPP
#define BARRY_HPP 

namespace barry {
  
  //! Tree class and TreeIterator class
  #include "typedefs.hpp"
  #include "cell-bones.hpp"
  #include "barray-bones.hpp"
  #include "barray-meat.hpp"
  
  #include "counters-bones.hpp"
  #include "statscounter.hpp"
  #include "support.hpp"
  #include "block-bones.hpp"
  #include "powerset.hpp"
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
  inline double (a) (const Array_Type * Array, uint i, uint j, Data_Type * data)\

#define COUNTER_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  Counter_fun_type<Array_Type, Data_Type> a = \
  [](const Array_Type * Array, uint i, uint j, Data_Type * data)

#define RULE_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  inline bool (a) (const Array_Type * Array, uint i, uint j, Data_Type * data)\

#define RULE_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
  Rule_fun_type<Array_Type, Data_Type> a = \
  [](const Array_Type * Array, uint i, uint j, Data_Type * data)

#endif