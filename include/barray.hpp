#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <map>

#ifndef BARRAY_HPP
#define BARRAY_HPP 

namespace barray {
  
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
  
  namespace counters {
    namespace network {
      #include "counters/network.hpp"
    }
    namespace phylo {
      #include "counters/phylo.hpp"
    }
  }
  
}

namespace netcounters = barray::counters::network;
namespace phylocounters = barray::counters::phylo;
#endif