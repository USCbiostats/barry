// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_ITERATOR_HPP 
#define BARRAY_ITERATOR_HPP 1


class ConstBArrayRowIter {
public:
  
  uint current_row, current_col;
  Row_type::const_iterator iter;
  const BArray * Array;
  
  ConstBArrayRowIter(const BArray * Array_) : Array(Array_) {
    
    // Finding the first entry of the iterator
    for (uint i = 0u; i < Array->N; ++i)
      if (A_ROW(i).size() != 0u) {
        iter = A_ROW(i).begin();
        break;
      }
    
    return;
  };
  ~ConstBArrayRowIter() {};
  
  // operat
    
};

#endif
