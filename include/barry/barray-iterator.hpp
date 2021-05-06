// #include <vector>
// #include <unordered_map>
// #include "typedefs.hpp"
// #include "barray-bones.hpp"

#ifndef BARRAY_ITERATOR_HPP 
#define BARRAY_ITERATOR_HPP 1

template<typename Cell_Type, typename Data_Type>
class ConstBArrayRowIter {
public:
    
    uint current_row, current_col;
    typename Row_type<Cell_Type>::const_iterator iter;
    const BArray<Cell_Type,Data_Type> * Array;
    
    ConstBArrayRowIter(const BArray<Cell_Type,Data_Type> * Array_) : Array(Array_) {
        
        // Finding the first entry of the iterator
        for (uint i = 0u; i < Array->nrow(); ++i)
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
