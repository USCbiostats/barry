#include "barray.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H

struct NodeData;
typedef barray::BArray<bool, NodeData> PhyloArray;

typedef unsigned int uint;

template <typename T>
using Vec = std::vector< T >;

struct NodeData {
  Vec< double > blengths;
  Vec< bool > states;
};


/**@brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */



namespace phylo_counters {

   
  // Functional gains
  COUNTER_FUNCTION(count_gains) {
    
    if (i == *data)
      return 1.0;
    else
      return 0.0;
    
  }
  
  COUNTER_FUNCTION(init_count_gains) {
    return 0.0;
  }
  
  barray::Counter<PhyloArray, uint> gains(
      count_gains<PhyloArray, uint>,
      init_count_gains<PhyloArray, uint>
  );

}


#endif
