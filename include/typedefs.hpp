#include <unordered_map>

#ifndef BARRAY_TYPEDEFS_HPP
#define BARRAY_TYPEDEFS_HPP 

// Constants
namespace CHECK {
  const int BOTH = -1;
  const int NONE = 0;
  const int ONE  = 1;
  const int TWO  = 2;
  
}

namespace EXISTS {
  const int BOTH = -1;
  const int NONE = 0;
  const int ONE  = 1;
  const int TWO  = 1;
  
  const int UKNOWN  = -1;
  const int AS_ZERO = 0;
  const int AS_ONE  = 1;
}




// Definition of the class structure

// Edgelist
typedef unsigned int uint;
class Cell;
typedef std::unordered_map< uint, Cell > umap_int_cell;
typedef std::unordered_map< uint, Cell* > umap_int_cell_ptr;

struct edge_list {
  std::vector< uint > source;
  std::vector< uint > target;
  std::vector< double > val;
};

#endif