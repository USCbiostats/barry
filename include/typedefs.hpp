#include <unordered_map>

#ifndef BARRAY_TYPEDEFS_HPP
#define BARRAY_TYPEDEFS_HPP 

// Mostly relevant for the BArray definition -----------------------------------
#define ROW(a) this->el_ij.at(a)
#define COL(a) this->el_ji.at(a)

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
template <class Type_A > class Cell;

template<typename Cell_Type>
using Row_type = std::unordered_map< uint, Cell<Cell_Type> >;

template<typename Cell_Type>
using Col_type = std::unordered_map< uint, Cell<Cell_Type>* >;

template<typename Cell_Type>
class Entries {
public:
  std::vector< uint > source;
  std::vector< uint > target;
  std::vector< Cell_Type > val;
  
  Entries() : source(0u), target(0u), val(0u) {};
  Entries(uint n) {
    source.reserve(n);
    target.reserve(n);
    val.reserve(n);
    return;
  };
  
  ~Entries() {};
  
  void resize(uint n) {
    source.resize(n);
    target.resize(n);
    val.resize(n);
    return;
  }
  
};
 
// Mostly relevant in the case of the stats count functions -------------------
template <typename Cell_Type> class BArray;
#define A_ROW(a) Array->el_ij.at(a)
#define A_COL(a) Array->el_ji.at(a)

template <typename Cell_Type>
using Counter_type = std::function<double(const BArray<Cell_Type> *, uint, uint)>;


#endif