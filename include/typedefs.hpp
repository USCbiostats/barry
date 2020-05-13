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
class Cell;
typedef std::unordered_map< uint, Cell > Row_type;
typedef std::unordered_map< uint, Cell* > Col_type;

class Entries {
public:
  std::vector< uint > source;
  std::vector< uint > target;
  std::vector< double > val;
  
  Entries() : source(0u), target(0u), val(0) {};
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
class BArray;
#define A_ROW(a) Array->el_ij.at(a)
#define A_COL(a) Array->el_ji.at(a)

typedef std::function<double(const BArray *, uint, uint)> Counter_type;


#endif