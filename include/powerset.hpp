#include <vector>
#include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef LBARRAY_POWERSET_HPP 
#define LBARRAY_POWERSET_HPP 1

/* A cell is a fundamental type that holds information about a cell
 * for now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
template <typename Cell_Type > 
class PowerSet {
  
public:
  std::vector< BArray<Cell_Type> > data;
  uint size, N, M;
  uint counter;
  
  /***
   * ! Creator and destructor functions
   */
  // PSet() : data(0u), size(0u), M(0u), N(0u) {};
  PowerSet(uint N_, uint M_) : data(0u), size(0u), M(M_), N(N_), counter(0u) {};
  ~PowerSet() {};
  
  void calc(uint i = 0u, BArray<Cell_Type> * a0 = nullptr);
  
  
};

/***
 * Function to generate the powerset of the 
 */
template <typename Cell_Type>
inline void PowerSet<Cell_Type>::calc(uint i, BArray<Cell_Type> * a0) {
  
  // If it is size 0, then need to recalculate the size, and initialize
  // the first datum
  if (a0 == nullptr) {
    
    this->data.reserve(pow(2.0, this->N * this->M));
    this->data.push_back(BArray<Cell_Type>(this->N, this->M));
    a0 = &this->data.at(0u);
    
  }
  
  // Here is the deal
  this->data.push_back(*a0); // Making it one
  this->data.at(++this->counter).insert_cell(i, false, false);
  
  // If the sum is even, we increase i, otherwise j
  if ((i + 1) < (this->N * this->M)) {
    pset(i + 1u, &this->data.at(this->counter));
    pset(i + 1u, a0);
  }
  
  return;
  
}

#endif
