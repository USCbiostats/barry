#include <vector>
#include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef LBARRAY_BONES_HPP 
#define LBARRAY_BONES_HPP 1

/* A cell is a fundamental type that holds information about a cell
 * for now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */

class LBArray {
  
public:
  std::vector< BArray > data;
  uint size, N, M;
  uint counter;
  
  /***
   * ! Creator and destructor functions
   */
  // LBArray() : data(0u), size(0u), M(0u), N(0u) {};
  LBArray(uint N_, uint M_) : data(0u), size(0u), M(M_), N(N_), counter(0u) {};
  ~LBArray() {};
  
  /***
   * Function to generate the powerset of the 
   */
  void pset(uint i = 0u, BArray * a0 = nullptr) {
    
    // If it is size 0, then need to recalculate the size, and initialize
    // the first datum
    if (a0 == nullptr) {
      
      this->data.reserve(pow(2.0, this->N * this->M));
      this->data.push_back(BArray(this->N, this->M));
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
  
};

#endif
