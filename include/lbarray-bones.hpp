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
  void pset(uint i = 0u, uint j = 0u, BArray * a0 = nullptr) {
    
    // If it is size 0, then need to recalculate the size, and initialize
    // the first datum
    if (a0 == nullptr) {
      this->data.reserve(pow(N * M, 2.0));
      BArray a_empty(this->N, this->M);
      a0 = &a_empty;
    }
    
    
    
    // Here is the deal
    this->data.push_back(*a0);
    this->data.push_back(*a0); // Making it one
    this->counter += 2;
    std::cout << "Coords ("<< i << ", " << j << "), current size: " << data.size() << std::endl;
    std::cout << "Size of [0] (" << this->data.at(0).N << ", " <<  this->data.at(0).M << ")" << std::endl;
    std::cout << "Size of [1] (" << this->data.at(1).N << ", " <<  this->data.at(1).M << ")" << std::endl;
    
    this->data.at(counter - 1).toggle_cell(i, j);
    
    a0          = &this->data.at(this->counter - 2);
    BArray * a1 = &this->data.at(this->counter - 1);
    
    std::cout << "Entering the recursion" << std::endl;
    
    // If the sum is even, we increase i, otherwise j
    if ((i + 1) < this->N)
      pset(i + 1, j, a1);
    
    if ((j + 1) < this->M)
      pset(i, j + 1, a1);
    
    if ((i + 1) < this->N)
      pset(i + 1, j, a0);
    
    if ((j + 1) < this->M)
      pset(i, j + 1, a0);
    
    return;
    
  }
};

#endif
