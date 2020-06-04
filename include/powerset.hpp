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
template <typename Cell_Type, typename Data_Type > 
class PowerSet {
  
public:
  std::vector< BArray<Cell_Type, Data_Type> > data;
  uint N, M;
  uint counter;
  
  /***
   * ! Creator and destructor functions
   */
  PowerSet() : data(0u), M(0u), N(0u) {};
  PowerSet(uint N_, uint M_) : data(0u), M(M_), N(N_), counter(0u) {};
  ~PowerSet() {};
  
  void calc(uint i = 0u, BArray<Cell_Type, Data_Type> * a0 = nullptr);
  void reset(uint N_, uint M_);
  
  /**
   * Getter functions
   */
  const std::vector< BArray<Cell_Type, Data_Type> > * get_data_ptr() const {return &data;};
  std::vector< BArray<Cell_Type, Data_Type> > get_data() const {return data;};
  
};

/***
 * Function to generate the powerset of the 
 */
template <typename Cell_Type, typename Data_Type>
inline void PowerSet<Cell_Type, Data_Type>::calc(
    uint i,
    BArray<Cell_Type,Data_Type> * a0
  ) {
  
  // If it is size 0, then need to recalculate the size, and initialize
  // the first datum
  if (a0 == nullptr) {
    
    this->data.reserve(pow(2.0, this->N * this->M));
    this->data.push_back(BArray<Cell_Type,Data_Type>(this->N, this->M));
    a0 = &this->data.at(0u);
    
  }
  
  // Here is the deal
  this->data.push_back(*a0); // Making it one
  this->data.at(++this->counter).insert_cell(i, true, false, false);
  
  // If the sum is even, we increase i, otherwise j
  if ((i + 1) < (this->N * this->M)) {
    calc(i + 1u, &this->data.at(this->counter));
    calc(i + 1u, a0);
  }
  
  return;
  
}

template <typename Cell_Type, typename Data_Type>
inline void PowerSet<Cell_Type,Data_Type>::reset(
    uint N_,
    uint M_
  ) {
  
  data.empty();
  N = N_, M = M_, counter = 0u;
  
  return;
}

#endif
