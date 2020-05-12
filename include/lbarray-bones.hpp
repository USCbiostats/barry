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



template <typename T>
struct vecHasher {
  std::size_t operator()(std::vector< T > const&  dat) const noexcept {
    
    std::hash< T > hasher;
    std::size_t hash = hasher(dat.at(0u));
    
    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    if (dat.size() > 1u)
      for (unsigned int i = 1u; i < dat.size(); ++i)
        hash ^= hasher(dat.at(i)) + 0x9e3779b9 + (hash<<6) + (hash>>2);
    
    return hash;
    
  }
};

typedef std::vector< std::pair< std::vector<double>, uint > > vec_pair_dbl_uint;

class SuffStats {
public:
  std::unordered_map< std::vector< double >, uint, vecHasher< double > > stats;
  SuffStats() {};
  ~SuffStats() {};
  
  void add(const std::vector< double > & x) {
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    if (stats.find(x) == stats.end()) {
      stats[x] = 1u;
    } else // We increment the counter
      stats.at(x)++;
   
    return; 
  };
  
  
  vec_pair_dbl_uint get_entries() const {
    
    vec_pair_dbl_uint ans;
    ans.reserve(stats.size());
    for (auto iter = stats.begin(); iter != stats.end(); ++iter)
      ans.push_back(*iter);
    
    
    return ans;
  };
  
  
};

#endif
