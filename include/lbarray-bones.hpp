#include <vector>
#include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "barray-counters.hpp"

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
  uint ncols;
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


/***
 * ! Users can a list of functions that can be used with this.
 * ! The baseline set of arguments is a pointer to a binary array.
 */ 
class StatsCounter {
public:
  
  // Should receive an array
  const BArray * Array;
  BArray EmptyArray;
  std::vector< double > current_stats;
  std::vector< double > change_stats;
   
  // We will save the data here
  SuffStats * stats;
  
  std::vector< Counter > counters;
  StatsCounter(const BArray * data, SuffStats * stats_) :
    counters(0u), Array(data), EmptyArray(data->N, data->M), stats(stats_) {}
  
  void add_counter(Counter f_);
  
  /***
   * ! This function recurses through the entries of Array and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_init(uint i, uint j);
  void count_current(uint i, uint j);
  void count_all();
  
};

inline void StatsCounter::add_counter(Counter f_){
  counters.push_back(f_);
  return;
}

inline void StatsCounter::count_init(uint i, uint j) {
  // Iterating through the functions, and updating the set of
  // statistics.
  current_stats.resize(counters.size(), 0.0);
  change_stats.resize(counters.size(), 0.0);
  for (uint n = 0u; n < counters.size(); ++n) 
    current_stats.at(n) = counters.at(n).init(Array, i, j);
  
  return;
}

inline void StatsCounter::count_current(uint i, uint j) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  for (uint n = 0u; n < counters.size(); ++n) {
    change_stats.at(n)   = counters.at(n).count(Array, i, j);
    current_stats.at(n) += change_stats.at(n);
  }

  return;
  
}

inline void StatsCounter::count_all() {
  
  // Initializing the counter on the empty array
  count_init(0u, 0u);
  
  // Setting it to zero.
  EmptyArray.clear();
  
  // Start iterating through the data
  uint N = Array->N;
  uint M = Array->M;
  for (uint row = 0; row < N; ++row) {
    
    // Any element?
    if (A_ROW(row).size() == 0u)
      continue;
    
    // If there's one, then update the statistic, by iterating
    for (auto col = A_ROW(row).begin(); col != A_ROW(row).end(); ++col) {
      
      // Adding a cell
      EmptyArray.insert_cell(row, col->first);
      
      // Computing the change statistics
      count_current(row, col->first);
     
    } 
    
  }
  
  // Adding to the sufficient statistics
  stats->add(current_stats);
  return;
  
}

/***
 * Some example functions
 */

#endif
