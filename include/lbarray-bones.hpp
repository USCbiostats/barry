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

typedef std::function<double(const BArray *, uint, uint)> f_array_uint_uint;
 

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
  
  std::vector< f_array_uint_uint > f;
  StatsCounter(const BArray * data, SuffStats * stats_) :
    f(0u), Array(data), EmptyArray(data->N, data->M), stats(stats_) {}
  
  void add_fun(f_array_uint_uint f_) {
    f.push_back(f_);
    return;
  }
  
  /***
   * ! This function recurses through the entries of Array and at each step of
   * ! adding a new cell it uses the functions to list the statistics.
   */
  void count_current(uint i, uint j);
  void count_all();
  
};

inline void StatsCounter::count_current(uint i, uint j) {
  
  // Iterating through the functions, and updating the set of
  // statistics.
  for (uint n = 0u; n < f.size(); ++n) {
    this->change_stats.at(n) = (f.at(n))(this->Array, i, j);
    this->current_stats.at(n) += this->change_stats.at(n);
  }

  return;
  
}

inline void StatsCounter::count_all() {
  
  // Initializing the counter on the empty array
  this->current_stats.resize(this->f.size(), 0.0);
  this->change_stats.resize(this->f.size());
  
  // Setting it to zero.
  this->EmptyArray.clear();
  
  // Start iterating through the data
  uint N = this->Array->N;
  uint M = this->Array->M;
  for (uint i = 0; i < N; ++i) {
    
    // Any element?
    if (Array->el_ij.at(i).size() == 0u)
      continue;
    
    // If there's one, then update the statistic, by iterating
    for (auto iter = Array->el_ij.at(i).begin(); iter != Array->el_ij.at(i).end(); ++iter) {
      
      // Adding a cell
      EmptyArray.insert_cell(i, iter->first);
      
      // Computing the change statistics
      this->count_current(i, iter->first);
     
    }
    
  }
  
  // Adding to the sufficient statistics
  this->stats->add(this->current_stats);
  return;
  
}

/***
 * Some example functions
 */

// Always adding, so this is trivially equal to +1
double counter_edges(const BArray * Array, uint i, uint j) {
  return 1.0;
}

double counter_mutual(const BArray * Array, uint i, uint j) {
  
  // We only count one direction
  if (i > j)
    return 0.0;
  
  // Is there any tie at ji? If not, then we have a new mutual!
  // but this only makes sence if the jth row and ith column exists
  if ((Array->N > j) && (Array->M > i) && !Array->is_empty(j, i)) {
    // std::cout << "Yep, Isee\n";
    return 1.0;
  }
  
  return 0.0;
  
}


#endif
