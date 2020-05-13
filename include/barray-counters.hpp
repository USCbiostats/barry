#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_COUNTERS_H
#define BARRAY_COUNTERS_H

/***
 * ! Counter class of function. Right now it has two members, 
 * ! the initializer and the actual counter. By default all counts are initialized
 * ! as zero.
 *  
 */
class Counter {
public:
  Counter_type count_fun;
  Counter_type init_fun;
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  Counter(Counter_type count_fun_) : count_fun(count_fun_), init_fun(nullptr) {};
  Counter(Counter_type count_fun_, Counter_type init_fun_):
    count_fun(count_fun_), init_fun(init_fun_) {};
  
  double count(const BArray * Array, uint i, uint j);
  double init(const BArray * Array, uint i, uint j);
};

inline double Counter::count(const BArray * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j);
}

inline double Counter::init(const BArray * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j);
}


namespace counters {

  // Edges counter
  inline double count_edges(const BArray * Array, uint i, uint j) {
    return 1.0;
  }

  Counter edges(count_edges);
  
  // Isolates counter
  inline double count_isolates(const BArray * Array, uint i, uint j) {
    
    double res = 0.0;
    // If it was previously empty, then he is to worry
    if (A_ROW(i).size() == 1u && A_COL(i).size())
    return -1.0;
    
  }
  inline double init_isolates(const BArray * Array, uint i, uint j) {
    return (double) (Array->N);
  }
  
  Counter isolates(count_isolates, init_isolates);
  
  // Mutuals
  inline double count_mutual(const BArray * Array, uint i, uint j) {
    
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
  
  Counter mutual(count_mutual);

}

#endif