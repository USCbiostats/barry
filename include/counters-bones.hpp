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
template <typename Cell_Type>
class Counter {
public:
  Counter_type<Cell_Type> count_fun;
  Counter_type<Cell_Type> init_fun;
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  Counter(Counter_type<Cell_Type> count_fun_) :
    count_fun(count_fun_), init_fun(nullptr) {};
  Counter(Counter_type<Cell_Type> count_fun_, Counter_type<Cell_Type> init_fun_):
    count_fun(count_fun_), init_fun(init_fun_) {};
  
  double count(const BArray<Cell_Type> * Array, uint i, uint j);
  double init(const BArray<Cell_Type> * Array, uint i, uint j);
};

template <typename Cell_Type>
inline double Counter<Cell_Type>::count(
    const BArray<Cell_Type> * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j);
}

template <typename Cell_Type>
inline double Counter<Cell_Type>::init(
    const BArray<Cell_Type> * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j);
}


namespace counters {

  // Edges counter
  template <typename Cell_Type>
  inline double count_edges(const BArray<Cell_Type> * Array, uint i, uint j) {
    return 1.0;
  }

  Counter<bool> edges(count_edges<bool>);
  
  // Isolates counter
  template <typename Cell_Type>
  inline double count_isolates(const BArray<Cell_Type> * Array, uint i, uint j) {
    
    if (i == j)
      return 0.0;
    
    double res = 0.0;
    
    // i is sending its first tie
    if (A_ROW(i).size() == 1u && A_COL(i).size() == 0u)
      res -= 1.0;
    
    // j is receiving its first tie, meaning that he
    // has no other tie but i's?
    if (A_ROW(j).size() == 0u && A_COL(j).size() == 1u)
      res -= 1.0;
    
    return res;
    
  }
  
  template <typename Cell_Type>
  inline double init_isolates(const BArray<Cell_Type> * Array, uint i, uint j) {
    return (double) (Array->N);
  }
  
  Counter<bool> isolates(count_isolates<bool>, init_isolates<bool>);
  
  // Mutuals
  template <typename Cell_Type>
  inline double count_mutual(const BArray<Cell_Type> * Array, uint i, uint j) {

    // Is there any tie at ji? If not, then we have a new mutual!
    // but this only makes sence if the jth row and ith column exists
    // if ((Array->N > j) && (Array->M > i)) 
    if (i == j)
      return 0.0;
    
    // printf("Checking if it is empty or not at (%i, %i)... ", i, j);
    if (!Array->is_empty(j, i, false)) {
      // printf("Yes, mutual.\n");
        return 1.0;
    }
    // printf("No, no mutual.\n");
    
    return 0.0;
    
  }
  
  Counter<bool> mutual(count_mutual<bool>);
  
  // 2-istars
  template<typename Cell_Type>
  inline double count_istar2(const BArray<Cell_Type> * Array, uint i, uint j) {
   
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
  
    if (A_COL(j).size() == 1u)
      return 0.0;
    
    return ((double) A_COL(j).size() - 1.0);
   
    // return 0.0; 
  }

  Counter<bool> istar2(count_istar2<bool>);
  
  // 2-ostars
  template<typename Cell_Type>
  inline double count_ostar2(const BArray<Cell_Type> * Array, uint i, uint j) {
    
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
    
    if (A_ROW(i).size() == 1u)
      return 0.0;
    
    return ((double) A_ROW(i).size() - 1.0);
    
    // return 0.0; 
  }
  
  Counter<bool> ostar2(count_ostar2<bool>);
  
  // ttriads
  template <typename Cell_Type>
  inline double count_ttriads(const BArray<Cell_Type> * Array, uint i, uint j) {
    
    // i-j, j-k, i-k
    double ans = 0.0;
    if (A_ROW(i).size() > A_ROW(j).size()) {
      
      for (auto col = A_ROW(j).begin(); col != A_ROW(j).end(); ++col)
        if (!Array->is_empty(i, col->first))
          ans += 1.0;
      
    } else {
      
      for (auto col = A_ROW(i).begin(); col != A_ROW(i).end(); ++col)
        if (!Array->is_empty(j, col->first))
          ans += 1.0;
      
    }
    
    return ans;

  }
  
  Counter<bool> ttriads(count_ttriads<bool>);
}

#endif