#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_COUNTERS_H
#define BARRAY_COUNTERS_H

/**
 * @brief A counter function based on change statistics.
 * 
 * This class is used by `CountStats` and `StatsCounter` as a way to count
 * statistics using change statistics.
 */
template <typename Cell_Type>
class Counter {
public:
  
  Counter_type<Cell_Type> count_fun;
  Counter_type<Cell_Type> init_fun;

  /***
   * ! Initializers
   */
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  
  /**
   * @brief Creator passing only a counter function
   * @param count_fun The main counter function.
   */
  Counter(Counter_type<Cell_Type> count_fun_) :
    count_fun(count_fun_), init_fun(nullptr) {};
  /**
   * @brief Creator passing a counter and an initializer
   * 
   * @param count_fun The main counter function.
   * @param init_fun The initializer function can also be used to check if the
   *  `BArray` has the required metadata (e.g. is it a directed graph?).
   */
  Counter(Counter_type<Cell_Type> count_fun_, Counter_type<Cell_Type> init_fun_):
    count_fun(count_fun_), init_fun(init_fun_) {};
  
  ~Counter() {};
  
  /***
   * ! Main functions.
   */
  double count(BArray<Cell_Type> * Array, uint i, uint j);
  double init(BArray<Cell_Type> * Array, uint i, uint j);
};

template <typename Cell_Type>
inline double Counter<Cell_Type>::count(
    BArray<Cell_Type> * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j, &this);
}

template <typename Cell_Type>
inline double Counter<Cell_Type>::init(
    BArray<Cell_Type> * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j, &this);
}


namespace counters {

  // Edges counter
  template <typename Cell_Type>
  inline double count_edges(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
    ) {
    return 1.0;
  }

  Counter<bool> edges(count_edges<bool>);
   
  // Isolates counter
  template <typename Cell_Type>
  inline double count_isolates(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
    ) {
    
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
  inline double init_isolates(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    return (double) (Array->N);
  }
  
  Counter<bool> isolates(count_isolates<bool>, init_isolates<bool>);
  
  // Mutuals -------------------------------------------------------------------
  template <typename Cell_Type>
  inline double init_mutual(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
    if (Array->N != Array->M)
      throw std::logic_error("The -mutual- counter only works on square arrays.");

    if (Array->meta.is("symmetric") || Array->meta.is("undirected"))
      throw std::logic_error("The -mutual- counter only works on directed (non-symmetric) arrays.");
    
    return 0.0;
  }
  
  template <typename Cell_Type>
  inline double count_mutual(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
    ) {

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
  
  Counter<bool> mutual(count_mutual<bool>, init_mutual<bool>);
  
  // 2-istars
  template<typename Cell_Type>
  inline double count_istar2(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
   
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
  inline double count_ostar2(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
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
  inline double count_ttriads(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
    // Self ties do not count
    if (i == j)
      return 0.0;
    
    double ans = 0.0;
    
    // Case 1: i-j, i-k, j-k
    if (A_ROW(j).size() < A_ROW(i).size()) {
      
      for (auto j_row = A_ROW(j).begin(); j_row != A_ROW(j).end(); ++j_row) 
        if ((j != j_row->first) && (i != j_row->first) && !Array->is_empty(i, j_row->first, false))
          ans += 1.0;

    } else {
      
      for (auto i_row = A_ROW(i).begin(); i_row != A_ROW(i).end(); ++i_row) 
        if ((i != i_row->first) && (i_row->first != j) && !Array->is_empty(j, i_row->first, false))
          ans += 1.0;
        
    }
    
    // Case 2: i-j, i-k, k-j  
    if (A_ROW(i).size() > A_COL(j).size()) {
      
      for (auto j_col = A_COL(j).begin(); j_col != A_COL(j).end(); ++j_col)
        if ((j != j_col->first) && (i != j_col->first) && !Array->is_empty(i, j_col->first, false))
          ans += 1.0;
        
    } else {
      
      for (auto i_row = A_ROW(i).begin(); i_row != A_ROW(i).end(); ++i_row) 
        if ((i != i_row->first) && (j != i_row->first) && !Array->is_empty(i_row->first, j, false))
          ans += 1.0;

    }
    
    // Case 3: 
    if (A_COL(i).size() > A_COL(j).size()) {
      
      for (auto j_col = A_COL(j).begin(); j_col != A_COL(j).end(); ++j_col)
        if ((j != j_col->first) && (i != j_col->first) && !Array->is_empty(j_col->first, i, false))
          ans += 1.0;
        
    } else {
      
      for (auto i_col = A_COL(i).begin(); i_col != A_COL(i).end(); ++i_col) 
        if ((i != i_col->first) && (j != i_col->first) && !Array->is_empty(i_col->first, j, false))
          ans += 1.0;
        
    }
   
    
    // The regular counter double counts
    return ans;

  }
  
  Counter<bool> ttriads(count_ttriads<bool>);
  
  // Cycle triads --------------------------------------------------------------
  template <typename Cell_Type>
  inline double count_ctriads(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
    if (i == j)
      return 0.0;
    
    double ans = 0.0;
    if (A_COL(i).size() < A_ROW(j).size()) {
      
      for (auto i_col = A_COL(i).begin(); i_col != A_COL(i).end(); ++i_col) 
        if ((i != i_col->first) && (j != i_col->first) && !Array->is_empty(j, i_col->first, false))
          ans += 1.0;
      
    } else {
      
      for (auto j_row = A_ROW(j).begin(); j_row != A_ROW(j).end(); ++j_row) 
        if ((j != j_row->first) && (i != j_row->first) && !Array->is_empty(j_row->first, i, false))
          ans += 1.0;
      
    }
    
    return ans;
    
  }
  
  Counter<bool> ctriads(count_ctriads<bool>);
  
  // Density --------------------------------------------------------------
  template <typename Cell_Type>
  inline double count_density(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    return 1.0/(Array->N * (Array->M - 1));
  }
  
  Counter<bool> density(count_density<bool>);
  
  // idegree1.5  -------------------------------------------------------------
  template <typename Cell_Type>
  inline double count_idegree15(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
    // In case of the first, we need to add
    if (A_COL(j).size() == 1u)
      return 1.0;
    
    return 
      pow((double) A_COL(j).size(), 1.5) - pow((double) A_COL(j).size() - 1, 1.5)
    ;
  }
  
  Counter<bool> idegree15(count_idegree15<bool>);
  
  // odegree1.5  -------------------------------------------------------------
  template <typename Cell_Type>
  inline double count_odegree15(
      BArray<Cell_Type> * Array, uint i, uint j,
      Counter<Cell_Type> * counter
  ) {
    
    // In case of the first, we need to add
    if (A_ROW(i).size() == 1u)
      return 1;
    
    return 
      pow((double) A_ROW(i).size(), 1.5) - pow((double) A_ROW(i).size() - 1, 1.5)
    ;
  }
  
  Counter<bool> odegree15(count_odegree15<bool>);
}

#endif
