#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_COUNTERS_H
#define BARRAY_COUNTERS_H

#define COUNTER_FUNCTION(a) template <typename Array_Type, typename Data_Type> \
  inline double (a) (Array_Type * Array, uint i, uint j, Data_Type * counter)

/**
 * @brief A counter function based on change statistics.
 * 
 * This class is used by `CountStats` and `StatsCounter` as a way to count
 * statistics using change statistics.
 */
template <typename Array_Type, typename Data_Type>
class Counter {
public:
  
  Data_Type * data = nullptr;
  Counter_fun_type<Array_Type,Data_Type> count_fun;
  Counter_fun_type<Array_Type,Data_Type> init_fun;

  /***
   * ! Initializers
   */
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  
  /**
   * @brief Creator passing only a counter function
   * @param count_fun The main counter function.
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_
    ) :
    count_fun(count_fun_), init_fun(nullptr) {};
  /**
   * @brief Creator passing a counter and an initializer
   * 
   * @param count_fun The main counter function.
   * @param init_fun The initializer function can also be used to check if the
   *  `BArray` has the required metadata (e.g. is it a directed graph?).
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_
    ): count_fun(count_fun_), init_fun(init_fun_) {};
  
  ~Counter() {};
  
  /***
   * ! Main functions.
   */
  double count(Array_Type * Array, uint i, uint j);
  double init(Array_Type * Array, uint i, uint j);
};

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::count(
    Array_Type * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j, data);
}

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::init(
    Array_Type * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j, data);
}


namespace counters {

  // Edges counter
  COUNTER_FUNCTION(count_edges) {
    return 1.0;
  } 

  Counter<BArray<bool,bool>, bool> edges(count_edges<BArray<bool,bool>, bool>);
   
  // Isolates counter
  COUNTER_FUNCTION(count_isolates) {
    
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
  
  COUNTER_FUNCTION(init_isolates) {
    return (double) (Array->N);
  }
  
  Counter<BArray<bool,bool>, bool> isolates(
      count_isolates<BArray<bool,bool>, bool>,
      init_isolates<BArray<bool,bool>, bool>
  );
  
  // Mutuals -------------------------------------------------------------------
  COUNTER_FUNCTION(init_mutual) {
    
    if (Array->N != Array->M)
      throw std::logic_error("The -mutual- counter only works on square arrays.");

    if (Array->meta.is("symmetric") || Array->meta.is("undirected"))
      throw std::logic_error("The -mutual- counter only works on directed (non-symmetric) arrays.");
    
    return 0.0;
  }
  
  template <typename Array_Type, typename Data_Type>
  inline double count_mutual(
      Array_Type * Array, uint i, uint j,
      Data_Type * counter
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
  
  Counter<BArray<bool,bool>, bool> mutual(
      count_mutual<BArray<bool,bool>, bool>,
      init_mutual<BArray<bool,bool>, bool>
  );
  
  // 2-istars
  COUNTER_FUNCTION(count_istar2) {
   
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
  
    if (A_COL(j).size() == 1u)
      return 0.0;
    
    return ((double) A_COL(j).size() - 1.0);
   
    // return 0.0; 
  }

  Counter<BArray<bool,bool>, bool> istar2(count_istar2<BArray<bool,bool>, bool>);
  
  // 2-ostars
  COUNTER_FUNCTION(count_ostar2) {
    
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
    
    if (A_ROW(i).size() == 1u)
      return 0.0;
    
    return ((double) A_ROW(i).size() - 1.0);
    
    // return 0.0; 
  }
  
  Counter<BArray<bool,bool>, bool> ostar2(count_ostar2<BArray<bool,bool>, bool>);
  
  // ttriads
  COUNTER_FUNCTION(count_ttriads) {
    
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
  
  Counter<BArray<bool,bool>, bool> ttriads(count_ttriads<BArray<bool,bool>, bool>);
  
  // Cycle triads --------------------------------------------------------------
  COUNTER_FUNCTION(count_ctriads) {
    
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
  
  Counter<BArray<bool,bool>, bool> ctriads(count_ctriads<BArray<bool,bool>, bool>);
  
  // Density --------------------------------------------------------------
  COUNTER_FUNCTION (count_density) {
    return 1.0/(Array->N * (Array->M - 1));
  }
  
  Counter<BArray<bool,bool>, bool> density(count_density<BArray<bool,bool>, bool>);
  
  // idegree1.5  -------------------------------------------------------------
  COUNTER_FUNCTION(count_idegree15) {
    
    // In case of the first, we need to add
    if (A_COL(j).size() == 1u)
      return 1.0;
    
    return 
      pow((double) A_COL(j).size(), 1.5) - pow((double) A_COL(j).size() - 1, 1.5)
    ;
  }
  
  Counter<BArray<bool,bool>, bool> idegree15(count_idegree15<BArray<bool,bool>, bool>);
  
  // odegree1.5  -------------------------------------------------------------
  COUNTER_FUNCTION(count_odegree15) {
    
    // In case of the first, we need to add
    if (A_ROW(i).size() == 1u)
      return 1;
    
    return 
      pow((double) A_ROW(i).size(), 1.5) - pow((double) A_ROW(i).size() - 1, 1.5)
    ;
  }
  
  Counter<BArray<bool,bool>, bool> odegree15(count_odegree15<BArray<bool,bool>,bool>);
}

#endif
