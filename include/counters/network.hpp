#include "../counters-bones.hpp"

#ifndef BARRAY_NETWORK_H
#define BARRAY_NETWORK_H

namespace network {

  // Edges counter
  COUNTER_FUNCTION_DEFAULT(count_edges) {
    return 1.0;
  } 

  Counter<> edges(count_edges<>);
   
  // Isolates counter
  COUNTER_FUNCTION_DEFAULT(count_isolates) {
    
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
  
  COUNTER_FUNCTION_DEFAULT(init_isolates) {
    return (double) (Array->N);
  }
  
  Counter<> isolates(count_isolates<>, init_isolates<>);
  
  // Mutuals -------------------------------------------------------------------
  COUNTER_FUNCTION_DEFAULT(init_mutual) {
    
    if (Array->N != Array->M)
      throw std::logic_error("The -mutual- counter only works on square arrays.");

    if (Array->meta.is("symmetric") || Array->meta.is("undirected"))
      throw std::logic_error("The -mutual- counter only works on directed (non-symmetric) arrays.");
    
    return 0.0;
  }
  
  COUNTER_FUNCTION_DEFAULT(count_mutual) {

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
  
  Counter<> mutual(count_mutual<>, init_mutual<>);
  
  // 2-istars
  COUNTER_FUNCTION_DEFAULT(count_istar2) {
   
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
  
    if (A_COL(j).size() == 1u)
      return 0.0;
    
    return ((double) A_COL(j).size() - 1.0);
   
    // return 0.0; 
  }

  Counter<> istar2(count_istar2<>);
  
  // 2-ostars
  COUNTER_FUNCTION_DEFAULT(count_ostar2) {
    
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
    
    if (A_ROW(i).size() == 1u)
      return 0.0;
    
    return ((double) A_ROW(i).size() - 1.0);
    
    // return 0.0; 
  }
  
  Counter<> ostar2(count_ostar2<>);
  
  // ttriads
  COUNTER_FUNCTION_DEFAULT(count_ttriads) {
    
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
  
  Counter<> ttriads(count_ttriads<>);
  
  // Cycle triads --------------------------------------------------------------
  COUNTER_FUNCTION_DEFAULT(count_ctriads) {
    
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
  
  Counter<> ctriads(count_ctriads<>);
  
  // Density --------------------------------------------------------------
  COUNTER_FUNCTION_DEFAULT (count_density) {
    return 1.0/(Array->N * (Array->M - 1));
  }
  
  Counter<> density(count_density<>);
  
  // idegree1.5  -------------------------------------------------------------
  COUNTER_FUNCTION_DEFAULT(count_idegree15) {
    
    // In case of the first, we need to add
    if (A_COL(j).size() == 1u)
      return 1.0;
    
    return 
      pow((double) A_COL(j).size(), 1.5) - pow((double) A_COL(j).size() - 1, 1.5)
    ;
  }
  
  Counter<> idegree15(count_idegree15<>);
  
  // odegree1.5  -------------------------------------------------------------
  COUNTER_FUNCTION_DEFAULT(count_odegree15) {
    
    // In case of the first, we need to add
    if (A_ROW(i).size() == 1u)
      return 1;
    
    return 
      pow((double) A_ROW(i).size(), 1.5) - pow((double) A_ROW(i).size() - 1, 1.5)
    ;
  }
  
  Counter<> odegree15(count_odegree15<>);
}

#endif
