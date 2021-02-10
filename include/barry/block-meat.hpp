#include "block-bones.hpp"

#ifndef BARRY_BLOCK_MEAT_HPP
#define BARRY_BLOCK_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline bool CellSeq<Array_Type,Data_Type>::is_blocked(
    const Array_Type * A,
    uint i, uint j
  ) {

  return rule(A, i, j, rule_data);
  
}

template <typename Array_Type, typename Data_Type>
inline void CellSeq<Array_Type, Data_Type>::get_seq(
    std::vector< std::pair< uint, uint> > * seq,
    const Array_Type * A,
    bool invert
) {
  
  uint N = A->nrow();
  uint K = A->ncol();
  
  if (!invert) {
    for (uint i = 0u; i < N; ++i) {
      for (uint j = 0u; j < K; ++j) {
        if (!this->is_blocked(A, i, j))
          seq->push_back({i, j});
      }
    }
  } else {
    for (uint i = 0u; i < N; ++i) {
      for (uint j = 0u; j < K; ++j) {
        if (this->is_blocked(A, i, j))
          seq->push_back({i, j});
      }
    }
  }
  
  
  return;
}



// Common rules

// Zero diag
template <typename Array_Type, typename Data_Type>
inline bool sequence_rule_zero_diag(const Array_Type * A, uint i, uint j, Data_Type *) {
  
  return i == j;
}

template <typename Array_Type, typename Data_Type>
inline bool sequence_rule_lower_tri(const Array_Type * A, uint i, uint j, Data_Type *) {
  
  return i < j;
}

template <typename Array_Type, typename Data_Type>
inline bool sequence_rule_upper_tri(const Array_Type * A, uint i, uint j, Data_Type *) {
  
  return i > j;
  
}

/*

// Decendent class: Range ------------------------------------------------------
class CellSeqRange : public CellSeq {
  
  uint start, end;
  uint start_i, start_j, end_i, end_j;
  
  CellSeqRange() {};
  ~CellSeqRange() {};
  
  uint next(uint i) const;
  
};

uint CellSeqRange::next(uint i) const {
  
  uint i_next = 1 + 1u;
  // Is it out of range?
  if ((i_next < start) || (i_next > end)) {
    
    // Checking range again
    if (i_next < (N * M))
      return i_next;
    else
      return 0u;
    
  } 
  
  // Need to skip to the next step
  if (++i_next < (N * M))
    return i_next;
  
  return 0u;
}

// Decendent class: Row --------------------------------------------------------
class CellSeqRow: public CellSeq {
  
  uint row;
  
  CellSeqRow() {};
  ~CellSeqRow() {};
  
};

// Decendent class: Col --------------------------------------------------------
class CellSeqCol: public CellSeq {
  
  uint col;
  
  CellSeqCol() {};
  ~CellSeqCol() {};
  
};

// Decendent class: Diag -------------------------------------------------------
class CellSeqDiag : public CellSeq {
  
  uint row0,col0;
  uint row1,col1;
  
  CellSeqDiag() {};
  ~CellSeqDiag() {};
  
};

*/
#endif