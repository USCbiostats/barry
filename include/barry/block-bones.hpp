// #include <vector>
// #include <stdexcept>

#include "typedefs.hpp"

#ifndef BARRAY_BLOCK_BONES_HPP
#define BARRAY_BLOCK_BONES_HPP 1




/***
 * Various ways to describe blocks of cells that are static. We have the
 * following ways:
 * - Sequence of coordinates
 * - Ranges
 * - Blocks
 * - Diagonal
 * - Upper or lower diagonal
 * - Rows
 * - Columns
 * 
 * The base class is call CellSeq and it should generally be able to have
 * the following functions.
 * 
 * - is_blocked(i, j), or is_blocked(i)
 * - operator-(): To substract from a full sequence.
 *
 * These should be encapsulated in a larger class called Constrains, which
 * is a collection of CellSeqs
 * 
 * - next_coord(i, j), or next_coord(i)
 * 
 * With the following common members
 * 
 * - N,
 * - M
 * 
 * 
 */

template <typename Array_Type, typename Data_Type>
inline bool default_rule(const Array_Type * A, uint i, uint j, Data_Type *) {
  return false;
}

/**@brief Sequence of cells indices.
 * 
 */
template<typename Array_Type, typename Data_Type>
class CellSeq {
protected:
  // Input data
  Rule_fun_type<Array_Type, Data_Type> rule = nullptr;
  Data_Type * rule_data                     = nullptr;

public:
  CellSeq() : rule(default_rule<Array_Type,Data_Type>) {};
  CellSeq(Rule_fun_type<Array_Type,Data_Type> rule_): rule(rule_) {};
  ~CellSeq() {};
  
  /***
   * Query functions
   */
  // bool is_blocked(const Array_Type * dat, uint & i) ;
  bool is_blocked(const Array_Type * A, uint i, uint j) ;
  
  /***
   * Getter function
   */
  void get_seq(
      std::vector< std::pair< uint , uint > > * seq,
      const Array_Type * A,
      bool invert = false
    ) ; 
  
};

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
