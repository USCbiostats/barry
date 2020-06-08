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

/**@brief Sequence of cells indices.
 * 
 */
class CellSeq {
private:
  // Since we are using column-major order, the first is the column.
  std::map< std::pair<uint,uint> , bool > data; 
  
protected:
  uint N, M;
  std::pair<uint,uint> coords;
  std::vector< std::pair< uint , uint > > Seq; 
  bool seq_computed = false;
  
public:
  CellSeq(std::vector< std::pair< uint, uint> > & data_, uint N_, uint M_);
  CellSeq() : N(0u), M(0u) {};
  ~CellSeq() {};
  
  /***
   * Query functions
   */
  bool is_blocked(uint & i) ;
  bool is_blocked(uint & i, uint & j) ;
  
  /***
   * These functions compute the next or the previous one.
   */
  bool next(uint & i, uint & j) const;
  bool next(uint & i) const;
  // virtual bool prev(uint & i, uint & j) const;
  // virtual bool prev(uint & i) const;
  
  /***
   * Update functions
   */
  void compute_seq();
  
  /***
   * Utility functions
   */
  void calc_coords(uint & i) {
    
    if (i == 0) coords = {0u,0u};
    else if ((i - 1u) == N*M) coords = {N-1u, M-1u};
    else if (i >= N*M) {
      
      throw std::range_error("The coordinate is out of range.");
      
    } else {
      
      coords = {
        i % N,
        floor((int) i / (int) N)
      };
      
    }
    
    return;
  }
  
  /***
   * Getter function
   */
  const std::vector< std::pair< uint , uint > > * get_seq() ; 
  
};

inline CellSeq::CellSeq(
    std::vector< std::pair< uint, uint> > & data_, uint N_, uint M_
) : N(N_), M(M_){
  
  
  if (N < 1u || M < 1u)
    throw std::logic_error("Either zero rows or columns.");
  
  for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
    
    // Checking the ranges
    if (iter->first >= N)
      throw std::range_error("One element is out of range (row).");
    if (iter->second >= M)
      throw std::range_error("One element is out of range (col).");
    
    if (data.find(*iter) == data.end())
      data[*iter] = true;
  }
  
  return;
}

inline bool CellSeq::is_blocked(uint & i) {
  
  if (data.size() == 0u)
    return false;
  
  calc_coords(i);
  return data.find(coords) != data.end();
  
}

inline bool CellSeq::is_blocked(uint & i, uint & j) {
  
  if (data.size() == 0u)
    return false;
  
  return data.find({i, j}) != data.end();
  
}



inline bool CellSeq::next(uint & i, uint & j) const {
  
  // Have we reached the end already?
  
  // Iterating until this is out of the range
  auto iter = data.begin();
  uint safeward = N*M, step = 0u;
  while (true) {
    
    // Basic check
    if (safeward <= step++)
      throw std::logic_error("Upsi, the next function went overflown (shouldn't happen)!.");
    
    // Increasing the row, if possible
    if (++i >= N) {
      
      i = 0u; // Restarts the list, so we need to increase the column.
      
      // Need to check on j?
      if (++j >= M) { // If we reached this far, then we reached the end
        j--;
        i = N - 1u;
        return true;
      }
      
    }
    
    // Is it on the list?
    if (data.find({i, j}) == data.end()) 
      break;
    
    continue;
  }
  
  return false;
  
}


inline void CellSeq::compute_seq() {
  
  if (seq_computed)
    return;
  
  // Starting from zero
  coords = {0u, 0u};
  bool reached_end = false;
  while (is_blocked(coords.first, coords.second)) 
    reached_end = next(coords.first, coords.second);
  
  if (reached_end)
    return;
  
  // Adding the first to the sequence
  Seq.push_back(coords);
  
  while (!next(coords.first, coords.second)) 
    Seq.push_back(coords);
  
  return;
  
}

inline const std::vector< std::pair< uint, uint> > * CellSeq::get_seq() {
  
  // Updating (if needed)
  compute_seq();
  
  return &Seq;
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
