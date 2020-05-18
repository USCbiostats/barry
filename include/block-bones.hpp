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

class CellSeq;
class Constraints {
public:
  friend class CellSeq;
  uint N, M;
  std::vector< CellSeq > S;
  
  Constraints() : S(0u) {};
  
  // Key functions
  void add(const CellSeq & x) {S.push_back(x); return;};
  std::vector< std::pair< uint, uint > > get_seq() const;
  std::pair< uint, uint > next_coord(uint i, uint j) const;
  
  
  
};

class CellSeq {
  
};

class CellRange : public CellSeq {
  
};

class CellDiag : public CellSeq {
  
};



#endif