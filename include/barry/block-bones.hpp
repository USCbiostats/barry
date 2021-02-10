// #include <vector>
// #include <stdexcept>

#include "typedefs.hpp"

#ifndef BARRY_BLOCK_BONES_HPP
#define BARRY_BLOCK_BONES_HPP 1

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


#endif
