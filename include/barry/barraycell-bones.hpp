#include "typedefs.hpp"

#ifndef BARRY_BARRAYCELL_BONES_HPP
#define BARRY_BARRAYCELL_BONES_HPP 1

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayCell {
private:
  
  BArray<Cell_Type,Data_Type> * Array;
  uint row;
  uint col;
  
public:
  
  BArrayCell(BArray<Cell_Type,Data_Type> * Array_, uint row_, uint col_) : 
  Array(Array_), row(row_), col(col_) {};
  ~BArrayCell(){};
  void operator=(const Cell_Type & val);
  void operator+=(const Cell_Type & val);
  void operator-=(const Cell_Type & val);
  void operator*=(const Cell_Type & val);
  void operator/=(const Cell_Type & val);
  
};

#endif