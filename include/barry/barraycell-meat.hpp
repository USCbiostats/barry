#include "barraycell-bones.hpp"

#ifndef BARRY_BARRAYCELL_MEAT_HPP
#define BARRY_BARRAYCELL_MEAT_HPP 1

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {
  
  if (Array->is_empty(row, col, true)) {
    Array->insert_cell(row, col, val, false, false);
  } else {
    Array->el_ij.at(row).at(col).value = val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
  
  if (Array->is_empty(row, col, true)) {
    Array->insert_cell(row, col, val, false, false);
  } else {
    Array->el_ij.at(row).at(col).value += val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
  
  if (Array->is_empty(row, col, true)) {
    Array->insert_cell(row, col, -val, false, false);
  } else {
    Array->el_ij.at(row).at(col).value -= val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
  
  if (!Array->is_empty(row, col, true)) {
    Array->el_ij.at(row).at(col).value *= val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
  
  if (!Array->is_empty(row, col, true)) {
    Array->el_ij.at(row).at(col).value /= val;
  }

}

#endif