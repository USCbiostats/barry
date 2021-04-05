#include "barraycell-bones.hpp"

#ifndef BARRY_BARRAYCELL_MEAT_HPP
#define BARRY_BARRAYCELL_MEAT_HPP 1

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {
  
  if (Array->is_empty(i, j, false)) {
    Array->insert_cell(i, j, val, false, false);
  } else {
    Array->el_ij.at(i).at(j).value = val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
  
  if (Array->is_empty(i, j, false)) {
    Array->insert_cell(i, j, val, false, false);
  } else {
    Array->el_ij.at(i).at(j).value += val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
  
  if (Array->is_empty(i, j, false)) {
    Array->insert_cell(i, j, -val, false, false);
  } else {
    Array->el_ij.at(i).at(j).value -= val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
  
  if (!Array->is_empty(i, j, false)) {
    Array->el_ij.at(i).at(j).value *= val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
  
  if (!Array->is_empty(i, j, false)) {
    Array->el_ij.at(i).at(j).value /= val;
  }

}

template<typename Cell_Type,typename Data_Type>
inline BArrayCell<Cell_Type,Data_Type>::operator Cell_Type() const {
    return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
  return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);  
}

template<typename Cell_Type,typename Data_Type>
inline BArrayCell_const<Cell_Type,Data_Type>::operator Cell_Type() const {
    return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
  return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator!=(const Cell_Type & val) const {
  return !(this->operator==(val));
}

#endif