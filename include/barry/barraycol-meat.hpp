#include "barraycol-bones.hpp"

#ifndef BARRY_BARRAYCOL_MEAT_HPP
#define BARRY_BARRAYCOL_MEAT_HPP 1

template<typename Cell_Type,typename Data_Type>
inline void BArrayCol<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value = val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCol<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value += val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCol<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, -val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value -= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCol<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    if (!Array->is_empty(i, j, false)) {
        Array->el_ij.at(i).at(j).value *= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCol<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    if (!Array->is_empty(i, j, false)) {
        Array->el_ij.at(i).at(j).value /= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline BArrayCol<Cell_Type,Data_Type>::operator Cell_Type() const {
        return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);  
}

template<typename Cell_Type,typename Data_Type>
inline Col_type<Cell_Type>::iterator BArrayCol<Cell_Type,Data_Type>::begin() noexcept {
    return this->Array->el_ji[this->i].begin();
}

template<typename Cell_Type,typename Data_Type>
inline Col_type<Cell_Type>::iterator BArrayCol<Cell_Type,Data_Type>::end() noexcept {
    return this->Array->el_ji[this->i].end();
}

/*******************************************************************************
 * Const Col 
 * ****************************************************************************/

template<typename Cell_Type,typename Data_Type>
inline BArrayCol_const<Cell_Type,Data_Type>::operator Cell_Type() const {
        return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator!=(const Cell_Type & val) const {
    return !(this->operator==(val));
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator<(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) < static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator>(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) > static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator<=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) <= static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCol_const<Cell_Type,Data_Type>::operator>=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) >= static_cast<Cell_Type>(val);    
}

#endif