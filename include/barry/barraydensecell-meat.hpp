#include "barraydensecell-bones.hpp"

#ifndef BARRY_BARRAYDENSECELL_MEAT_HPP
#define BARRY_BARRAYDENSECELL_MEAT_HPP 1

#define POS(a, b) (a) + (b) * Array->N 

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {

    (*dat) = val;
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    (*dat) += val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    (*dat) -= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    (*dat) *= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    (*dat) /= val;

}

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type>::operator Cell_Type() const {
        return (*dat);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return (*dat) == val;  
}

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell_const<Cell_Type,Data_Type>::operator Cell_Type() const {
        return dat;
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return dat == val;
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator!=(const Cell_Type & val) const {
    return !(this->operator==(val));
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator<(const Cell_Type & val) const {
    return dat < val;    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator>(const Cell_Type & val) const {
    return dat > val;    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator<=(const Cell_Type & val) const {
    return dat <= val;    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell_const<Cell_Type,Data_Type>::operator>=(const Cell_Type & val) const {
    return dat >= val;    
}

class NetworkData;


template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator==(const int & val) const {
    return dat == static_cast<int>(val);
}

template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator!=(const int & val) const {
    return !(this->operator==(static_cast<int>(val)));
}

template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator<(const int & val) const {
    return dat < static_cast<int>(val);    
}

template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator>(const int & val) const {
    return dat > static_cast<int>(val);    
}

template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator<=(const int & val) const {
    return dat <= static_cast<int>(val);    
}

template<>
inline bool BArrayDenseCell_const<int,NetworkData>::operator>=(const int & val) const {
    return dat >= static_cast<int>(val);    
}


#undef POS

#endif