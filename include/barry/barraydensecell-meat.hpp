#include "barraydensecell-bones.hpp"

#ifndef BARRY_BARRAYDENSECELL_MEAT_HPP
#define BARRY_BARRAYDENSECELL_MEAT_HPP 1

#define POS(a, b) (a) + (b) * dat->N 

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {

    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    dat->el[POS(i,j)]  += val;
    dat->el_rowsums[i] += val;
    dat->el_colsums[j] += val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    dat->el[POS(i,j)]  -= val;
    dat->el_rowsums[i] -= val;
    dat->el_colsums[j] -= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    dat->el_colsums[j] += (dat->el[POS(i,j)] * val - dat->el[POS(i,j)]);
    dat->el[POS(i,j)]  += (dat->el[POS(i,j)] * val - dat->el[POS(i,j)]);
    dat->el_rowsums[i] *= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_rowsums[i] += (old/val - old);
    dat->el_colsums[j] += (old/val - old);
    dat->el[POS(i,j)]  /= val;

}

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type>::operator Cell_Type() const {
        return dat->el[POS(i,j)];
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return dat->el[POS(i,j)] == val;  
}

#undef POS

#endif