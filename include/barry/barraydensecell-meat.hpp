// #include "barraydensecell-bones.hpp"

#ifndef BARRY_BARRAYDENSECELL_MEAT_HPP
#define BARRY_BARRAYDENSECELL_MEAT_HPP 1

#define POS(a, b) (a) + (b) * dat->N 

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type>& BArrayDenseCell<Cell_Type,Data_Type>::operator=(
    const BArrayDenseCell<Cell_Type,Data_Type> & other
    ) {
    
    Cell_Type val = static_cast<Cell_Type>(other);
    #ifdef BARRY_DEBUG
    Cell_Type old      =  dat->el.at(POS(i,j));
    dat->el.at(POS(i,j))  =  val;
    dat->el_rowsums.at(i) += (val - old);
    dat->el_colsums.at(j) += (val - old);
    #else
    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    #endif

    return *this;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {

    #ifdef BARRY_DEBUG
    Cell_Type old      =  dat->el.at(POS(i,j));
    dat->el.at(POS(i,j))  =  val;
    dat->el_rowsums.at(i) += (val - old);
    dat->el_colsums.at(j) += (val - old);
    #else
    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    #endif
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    dat->el.at(POS(i,j))  += val;
    dat->el_rowsums.at(i) += val;
    dat->el_colsums.at(j) += val;
    #else
    dat->el[POS(i,j)]  += val;
    dat->el_rowsums[i] += val;
    dat->el_colsums[j] += val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    dat->el.at(POS(i,j))  -= val;
    dat->el_rowsums.at(i) -= val;
    dat->el_colsums.at(j) -= val;
    #else
    dat->el[POS(i,j)]  -= val;
    dat->el_rowsums[i] -= val;
    dat->el_colsums[j] -= val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    Cell_Type old = dat->el.at(POS(i,j));
    dat->el_colsums.at(j) += (old * val - old);
    dat->el_rowsums.at(i) += (old * val - old);
    dat->el.at(POS(i,j)) *= val;
    #else
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_colsums[j] += (old * val - old);
    dat->el_rowsums[i] += (old * val - old);
    dat->el[POS(i,j)] *= val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    Cell_Type old = dat->el.at(POS(i,j));
    dat->el_rowsums.at(i) += (old/val - old);
    dat->el_colsums.at(j) += (old/val - old);
    dat->el.at(POS(i,j))  /= val;
    #else
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_rowsums[i] += (old/val - old);
    dat->el_colsums[j] += (old/val - old);
    dat->el[POS(i,j)]  /= val;
    #endif

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