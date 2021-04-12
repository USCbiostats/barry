// #include <stdexcept>
#include "barray-bones.hpp"

#ifndef BARRY_BARRAY_MEAT_OPERATORS_HPP
#define BARRY_BARRAY_MEAT_OPERATORS_HPP 1

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]

template<typename Cell_Type, typename Data_Type>
inline void checkdim_(
    const BArray<Cell_Type, Data_Type>& lhs,
    const BArray<Cell_Type, Data_Type>& rhs
) {

    if (lhs.ncol() != rhs.ncol())
        throw std::length_error("Number of columns do not match.");

    if (lhs.nrow() != rhs.nrow())
        throw std::length_error("Number of rows do not match.");

    return;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator+=(
    const BArray<Cell_Type, Data_Type>& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs.get_cell(i, j);
        }
    }

    return *this;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator+=(
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs;
        }
    }

    return *this;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator-=(
    const BArray<Cell_Type, Data_Type>& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs.get_cell(i, j);
        }
    }

    return *this;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator-=(
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs;
        }
    }

    return *this;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator*=(
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {

        if (ROW(i).size() == 0u)
            continue;

        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            this->operator()(i, col->first) *= rhs;
        }
    }

    return *this;
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>& BArray<Cell_Type, Data_Type>::operator/=(
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {

        if (ROW(i).size() == 0u)
            continue;

        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            this->operator()(i, col->first) /= rhs;
        }
    }

    return *this;
}

#undef ROW
#undef COL

#endif