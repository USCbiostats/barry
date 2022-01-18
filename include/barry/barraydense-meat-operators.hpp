// #include <stdexcept>
#include "barraydense-bones.hpp"

#ifndef BARRY_BARRAYDENSE_MEAT_OPERATORS_HPP
#define BARRY_BARRAYDENSE_MEAT_OPERATORS_HPP 1

#define BDENSE_TYPE() BArrayDense<Cell_Type, Data_Type>

#define BDENSE_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BDENSE_TEMPLATE(a,b) \
    template BDENSE_TEMPLATE_ARGS() inline a BDENSE_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template BDENSE_TEMPLATE_ARGS()
inline void checkdim_(
    const BDENSE_TYPE()& lhs,
    const BDENSE_TYPE()& rhs
) {

    if (lhs.ncol() != rhs.ncol())
        throw std::length_error("Number of columns do not match.");

    if (lhs.nrow() != rhs.nrow())
        throw std::length_error("Number of rows do not match.");

    return;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator+=) (
    const BDENSE_TYPE()& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (uint i = 0u; i < nrow(); ++i)
        for (uint j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs;
        }
    }

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator-=) (
    const BDENSE_TYPE()& rhs
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

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < ncol(); ++j) 
            this->operator()(i, j) -= rhs;
        
    

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator*=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < nrow(); ++j)
            el[POS(i, j)] *= rhs;

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator/=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < nrow(); ++j)
            el[POS(i, j)] /= rhs;

    return *this;
}

#undef BDENSE_TYPE
#undef BDENSE_TEMPLATE_ARGS
#undef BDENSE_TEMPLATE

#undef ROW
#undef COL
#undef POS
#undef POS_N

#endif