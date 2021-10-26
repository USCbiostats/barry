// #include <stdexcept>
#include "barray-bones.hpp"

#ifndef BARRY_BARRAY_MEAT_OPERATORS_HPP
#define BARRY_BARRAY_MEAT_OPERATORS_HPP 1

#define BARRAY_TYPE() BArray<Cell_Type, Data_Type>

#define BARRAY_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BARRAY_TEMPLATE(a,b) \
    template BARRAY_TEMPLATE_ARGS() inline a BARRAY_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]

template BARRAY_TEMPLATE_ARGS()
inline void checkdim_(
    const BARRAY_TYPE()& lhs,
    const BARRAY_TYPE()& rhs
) {

    if (lhs.ncol() != rhs.ncol())
        throw std::length_error("Number of columns do not match.");

    if (lhs.nrow() != rhs.nrow())
        throw std::length_error("Number of rows do not match.");

    return;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator+=) (
    const BArray<Cell_Type, Data_Type>& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (uint i = 0u; i < nrow(); ++i)
        for (uint j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator-=) (
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

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator*=) (
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

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator/=) (
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

#undef BARRAY_TYPE
#undef BARRAY_TEMPLATE_ARGS
#undef BARRAY_TEMPLATE

#undef ROW
#undef COL

#endif