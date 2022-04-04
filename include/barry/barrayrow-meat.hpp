#ifndef BARRY_BARRAYROW_MEAT_HPP
#define BARRY_BARRAYROW_MEAT_HPP 1

#define BROW_TYPE() BArrayRow<Cell_Type, Data_Type>

#define BROW_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BROW_TEMPLATE(a,b) \
    template BROW_TEMPLATE_ARGS() inline a BROW_TYPE()::b

BROW_TEMPLATE(void, operator=) (const BROW_TYPE() & val) {
    
    // First, zeroout the row
    this->Array->zero_row(j);

    // Then, iterate throught the values of val and add it
    for (auto & v: val)
        Array->inser_cell(i, v.first, v.second);

    // Return
    return;

}

BROW_TEMPLATE(void, operator+=) (const BROW_TYPE() & val) {
    
    for (auto & v : val)
        this->Array->operator(i, v.first) += v.second;
    
    return;

}

BROW_TEMPLATE(void, operator-=) (
    const BROW_TYPE() & val
) {
    
    for (auto & v : val)
        this->Array->operator(i, v.first) -= v.second;
    
    return;

}

BROW_TEMPLATE(void, operator*=) (
    const BROW_TYPE() & val
) {
    
    if (!Array->is_empty(i, j, false)) {
        Array->el_ij.at(i).at(j).value *= val;
    }

}

BROW_TEMPLATE(void, operator/=) (
    const BROW_TYPE() & val
) {
    
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

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator<(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) < static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator>(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) > static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator<=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) <= static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator>=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) >= static_cast<Cell_Type>(val);    
}

#undef BROW_TYPE
#undef BROW_TEMPLATE_ARGS
#undef BROW_TEMPLATE

#endif