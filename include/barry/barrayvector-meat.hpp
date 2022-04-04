#ifndef BARRY_BARRAYVECTOR_MEAT_HPP
#define BARRY_BARRAYVECTOR_MEAT_HPP 1

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::init_vec() {

    if (vec_initialized)
        return;

    if (dim == 0u)
    {

        for (const auto& a : Array->el_ij[i])
            vec.push_back(a);
            
    } else {

        for (const auto& a : Array->el_ji[i])
            vec.push_back(std::make_pair<uint, Cell<Cell_Type>>(a.first, *(a.second)));

    }

    vec_initialized = true;

    return;
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector<Cell_Type,Data_Type>::is_row() const noexcept {
    return dim == 0;
}

template<typename Cell_Type, typename Data_Type>
inline bool BArrayVector<Cell_Type,Data_Type>::is_col() const noexcept {
    return dim == 1;
}

template<typename Cell_Type, typename Data_Type>
inline uint BArrayVector<Cell_Type,Data_Type>::size() const noexcept {

    if (dim == 0u)
        return Array->el_ij[i].size();
    else
        return Array->el_ji[i].size();
    

}

template<typename Cell_Type, typename Data_Type>
inline std::vector< Cell_Type >::const_iterator BArrayVector<Cell_Type,Data_Type>::begin() noexcept {
    
    // For this, we will need the iterator
    init_vec();

    if (dim = 0u)
    {

    } else {

    }
}

template<typename Cell_Type, typename Data_Type>
inline std::vector< Cell_Type >::const_iterator BArrayVector<Cell_Type,Data_Type>::end() noexcept {

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {
    
    uint k = 0u;
    uint N_ = (dim == 0u) ? Array->nrow() : Array->ncol();
    
    if (dim == 0u)
    {

        for (auto j = 0u; j < N_; ++j)
            Array(i, j) = val;

    } else {

        for (auto j = 0u; j < N_; ++j)
            Array(j, i) = val;

    }


}

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    uint k = 0u;
    uint N_ = (dim == 0u) ? Array->nrow() : Array->ncol();
    
    if (dim == 0u)
    {

        for (auto j = 0u; j < N_; ++j)
            Array(i, j) += val;

    } else {

        for (auto j = 0u; j < N_; ++j)
            Array(j, i) += val;

    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    uint k = 0u;
    uint N_ = (dim == 0u) ? Array->nrow() : Array->ncol();
    
    if (dim == 0u)
    {

        for (auto j = 0u; j < N_; ++j)
            Array(i, j) -= val;

    } else {

        for (auto j = 0u; j < N_; ++j)
            Array(j, i) -= val;

    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    uint k = 0u;
    uint N_ = (dim == 0u) ? Array->nrow() : Array->ncol();
    
    if (dim == 0u)
    {

        for (auto j = 0u; j < N_; ++j)
            Array(i, j) *= val;

    } else {

        for (auto j = 0u; j < N_; ++j)
            Array(j, i) *= val;

    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayVector<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    uint k = 0u;
    uint N_ = (dim == 0u) ? Array->nrow() : Array->ncol();
    
    if (dim == 0u)
    {

        for (auto j = 0u; j < N_; ++j)
            Array(i, j) /= val;

    } else {

        for (auto j = 0u; j < N_; ++j)
            Array(j, i) /= val;

    }

}

template<typename Cell_Type,typename Data_Type>
inline BArrayVector<Cell_Type,Data_Type>::operator std::vector< Cell_Type >() const {

    if (dim == 0u)
        return Array->get_row_vec(i, false);
    else
        return Array->get_col_vec(i, false);
        
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {

    if (dim == 0u)
    {
        for (uint j = 0u; j < Array->ncol(); ++j)
        {
            if (Array(i, j) != val)
                return false;

        }

    } else {

        for (uint j = 0u; j < Array->nrow(); ++j)
        {
            if (Array(j, i) != val)
                return false;

        }

    }

    return true;
    
}

template<typename Cell_Type,typename Data_Type>
inline BArrayVector_const<Cell_Type,Data_Type>::operator std::vector< Cell_Type >() const {

    if (dim == 0u)
        return Array->get_row_vec(i, false);
    else
        return Array->get_col_vec(i, false);
        
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    
    if (dim == 0u)
    {
        for (uint j = 0u; j < Array->ncol(); ++j)
        {
            if (Array(i, j) != val)
                return false;

        }

    } else {

        for (uint j = 0u; j < Array->nrow(); ++j)
        {
            if (Array(j, i) != val)
                return false;

        }

    }

    return true;

}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator!=(const Cell_Type & val) const {
    return !(this->operator==(val));
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator<(const Cell_Type & val) const {
    
    if (dim == 0u)
    {
        for (uint j = 0u; j < Array->ncol(); ++j)
        {
            if (Array(i, j) >= val)
                return false;

        }

    } else {

        for (uint j = 0u; j < Array->nrow(); ++j)
        {
            if (Array(j, i) >= val)
                return false;

        }

    }

    return true;

}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator<=(const Cell_Type & val) const {
    
    if (dim == 0u)
    {
        for (uint j = 0u; j < Array->ncol(); ++j)
        {
            if (Array(i, j) > val)
                return false;

        }

    } else {

        for (uint j = 0u; j < Array->nrow(); ++j)
        {
            if (Array(j, i) > val)
                return false;

        }

    }

    return true;

}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator>(const Cell_Type & val) const {
    return !(this->operator<=(val));
}



template<typename Cell_Type,typename Data_Type>
inline bool BArrayVector_const<Cell_Type,Data_Type>::operator>=(const Cell_Type & val) const {
    return !(this->operator<(val));    
} 

#endif