// #include<vector>
// #include <stdexcept>
// #include "typedefs.hpp"

/**
 * @brief Dense bi-dimensional array
 * 
 * @details elements is stored in a std::vector, in col-major order.
 * 
 * @tparam Cell_Type 
 */
template<typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDense {
private:
    unsigned int N; ///< Number of rows
    unsigned int M; ///< Number of columns
    std::vector< Cell_Type > elements;
    Data_Type * data = nullptr;

public:
    BArrayDense();
    BArrayDense(
        uint N_, uint M_,
        std::vector< Cell_Type > elements_ = {}
        );
    ~BArrayDense() {};

    void fill(const Cell_Type & d);
    
    const std::vector< Cell_Type > & elements_raw() const noexcept;
    const std::vector< Cell_Type > * elements_ptr() const noexcept;
    uint nrow() const noexcept;
    uint ncol() const noexcept;

    Cell_Type operator()(uint i, uint j, bool check_bounds = true) const;
    Cell_Type operator[](uint i) const;

    void print() const;

};

template<typename Cell_Type,typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type>::BArrayDense()
{

    N = 0u;
    M = 0u;
    elements.resize(0u);

    elements.size()

}

template<typename Cell_Type,typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type>::BArrayDense(
    uint N_, uint M_,
    std::vector< Cell_Type > elements_
) : N(N_), M(M_)
{

    if (elements_->size() != static_cast<uint>(N*M))
        throw std::logic_error("N*M don't match thedim of the elements.");

    elements.resize(elements_.size());
    std::copy(elements_.begin(), elements_.end(), elements.begin());
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDense<Cell_Type,Data_Type>::fill(const Cell_Type & d)
{
    std::fill(elements.begin(), elements.end(), d);
}


template<typename Cell_Type,typename Data_Type>
inline const std::vector< Cell_Type > &
    BArrayDense<Cell_Type,Data_Type>::elements_raw() const noexcept
{
    return elements;
}
 
template<typename Cell_Type,typename Data_Type>
inline const std::vector< Cell_Type > *
    BArrayDense<Cell_Type,Data_Type>::elements_ptr() const noexcept
{
    return *elements;
}

template<typename Cell_Type,typename Data_Type>
inline uint BArrayDense<Cell_Type,Data_Type>::nrow() const noexcept
{
    return N;
}

template<typename Cell_Type,typename Data_Type>
inline uint BArrayDense<Cell_Type,Data_Type>::ncol() const noexcept
{
    return M;
}

template<typename Cell_Type,typename Data_Type>
inline Cell_Type BArrayDense<Cell_Type,Data_Type>::operator()(
    uint i, uint j, bool check_bounds
) const 
{
    if (check_bounds)
    {

        if (i >= N)
            throw std::range_error("Row index i out of range.");
        else if (j >= M)
            throw std::range_error("Col index j out of range.");

    }

    return elements[j * N + i];

}

template<typename Cell_Type,typename Data_Type>
inline Cell_Type BArrayDense<Cell_Type,Data_Type>::operator[](uint i) const 
{           
    return elements[i];
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDense<Cell_Type,Data_Type>::print() const {

    if (N*M == 0u)
        printf_barry("< empty dense array >\n");

    for (uint i = 0u; i < N; ++i)
    {
        printf_barry("[%3i,]", i);
        for (uint j = 0u; j < M; ++j)
            printf_barry(" %.2f ", elements[j*N + i]);

        printf_barry("\n");
    }

}