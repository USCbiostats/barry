#include "typedefs.hpp"

#ifndef BARRY_BARRAYVECTOR_BONES_HPP
#define BARRY_BARRAYVECTOR_BONES_HPP 1

/**
 * @brief Row or column of a `BArray`
 * 
 * @tparam Cell_Type 
 * @tparam Data_Type 
 */
template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayVector {
private:
  
    BArray<Cell_Type,Data_Type> * Array;
    std::vector< std::pair< uint, Cell_Type > > vec;
    uint dim;
    uint i;

    void init_vec();
    bool vec_initialized = false;
  
public:
  
    /**
     * @brief Construct a new BArrayVector object
     * 
     * @param Array_ Pointer to a `BArray` object
     * @param dim_ Dimension. 0 means row and 1 means column.
     * @param i_ Element to point.
     * @param check_bounds When `true`, check boundaries.
     */
    BArrayVector(
        BArray<Cell_Type,Data_Type> * Array_,
        uint & dim_
        uint & i_,
        bool check_bounds = true
        ) : 
    Array(Array_), vec(0u), dim(dim_), i(i_) {

        if (dim > 1u)
            throw std::range_error("-dim_- should be either 0 (row) or 1 (col).");

        if (check_bounds) {

            if ((dim == 0u) && (i >= Array->nrow()))
                throw std::length_error("Row out of range.");
            if ((dim == 1u) && (j >= Array->ncol()))
                throw std::length_error("Col out of range.");

        }
    };

    ~BArrayVector() {};

    bool is_row() const noexcept;
    bool is_col() const noexcept;
    uint size() const noexcept;
    std::vector< Cell_Type >::const_iterator begin() noexcept;
    std::vector< Cell_Type >::const_iterator end() noexcept;

    void operator=(const Cell_Type & val);
    void operator+=(const Cell_Type & val);
    void operator-=(const Cell_Type & val);
    void operator*=(const Cell_Type & val);
    void operator/=(const Cell_Type & val);

    operator std::vector< Cell_Type >() const;
    bool operator==(const Cell_Type & val) const;
  
};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayVector_const {
private:
    
    const BArray<Cell_Type,Data_Type> * Array;
    std::vector< std::pair< uint, Cell_Type > > vec;
    uint dim;
    uint i;

    void init_vec();
    bool vec_initialized = false;
    
public:
  
    BArrayVector_const(
        const BArray<Cell_Type,Data_Type> * Array_,
        uint & dim_
        uint & i_,
        bool check_bounds = true
        ) : 
    Array(Array_), vec(0u), dim(dim_), i(i_) {

        if (dim > 1u)
            throw std::range_error("-dim_- should be either 0 (row) or 1 (col).");

        if (check_bounds) {

            if ((dim == 0u) && (i >= Array->nrow()))
                throw std::length_error("Row out of range.");
            if ((dim == 1u) && (j >= Array->ncol()))
                throw std::length_error("Col out of range.");

        }

    };
    
    ~BArrayVector_const() {};

    bool is_row() const noexcept;
    bool is_col() const noexcept;
    uint size() const noexcept;
    std::vector< Cell_Type >::const_iterator begin() noexcept;
    std::vector< Cell_Type >::const_iterator end() noexcept;
    
    operator std::vector<Cell_Type>() const;
    bool operator==(const Cell_Type & val) const;
    bool operator!=(const Cell_Type & val) const;
    bool operator<(const Cell_Type & val) const;
    bool operator>(const Cell_Type & val) const;
    bool operator<=(const Cell_Type & val) const;
    bool operator>=(const Cell_Type & val) const;
  
};

#endif