#include "typedefs.hpp"

#ifndef BARRY_BARRAYDENSECELL_BONES_HPP
#define BARRY_BARRAYDENSECELL_BONES_HPP 1

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCell {
private:
  
    BArrayDense<Cell_Type,Data_Type> * Array;
    uint i;
    uint j;
  
public:
  
    BArrayDenseCell(BArrayDense<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {
        if (check_bounds) {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }
    };

    ~BArrayDenseCell(){};
    void operator=(const Cell_Type & val);
    void operator+=(const Cell_Type & val);
    void operator-=(const Cell_Type & val);
    void operator*=(const Cell_Type & val);
    void operator/=(const Cell_Type & val);

    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
  
};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCell_const {
private:
    
    const BArrayDense<Cell_Type,Data_Type> * Array;
    uint i;
    uint j;
    
public:
  
    BArrayDenseCell_const(const BArrayDense<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {
        if (check_bounds) {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }
    };
    
    ~BArrayDenseCell_const(){};
    
    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
    bool operator!=(const Cell_Type & val) const;
    bool operator<(const Cell_Type & val) const;
    bool operator>(const Cell_Type & val) const;
    bool operator<=(const Cell_Type & val) const;
    bool operator>=(const Cell_Type & val) const;
  
};

#endif