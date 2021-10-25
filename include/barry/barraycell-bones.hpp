#include "typedefs.hpp"

#ifndef BARRY_BARRAYCELL_BONES_HPP
#define BARRY_BARRAYCELL_BONES_HPP 1

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayCell {
private:
  
    BArray<Cell_Type,Data_Type> * Array;
    uint i;
    uint j;
  
public:
  
    BArrayCell(BArray<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {

        if (check_bounds)
        {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }

    };

    ~BArrayCell(){};
    void operator=(const Cell_Type & val);
    void operator+=(const Cell_Type & val);
    void operator-=(const Cell_Type & val);
    void operator*=(const Cell_Type & val);
    void operator/=(const Cell_Type & val);

    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
  
};



template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayCell_const {
private:
    
    const BArray<Cell_Type,Data_Type> * Array;
    uint i;
    uint j;
    
public:
  
    BArrayCell_const(const BArray<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {
        if (check_bounds) {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }
    };
    
    ~BArrayCell_const(){};
    
    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
    bool operator!=(const Cell_Type & val) const;
    bool operator<(const Cell_Type & val) const;
    bool operator>(const Cell_Type & val) const;
    bool operator<=(const Cell_Type & val) const;
    bool operator>=(const Cell_Type & val) const;
  
};

#endif