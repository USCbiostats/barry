#include "typedefs.hpp"

#ifndef BARRY_BARRAYDENSECELL_BONES_HPP
#define BARRY_BARRAYDENSECELL_BONES_HPP 1

#define POS(a, b) (a) + (b) * N

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCell {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCol<Cell_Type,Data_Type>;
    friend class BArrayDenseCol_const<Cell_Type,Data_Type>;
private:
  
    BArrayDense<Cell_Type,Data_Type> * dat;
    uint i;
    uint j;
  
public:
  
    BArrayDenseCell(
        BArrayDense<Cell_Type,Data_Type> * Array_,
        uint i_,
        uint j_,
        bool check_bounds = true
        ) : 
    i(i_), j(j_)
    {

        if (check_bounds)
        {

            if (i >= Array_->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array_->ncol())
                throw std::length_error("Col out of range.");

        }
        dat = Array_;

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


#undef POS

#endif