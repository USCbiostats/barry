#ifndef BARRY_BARRAYROW_BONES_HPP
#define BARRY_BARRAYROW_BONES_HPP 1

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayRow {
private:
  
    BArray<Cell_Type,Data_Type> * Array;
    size_t i;
  
public:
  
    BArrayRow(BArray<Cell_Type,Data_Type> * Array_, size_t i_,, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {

        if (check_bounds)
        {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");

        }

    };

    ~BArrayRow(){};
    void operator=(const BArrayRow<Cell_Type,Data_Type> & val);
    void operator+=(const BArrayRow<Cell_Type,Data_Type> & val);
    void operator-=(const BArrayRow<Cell_Type,Data_Type> & val);
    void operator*=(const BArrayRow<Cell_Type,Data_Type> & val);
    void operator/=(const BArrayRow<Cell_Type,Data_Type> & val);

    operator BArrayRow<Cell_Type,Data_Type>() const;
    bool operator==(const BArrayRow<Cell_Type,Data_Type> & val) const;
  
};



template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayRow_const {
private:
    
    const BArray<Cell_Type,Data_Type> * Array;
    size_t i;
    
public:
  
    BArrayRow_const(const BArray<Cell_Type,Data_Type> * Array_, size_t i_, bool check_bounds = true) : 
    Array(Array_), i(i_), {
        if (check_bounds) {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");

        }
    };
    
    ~BArrayRow_const(){};
    
    operator BArrayRow_const<Cell_Type,Data_Type>() const;
    bool operator==(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
    bool operator!=(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
    bool operator<(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
    bool operator>(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
    bool operator<=(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
    bool operator>=(const BArrayRow_const<Cell_Type,Data_Type> & val) const;
  
};

#endif