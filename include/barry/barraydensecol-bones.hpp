#ifndef BARRY_BARRAYDENSECOL_BONES 
#define BARRY_BARRAYDENSECOL_BONES

#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)
#define ZERO_CELL static_cast<Cell_Type>(0.0)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    Col_type<Cell_Type> col;
    size_t index;
    bool col_filled = false;

    void fill_if_needed()
    {
        if (!col_filled)
        {

            for (size_t i = 0u; i < array->N; ++i)
            {
                
                if (array->el[POS_N(i, index, array->N)] != ZERO_CELL)
                    col[i] = col[POS_N(i, index, array->N)];
                    
            }

            col_filled = true;
            
        }
    }

public:
    BArrayDenseCol(
        BArrayDense< Cell_Type,Data_Type > & array_,
        size_t j
    ) : array(&array_), index(j) {};


    typename Col_type<Cell_Type>::iterator & begin()
    {
        fill_if_needed();
        return col.begin();
    };

    typename Col_type<Cell_Type>::iterator & end()
    {
        fill_if_needed();
        return col.end();
    };

    size_t size() const noexcept
    {
        fill_if_needed();
        return col.size();
    };

    std::pair<size_t,Cell_Type*> & operator()(size_t i)
    {
        fill_if_needed();
        return col[i];
    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    size_t index;
    Col_type<Cell_Type> col;

public:
    BArrayDenseCol_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        size_t j
    ) : array(&array_), index(j)
    {

        for (size_t i = 0u; i < array->N; ++i)
        {
            
            if (array->el[POS_N(i, index, array->N)] != ZERO_CELL)
                col[i] = col[POS_N(i, index, array->N)];
                
        }

    };

    typename Col_type<Cell_Type>::iterator begin()
    {
        return col.begin();
    };

    typename Col_type<Cell_Type>::iterator end()
    {
        return col.end();
    };


    size_t size() const noexcept
    {
        return col.size();
    };

    const std::pair<size_t,Cell_Type*> operator()(size_t i) const
    {
        return col[i];
    }

};

#undef POS
#undef POS_N
#undef ZERO_CELL

#endif