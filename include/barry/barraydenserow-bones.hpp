#ifndef BARRY_BARRAYDENSEROW_BONES_HPP
#define BARRY_BARRAYDENSEROW_BONES_HPP

#define POS(a,b) (b) * N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)
#define ZERO_CELL static_cast< Cell_Type >(0.0)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    Row_type< Cell_Type > row;
    unsigned int index;
    bool row_filled = false; // after row is filled

    void fill_if_needed()
    {
        if (!row_filled)
        {

            for (unsigned int j = 0u; j < array->M; ++j)
            {
                
                if (array->el[POS_N(index, j, array->N)] != ZERO_CELL)
                    row[j] = row[POS_N(index, j, array->N)];
                    
            }

            row_filled = true;
            
        }
    }


public:

    BArrayDenseRow(
        BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int i
    ) : array(&array_), index(i) {};

    typename Row_type<Cell_Type>::iterator & begin()
    {

        fill_if_needed();
        return row.begin();

    };

    typename Row_type<Cell_Type>::iterator & end()
    {

        fill_if_needed();
        return row.end();

    };

    size_t size() const noexcept
    {

        fill_if_needed();
        return row.size();

    };

    std::pair<unsigned int,Cell<Cell_Type>> & operator()(unsigned int i)
    {

        fill_if_needed();
        return row[i];

    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    Row_type< Cell_Type > row;
    unsigned int index;

public:
    BArrayDenseRow_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int i
    ) : array(&array_), index(i)
    {

        for (unsigned int j = 0u; j < array->M; ++j)
        {
            
            if (array->el[POS_N(index, j, array->N)] != ZERO_CELL)
                row[j] = row[POS_N(index, j, array->N)];
                
        }

        return;


    };

    typename Row_type< Cell_Type >::const_iterator begin() const
    {
        return row.begin();
    };

    typename Row_type< Cell_Type >::const_iterator end() const
    {
        return row.end();
    };

    size_t size() const noexcept
    {
        return row.size();
    };

    const std::pair<unsigned int,Cell<Cell_Type>> operator()(unsigned int i) const
    {
        return row[i];
    }

};

#undef POS
#undef POS_N
#undef ZERO_CELL

#endif