#ifndef BARRY_BARRAYDENSEROW_BONES_HPP
#define BARRY_BARRAYDENSEROW_BONES_HPP

#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    std::vector< std::pair<unsigned int, Cell<Cell_Type> > > dat;
    unsigned int index;

public:
    BArrayDenseRow(
        BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int i
    ) : array(&array_), index(i)
    {

        const unsigned int N = array->N;

        dat.resize(array->M);
        for (unsigned int j = 0u; j < dat.size(); ++j)
            dat[i] = std::pair<unsigned int,Cell<Cell_Type>>(
                j,
                Cell<Cell_Type>(array->el[POS_N(i, j, N)]))
                ;

        return;

    };

    typename std::vector< std::pair<unsigned int,Cell<Cell_Type>>  >::iterator & begin()
    {
        return dat.begin();
    };

    typename std::vector< std::pair<unsigned int,Cell<Cell_Type>>  >::iterator & end()
    {
        return dat.end();
    };

    size_t size() const noexcept
    {
        return dat.size();
    };

    std::pair<unsigned int,Cell<Cell_Type>> & operator()(unsigned int i)
    {
        return dat[i];
    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    std::vector< std::pair<unsigned int, Cell<Cell_Type> > > dat;
    unsigned int index;

public:
    BArrayDenseRow_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int i
    ) : array(&array_), index(i)
    {

        const unsigned int N = array->N;

        dat.resize(array->M);
        for (unsigned int j = 0u; j < dat.size(); ++j)
            dat[i] = std::pair<unsigned int, Cell_Type>(
                j,
                Cell<Cell_Type>(array->el[POS_N(i, j, N)])
                );

        return;


    };

    typename std::vector< std::pair<unsigned int,Cell<Cell_Type>>  >::const_iterator begin()
    {
        return dat.begin();
    };

    typename std::vector< std::pair<unsigned int,Cell<Cell_Type>>  >::const_iterator end()
    {
        return dat.end();
    };

    size_t size() const noexcept
    {
        return dat.size();
    };

    const std::pair<unsigned int,Cell<Cell_Type>> operator()(unsigned int i) const
    {
        return dat[i];
    }

};

#undef POS
#undef POS_N

#endif