#ifndef BARRY_BARRAYDENSECOL_BONES 
#define BARRY_BARRAYDENSECOL_BONES

#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    std::vector< std::pair<unsigned int, Cell_Type* > > dat;
    unsigned int index;

public:
    BArrayDenseCol(
        BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int j
    ) : array(&array_), index(j)
    {

        const unsigned int N = array->N;

        dat.resize(N);
        for (unsigned int i = 0u; i < dat.size(); ++i)
            dat[i] = std::pair<unsigned int, Cell_Type*>(j, &(array->el[POS_N(i, j, N)]));

        return;

    };


    typename std::vector< std::pair<unsigned int,Cell_Type*>  >::iterator & begin()
    {
        return dat.begin();
    };

    typename std::vector< std::pair<unsigned int,Cell_Type*>  >::iterator & end()
    {
        return dat.end();
    };

    size_t size() const noexcept
    {
        return dat.size();
    };

    std::pair<unsigned int,Cell_Type*> & operator()(unsigned int i)
    {
        return dat[i];
    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    std::vector< std::pair<unsigned int, Cell_Type* > > dat;
    std::vector< Cell_Type > cell_dat;
    unsigned int index;

public:
    BArrayDenseCol_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int j
    ) : array(&array_), index(j)
    {

        const unsigned int N = array->N;

        dat.resize(N);
        cell_dat.resize(N);
        for (unsigned int i = 0u; i < dat.size(); ++i)
        {
            cell_dat[i] = array->el[POS_N(i, j, N)];
            dat[i]      = std::pair<unsigned int, Cell_Type*>(j, &cell_dat[i]);
        }

        return;


    };

    typename std::vector< std::pair<unsigned int,Cell_Type*>  >::iterator begin()
    {
        return dat.begin();
    };

    typename std::vector< std::pair<unsigned int,Cell_Type*>  >::iterator end()
    {
        return dat.end();
    };


    size_t size() const noexcept
    {
        return dat.size();
    };

    const std::pair<unsigned int,Cell_Type*> operator()(unsigned int i) const
    {
        return dat[i];
    }

};

#undef POS
#undef POS_N

#endif