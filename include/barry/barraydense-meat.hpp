// #include <stdexcept>
#include "barraydense-bones.hpp"

#ifndef BARRY_BARRAYDENSE_MEAT_HPP
#define BARRY_BARRAYDENSE_MEAT_HPP 

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell_const;

#define BDENSE_TYPE() BArrayDense<Cell_Type, Data_Type>

#define BDENSE_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BDENSE_TEMPLATE(a,b) \
    template BDENSE_TEMPLATE_ARGS() inline a BDENSE_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template<typename Cell_Type, typename Data_Type>
Cell<Cell_Type> BArrayDense<Cell_Type,Data_Type>::Cell_default = Cell<Cell_Type>(); 

#define ZERO_CELL Cell< Cell_Type >( (Cell_Type) 0.0, false, false)

// Edgelist with data
BDENSE_TEMPLATE(,BArrayDense)(
    uint N_,
    uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< Cell_Type > & value,
    bool add
) {
  
    if (source.size() != target.size())
        throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
        throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;

    el.resize(N, M);
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        bool empty = is_empty(source[i], target[i], true);
        if (add && !empty)
        {

            el[POS(source[i], target[i])].add(value[i]);
            continue;

        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        this->insert_cell(source[i], target[i], value[i], false, false);

    }
    
    return;
  
}

// Edgelist with data
BDENSE_TEMPLATE(, BArrayDense)(
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    bool add
) {
  
    std::vector< Cell_Type > value(source.size(), static_cast<Cell_Type>(1.0));

    if (source.size() != target.size())
      throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
      throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;
    
    el.resize(N * M);
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        if ((source[i] >= N_) | (target[i] >= M_))
            throw std::range_error("Either source or target point to an element outside of the range by (N,M).");
        
        // Checking if it exists
        if (el[POS(source[i], target[i])].active)
        {

            if (!add)
                throw std::logic_error("The value already exists. Use 'add = true'.");
          
            // Increasing the value (this will automatically update the
            // other value)
            el[POS(source[i], target[i])].add(value[i]);
            continue;

        }
        
        // Adding the value and creating a pointer to it
        el[POS(source[i], target[i])].value   = value[i];
        el[POS(source[i], target[i])].visited = visited;
        el[POS(source[i], target[i])].active  = true;
        
        NCells++;

    }
    
    return;
  
}

BDENSE_TEMPLATE(, BArrayDense)(
    const BDENSE_TYPE() & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M){
  
    // Dimensions
    el.resize(N * M);
    
    std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));

    this->NCells  = Array_.NCells;
    this->visited = Array_.visited;
    
    // Data
    if (Array_.data != nullptr)
    {

        if (copy_data)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        } else {

            data = Array_.data;
            delete_data = false;

        }

    }
    
    return;
  
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator=) (
    const BDENSE_TYPE() & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        clear(true);
        resize(Array_.N, Array_.M);
        
        // Entries
        for (uint i = 0u; i < N; ++i)
        {
          
            if (Array_.nnozero() == nnozero())
                break;
            
            for (auto& r : Array_.row(i, false)) 
                this->insert_cell(i, r.first, r.second.value, false, false);
          
        }
      
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;

        }

        if (Array_.data != nullptr)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        }
      
    }
      
    return *this;
  
}

BDENSE_TEMPLATE(, BArrayDense)(
    BDENSE_TYPE() && x
    ) noexcept :
    N(std::move(x.N)), M(std::move(x.M)),
    NCells(std::move(x.NCells)),
    el(std::move(x.el)),
    data(std::move(x.data)),
    delete_data(std::move(x.delete_data))
{

      x.data        = nullptr;
      x.delete_data = false;

}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator=)(
    BDENSE_TYPE() && x
) noexcept {
  
    // Clearing
    if (this != &x)
    {
      
        N      = x.N;
        M      = x.M;
        NCells = x.NCells;
        
        std::swap(el, x.el);
              
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;

        }

        if (x.data != nullptr)
        {

            data        = std::move(x.data);
            delete_data = x.delete_data;

            x.delete_data = false;
            x.data = nullptr;

        }
      
    }
      
    return *this;
  
}

BDENSE_TEMPLATE(bool, operator==) (
    const BDENSE_TYPE() & Array_
) {
    
    // Dimension and number of cells used
    if ((N != Array_.nrow()) | (M != Array_.ncol()) | (NCells != Array_.nnozero()))
        return false;
    
    // One holds, and the other doesn't.
    if ((!data & Array_.data) | (data & !Array_.data))
        return false;
    
    if (this->el != Array_.el)
        return false;
    
    return true;
}

BDENSE_TEMPLATE(, ~BArrayDense) () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

BDENSE_TEMPLATE(void, set_data) (
    Data_Type * data_,
    bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

BDENSE_TEMPLATE(Data_Type *, D) () {
    return this->data;
}

BDENSE_TEMPLATE(const Data_Type *, D) () const {
    return this->data;
}

BDENSE_TEMPLATE(void, out_of_range) (
    uint i,
    uint j
) const {

    if (i >= N)
        throw std::range_error("The row is out of range.");
    else if (j >= M)
        throw std::range_error("The column is out of range.");

    return;

}
    
BDENSE_TEMPLATE(Cell_Type, get_cell) (
    uint i,
    uint j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    return el[POS(i, j)].value;
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_row_vec) (
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans(ncol(), static_cast< Cell_Type >(false));
    for (uint j = 0u; j < M; ++j) 
        ans[j] = el[POS(i, j)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_row_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (uint j = 0u; j < M; ++j) 
        x->at(j) = el[POS(i, j)];
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_col_vec)(
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans(nrow(), static_cast< Cell_Type >(false));
    for (uint j = 0u; j < N; ++j) 
        ans[j] = el[POS(j, i)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_col_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (uint j = 0u; j < N; ++j) 
        x->at(j) = el[POS(j, i)];//this->get_cell(iter->first, i, false);
    
}
template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseRow_const<Cell_Type,Data_Type> BDENSE_TYPE()::row(
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow_const<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseRow<Cell_Type,Data_Type> & BDENSE_TYPE()::row(
    uint i,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseCol_const<Cell_Type,Data_Type> 
BArrayDense<Cell_Type,Data_Type>::col(
    uint j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol_const<Cell_Type,Data_Type>(*this, j);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseCol<Cell_Type,Data_Type> & 
BArrayDense<Cell_Type,Data_Type>::col(
    uint j,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol<Cell_Type,Data_Type>(*this, j);

}

BDENSE_TEMPLATE(Entries< Cell_Type >, get_entries)() const {
    
    Entries<Cell_Type> res(NCells);
    
    for (uint i = 0u; i < N; ++i)
    {
        for (uint j = 0u; col < M; ++j)
        {

            if (el[POS(i, j)].active)
            res.source.push_back(i),
            res.target.push_back(j),
            res.val.push_back(el[POS(i, j)].value);

        }

    }
    
    return res;

}

BDENSE_TEMPLATE(bool, is_empty)(
    uint i,
    uint j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i, j)].active;
    
}

BDENSE_TEMPLATE(unsigned int, nrow)() const noexcept {
    return N;
}

BDENSE_TEMPLATE(unsigned int, ncol)() const noexcept {
    return M;
}

BDENSE_TEMPLATE(unsigned int, nnozero)() const noexcept {
    return NCells;
}

BDENSE_TEMPLATE(Cell< Cell_Type>, default_val)() const {
    return this->Cell_default;
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator+=)(
    const std::pair<uint,uint> & coords
) {
    
    this->insert_cell(
        coords.first,
        coords.second,
        this->Cell_default,
        true,
        true
    );
    
    return *this;
    
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator-=)(
    const std::pair<uint,uint> & coords
) {
    
    this->rm_cell(
        coords.first,
        coords.second,
        true,
        true
    );
    
    return *this;
    
}

template BDENSE_TEMPLATE_ARGS()
inline BArrayDenseCell<Cell_Type,Data_Type> BDENSE_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) {
    
    return BArrayDenseCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template BDENSE_TEMPLATE_ARGS()
inline const BArrayDenseCell_const<Cell_Type,Data_Type> BDENSE_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) const {
    
    return BArrayDenseCell_const<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

BDENSE_TEMPLATE(void, rm_cell) (
    uint i,
    uint j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);
        
    // Remove the pointer first (so it wont point to empty)
    el[POS(i, j)].active = false;
    
    NCells--;
    
    return;

}

BDENSE_TEMPLATE(void, insert_cell) (
    uint i,
    uint j,
    const Cell< Cell_Type> & v,
    bool check_bounds,
    bool check_exists
) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists)
    {
        
        // Checking if nothing here, then we move along
        if (NCells == 0u)
        {
            Cell< Cell_Type > tmp(v);
            el[POS(i, j)] = tmp;
            NCells++;
            return;
            
        } else
            throw std::logic_error("The cell already exists.");        
        
    } else {
        
        Cell< Cell_Type > tmp(v);
        el[POS(i, j)] = tmp;
        NCells++;
        
    }
    
    return;
    
}

BDENSE_TEMPLATE(void, insert_cell) (
        uint i,
        uint j,
        Cell< Cell_Type> && v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists)
    {
        
        // Checking if nothing here, then we move along
        if (NCells == 0u)
        {
            
            el[POS(i, j)] = std::move(v);
            NCells++;
            return;
            
        } else 
            throw std::logic_error("The cell already exists.");
        
        
    } else {
        
        el[POS(i, j)] = std::move(v);
        NCells++;
        return;
        
    }
    
    return;
    
}

BDENSE_TEMPLATE(void, insert_cell)(
    uint i,
    uint j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
        
    return insert_cell(i, j, Cell<Cell_Type>(v, visited), check_bounds, check_exists);

}

BDENSE_TEMPLATE(void, swap_cells) (
        uint i0, uint j0,
        uint i1, uint j1,
        bool check_bounds,
        int check_exists,
        int * report
) {
    
    if (check_bounds) {
        out_of_range(i0,j0);
        out_of_range(i1,j1);
    }
    
    // Simplest case, we know both exists, so we don't need to check anything
    if (check_exists == CHECK::NONE)
    {
        
        // Just in case, if this was passed
        if (report != nullptr)
            (*report) = EXISTS::BOTH;
        
        // If source and target coincide, we do nothing
        if ((i0 == i1) && (j0 == j1)) 
            return;
        
        // Using the initializing by move, after this, the cell becomes
        // invalid. We use pointers instead as this way we access the Heap memory,
        // which should be faster to access.
        Cell<Cell_Type> c0(std::move(el[POS(i0, j0)]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(el[POS(i1, j1)]));
        rm_cell(i1, j1, false, false);
        
        // Inserting the cells by reference, these will be deleted afterwards
        insert_cell(i0, j0, c1, false, false);
        insert_cell(i1, j1, c0, false, false);
        
        return;
        
    }
    
    bool check0, check1;
    if (check_exists == CHECK::BOTH)
    {
        
        check0 = !is_empty(i0, j0, false);
        check1 = !is_empty(i1, j1, false);
        
    } else if (check_exists == CHECK::ONE)
    {
        
        check0 = !is_empty(i0, j0, false);
        check1 = true;
        
    } else if (check_exists == CHECK::TWO)
    {
        
        check0 = true;
        check1 = !is_empty(i1, j1, false);
        
    }
    
    if (report != nullptr) 
        (*report) = EXISTS::NONE;
    
    // If both cells exists
    if (check0 & check1)
    {
        
        if (report != nullptr) 
            (*report) = EXISTS::BOTH;
        
        // If source and target coincide, we do nothing
        if ((i0 == i1) && (j0 == j1)) 
            return;
        
        Cell<Cell_Type> c0(std::move(el[POS(i0, j0)]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(el[POS(i1, j1)]));
        rm_cell(i1, j1, false, false);
        
        insert_cell(i0, j0, c1, false, false);
        insert_cell(i1, j1, c0, false, false);
        
    } else if (!check0 & check1) { // If only the second exists
        
        if (report != nullptr) 
            (*report) = EXISTS::TWO;
        
        insert_cell(i0, j0, el[POS(i1, j1)], false, false);
        rm_cell(i1, j1, false, false);
        
    } else if (check0 & !check1) {
        
        if (report != nullptr) 
            (*report) = EXISTS::ONE;
        
        insert_cell(i1, j1, el[POS(i0, j0)], false, false);
        rm_cell(i0, j0, false, false);
        
    }
    
    return;

}

BDENSE_TEMPLATE(void, toggle_cell) (
    uint i,
    uint j,
    bool check_bounds,
    int check_exists
) {
    
    if (check_bounds)
        out_of_range(i, j);
    
    if (check_exists == EXISTS::UKNOWN)
    {
        
        if (is_empty(i, j, false))
        {

            insert_cell(i, j, BArrayDense<Cell_Type, Data_Type>::Cell_default, false, false);
            el[POS(i, j)].visited = visited;

        } else
            rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ONE) {
        
        rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ZERO) {
        
        insert_cell(i, j, BArrayDense<Cell_Type,Data_Type>::Cell_default, false, false);
        el[POS(i, j)].visited = visited;
        
    }
    
    return;
    
}

BDENSE_TEMPLATE(void, swap_rows) (
    uint i0,
    uint i1,
    bool check_bounds
) {
  
    if (check_bounds)
    {

        out_of_range(i0,0u);
        out_of_range(i1,0u);

    }
     
    if (NCells == 0u)
        return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    for (uint j = 0u; j < M; ++j)
        std::swap(el[POS(i0, j)], el[POS(i1, j)]);
    
    return;
}

// This swapping is more expensive overall
BDENSE_TEMPLATE(void, swap_cols) (
    uint j0,
    uint j1,
    bool check_bounds
) {

    if (check_bounds)
    {

        out_of_range(0u, j0);
        out_of_range(0u, j1);

    }
    
    if (NCells == 0u)
        return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    for (uint i = 0u; i < N; ++i)
        std::swap(el[POS(i, j0)], el[POS(i, j1)]);
    
    return;
}

BDENSE_TEMPLATE(void, zero_row) (
    uint i,
    bool check_bounds
    ) {
  
    if (check_bounds)
        out_of_range(i, 0u);
    
    // If already empty, nothing to do
    if (NCells == 0u)
        return;

    // Else, remove all elements
    for (unsigned int col = 0u; col < M; col++) 
        rm_cell(i, col, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, zero_col) (
    uint j,
    bool check_bounds
  ) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    // Nothing to do
    if (NCells == 0u)
        return;
    
    // Else, remove all elements
    for (unsigned int row = 0u; row < N; row++) 
        rm_cell(row, j, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, transpose) () {
  
    if (NCells == 0u)
    {

        std::swap(N, M);
        return;

    }

    // Start by flipping the switch 
    visited = !visited;

    // uint N0 = N, M0 = M;
    std::vector< Cell< Cell_Type > > tmp_el(std::move(el));
    el.resize(N * M, ZERO_CELL);
    for (uint i = 0u; i < N; ++i) 
        for (uint j = 0u; j < M; ++j)
            std::swap(tmp_el[POS(i, j)], el[POS_N(j, i, M)]);
    
    // Swapping the values
    std::swap(N, M);
    
    return;

}

BDENSE_TEMPLATE(void, clear) (
    bool hard
) {
    
    if (hard)
    {
      
        el.clear();
        el.resize(N * M);
        NCells = 0u;
      
    } else {
        
        for (auto & i : el)
            i.active = false;
        
    }
    
    return;
    
}

BDENSE_TEMPLATE(void, resize) (
    uint N_,
    uint M_
) {

    // If already empty, it is very simple
    if (NCells == 0u)
    {
        
        el.resize(N_ * M_);
        N = N_;
        M = M_;
        return;

    }
  
    // Moving stuff around
    std::vector< Cell< Cell_Type > > el_tmp(std::move(el));
    el.resize(N_ * M_, ZERO_CELL);
    for (unsigned int i = 0u; i < N; ++i)
    {
        // If reached the end
        if (i >= N_)
            break;

        for (unsigned int j = 0u; j < M; ++j)
        {

            if (j >= M_)
                break;

            std::swap(el_tmp[POS_N(i, j, N_)], el[POS(i, j)]);

        }

    }

    N = N_;
    M = M_;
    
    return;

}

BDENSE_TEMPLATE(void, reserve) () {

#ifdef BARRAY_USE_UNORDERED_MAP
    tmp_row.reserve(M);
    tmp_col.reserve(N);
#endif
    return;
  
}

BDENSE_TEMPLATE(void, print) () const {
  
    for (uint i = 0u; i < N; ++i)
    {

        printf_barry("[%3i,] ", i);

        for (uint j = 0u; j < M; ++j)
        {

            if (this->is_empty(i, j, false))
                printf_barry("    . ");
            else 
                printf_barry(" %.2f ", static_cast<double>(this->get_cell(i, j, false)));
            
        }

        printf_barry("\n");

    }
    
    return;
    
}

#undef ROW
#undef COL
#undef POS
#undef POS_N

#undef BDENSE_TYPE
#undef BDENSE_TEMPLATE_ARGS
#undef BDENSE_TEMPLATE

#endif

