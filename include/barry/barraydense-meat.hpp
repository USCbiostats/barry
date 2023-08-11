// #include <stdexcept>
// #include "barraydense-bones.hpp"

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


#define BDENSE_TYPE() BArrayDense<Cell_Type, Data_Type>

#define BDENSE_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BDENSE_TEMPLATE(a,b) \
    template BDENSE_TEMPLATE_ARGS() inline a BDENSE_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template<typename Cell_Type, typename Data_Type>
Cell_Type BArrayDense<Cell_Type,Data_Type>::Cell_default = static_cast< Cell_Type >(1.0); 

#define ZERO_CELL static_cast<Cell_Type>(0.0)

// Edgelist with data
BDENSE_TEMPLATE(,BArrayDense)(
    size_t N_,
    size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
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

    el.resize(N * M, ZERO_CELL);
    el_rowsums.resize(N, ZERO_CELL);
    el_colsums.resize(M, ZERO_CELL);
    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        bool empty = is_empty(source[i], target[i], true);
        if (add && !empty)
        {

            Cell_Type tmp = el[POS(source[i], target[i])];
            
            el_rowsums[source[i]] += (value[i] - tmp);
            el_colsums[target[i]] += (value[i] - tmp);

            el[POS(source[i], target[i])] += value[i];
            
            continue;

        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        el[POS(source[i], target[i])] = value[i];

        el_rowsums[source[i]] += value[i];
        el_colsums[target[i]] += value[i];
        

    }
    
    return;
  
}

// Edgelist without data
BDENSE_TEMPLATE(, BArrayDense)(
    size_t N_, size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
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

    el.resize(N * M, ZERO_CELL);
    el_rowsums.resize(N, ZERO_CELL);
    el_colsums.resize(M, ZERO_CELL);
    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        bool empty = is_empty(source[i], target[i], true);
        if (add && !empty)
        {

            Cell_Type tmp = el[POS(source[i], target[i])];
            
            el_rowsums[source[i]] += (value[i] - tmp);
            el_colsums[target[i]] += (value[i] - tmp);

            el[POS(source[i], target[i])] += value[i];
            
            continue;

        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        el[POS(source[i], target[i])] = value[i];

        el_rowsums[source[i]] += value[i];
        el_colsums[target[i]] += value[i];
        

    }
  
}

BDENSE_TEMPLATE(, BArrayDense)(
    const BDENSE_TYPE() & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M){
  
    // Dimensions
    el.resize(0u);
    el_rowsums.resize(0u);
    el_colsums.resize(0u);
    
    std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
    std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
    std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));

    // this->NCells  = Array_.NCells;
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
      
        el.resize(0u);
        el_rowsums.resize(0u);
        el_colsums.resize(0u);
        
        // Entries
        std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
        std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
        std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));


        // this->NCells = Array_.NCells;
        this->N      = Array_.N;
        this->M      = Array_.M;
      
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
    // NCells(std::move(x.NCells)),
    el(std::move(x.el)),
    el_rowsums(std::move(x.el_rowsums)),
    el_colsums(std::move(x.el_colsums)),
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
        // NCells = x.NCells;
        
        std::swap(el, x.el);
        std::swap(el_rowsums, x.el_rowsums);
        std::swap(el_colsums, x.el_colsums);
              
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
    if ( (N != Array_.nrow()) | (M != Array_.ncol()) )
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

BDENSE_TEMPLATE(Data_Type *, D_ptr) () {
    return this->data;
}

BDENSE_TEMPLATE(const Data_Type *, D_ptr) () const {
    return this->data;
}

BDENSE_TEMPLATE(Data_Type &, D) () {
    return *this->data;
}

BDENSE_TEMPLATE(const Data_Type &, D) () const {
    return *this->data;
}

BDENSE_TEMPLATE(void, out_of_range) (
    size_t i,
    size_t j
) const {

    if (i >= N)
    {
        std::string err_msg = "The row is out of range: " + std::to_string(i) + " >= " + std::to_string(N);
        throw std::range_error(err_msg);

    } else if (j >= M)
    {
        std::string err_msg = "The column is out of range: " + std::to_string(j) + " >= " + std::to_string(M);
        throw std::range_error(err_msg);
    }

    return;

}
    
BDENSE_TEMPLATE(Cell_Type, get_cell) (
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    return el[POS(i, j)];
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_row_vec) (
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans(ncol(), static_cast< Cell_Type >(false));
    for (size_t j = 0u; j < M; ++j) 
        ans[j] = el[POS(i, j)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_row_vec) (
    std::vector<Cell_Type> * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (size_t j = 0u; j < M; ++j) 
        x->at(j) = el[POS(i, j)];
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_col_vec)(
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans(nrow(), static_cast< Cell_Type >(false));
    for (size_t j = 0u; j < N; ++j) 
        ans[j] = el[POS(j, i)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_col_vec) (
    std::vector<Cell_Type> * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (size_t j = 0u; j < N; ++j) 
        x->at(j) = el[POS(j, i)];//this->get_cell(iter->first, i, false);
    
}
template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseRow_const<Cell_Type,Data_Type> BDENSE_TYPE()::row(
    size_t i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow_const<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseRow<Cell_Type,Data_Type> & BDENSE_TYPE()::row(
    size_t i,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseCol_const<Cell_Type,Data_Type> 
BArrayDense<Cell_Type,Data_Type>::col(
    size_t j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol_const<Cell_Type,Data_Type>(*this, j);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseCol<Cell_Type,Data_Type> & 
BArrayDense<Cell_Type,Data_Type>::col(
    size_t j,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol<Cell_Type,Data_Type>(*this, j);

}

BDENSE_TEMPLATE(Entries< Cell_Type >, get_entries)() const {
    
    size_t nzero = this->nnozero();

    Entries<Cell_Type> res(nzero);
    
    for (size_t i = 0u; i < N; ++i)
    {
        for (size_t j = 0u; j < M; ++j)
        {

            if (el[POS(i, j)] != BARRY_ZERO_DENSE)
            {

                res.source.push_back(i),
                res.target.push_back(j),
                res.val.push_back(el[POS(i, j)]);

            }
            

        }

    }
    
    return res;

}

BDENSE_TEMPLATE(bool, is_empty)(
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i, j)] == ZERO_CELL;
    
}

BDENSE_TEMPLATE(size_t, nrow)() const noexcept {
    return N;
}

BDENSE_TEMPLATE(size_t, ncol)() const noexcept {
    return M;
}

BDENSE_TEMPLATE(size_t, nnozero)() const noexcept {

    size_t nzero = 0u;
    for (auto & v : el)
        if (v != BARRY_ZERO_DENSE)
            nzero++;

    return nzero;
}

BDENSE_TEMPLATE(Cell< Cell_Type>, default_val)() const {
    return this->Cell_default;
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator+=)(
    const std::pair<size_t,size_t> & coords
) {
    

    size_t i = coords.first;
    size_t j = coords.second;

    out_of_range(i, j);

    el[POS(i,j)]  += 1;
    el_rowsums[i] += 1;
    el_colsums[j] += 1;
    
    return *this;
    
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator-=)(
    const std::pair<size_t,size_t> & coords
) {
    
    size_t i = coords.first;
    size_t j = coords.second;

    out_of_range(i, j);

    Cell_Type old = el[POS(i,j)];

    el[POS(i,j)]   = ZERO_CELL;
    el_rowsums[i] -= old;
    el_colsums[j] -= old;
    
    return *this;
    
}

template BDENSE_TEMPLATE_ARGS()
inline BArrayDenseCell<Cell_Type,Data_Type> BDENSE_TYPE()::operator()(  
    size_t i,
    size_t j,
    bool check_bounds
) {
    
    return BArrayDenseCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template BDENSE_TEMPLATE_ARGS()
inline const Cell_Type BDENSE_TYPE()::operator()(  
    size_t i,
    size_t j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i,j)];
    
}

BDENSE_TEMPLATE(void, rm_cell) (
    size_t i,
    size_t j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);

    BARRY_UNUSED(check_exists)
        
    // Remove the pointer first (so it wont point to empty)
    el_rowsums[i] -= el[POS(i, j)];
    el_colsums[j] -= el[POS(i, j)];    
    el[POS(i, j)] = BARRY_ZERO_DENSE;
    
    return;

}

BDENSE_TEMPLATE(void, insert_cell) (
    size_t i,
    size_t j,
    const Cell< Cell_Type> & v,
    bool check_bounds,
    bool check_exists
) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    BARRY_UNUSED(check_exists)

    if (el[POS(i,j)] == BARRY_ZERO_DENSE)
    {

        el_rowsums[i] += v.value;
        el_colsums[j] += v.value;
        
    } 
    else
    {

        Cell_Type old = el[POS(i,j)];
        el_rowsums[i] += (v.value - old);
        el_colsums[j] += (v.value - old);

    }

    el[POS(i, j)] = v.value;

    return;

    
}

BDENSE_TEMPLATE(void, insert_cell)(
    size_t i,
    size_t j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
    
    if (check_bounds)
        out_of_range(i,j);

    BARRY_UNUSED(check_exists)
        
    if (el[POS(i,j)] == BARRY_ZERO_DENSE)
    {

        el_rowsums[i] += v;
        el_colsums[j] += v;
        
    } 
    else
    {

        Cell_Type old = el[POS(i,j)];
        el_rowsums[i] += (v - old);
        el_colsums[j] += (v - old);

    }

    el[POS(i, j)] = v;

}

BDENSE_TEMPLATE(void, swap_cells) (
        size_t i0, size_t j0,
        size_t i1, size_t j1,
        bool check_bounds,
        int check_exists,
        int * report
) {
    
    if (check_bounds) {
        out_of_range(i0,j0);
        out_of_range(i1,j1);
    }
    
        
    // Just in case, if this was passed
    if (report != nullptr)
        (*report) = EXISTS::BOTH;
    
    // If source and target coincide, we do nothing
    if ((i0 == i1) && (j0 == j1)) 
        return;

    // Updating rowand col sumns    
    Cell_Type val0 = el[POS(i0,j0)];
    Cell_Type val1 = el[POS(i1,j1)];

    rm_cell(i0, j0, false, false);
    rm_cell(i1, j1, false, false);
    
    // Inserting the cells by reference, these will be deleted afterwards
    insert_cell(i0, j0, val1, false, false);
    insert_cell(i1, j1, val0, false, false);
    
    return;

}

BDENSE_TEMPLATE(void, toggle_cell) (
    size_t i,
    size_t j,
    bool check_bounds,
    int check_exists
) {

    if (check_bounds)
        out_of_range(i, j);

    if (el[POS(i,j)] == ZERO_CELL)
        insert_cell(i,j,1,false,false);
    else
        rm_cell(i,j,false,false);
    
    return;
    
}

BDENSE_TEMPLATE(void, swap_rows) (
    size_t i0,
    size_t i1,
    bool check_bounds
) {
  
    if (check_bounds)
    {

        out_of_range(i0,0u);
        out_of_range(i1,0u);

    }
     
    // if (NCells == 0u)
    //     return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    for (size_t j = 0u; j < M; ++j)
        std::swap(el[POS(i0, j)], el[POS(i1, j)]);

    std::swap(el_rowsums[i0], el_rowsums[i1]);
    
    return;
}

// This swapping is more expensive overall
BDENSE_TEMPLATE(void, swap_cols) (
    size_t j0,
    size_t j1,
    bool check_bounds
) {

    if (check_bounds)
    {

        out_of_range(0u, j0);
        out_of_range(0u, j1);

    }
    
    if ((el_colsums[j0] == ZERO_CELL) && el_colsums[j1] == ZERO_CELL)
        return;

    // Swapping happens naturally, need to take care of the pointers
    // though
    for (size_t i = 0u; i < N; ++i)
        std::swap(el[POS(i, j0)], el[POS(i, j1)]);

    std::swap(el_colsums[j0], el_colsums[j1]);
    
    return;
}

BDENSE_TEMPLATE(void, zero_row) (
    size_t i,
    bool check_bounds
    ) {
  
    if (check_bounds)
        out_of_range(i, 0u);

    if (el_rowsums[i] == ZERO_CELL)
        return;

    // Else, remove all elements
    for (size_t col = 0u; col < M; col++) 
        rm_cell(i, col, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, zero_col) (
    size_t j,
    bool check_bounds
  ) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    if (el_colsums[j] == ZERO_CELL)
        return;
    
    // Else, remove all elements
    for (size_t row = 0u; row < N; row++) 
        rm_cell(row, j, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, transpose) () {
  
    // if (NCells == 0u)
    // {

    //     std::swap(N, M);
    //     return;

    // }

    // Start by flipping the switch 
    visited = !visited;

    // size_t N0 = N, M0 = M;
    std::vector< Cell< Cell_Type > > tmp_el(std::move(el));
    el.resize(N * M, ZERO_CELL);
    for (size_t i = 0u; i < N; ++i) 
        for (size_t j = 0u; j < M; ++j)
            std::swap(tmp_el[POS(i, j)], el[POS_N(j, i, M)]);
    
    // Swapping the values
    std::swap(N, M);
    std::swap(el_rowsums, el_colsums);
    
    return;

}

BDENSE_TEMPLATE(void, clear) (
    bool hard
) {
    
    BARRY_UNUSED(hard)
    
    std::fill(el.begin(), el.end(), ZERO_CELL);
    std::fill(el_rowsums.begin(), el_rowsums.end(), ZERO_CELL);
    std::fill(el_colsums.begin(), el_colsums.end(), ZERO_CELL);
    
    return;
    
}

BDENSE_TEMPLATE(void, resize) (
    size_t N_,
    size_t M_
) {

    // Moving stuff around
    std::vector< Cell_Type > el_tmp(el);
    el.resize(N_ * M_, ZERO_CELL);
    el_rowsums.resize(N_, ZERO_CELL);
    el_colsums.resize(M_, ZERO_CELL);

    for (size_t i = 0u; i < N; ++i)
    {
        // If reached the end
        if (i >= N_)
            break;

        for (size_t j = 0u; j < M; ++j)
        {

            if (j >= M_)
                break;

            insert_cell(i, j, el_tmp[POS_N(i, j, N_)], false, false);

        }

    }

    N = N_;
    M = M_;
    
    return;

}

BDENSE_TEMPLATE(void, reserve) () {

    el.reserve(N * M);
    el_rowsums.reserve(N);
    el_colsums.reserve(M);
    return;
  
}

BDENSE_TEMPLATE(void, print) (
    const char * fmt,
    ...
) const
{
  
    std::va_list args;
    va_start(args, fmt);
    printf_barry(fmt, args);
    va_end(args);

    for (size_t i = 0u; i < N; ++i)
    {

        printf_barry("[%3li,] ", i);

        for (size_t j = 0u; j < M; ++j)
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

BDENSE_TEMPLATE(const std::vector< Cell_Type > &, get_data)() const
{
    return el;
}

BDENSE_TEMPLATE(const Cell_Type, rowsum)(size_t i) const
{
    return el_rowsums[i];
}

BDENSE_TEMPLATE(const Cell_Type, colsum)(size_t j) const
{
    return el_colsums[j];
}

#undef ROW
#undef COL
#undef POS
#undef POS_N

#undef BDENSE_TYPE
#undef BDENSE_TEMPLATE_ARGS
#undef BDENSE_TEMPLATE
#undef ZERO_CELL

#endif

