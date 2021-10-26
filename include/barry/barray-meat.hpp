// #include <stdexcept>
#include "barray-bones.hpp"

#ifndef BARRY_BARRAY_MEAT_HPP
#define BARRY_BARRAY_MEAT_HPP 

#define BARRAY_TYPE() BArray<Cell_Type, Data_Type>

#define BARRAY_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BARRAY_TEMPLATE(a,b) \
    template BARRAY_TEMPLATE_ARGS() inline a BARRAY_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]


template<typename Cell_Type, typename Data_Type>
Cell<Cell_Type> BArray<Cell_Type,Data_Type>::Cell_default = Cell<Cell_Type>(); 


// Edgelist with data
BARRAY_TEMPLATE(,BArray) (
    uint N_, uint M_,
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

    el_ij.resize(N);
    el_ji.resize(M);
    
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i) {
      
        // Checking range
        bool empty = this->is_empty(source[i], target[i], true);
        if (add && !empty) {
            ROW(source[i])[target[i]].add(value[i]);
            continue;
        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        this->insert_cell(source[i], target[i], value[i], false, false);
    }
    
    return;
  
}

// Edgelist with data
BARRAY_TEMPLATE(,BArray) (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    bool add
) {
  
    std::vector< Cell_Type > value(source.size(), (Cell_Type) 1.0);

    if (source.size() != target.size())
      throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
      throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;
    
    el_ij.resize(N);
    el_ji.resize(M);
    
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i) {
      
        // Checking range
        if ((source[i] >= N_) | (target[i] >= M_))
            throw std::range_error("Either source or target point to an element outside of the range by (N,M).");
        
        // Checking if it exists
        auto search = ROW(source[i]).find(target[i]);
        if (search != ROW(source[i]).end()) {
            if (!add)
                throw std::logic_error("The value already exists. Use 'add = true'.");
          
            // Increasing the value (this will automatically update the
            // other value)
            ROW(source[i])[target[i]].add(value[i]);
            continue;
        }
        
        // Adding the value and creating a pointer to it
        ROW(source[i]).emplace(
            std::pair<uint, Cell< Cell_Type> >(
                target[i],
                Cell< Cell_Type >(value[i], visited)
            )
        );
        
        COL(target[i]).emplace(
            source[i],
            &ROW(source[i])[target[i]]
        );

        NCells++;

    }
    
    return;
  
}

BARRAY_TEMPLATE(,BArray) (
    const BArray<Cell_Type,Data_Type> & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M)
{
  
    // Dimensions
    // el_ij.resize(N);
    // el_ji.resize(M);
    
    std::copy(Array_.el_ij.begin(), Array_.el_ij.end(), std::back_inserter(el_ij));
    std::copy(Array_.el_ji.begin(), Array_.el_ji.end(), std::back_inserter(el_ji));

    // Taking care of the pointers
    for (uint i = 0u; i < N; ++i)
    {

        for (auto& r: row(i, false))
            COL(r.first)[i] = &ROW(i)[r.first];

    }

    this->NCells  = Array_.NCells;
    this->visited = Array_.visited;
    
    // Data
    if (Array_.data != nullptr)
    {

        if (copy_data)
        {

            data = new Data_Type(* Array_.data );
            delete_data = true;

        } else {

            data = Array_.data;
            delete_data = false;

        }

    }
    
    return;
  
}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator=) (
    const BArray<Cell_Type,Data_Type> & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        this->clear(true);
        this->resize(Array_.N, Array_.M);
        
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
            delete_data = false;

        }

        if (Array_.data != nullptr)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        }
      
    }
      
    return *this;
  
}

BARRAY_TEMPLATE(,BArray) (
    BARRAY_TYPE() && x
  ) noexcept :
  N(0u), M(0u), NCells(0u),
  data(nullptr),
  delete_data(x.delete_data)
  {

    this->clear(true);
    this->resize(x.N, x.M);
    
    // Entries
    for (uint i = 0u; i < N; ++i) {
      
        if (x.nnozero() == nnozero())
            break;
        
        for (auto& r : x.row(i, false)) 
            this->insert_cell(i, r.first, r.second.value, false, false);
      
    }

    // Managing data
    if (x.data != nullptr)
    {

        if (x.delete_data)
        {

            data = new Data_Type(*x.data);
            delete_data = true;

        } else {
            data = x.data;
            delete_data = false;
        }


    }

}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator=) (
    BARRAY_TYPE() && x
) noexcept {
  
    // Clearing
    if (this != &x) {
      
        this->clear(true);
        this->resize(x.N, x.M);
        
        // Entries
        for (uint i = 0u; i < N; ++i) {
          
            if (x.nnozero() == nnozero())
                break;
            
            for (auto& r : x.row(i, false)) 
                this->insert_cell(i, r.first, r.second.value, false, false);
          
        }
      
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;
            delete_data = false;

        }

        if (x.data != nullptr)
        {

            data = new Data_Type( *x.data );
            delete_data = true;

        }

        // x.data = nullptr;
        // x.delete_data = false;
      
    }
      
    return *this;
  
}

BARRAY_TEMPLATE(bool, operator==) (
    const BARRAY_TYPE() & Array_
) {
    
    // Dimension and number of cells used
    if ((N != Array_.nrow()) | (M != Array_.ncol()) | (NCells != Array_.nnozero()))
        return false;
    
    // One holds, and the other doesn't.
    if ((!data & Array_.data) | (data & !Array_.data))
        return false;
    
    if (this->el_ij != Array_.el_ij)
        return false;
    
    return true;
}

BARRAY_TEMPLATE(,~BArray) () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

BARRAY_TEMPLATE(void, set_data) (
    Data_Type * data_, bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

BARRAY_TEMPLATE(Data_Type *, D) ()
{
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type * BArray<Cell_Type,Data_Type>::D() const
{
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
inline void BArray<Cell_Type,Data_Type>::flush_data()
{

    if (delete_data)
    {
        delete data;
        delete_data = false;
    }

    data = nullptr;

    return;

}

BARRAY_TEMPLATE(void, out_of_range) (
    uint i,
    uint j
) const {

    if (i >= N)
        throw std::range_error("The row is out of range.");
    else if (j >= M)
        throw std::range_error("The column is out of range.");
    return;

}
    
BARRAY_TEMPLATE(Cell_Type, get_cell) (
    uint i,
    uint j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    if (ROW(i).size() == 0u)
        return (Cell_Type) 0.0;
    
    // If it is not empty, then find and return
    auto search = ROW(i).find(j);
    if (search != ROW(i).end())
        return search->second.value;
    
    // This is if it is empty
    return (Cell_Type) 0.0;
    
}

BARRAY_TEMPLATE(std::vector< Cell_Type >, get_row_vec) (
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans(ncol(), (Cell_Type) false);
    for (const auto & iter : row(i, false)) 
        ans[iter.first] = iter.second.value; //this->get_cell(i, iter->first, false);
    

    return ans;
}

BARRAY_TEMPLATE(void, get_row_vec) (
    std::vector< Cell_Type > * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (const auto & iter : row(i, false)) 
        x->at(iter.first) = iter.second.value; // this->get_cell(i, iter->first, false);
    
}

BARRAY_TEMPLATE(std::vector< Cell_Type >, get_col_vec) (
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans(nrow(), (Cell_Type) false);
    for (const auto iter : col(i, false)) 
        ans[iter.first] = iter.second->value;//this->get_cell(iter->first, i, false);
    
    return ans;

}

BARRAY_TEMPLATE(void, get_col_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (const auto & iter : col(i, false)) 
        x->at(iter.first) = iter.second->value;//this->get_cell(iter->first, i, false);
    
}

BARRAY_TEMPLATE(const Row_type< Cell_Type > &, row) (
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return this->el_ij[i];

}

BARRAY_TEMPLATE(const Col_type< Cell_Type > &, col) (
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, i);

    return this->el_ji[i];
    
}

BARRAY_TEMPLATE(Entries< Cell_Type >, get_entries) () const {
    
    Entries<Cell_Type> res(NCells);
    
    for (uint i = 0u; i < N; ++i) {
        
        if (ROW(i).size() == 0u)
            continue;
        
        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            res.source.push_back(i),
            res.target.push_back(col->first),
            res.val.push_back(col->second.value);
        }
    }
    
    return res;
}

BARRAY_TEMPLATE(bool, is_empty) (
    uint i,
    uint j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    if (ROW(i).size() == 0u)
        return true;
    else if (COL(j).size() == 0u)
        return true;
    
    if (ROW(i).find(j) == ROW(i).end())
        return true;
    
    return false;
    
}


BARRAY_TEMPLATE(unsigned int, nrow) () const noexcept {
    return N;
}


BARRAY_TEMPLATE(unsigned int, ncol) () const noexcept {
    return M;
}


BARRAY_TEMPLATE(unsigned int, nnozero) () const noexcept {
    return NCells;
}

BARRAY_TEMPLATE(Cell< Cell_Type >, default_val) () const {
    return this->Cell_default;
}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator+=) (
    const std::pair<uint,uint> & coords
) {
    
    this->insert_cell(
            coords.first,
            coords.second,
            this->Cell_default,
            true, true
    );
    
    return *this;
    
}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator-=) (
    const std::pair<uint,uint> & coords
) {
    
    this->rm_cell(
            coords.first,
            coords.second,
            true, true
    );
    
    return *this;
    
}

template BARRAY_TEMPLATE_ARGS()
inline BArrayCell<Cell_Type,Data_Type> BARRAY_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) {
    
    return BArrayCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template BARRAY_TEMPLATE_ARGS()
inline const BArrayCell_const<Cell_Type,Data_Type> BARRAY_TYPE()::operator() (  
    uint i,
    uint j,
    bool check_bounds
) const {
    
    return BArrayCell_const<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

BARRAY_TEMPLATE(void, rm_cell) (
    uint i,
    uint j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);
    
    if (check_exists) {
        // Nothing to do
        if (ROW(i).size() == 0u)
            return;
        
        // Checking the counter part
        if (COL(j).size() == 0u)
            return;
        
        // Hard work, need to remove it from both, if it exist
        if (ROW(i).find(j) == ROW(i).end())
            return;
    }
    
    // Remove the pointer first (so it wont point to empty)
    COL(j).erase(i);
    ROW(i).erase(j);
    
    NCells--;
    
    return;
    
}

BARRAY_TEMPLATE(void, insert_cell) (
        uint i,
        uint j,
        const Cell< Cell_Type> & v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, insert_cell) (
        uint i,
        uint j,
        Cell< Cell_Type> && v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, insert_cell) (
    uint i,
    uint j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
        
    return insert_cell(i, j, Cell<Cell_Type>(v, visited), check_bounds, check_exists);

}

BARRAY_TEMPLATE(void, swap_cells) (
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
        Cell<Cell_Type> c0(std::move(ROW(i0)[j0]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(ROW(i1)[j1]));
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
        
    } else if (check_exists == CHECK::ONE) {
        
        check0 = !is_empty(i0, j0, false);
        check1 = true;
        
    } else if (check_exists == CHECK::TWO) {
        
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
        
        Cell<Cell_Type> c0(std::move(ROW(i0)[j0]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(ROW(i1)[j1]));
        rm_cell(i1, j1, false, false);
        
        insert_cell(i0, j0, c1, false, false);
        insert_cell(i1, j1, c0, false, false);
        
    } else if (!check0 & check1) { // If only the second exists
        
        if (report != nullptr) 
            (*report) = EXISTS::TWO;
        
        insert_cell(i0, j0, ROW(i1)[j1], false, false);
        rm_cell(i1, j1, false, false);
        
    } else if (check0 & !check1) {
        
        if (report != nullptr) 
            (*report) = EXISTS::ONE;
        
        insert_cell(i1, j1, ROW(i0)[j0], false, false);
        rm_cell(i0, j0, false, false);
        
    }
    
    return;
}

BARRAY_TEMPLATE(void, toggle_cell) (
    uint i,
    uint j,
    bool check_bounds,
    int check_exists
) {
    
    if (check_bounds)
        out_of_range(i, j);
    
    if (check_exists == EXISTS::UKNOWN) {
        
        if (is_empty(i, j, false)) {
            insert_cell(i, j, BArray<Cell_Type, Data_Type>::Cell_default, false, false);
            ROW(i)[j].visited = visited;
        } else
            rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ONE) {
        
        rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ZERO) {
        
        insert_cell(i, j, BArray<Cell_Type,Data_Type>::Cell_default, false, false);
        ROW(i)[j].visited = visited;
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, swap_rows) (
    uint i0,
    uint i1,
    bool check_bounds
) {
  
    if (check_bounds) {
        out_of_range(i0,0u);
        out_of_range(i1,0u);
    }
    
    bool move0=true, move1=true;
    if (ROW(i0).size() == 0u) move0 = false;
    if (ROW(i1).size() == 0u) move1 = false;
    
    if (!move0 && !move1)
        return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    ROW(i0).swap(ROW(i1));
    
    // Delete the thing
    if (move0)
        for (auto& i: row(i1, false))
            COL(i.first).erase(i0);
    
    if (move1)
        for (auto& i: row(i0, false))
            COL(i.first).erase(i1);
    
    // Now, point to the thing, if it has something to point at. Recall that
    // the indices swapped.
    if (move1)
        for (auto& i: row(i0, false))
            COL(i.first)[i0] = &ROW(i0)[i.first];
    
    if (move0)
        for (auto& i: row(i1, false))
            COL(i.first)[i1] = &ROW(i1)[i.first];
    
    return;

}

// This swapping is more expensive overall
BARRAY_TEMPLATE(void, swap_cols) (
    uint j0,
    uint j1,
    bool check_bounds
) {
  
    if (check_bounds) {
        out_of_range(0u, j0);
        out_of_range(0u, j1);
    }
    
    // Which ones need to be checked
    bool check0 = true, check1 = true;
    if (COL(j0).size() == 0u) check0 = false;
    if (COL(j1).size() == 0u) check1 = false;
    
    if (check0 && check1) {
      
        // Just swapping one at a time
        int status;
        Col_type<Cell_Type> col_tmp = COL(j1);
        Col_type<Cell_Type> col1 = COL(j0);
        for (auto iter = col1.begin(); iter != col1.end(); ++iter) {
            
            // Swapping values (col-wise)
            swap_cells(iter->first, j0, iter->first, j1, false, CHECK::TWO, &status);
            
            // Need to remove it, so we don't swap that as well
            if (status == EXISTS::BOTH)
                col_tmp.erase(iter->first);
        }
        
        // If there's anything left to move, we start moving it, otherwise, we just
        // skip it
        if (col_tmp.size() != 0u) {
          
            for (auto iter = col_tmp.begin(); iter != col_tmp.end(); ++iter) {
                insert_cell(iter->first, j0, *iter->second, false, false);
                rm_cell(iter->first, j1);
            }
          
        }
      
    } else if (check0 && !check1) {
      
        // 1 is empty, so we just add new cells and remove the other ones
        for (auto iter = COL(j0).begin(); iter != COL(j0).begin(); ++iter)
            insert_cell(iter->first, j1, *iter->second, false, false);
        
        // Setting the column to be zero
        COL(j0).empty();
      
    } else if (!check0 && check1) {
      
        // 1 is empty, so we just add new cells and remove the other ones
        for (auto iter = COL(j1).begin(); iter != COL(j1).begin(); ++iter) {
          
            // Swapping values (col-wise)
            insert_cell(iter->first, j0, *iter->second, false, false);

        }
        
        // Setting the column to be zero
        COL(j1).empty();
      
    }
    
    
    return;
}

BARRAY_TEMPLATE(void, zero_row) (
    uint i,
    bool check_bounds
) {
  
    if (check_bounds)
        out_of_range(i, 0u);
    
    // Nothing to do
    if (ROW(i).size() == 0u)
        return;
    
    // Else, remove all elements
    auto row0 = ROW(i);
    for (auto row = row0.begin(); row != row0.end(); ++row) 
        rm_cell(i, row->first, false, false);
    
    return;
  
}

BARRAY_TEMPLATE(void, zero_col) (
    uint j,
    bool check_bounds
) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    // Nothing to do
    if (COL(j).size() == 0u)
        return;
    
    // Else, remove all elements
    auto col0 = COL(j);
    for (auto col = col0.begin(); col != col0.end(); ++col) 
        rm_cell(col->first, j, false, false);
    
    return;
  
}

BARRAY_TEMPLATE(void, transpose) () {
  
    // Start by flipping the switch 
    visited = !visited;
    
    // Do we need to resize (increase) either?
    if      (N > M) el_ji.resize(N);
    else if (N < M) el_ij.resize(M);
    
    // uint N0 = N, M0 = M;
    int status;
    for (uint i = 0u; i < N; ++i)
    {
        
        // Do we need to move anything?
        if (ROW(i).size() == 0u)
            continue;
        
        // We now iterate changing rows
        Row_type<Cell_Type> row = ROW(i);
        for (auto col = row.begin(); col != row.end(); ++col)
        {
          
            // Skip if in the diagoal
            if (i == col->first)
            {
                ROW(i)[i].visited = visited;
                continue;
            }
            
            // We have not visited this yet, we need to change that
            if (ROW(i)[col->first].visited != visited)
            {
                
                // First, swap the contents
                swap_cells(i, col->first, col->first, i, false, CHECK::TWO, &status);
                
                // Changing the switch
                if (status == EXISTS::BOTH)
                    ROW(i)[col->first].visited = visited;
                
                ROW(col->first)[i].visited = visited;
              
            }
          
        }
      
    }
    
    // Shreding. Note that no information should have been lost since, hence, no
    // change in NCells.
    if (N > M) el_ij.resize(M);
    else if (N < M) el_ji.resize(N);
    
    // Swapping the values
    std::swap(N, M);
    
    return;

}

BARRAY_TEMPLATE(void, clear) (
    bool hard
) {
    
    if (hard)
    {
      
        el_ji.clear();
        el_ij.clear();
        
        el_ij.resize(N);
        el_ji.resize(M);
        NCells = 0u;
      
    } else {
        
        for (unsigned int i = 0u; i < N; ++i)
            zero_row(i, false);
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, resize) (
    uint N_,
    uint M_
) {
  
    // Removing rows
    if (N_ < N)
        for (uint i = N_; i < N; ++i)
            zero_row(i, false);
    
    // Removing cols
    if (M_ < M)
        for (uint j = M_; j < M; ++j)
            zero_col(j, false);
    
    // Resizing will invalidate pointers and values out of range
    if (M_ != M) {
        el_ji.resize(M_);
        M = M_;
    }
    
    if (N_ != N) {
        el_ij.resize(N_);
        N = N_;
    }
    
    
    return;

}

BARRAY_TEMPLATE(void, reserve) () {
#ifdef BARRAY_USE_UNORDERED_MAP
    for (uint i = 0u; i < N; i++)
        ROW(i).reserve(M);
    
    for (uint i = 0u; i < M; i++)
        COL(i).reserve(N);
#endif
    return;
  
}

BARRAY_TEMPLATE(void, print) (
    const char * fmt,
    ...
) const {
  
    std::va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    for (uint i = 0u; i < N; ++i)
    {

        #ifdef BARRY_DEBUG_LEVEL
            #if BARRY_DEBUG_LEVEL > 1
                printf_barry("%s [%3i,]", BARRY_DEBUG_HEADER, i);
            #endif
        #else
        printf_barry("[%3i,] ", i);
        #endif
        for (uint j = 0u; j < M; ++j) {
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

#undef BARRAY_TYPE
#undef BARRAY_TEMPLATE_ARGS
#undef BARRAY_TEMPLATE

#endif

