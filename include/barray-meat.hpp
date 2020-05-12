#include <stdexcept>
#include "barray-bones.hpp"

#ifndef BARRAY_MEAT_HPP
#define BARRAY_MEAT_HPP 

// Edgelist with data
BArray::BArray (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value,
    bool add
) {
  
  if (source.size() != target.size())
    throw std::length_error("-source- and -target- don't match on length.");
  if (source.size() != value.size())
    throw std::length_error("-sorce- and -value- don't match on length.");
  
  // Initializing
  this->N = N_;
  this->M = M_;
  this->visited = false;
  this->NCells = 0u;
  
  el_ij.resize(N);
  el_ji.resize(M);
  
  
  // Writing the data
  for (uint i = 0u; i < source.size(); ++i) {
    
    // Checking range
    if (source.at(i) >= N_ | target.at(i) >= M_)
      throw std::range_error("Either source or target point to an element outside of the range by (N,M).");
    
    // Checking if it exists
    auto search = el_ij.at(source.at(i)).find(target.at(i));
    if (search != el_ij.at(source.at(i)).end()) {
      if (!add)
        throw std::logic_error("The value already exists. Use 'add = true'.");
      
      // Increasing the value (this will automatically update the
      // other value)
      el_ij.at(source.at(i))[target.at(i)].add(value.at(i));
      continue;
    }
    
    // Adding the value and creating a pointer to it
    el_ij.at(source.at(i)).emplace(target.at(i), Cell(value.at(i), this->visited));
    el_ji.at(target.at(i)).emplace(source.at(i), &el_ij.at(source.at(i))[target.at(i)]);
    this->NCells++;
  }
  
  return;
  
}

void BArray::out_of_range(uint i, uint j) const {
  if (i >= this->N)
    throw std::range_error("The row is out of range.");
  else if (j >= this->M)
    throw std::range_error("The column is out of range.");
  return;
}
  

double BArray::get_cell(uint i, uint j, bool check_bounds) const {
  
  // Checking boundaries  
  if (check_bounds)
    out_of_range(i,j);
  

  if (this->el_ij.at(i).size() == 0u)
    return 0.0;
  
  // If it is not empty, then find and return
  auto search = el_ij.at(i).find(j);
  if (search != el_ij.at(i).end())
    return search->second.value;
  
  // This is if it is empty
  return 0.0;
  
}

const umap_int_cell * BArray::get_row(uint i, bool check_bounds) const {

  // Checking boundaries  
  if (check_bounds) 
    out_of_range(i, 0u);

  return &(el_ij.at(i));
}

const umap_int_cell_ptr * BArray::get_col(uint i, bool check_bounds) const {
  
  // Checking boundaries  
  if (check_bounds) 
    out_of_range(0u, i);
  
  return &(el_ji.at(i));
}


bool BArray::is_empty(uint i, uint j, bool check_bounds) const {
  
  if (check_bounds)
    out_of_range(i, j);
  
  if (this->el_ij.at(i).size() == 0u)
    return true;
  else if (this->el_ji.at(j).size() == 0u)
    return true;
  
  if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end())
    return true;
  
  return false;
  
}

void BArray::rm_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  
  // Checking the boundaries
  if (check_bounds)
    out_of_range(i,j);
  
  if (check_exists) {
    // Nothing to do
    if (this->el_ij.at(i).size() == 0u)
      return;
    
    // Checking the counter part
    if (this->el_ji.at(j).size() == 0u)
      return;
    
    // Hard work, need to remove it from both, if it exist
    if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end())
      return;
  }
  
  // Remove the pointer first (so it wont point to empty)
  this->el_ji.at(j).erase(i);
  this->el_ij.at(i).erase(j);
  
  this->NCells--;
  
  return;
}

void BArray::insert_cell(uint i, uint j, Cell v, bool check_bounds, bool check_exists) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  if (check_exists) {
    
    // Checking if nothing here, then we move along
    if (this->el_ij.at(i).size() == 0u) {
      
      this->el_ij.at(i).emplace(j, v);
      this->el_ji.at(j).emplace(i, &this->el_ij.at(i).at(j));
      this->NCells++;
      return;
      
    }
    
    // In this case, the row exists, but we are checking that the value is empty  
    if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end()) {
      this->el_ij.at(i).emplace(j, v);
      this->el_ji.at(j).emplace(i, &this->el_ij.at(i).at(j));
      this->NCells++;
    } else {
      throw std::logic_error("The cell already exists.");
    }
    
  } else {
    
    this->el_ij.at(i).emplace(j, v);
    this->el_ji.at(j).emplace(i, &this->el_ij.at(i).at(j));
    
  }
  
  return;
  
}

void BArray::insert_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  return this->insert_cell(i, j, Cell(1.0, this->visited), check_bounds, check_exists);
}

void BArray::insert_cell(uint i, uint j, double v, bool check_bounds, bool check_exists) {
  return this->insert_cell(i, j, Cell(v, this->visited), check_bounds, check_exists);
}

void BArray::insert_cell(uint i, bool check_bounds, bool check_exists) {
  // std::cout << "inserting by single i: " << i <<" at(" << i0 <<", " << j0 << " )"  << std::endl;
  return this->insert_cell(
      (int) i % (int) this->N,
      floor((int) i / (int) this->N),
      Cell(1.0, this->visited),
      check_bounds,
      check_exists
  );
}

void BArray::swap_cells(
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
  if (check_exists == CHECK::NONE) {
    
    // Just in case, if this was passed
    if (report != nullptr)
      (*report) = EXISTS::BOTH;
    
    // If source and target coincide, we do nothing
    if ((i0 == i1) && (j0 == j1)) 
      return;
    
    Cell c0 = this->el_ij.at(i0).at(j0);
    this->rm_cell(i0, j0, false, false);
    Cell c1 = this->el_ij.at(i1).at(j1);
    this->rm_cell(i1, j1, false, false);
    
    this->insert_cell(i0, j0, c1, false, false);
    this->insert_cell(i1, j1, c0, false, false);
    
    return;
    
  }
  
  bool check0, check1;
  if (check_exists == CHECK::BOTH) {
    
    check0 = !this->is_empty(i0, j0, false);
    check1 = !this->is_empty(i1, j1, false);
    
  } else if (check_exists == CHECK::ONE) {
    
    check0 = !this->is_empty(i0, j0, false);
    check1 = true;
    
  } else if (check_exists == CHECK::TWO) {
    
    check0 = true;
    check1 = !this->is_empty(i1, j1, false);
    
  }
  
  if (report != nullptr) 
    (*report) = EXISTS::NONE;
  
  // If both cells exists
  if (check0 & check1) {
    
    if (report != nullptr) 
      (*report) = EXISTS::BOTH;
    
    // If source and target coincide, we do nothing
    if ((i0 == i1) && (j0 == j1)) 
      return;
    
    Cell c0 = this->el_ij.at(i0).at(j0);
    this->rm_cell(i0, j0, false, false);
    Cell c1 = this->el_ij.at(i1).at(j1);
    this->rm_cell(i1, j1, false, false);
    
    this->insert_cell(i0, j0, c1, false, false);
    this->insert_cell(i1, j1, c0, false, false);
    
  } else if (!check0 & check1) { // If only the second exists
    
    if (report != nullptr) 
      (*report) = EXISTS::TWO;
    
    this->insert_cell(i0, j0, this->el_ij.at(i1).at(j1), false, false);
    this->rm_cell(i1, j1, false, false);
    
  } else if (check0 & !check1) {
    
    if (report != nullptr) 
      (*report) = EXISTS::ONE;
    
    this->insert_cell(i1, j1, this->el_ij.at(i0).at(j0), false, false);
    this->rm_cell(i0, j0, false, false);
    
  }
  
  return;
}

void BArray::toggle_cell(uint i, uint j, bool check_bounds, int check_exists) {
  
  if (check_bounds)
    out_of_range(i, j);
  
  if (check_exists == EXISTS::UKNOWN) {
    
    if (this->is_empty(i, j, false))
      this->insert_cell(i, j, 1.0, false, false);
    else
      this->rm_cell(i, j, false, false);
    
  } else if (check_exists == EXISTS::AS_ONE) {
    
    this->rm_cell(i, j, false, false);
    
  } else if (check_exists == EXISTS::AS_ZERO) {
    
    this->insert_cell(i, j, 1.0, false, false);
    
  }
  
  return;
  
}

void BArray::swap_rows(uint i0, uint i1, bool check_bounds) {
  
  if (check_bounds) {
    out_of_range(i0,0u);
    out_of_range(i1,0u);
  }
  
  bool move0=true, move1=true;
  if (this->el_ij.at(i0).size() == 0u) move0 = false;
  if (this->el_ij.at(i1).size() == 0u) move1 = false;
  
  if (!move0 && !move1)
    return;
  
  // Swapping happens naturally, need to take care of the pointers
  // though
  this->el_ij.at(i0).swap(this->el_ij.at(i1));
  
  // Delete the thing
  if (move0)
    for (auto iter = el_ij.at(i1).begin(); iter != el_ij.at(i1).end(); ++iter)
      el_ji.at(iter->first).erase(i0);
  
  if (move1)
    for (auto iter = el_ij.at(i0).begin(); iter != el_ij.at(i0).end(); ++iter)
      el_ji.at(iter->first).erase(i1);
  
  // Now, point to the thing, if it has something to point at. Recall that
  // the indices swapped.
  if (move1)
    for (auto iter = el_ij.at(i0).begin(); iter != el_ij.at(i0).end(); ++iter)
      el_ji.at(iter->first)[i0] = &el_ij.at(i0).at(iter->first);
  
  if (move0)
    for (auto iter = el_ij.at(i1).begin(); iter != el_ij.at(i1).end(); ++iter)
      el_ji.at(iter->first)[i1] = &el_ij.at(i1).at(iter->first);
  
  return;
}

// This swapping is more expensive overall
void BArray::swap_cols(uint j0, uint j1, bool check_bounds) {
  
  if (check_bounds) {
    out_of_range(0u, j0);
    out_of_range(0u, j1);
  }
  
  bool move0=true, move1=true;
  if (this->el_ji.at(j0).size() == 0u) move0 = false;
  if (this->el_ji.at(j1).size() == 0u) move1 = false;
  
  if (!move0 && !move1)
    return;
  
  // We now need to swap the contents internally
  if (move0 && move1) {
    
    // We get a copy that we use to check whether we already checked these
    // values or not
    umap_int_cell_ptr col0 = this->el_ji.at(j0);
    umap_int_cell_ptr col1 = this->el_ji.at(j1);
    int status;
    for (auto col = col0.begin(); col != col0.end(); ++col) {
      
      // We first swap
      this->swap_cells(col->first, j0, col->first, j1, false, CHECK::TWO, &status);
      
      // We tagged the ones that were already swapped as done
      if (status == EXISTS::BOTH)
        col1.erase(col->first);
      
    }
    
    for (auto col = col1.begin(); col != col1.end(); ++col) 
      this->swap_cells(col->first, j1, col->first, j0, false, CHECK::ONE);
    
  }
  
  
  
  return;
}

void BArray::zero_row(uint i, bool check_bounds) {
  
  if (check_bounds)
    this->out_of_range(i, 0u);
  
  // Nothing to do
  if (this->el_ij.at(i).size() == 0u)
    return;
  
  // Else, remove all elements
  auto row0 = this->el_ij.at(i);
  for (auto row = row0.begin(); row != row0.end(); ++row) 
    this->rm_cell(i, row->first, false, false);
  
  return;
  
}

void BArray::zero_col(uint j, bool check_bounds) {
  
  if (check_bounds)
    this->out_of_range(0u, j);
  
  // Nothing to do
  if (this->el_ji.at(j).size() == 0u)
    return;
  
  // Else, remove all elements
  auto col0 = this->el_ji.at(j);
  for (auto col = col0.begin(); col != col0.end(); ++col) 
    this->rm_cell(col->first, j, false, false);
  
  return;
  
}

void BArray::transpose() {
  
  // Start by flipping the switch 
  this->visited = !this->visited;
  
  // Do we need to resize (increase) either?
  if      (this->N > this->M) this->el_ji.resize(this->N);
  else if (this->N < this->M) this->el_ij.resize(this->M);
  
  // uint N0 = this->N, M0 = this->M;
  int status;
  for (uint i = 0u; i < this->N; ++i) {
    
    // Do we need to move anything?
    if (this->el_ij.at(i).size() == 0u)
      continue;
    
    // We now iterate changing rows
    umap_int_cell row = this->el_ij.at(i);
    for (auto col = row.begin(); col != row.end(); ++col) {
      
      // We have not visited this yet, we need to change that
      if (this->el_ij.at(i).at(col->first).visited != this->visited) {
        
        // First, swap the contents
        this->swap_cells(i, col->first, col->first, i, false, CHECK::TWO, &status);
        
        // Changing the switch
        if (status == EXISTS::BOTH)
          this->el_ij.at(i).at(col->first).visited = this->visited;
        
        this->el_ij.at(col->first).at(i).visited = this->visited;
        
      }
      
    }
    
  }
  
  // Shreding
  if (this->N > this->M) this->el_ij.resize(this->M);
  else if (this->N < this->M) this->el_ji.resize(this->N);
  
  // Swapping the values
  std::swap(this->N, this->M);
  
  return;
}

void BArray::clear() {
  this->el_ji.clear();
  this->el_ij.clear();
  
  this->el_ij.resize(this->N);
  this->el_ji.resize(this->M);
  this->NCells = 0u;
  
  return;
  
}

void BArray::resize(uint N_, uint M_) {
  
  // Removing rows
  if (N_ < this->N)
    for (uint i = N_; i < this->N; ++i)
      this->zero_row(i, false);
  
  // Removing cols
  if (M_ < this->M)
    for (uint j = M_; j < this->M; ++j)
      this->zero_col(j, false);
  
  // Resizing will invalidate pointers and values out of range
  if (N_ != this->N) {
    this->el_ji.resize(M_);
    this->N = N_;
  }
  
  if (N_ != this->M) {
    this->el_ij.resize(N_);
    this->M = M_;
  }
  
  return;

}

class LArray {
public:
  std::vector< BArray* > data;
  LArray() {};
  LArray(uint n) : data(n) {};
  ~LArray() {};
};

#endif

