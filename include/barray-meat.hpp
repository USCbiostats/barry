#include <stdexcept>
#include "barray-bones.hpp"

#ifndef BARRAY_MEAT_HPP
#define BARRAY_MEAT_HPP 

// Edgelist with data
inline BArray::BArray (
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

  el_ij.resize(N);
  el_ji.resize(M);
  
  
  // Writing the data
  for (uint i = 0u; i < source.size(); ++i) {
    
    // Checking range
    if (source.at(i) >= N_ | target.at(i) >= M_)
      throw std::range_error("Either source or target point to an element outside of the range by (N,M).");
    
    // Checking if it exists
    auto search = ROW(source.at(i)).find(target.at(i));
    if (search != ROW(source.at(i)).end()) {
      if (!add)
        throw std::logic_error("The value already exists. Use 'add = true'.");
      
      // Increasing the value (this will automatically update the
      // other value)
      ROW(source.at(i))[target.at(i)].add(value.at(i));
      continue;
    }
    
    // Adding the value and creating a pointer to it
    ROW(source.at(i)).emplace(std::pair<uint, Cell>(target.at(i), Cell(value.at(i), this->visited)));
    COL(target.at(i)).emplace(
        source.at(i),
        &ROW(source.at(i))[target.at(i)]
    );
    this->NCells++;
  }
  
  return;
  
}

inline void BArray::out_of_range(uint i, uint j) const {
  if (i >= this->N)
    throw std::range_error("The row is out of range.");
  else if (j >= this->M)
    throw std::range_error("The column is out of range.");
  return;
}
  

inline double BArray::get_cell(uint i, uint j, bool check_bounds) const {
  
  // Checking boundaries  
  if (check_bounds)
    this->out_of_range(i,j);
  

  if (ROW(i).size() == 0u)
    return 0.0;
  
  // If it is not empty, then find and return
  auto search = ROW(i).find(j);
  if (search != ROW(i).end())
    return search->second.value;
  
  // This is if it is empty
  return 0.0;
  
}

inline const umap_int_cell * BArray::get_row(uint i, bool check_bounds) const {

  // Checking boundaries  
  if (check_bounds) 
    this->out_of_range(i, 0u);

  return &ROW(i);
}

inline const umap_int_cell_ptr * BArray::get_col(uint i, bool check_bounds) const {
  
  // Checking boundaries  
  if (check_bounds) 
    this->out_of_range(0u, i);
  
  return &COL(i);
}

inline Entries BArray::get_entries() const {
  
  Entries res(this->NCells);
  
  for (uint i = 0u; i < this->N; ++i) {
    
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

inline bool BArray::is_empty(uint i, uint j, bool check_bounds) const {
  
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

inline void BArray::rm_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  
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
  
  this->NCells--;
  
  return;
}

inline void BArray::insert_cell(uint i, uint j, std::pair< double, bool> v, bool check_bounds, bool check_exists) { 
  
  if (check_bounds)
    out_of_range(i,j);
  
  if (check_exists) {
    
    // Checking if nothing here, then we move along
    if (ROW(i).size() == 0u) {
      
      ROW(i).emplace(j, std::move(Cell(v.first, v.second)));
      COL(j).emplace(i, &ROW(i).at(j));
      this->NCells++;
      return;
      
    }
    
    // In this case, the row exists, but we are checking that the value is empty  
    if (ROW(i).find(j) == ROW(i).end()) {
      ROW(i).emplace(j, std::move(Cell(v.first, v.second)));
      COL(j).emplace(i, &ROW(i).at(j));
      this->NCells++;
    } else {
      throw std::logic_error("The cell already exists.");
    }
    
  } else {
    
    ROW(i).emplace(j, std::move(Cell(v.first, v.second)));
    COL(j).emplace(i, &ROW(i).at(j));
    
  }
  
  return;
  
}

inline void BArray::insert_cell(uint i, uint j, Cell & v, bool check_bounds, bool check_exists) {
  return this->insert_cell(i, j, std::pair<double, bool>(v.value, v.visited), check_bounds, check_exists);
}

inline void BArray::insert_cell(uint i, uint j, double v, bool check_bounds, bool check_exists) {
  return this->insert_cell(i, j, std::pair<double, bool>(v, this->visited), check_bounds, check_exists);
}

inline void BArray::insert_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  return this->insert_cell(i, j, std::pair<double,bool>(1.0, this->visited), check_bounds, check_exists);
  
}

inline void BArray::insert_cell(uint i, std::pair< double, bool > v, bool check_bounds, bool check_exists) {
  
  return this->insert_cell(
      (int) i % (int) this->N,
      floor((int) i / (int) this->N),
      v,
      check_bounds,
      check_exists
  );
}

inline void BArray::insert_cell(uint i, Cell & v, bool check_bounds, bool check_exists) {
  
  return this->insert_cell(
      (int) i % (int) this->N,
      floor((int) i / (int) this->N),
      std::pair<double,bool>(v.value, v.visited), check_bounds, check_exists
  );
  
}

inline void BArray::insert_cell(uint i, double v, bool check_bounds, bool check_exists) {
  
  return this->insert_cell(
      (int) i % (int) this->N,
      floor((int) i / (int) this->N),
      std::pair<double,bool>(v, this->visited),
      check_bounds,
      check_exists
  );
}

inline void BArray::insert_cell(uint i, bool check_bounds, bool check_exists) {

  return this->insert_cell(
      (int) i % (int) this->N,
      floor((int) i / (int) this->N),
      std::pair<double,bool>(1.0, this->visited),
      check_bounds,
      check_exists
  );
}

inline void BArray::swap_cells(
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
    
    // Using the initializing by move, after this, the cell becomes
    // invalid. We use pointers instead as this way we access the Heap memory,
    // which should be faster to access.
    Cell c0(std::move(ROW(i0).at(j0)));
    this->rm_cell(i0, j0, false, false);
    Cell c1(std::move(ROW(i1).at(j1)));
    this->rm_cell(i1, j1, false, false);
    
    // Inserting the cells by reference, these will be deleted afterwards
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
    
    Cell c0(std::move(ROW(i0).at(j0)));
    this->rm_cell(i0, j0, false, false);
    Cell c1(std::move(ROW(i1).at(j1)));
    this->rm_cell(i1, j1, false, false);
    
    this->insert_cell(i0, j0, c1, false, false);
    this->insert_cell(i1, j1, c0, false, false);
    
  } else if (!check0 & check1) { // If only the second exists
    
    if (report != nullptr) 
      (*report) = EXISTS::TWO;
    
    this->insert_cell(i0, j0, ROW(i1).at(j1), false, false);
    this->rm_cell(i1, j1, false, false);
    
  } else if (check0 & !check1) {
    
    if (report != nullptr) 
      (*report) = EXISTS::ONE;
    
    this->insert_cell(i1, j1, ROW(i0).at(j0), false, false);
    this->rm_cell(i0, j0, false, false);
    
  }
  
  return;
}

inline void BArray::toggle_cell(uint i, uint j, bool check_bounds, int check_exists) {
  
  if (check_bounds)
    out_of_range(i, j);
  
  if (check_exists == EXISTS::UKNOWN) {
    
    if (this->is_empty(i, j, false))
      this->insert_cell(i, j, std::pair<double,bool>(1.0, this->valued), false, false);
    else
      this->rm_cell(i, j, false, false);
    
  } else if (check_exists == EXISTS::AS_ONE) {
    
    this->rm_cell(i, j, false, false);
    
  } else if (check_exists == EXISTS::AS_ZERO) {
    
    this->insert_cell(i, j, std::pair<double,bool>(1.0, this->valued), false, false);
    
  }
  
  return;
  
}

inline void BArray::swap_rows(uint i0, uint i1, bool check_bounds) {
  
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
inline void BArray::swap_cols(uint j0, uint j1, bool check_bounds) {
  
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
    umap_int_cell_ptr col_tmp = COL(j1);
    umap_int_cell_ptr col1 = COL(j0);
    for (auto iter = col1.begin(); iter != col1.end(); ++iter) {
      
      // Swapping values (col-wise)
      this->swap_cells(iter->first, j0, iter->first, j1, false, CHECK::TWO, &status);
      
      // Need to remove it, so we don't swap that as well
      if (status == EXISTS::BOTH)
        col_tmp.erase(iter->first);
    }
    
    // If there's anything left to move, we start moving it, otherwise, we just
    // skip it
    if (col_tmp.size() != 0u) {
      
      for (auto iter = col_tmp.begin(); iter != col_tmp.end(); ++iter) {
        this->insert_cell(iter->first, j0, *iter->second, false, false);
        this->rm_cell(iter->first, j1);
      }
      
    }
    
  } else if (check0 && !check1) {
    
    // 1 is empty, so we just add new cells and remove the other ones
    for (auto iter = COL(j0).begin(); iter != COL(j0).begin(); ++iter) 
      this->insert_cell(iter->first, j1, *iter->second, false, false);
    
    // Setting the column to be zero
    COL(j0).empty();
    
  } else if (!check0 && check1) {
    
    // 1 is empty, so we just add new cells and remove the other ones
    for (auto iter = COL(j1).begin(); iter != COL(j1).begin(); ++iter) {
      
      // Swapping values (col-wise)
      this->insert_cell(iter->first, j0, *iter->second, false, false);

    }
    
    // Setting the column to be zero
    COL(j1).empty();
    
  }
  
  
  return;
}

inline void BArray::zero_row(uint i, bool check_bounds) {
  
  if (check_bounds)
    this->out_of_range(i, 0u);
  
  // Nothing to do
  if (ROW(i).size() == 0u)
    return;
  
  // Else, remove all elements
  auto row0 = ROW(i);
  for (auto row = row0.begin(); row != row0.end(); ++row) 
    this->rm_cell(i, row->first, false, false);
  
  return;
  
}

inline void BArray::zero_col(uint j, bool check_bounds) {
  
  if (check_bounds)
    this->out_of_range(0u, j);
  
  // Nothing to do
  if (COL(j).size() == 0u)
    return;
  
  // Else, remove all elements
  auto col0 = COL(j);
  for (auto col = col0.begin(); col != col0.end(); ++col) 
    this->rm_cell(col->first, j, false, false);
  
  return;
  
}

inline void BArray::transpose() {
  
  // Start by flipping the switch 
  this->visited = !this->visited;
  
  // Do we need to resize (increase) either?
  if      (this->N > this->M) this->el_ji.resize(this->N);
  else if (this->N < this->M) this->el_ij.resize(this->M);
  
  // uint N0 = this->N, M0 = this->M;
  int status;
  for (uint i = 0u; i < this->N; ++i) {
    
    // Do we need to move anything?
    if (ROW(i).size() == 0u)
      continue;
    
    // We now iterate changing rows
    umap_int_cell row = ROW(i);
    for (auto col = row.begin(); col != row.end(); ++col) {
      
      // Skip if in the diagoal
      if (i == col->first) {
        ROW(i).at(i).visited = this->visited;
        continue;
      }
      
      // We have not visited this yet, we need to change that
      if (ROW(i).at(col->first).visited != this->visited) {
        
        // First, swap the contents
        this->swap_cells(i, col->first, col->first, i, false, CHECK::TWO, &status);
        
        // Changing the switch
        if (status == EXISTS::BOTH)
          ROW(i).at(col->first).visited = this->visited;
        
        ROW(col->first).at(i).visited = this->visited;
        
      }
      
    }
    
  }
  
  // Shreding. Note that no information should have been lost since, hence, no
  // change in NCells.
  if (this->N > this->M) this->el_ij.resize(this->M);
  else if (this->N < this->M) this->el_ji.resize(this->N);
  
  // Swapping the values
  std::swap(this->N, this->M);
  
  return;
}

inline void BArray::clear() {
  this->el_ji.clear();
  this->el_ij.clear();
  
  this->el_ij.resize(this->N);
  this->el_ji.resize(this->M);
  this->NCells = 0u;
  
  return;
  
}

inline void BArray::resize(uint N_, uint M_) {
  
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

#endif

