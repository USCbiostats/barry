#include <Rcpp.h>
using namespace Rcpp;

#ifndef ARRAY_H
#define ARRAY_H 1


// Constants
#define ARRAY_CHECK_BOTH -1
#define ARRAY_CHECK_NONE 0
#define ARRAY_CHECK_ONE 1
#define ARRAY_CHECK_TWO 2

// Definition of the class structure

// Edgelist
typedef unsigned int uint;
class Cell;
typedef std::unordered_map< uint, Cell > umap_int_cell;
typedef std::unordered_map< uint, Cell* > umap_int_cell_ptr;

/* A cell is a fundamental type that holds information about a cell
 * for now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
class Cell {
public:
  double value;
  bool visited;
  Cell() {};
  Cell(double value_) : value(value_), visited(false) {};
  Cell(double value_, bool visited_) : value(value_), visited(visited_) {};
  ~Cell() {};
  
  void add(double x) {this->value+=x;};
};

class Array {
public:
  uint N;
  uint M;
  std::vector< umap_int_cell >  el_ij;
  std::vector< umap_int_cell_ptr > el_ji;
  
  // This is as a reference, if we need to iterate through the cells and we need
  // to keep track which were visited, we use this as a reference. So that if
  // cell.visited = true and visited = true, it means that we haven't been here
  // yet. Ideally, any routine using this->visited should switch it at the
  // beginning of the routine.
  bool visited;
  
  // Empty datum
  Array (uint N_, uint M_) : N(N_), M(M_), el_ij(N_), el_ji(M_) {};
  
  // Edgelist with data
  Array (
      uint N_, uint M_,
      const std::vector< uint > & source,
      const std::vector< uint > & target,
      const std::vector< double > & value,
      bool add = true
  );
  
  // Function to access the elements
  // bool check_cell
  void out_of_range(uint i, uint j) const {
    bool ans = ((i >= this->N) || (j >= this->M)) ? true : false;
    if (ans)
      Rcpp::stop("Out of range!.");
    return;
  }
  double get_cell(uint i, uint j, bool check_bounds = true) const;
  const umap_int_cell * get_row(uint i, bool check_bounds = true) const {return &(el_ij.at(i));}
  const umap_int_cell_ptr * get_col(uint i, bool check_bounds = true) const {return &(el_ji.at(i));}
  bool is_empty(uint i, uint j, bool check_bounds = true) const;
  
  // Deletion addition operations
  void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void insert_cell(uint i, uint j, double v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, Cell v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void swap_cells(uint i0, uint j0, uint i1, uint j1, bool check_bounds = true, int check_exists = ARRAY_CHECK_BOTH);
  
  void swap_rows(uint i0, uint i1, bool check_bounds = true);
  void swap_cols(uint j0, uint j1, bool check_bounds = true);
  
  void transpose();
  
  
};

// Edgelist with data
Array::Array (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value,
    bool add
) {
  
  if (source.size() != target.size())
    Rcpp::stop("Must match the size.");
  if (source.size() != value.size())
    Rcpp::stop("Must match the size.");
  
  // Initializing
  this->N = N_;
  this->M = M_;
  this->visited = false;
  el_ij.resize(N);
  el_ji.resize(M);
  
  
  // Writing the data
  for (uint i = 0u; i < source.size(); ++i) {
    
    // Checking range
    if (source.at(i) >= N_ | target.at(i) >= M_)
      Rcpp::stop("Out of range.");
    
    // Checking if it exists
    auto search = el_ij.at(source.at(i)).find(target.at(i));
    if (search != el_ij.at(source.at(i)).end()) {
      if (!add)
        Rcpp::stop("The value already exists");
      
      // Increasing the value (this will automatically update the
      // other value)
      el_ij.at(source.at(i))[target.at(i)].add(value.at(i));
      continue;
    }
    
    // Adding the value and creating a pointer to it
    el_ij.at(source.at(i))[target.at(i)] = Cell(value.at(i), this->visited);
    el_ji.at(target.at(i))[source.at(i)] = &el_ij.at(source.at(i))[target.at(i)];
  }
  
  return;
  
}

double Array::get_cell(uint i, uint j, bool check_bounds) const {
  
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

bool Array::is_empty(uint i, uint j, bool check_bounds) const {
  
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

void Array::rm_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  
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
  
  return;
}

/***
 * This member adds a new object to the edgelist
 */
void Array::insert_cell(uint i, uint j, double v, bool check_bounds, bool check_exists) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  if (check_exists) {
    
    // Checking if nothing here, then we move along
    if (this->el_ij.at(i).size() == 0u) {
      
      this->el_ij.at(i)[j] = Cell(v);
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
      return;
      
    }
    
    // In this case, the row exists, but we are checking that the value is empty
    if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end()) {
      this->el_ij.at(i)[j] = Cell(v);
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    }
    
  } else 
    this->el_ij.at(i).at(j).add(v);
  
  return;
  
}

void Array::insert_cell(uint i, uint j, Cell v, bool check_bounds, bool check_exists) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  if (check_exists) {
    
    // Checking if nothing here, then we move along
    if (this->el_ij.at(i).size() == 0u) {
      
      this->el_ij.at(i)[j] = v;
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
      return;
      
    }
  
    // In this case, the row exists, but we are checking that the value is empty  
    if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end()) {
      this->el_ij.at(i)[j] = v;
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    } else {
      Rcpp::stop("Cells cannot be added together.");
    }
      
  } else {
    
    this->el_ij.at(i)[j] = v;
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    
  }
  
  return;
  
}

void Array::insert_cell(uint i, uint j, bool check_bounds, bool check_exists) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  
  if (check_exists) {
    // Checking if nothing here, then we move along
    if (this->el_ij.at(i).size() == 0u) {
      
      this->el_ij.at(i)[j] = Cell(1.0);
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
      return;
      
    }
  
    // In this case, the row exists, but we are checking that the value is empty
    if (this->el_ij.at(i).find(j) == this->el_ij.at(i).end()) {
      this->el_ij.at(i)[j] = Cell(1.0);
      this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    } else
      this->el_ij.at(i).at(j).add(1.0);
      
    
  } else {
    
    // No checks to see if it exists, whatsoever, we just add it.
    this->el_ij.at(i)[j] = Cell(1.0);
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    
  }
  
  return;
  
}

void Array::swap_cells(uint i0, uint j0, uint i1, uint j1, bool check_bounds, int check_exists) {
  
  if (check_bounds) {
    out_of_range(i0,j0);
    out_of_range(i1,j1);
  }
  
  // If source and target coincide, we do nothing
  if ((i0 == i1) && (j0 == j1))
    return;
  
  // Simplest case, we know both exists, so we don't need to check anything
  if (check_exists == ARRAY_CHECK_NONE) {
    
    Cell c0 = this->el_ij.at(i0).at(j0);
    this->rm_cell(i0, j0, false, false);
    Cell c1 = this->el_ij.at(i1).at(j1);
    this->rm_cell(i1, j1, false, false);
    
    this->insert_cell(i0, j0, c1, false, false);
    this->insert_cell(i1, j1, c0, false, false);
    
    return;

  }
    
  bool move0 = true, move1 = true;
  
  // Do we need to move 0?
  bool check0 = (check_exists == ARRAY_CHECK_BOTH) || (check_exists == ARRAY_CHECK_ONE);
  bool check1 = (check_exists == ARRAY_CHECK_BOTH) || (check_exists == ARRAY_CHECK_TWO);
  
  if (check0) {
    if      (this->el_ij.at(i0).size() == 0u) move0 = false;
    else if (this->el_ji.at(j0).size() == 0u) move0 = false;
  }
  
  // How about 1?
  if (check1) {
    if      (this->el_ij.at(i1).size() == 0u) move1 = false;
    else if (this->el_ji.at(j1).size() == 0u) move1 = false;
  }
  
  
  // Case 1: Both need to be moved:
  if (move0 && move1) {
    
    // Swapping elements, the pointers will work OK as the memory
    // location hasn't change, right?
    check0 = !check0 || (this->el_ij.at(i0).find(j0) != this->el_ij.at(i0).end());
    check1 = !check1 || (this->el_ij.at(i1).find(j1) != this->el_ij.at(i1).end());
    
    // If both cells exists
    if (check0 & check1) {
      
      Cell c0 = this->el_ij.at(i0).at(j0);
      this->rm_cell(i0, j0, false, false);
      Cell c1 = this->el_ij.at(i1).at(j1);
      this->rm_cell(i1, j1, false, false);
      
      this->insert_cell(i0, j0, c1, false, false);
      this->insert_cell(i1, j1, c0, false, false);
      
    } else if (!check0 & check1) { // If only the second exists
      
      this->insert_cell(i0, j0, this->el_ij.at(i1).at(j1), false, false);
      this->rm_cell(i1, j1, false, false);
      
    } else if (check0 & !check1) {
      
      this->insert_cell(i1, j1, this->el_ij.at(i0).at(j0), false, false);
      this->rm_cell(i0, j0, false, false);
      
    }
    
  } else if (move0 && !move1) {
    
    check0 = !check0 || (this->el_ij.at(i0).find(j0) != this->el_ij.at(i0).end());
    
    // Adding and removing at the same time
    if (check0) {
      this->insert_cell(i1, j1, this->el_ij.at(i0).at(j0), false, false);
      this->rm_cell(i0, j0, false, false);
    }
    
  } else if (!move0 && move1) {
    
    check1 = !check1 || (this->el_ij.at(i1).find(j1) != this->el_ij.at(i1).end());
    
    // Likewise, but the other element
    if (check1) {
      this->insert_cell(i0, j0, this->el_ij.at(i1).at(j1), false, false);
      this->rm_cell(i1, j1, false, false);
    }
    
  }
  
  return;
}


void Array::swap_rows(uint i0, uint i1, bool check_bounds) {
  
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
void Array::swap_cols(uint j0, uint j1, bool check_bounds) {
  
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
    for (auto col = col0.begin(); col != col0.end(); ++col) {
      
      // We first swap
      this->swap_cells(col->first, j0, col->first, j1, false, ARRAY_CHECK_TWO);
      
      // We tagged the ones that were already swapped as done
      if (col1.find(col->first) != col1.end())
        col1.erase(col->first);

    }
    
    for (auto col = col1.begin(); col != col1.end(); ++col) 
      this->swap_cells(col->first, j1, col->first, j0, false, ARRAY_CHECK_ONE);
    
  }
  

  
  return;
}

void Array::transpose() {
  
  // Start by flipping the switch 
  this->visited = !this->visited;
  
  // Do we need to resize (increase) either?
  if      (this->N > this->M) this->el_ji.resize(this->N);
  else if (this->N < this->M) this->el_ij.resize(this->M);
  
  // uint N0 = this->N, M0 = this->M;

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
        this->swap_cells(i, col->first, col->first, i, false, ARRAY_CHECK_TWO);
        
        // Changing the switch
        if (!this->is_empty(i, col->first, false))
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

class LArray {
public:
  std::vector< Array* > data;
  LArray() {};
  LArray(uint n) : data(n) {};
  ~LArray() {};
};

#endif