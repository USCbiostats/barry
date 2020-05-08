#ifndef ARRAY_H
#define ARRAY_H 1

#include <Rcpp.h>
using namespace Rcpp;


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
  ~Cell() {};
  
  void add(double x) {this->value+=x;};
};

class Array {
public:
  uint N;
  uint M;
  std::vector< umap_int_cell >  el_ij;
  std::vector< umap_int_cell_ptr > el_ji;
  
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
  
  // Deletion addition operations
  void rm_cell(uint i, uint j, bool check_bounds = true);
  void insert_cell(uint i, uint j, double v, bool check_bounds = true);
  void insert_cell(uint i, uint j, bool check_bounds = true);
  void swap_cells(uint i0, uint j0, uint i1, uint j1, bool check_bounds = true);
  void swap_rows(uint i0, uint i1, bool check_bounds = true);
  void swap_col(uint j0, uint j1, bool check_bounds = true);
  
  
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
  N = N_;
  M = M_;
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
    el_ij.at(source.at(i))[target.at(i)] = Cell(value.at(i));
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

void Array::rm_cell(uint i, uint j, bool check_bounds) {
  
  // Checking the boundaries
  if (check_bounds)
    out_of_range(i,j);

  // Nothing to do
  if (this->el_ij.at(i).size() == 0u)
    return;
  
  // Checking the counter part
  if (this->el_ji.at(j).size() == 0u)
    return;
  
  // Hard work, need to remove it from both, if it exist
  auto search = this->el_ij.at(i).find(j);
  if (search == this->el_ij.at(i).end())
    return;
  
  // Remove the pointer first (so it wont point to empty)
  this->el_ji.at(j).erase(i);
  this->el_ij.at(i).erase(j);
  
  return;
}

/***
 * This member adds a new object to the edgelist
 */
void Array::insert_cell(uint i, uint j, double v, bool check_bounds) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  // Checking if nothing here, then we move along
  if (this->el_ij.at(i).size() == 0u) {
    
    this->el_ij.at(i)[j] = Cell(v);
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    return;
    
  }
  
  // In this case, the row exists, but we are checking that the value is empty
  auto search = this->el_ij.at(i).find(j);
  if (search == this->el_ij.at(i).end()) {
    this->el_ij.at(i)[j] = Cell(v);
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
  } else {
    this->el_ij.at(i).at(j).add(v);
  }
  
  return;
  
}

void Array::insert_cell(uint i, uint j, bool check_bounds) {
  
  if (check_bounds)
    out_of_range(i,j);
  
  // Checking if nothing here, then we move along
  if (this->el_ij.at(i).size() == 0u) {
    
    this->el_ij.at(i)[j] = Cell(1.0);
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
    return;
    
  }
  
  // In this case, the row exists, but we are checking that the value is empty
  auto search = this->el_ij.at(i).find(j);
  if (search == this->el_ij.at(i).end()) {
    this->el_ij.at(i)[j] = Cell(1.0);
    this->el_ji.at(j)[i] = &this->el_ij.at(i).at(j);
  } else {
    this->el_ij.at(i).at(j).add(1.0);
  }
  
  return;
  
}

void Array::swap_cells(uint i0, uint j0, uint i1, uint j1, bool check_bounds) {
  
  if (check_bounds) {
    out_of_range(i0,j0);
    out_of_range(i1,j1);
  }
    
  bool move0 = true, move1 = true;
  
  // Do we need to move 0?
  if (this->el_ij.at(i0).size() == 0u) move0 = false;
  else if (this->el_ji.at(j0).size() == 0u) move0 = false;
  
  // How about 1?
  if (this->el_ij.at(i1).size() == 0u) move1 = false;
  else if (this->el_ji.at(j1).size() == 0u) move1 = false;
  
  // Case 1: Both need to be moved:
  if (move0 && move1) {
    
    // Swapping elements, the pointers will work OK as the memory
    // location hasn't change, right?
    Cell c0 = this->el_ij.at(i0).at(j0);
    this->el_ij.at(i0).at(j0) = this->el_ij.at(i1).at(j1);
    this->el_ij.at(i1).at(j1) = c0;
    
  } else if (move0 && !move1) {
    
    // Adding and removing at the same time
    this->insert_cell(i1, j1, this->el_ij.at(i0).at(j0).value, false);
    this->rm_cell(i0, j0, false);
    
  } else if (!move0 && move1) {
    
    // Likewise, but the other element
    this->insert_cell(i0, j0, this->el_ij.at(i1).at(j1).value, false);
    this->rm_cell(i1, j1, false);
    
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
void Array::swap_col(uint j0, uint j1, bool check_bounds) {
  
  // if (check_bounds) {
  //   out_of_range(j0,0u);
  //   out_of_range(j1,0u);
  // }
  // 
  // bool move0=true, move1=true;
  // if (this->el_ji.at(j0).size() == 0u) move0 = false;
  // if (this->el_ji.at(j1).size() == 0u) move1 = false;
  // 
  // this->swap_cells()
  
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