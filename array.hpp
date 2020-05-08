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

class EdgeList {
public:
  uint N;
  uint M;
  std::vector< umap_int_cell >  el_ij;
  std::vector< umap_int_cell_ptr > el_ji;
  
  // Empty datum
  EdgeList (uint N_, uint M_) : N(N_), M(M_), el_ij(N_), el_ji(M_) {};
  
  // Edgelist with data
  EdgeList (
      uint N_, uint M_,
      const std::vector< uint > & source,
      const std::vector< uint > & target,
      const std::vector< double > & value,
      bool add = true
  );
  
  // Function to access the elements
  double get_cell(uint i, uint j) const;
  const umap_int_cell * get_row(uint i) const {return &(el_ij.at(i));}
  const umap_int_cell_ptr * get_col(uint i) const {return &(el_ji.at(i));}
  
  // Deletion addition operations
  void rm_cell(uint i, uint j);
  void insert_cell(uint i, uint j, double v);
  
  
};

// Edgelist with data
EdgeList::EdgeList (
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

double EdgeList::get_cell(uint i, uint j) const {
  
  if (this->el_ij.at(i).size() == 0u)
    return 0.0;
  
  // If it is not empty, then find and return
  auto search = el_ij.at(i).find(j);
  if (search != el_ij.at(i).end())
    return search->second.value;
  
  // This is if it is empty
  return 0.0;
  
}

void EdgeList::rm_cell(uint i, uint j) {
  
  // Checking the boundaries
  if (i >= this->N | j >= this->M)
    Rcpp::stop("Removing out of range.");
  
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
void EdgeList::insert_cell(uint i, uint j, double v) {
  
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

#endif