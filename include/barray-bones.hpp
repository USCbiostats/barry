#include <vector>
#include <unordered_map>
#include "typedefs.hpp"

#ifndef BARRAY_BONES_HPP 
#define BARRAY_BONES_HPP 1

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

class BArray {
public:
  uint N;
  uint M;
  uint NCells;
  std::vector< umap_int_cell >  el_ij;
  std::vector< umap_int_cell_ptr > el_ji;
  
  // This is as a reference, if we need to iterate through the cells and we need
  // to keep track which were visited, we use this as a reference. So that if
  // cell.visited = true and visited = true, it means that we haven't been here
  // yet. Ideally, any routine using this->visited should switch it at the
  // beginning of the routine.
  bool visited;
  
  // Empty datum
  BArray() : N(0u), M(0u), el_ij(0u), el_ji(0u), NCells(0u) {};
  BArray (uint N_, uint M_) : N(N_), M(M_), el_ij(N_), el_ji(M_), NCells(0u) {};
  
  // Edgelist with data
  BArray (
      uint N_, uint M_,
      const std::vector< uint > & source,
      const std::vector< uint > & target,
      const std::vector< double > & value,
      bool add = true
  );
  
  // Function to access the elements
  // bool check_cell
  void out_of_range(uint i, uint j) const;
  double get_cell(uint i, uint j, bool check_bounds = true) const;
  const umap_int_cell * get_row(uint i, bool check_bounds = true) const;
  const umap_int_cell_ptr * get_col(uint i, bool check_bounds = true) const;
  bool is_empty(uint i, uint j, bool check_bounds = true) const;

  // Deletion addition operations
  void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void insert_cell(uint i, uint j, Cell v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, double v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, bool check_bounds = true, bool check_exists = true);
  
  void swap_cells(
      uint i0, uint j0, uint i1, uint j1, bool check_bounds = true,
      int check_exists = CHECK::BOTH,
      int * report     = nullptr
      );
  
  void toggle_cell(uint i, uint j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
  
  void swap_rows(uint i0, uint i1, bool check_bounds = true);
  void swap_cols(uint j0, uint j1, bool check_bounds = true);
  
  void zero_row(uint i, bool check_bounds = true);
  void zero_col(uint j, bool check_bounds = true);
  
  void transpose();
  void clear();
  void resize(uint N_, uint M_);
  
  // Advances operators
  // void toggle_iterator
  
    
};

#endif
