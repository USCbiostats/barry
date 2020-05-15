#include <vector>
#include <unordered_map>
#include "typedefs.hpp"
#include "cell-bones.hpp"

#ifndef BARRAY_BONES_HPP 
#define BARRAY_BONES_HPP 1

template <typename Cell_Type> class BArray {
public:
  uint N;
  uint M;
  uint NCells = 0u;
  std::vector< Row_type< Cell_Type > > el_ij;
  std::vector< Col_type< Cell_Type > > el_ji;
  
  // This is as a reference, if we need to iterate through the cells and we need
  // to keep track which were visited, we use this as a reference. So that if
  // cell.visited = true and visited = true, it means that we haven't been here
  // yet. Ideally, any routine using this->visited should switch it at the
  // beginning of the routine.
  bool visited = false;
  
  /***
   * ! Other information regarding the object
   */
  bool symetric;
  bool valued;

  // Empty datum
  BArray() : N(0u), M(0u), el_ij(0u), el_ji(0u), NCells(0u) {};
  BArray (uint N_, uint M_) : N(N_), M(M_), el_ij(N_), el_ji(M_), NCells(0u) {};
  
  // Edgelist with data
  BArray (
      uint N_, uint M_,
      const std::vector< uint > & source,
      const std::vector< uint > & target,
      const std::vector< Cell_Type > & value,
      bool add = true
  );
   
  // Edgelist with no data (simpler)
  BArray (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    bool add = true
  );
  
  // Function to access the elements
  // bool check_cell
  void out_of_range(uint i, uint j) const;
  Cell_Type get_cell(uint i, uint j, bool check_bounds = true) const; 
  const Row_type< Cell_Type > * get_row(uint i, bool check_bounds = true) const;
  const Col_type< Cell_Type > * get_col(uint i, bool check_bounds = true) const;
  Entries<Cell_Type> get_entries() const;
  
  
  // Queries
  bool is_empty(uint i, uint j, bool check_bounds = true) const;

  // Deletion addition operations
  void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void insert_cell(uint i, uint j, std::pair< Cell_Type, bool> v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, Cell< Cell_Type > & v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, Cell_Type v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void insert_cell(uint i, std::pair< Cell_Type, bool> v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, Cell< Cell_Type > & v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, Cell_Type v, bool check_bounds = true, bool check_exists = true);
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
