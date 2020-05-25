#include <vector>
#include <unordered_map>
#include "meta.hpp"
#include "typedefs.hpp"
#include "cell-bones.hpp"

#ifndef BARRAY_BONES_HPP 
#define BARRAY_BONES_HPP 1


//! Baseline class for binary arrays.
/** `BArray` class objects are arbitrary arrays
 *  in which non-empty cells hold data of type `Cell_Type`. The non-empty cells
 *  are stored by row and indexed using `unordered_map`s, i.e.
 *  `std::vector< std::unordered_map<unsigned int,Cell_Type> >`.
 */
template <typename Cell_Type> class BArray {
public:
  uint N;
  uint M;
  uint NCells = 0u;
  std::vector< Row_type< Cell_Type > > el_ij;
  std::vector< Col_type< Cell_Type > > el_ji;
  
  static Cell< Cell_Type > Cell_default;
  
  /** This is as a reference, if we need to iterate through the cells and we need
   *  to keep track which were visited, we use this as a reference. So that if
   *  cell.visited = true and visited = true, it means that we haven't been here
   *  yet. Ideally, any routine using this->visited should switch it at the
   *  beginning of the routine.
   */
  bool visited = false;
  
  //! Other information regarding the object.
  Meta meta;

  // Empty datum
  BArray() : N(0u), M(0u), el_ij(0u), el_ji(0u), NCells(0u) {};
  BArray (uint N_, uint M_) : N(N_), M(M_), el_ij(N_), el_ji(M_), NCells(0u) {};
  
  //! Edgelist with data
  /** This function takes the following arguments
   * @param N_ Number of rows
   * @param M_ Number of columns
   * @param source An unsigned vector ranging from 0 to N_
   * @param target An unsigned int vector ranging from 0 to M_
   * @param target When `true` tries to add repeated observations.
   */
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

  //! Check whether the cell is empty
  /**
   * @param i,j Coordinates
   * @param check_bounds If `false` avoids checking bounds.
   */
  bool is_empty(uint i, uint j, bool check_bounds = true) const;

  // Deletion addition operations
  void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
  void insert_cell(uint i, uint j, Cell< Cell_Type > & v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, Cell_Type v, bool check_bounds = true, bool check_exists = true);
  void insert_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
  
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
