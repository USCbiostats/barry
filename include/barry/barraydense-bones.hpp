// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "cell-bones.hpp"
#include "barraycell-bones.hpp"

#ifndef BARRY_BARRAYDENSE_BONES_HPP 
#define BARRY_BARRAYDENSE_BONES_HPP 1

/**
 * @brief Baseline class for binary arrays.
 * 
 * `BArrayDense` class objects are arbitrary arrays
 * in which non-empty cells hold data of type `Cell_Type`. The non-empty cells
 * are stored by row and indexed using `unordered_map`s, i.e.
 * `std::vector< std::unordered_map<unsigned int,Cell_Type> >`.
 *
 * @tparam Cell_Type Type of cell (any type).
 * @tparam Data_Type Data type of the array (bool default).
 */
template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDense {
    friend class BArrayCell<Cell_Type,Data_Type>;
    friend class BArrayCell_const<Cell_Type,Data_Type>;
private:
    uint N;
    uint M;
    uint NCells = 0u;
    std::vector< Cell< Cell_Type > > el;
    Data_Type * data = nullptr;
    bool delete_data = false;

    static Cell< Cell_Type > Cell_default;

public:
    
    /** 
     * This is as a reference, if we need to iterate through the cells and we need
     * to keep track which were visited, we use this as a reference. So that if
     * cell.visited = true and visited = true, it means that we haven't been here
     * yet. Ideally, any routine using this->visited should switch it at the
     * beginning of the routine.
     */
    bool visited = false;
    

    /**
     * @name Constructors
     * 
     * @param N_ Number of rows
     * @param M_ Number of columns
     * @param source An unsigned vector ranging from 0 to N_
     * @param target An unsigned int vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     */
    ///@{
    
    /** @brief Zero-size array */
    BArrayDense() : N(0u), M(0u), NCells(0u), el(0u) {};
    
    /** @brief Empty array */
    BArrayDense (uint N_, uint M_) : N(N_), M(M_), NCells(0u), el(N_ * M_) {};
    
    /** @brief Edgelist with data */
    BArrayDense (
        uint N_, uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArrayDense (
        uint N_, uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
        bool add = true
    );
    
    /** @brief Copy constructor */
    BArrayDense(const BArrayDense<Cell_Type,Data_Type> & Array_, bool copy_data = false);
    
    /** @brief Assignment constructor */
    BArrayDense<Cell_Type,Data_Type> & operator=(const BArrayDense<Cell_Type,Data_Type> & Array_);

    /** @brief Move operator */
    BArrayDense(BArrayDense<Cell_Type,Data_Type> && x) noexcept;

    /** @brief Move assignment */
    BArrayDense<Cell_Type,Data_Type> & operator=(BArrayDense<Cell_Type,Data_Type> && x) noexcept;
    ///@}
    
    bool operator==(const BArrayDense<Cell_Type,Data_Type> & Array_);

    ~BArrayDense();
    
    // In principle, copy can be faster by using openmp on the rows
    // since those are independent.
    // BArrayDense(BArrayDense & A);
    
    /**
     * @brief Set the data object
     * 
     * @param data_ 
     * @param delete_data_ 
     */
    ///@{
    void set_data(Data_Type * data_, bool delete_data_ = false);
    Data_Type * D();
    const Data_Type * D() const;
    ///@}
    
    // Function to access the elements
    // bool check_cell
    void out_of_range(uint i, uint j) const;
    Cell_Type get_cell(uint i, uint j, bool check_bounds = true) const; 
    const Col_type< Cell_Type > * get_col(uint i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_col_vec(uint i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(uint i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    const Row_type< Cell_Type > & row(uint i, bool check_bounds = true) const;
    const Col_type< Cell_Type > & col(uint i, bool check_bounds = true) const;

    /**
     * @brief Get the edgelist
     * 
     * `Entries` is a class with three objects: Two `std::vector` with the row and
     * column coordinates respectively, and one `std::vector` with the corresponding
     * value of the cell.
     * 
     * @return Entries<Cell_Type> 
     */
    Entries<Cell_Type> get_entries() const;

    /**
     * @name Queries
     * @details `is_empty` queries a single cell. `nrow`, `ncol`, and `nnozero`
     * return the number of rows, columns, and non-zero cells respectively.
     * @param i,j Coordinates
     * @param check_bounds If `false` avoids checking bounds.
     */
    ///@{
    bool is_empty(uint i, uint j, bool check_bounds = true) const;
    uint nrow() const noexcept;
    uint ncol() const noexcept;
    uint nnozero() const noexcept;
    Cell<Cell_Type> default_val() const;
    ///@}

    /**
     * @name Cell-wise insertion/deletion
     * @param i,j Row,column
     * @param check_bounds When `true` and out of range, the function throws an
     * error.
     * @param check_exists Wither check if the cell exists (before trying to
     * delete/add), or, in the case of `swap_cells`, check if either of both
     * cells exists/don't exist.
     */
    ///@{  
    BArrayDense<Cell_Type,Data_Type> & operator+=(const std::pair<uint, uint> & coords);
    BArrayDense<Cell_Type,Data_Type> & operator-=(const std::pair<uint, uint> & coords);
    BArrayCell<Cell_Type,Data_Type> operator()(uint i, uint j, bool check_bounds = true);
    const BArrayCell_const<Cell_Type,Data_Type> operator()(uint i, uint j, bool check_bounds = true) const;
    
    void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(uint i, uint j, const Cell< Cell_Type > & v, bool check_bounds, bool check_exists);
    void insert_cell(uint i, uint j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(uint i, uint j, Cell_Type v, bool check_bounds, bool check_exists);
    
    void swap_cells(
        uint i0, uint j0, uint i1, uint j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(uint i, uint j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(uint i, uint j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(uint i0, uint i1, bool check_bounds = true);
    void swap_cols(uint j0, uint j1, bool check_bounds = true);
    
    void zero_row(uint i, bool check_bounds = true);
    void zero_col(uint j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(uint N_, uint M_);
    void reserve();

    // Advances operators
    // void toggle_iterator
    
    // Misc
    void print() const;

    /**
     * @name Arithmetic operators
     * 
     */
    ///@{
    BArrayDense<Cell_Type,Data_Type>& operator+=(const BArrayDense<Cell_Type,Data_Type>& rhs);
    BArrayDense<Cell_Type,Data_Type>& operator+=(const Cell_Type & rhs);

    BArrayDense<Cell_Type,Data_Type>& operator-=(const BArrayDense<Cell_Type,Data_Type>& rhs);
    BArrayDense<Cell_Type,Data_Type>& operator-=(const Cell_Type & rhs);
    
    BArrayDense<Cell_Type,Data_Type>& operator/=(const Cell_Type & rhs);
    BArrayDense<Cell_Type,Data_Type>& operator*=(const Cell_Type & rhs);
    ///@}
    
    // /**
    //  * @name Casting between types
    //  */
    // ///@{
    // operator BArrayDense<double,bool>() const;
    // operator BArrayDense<int,bool>() const;
    // operator BArrayDense<uint,bool>() const;
    // operator BArrayDense<bool,bool>() const;
    // ///@}

};

#endif
