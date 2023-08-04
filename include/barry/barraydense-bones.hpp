#ifndef BARRY_BARRAYDENSE_BONES_HPP 
#define BARRY_BARRAYDENSE_BONES_HPP 1

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell_const;

/**
 * @brief Baseline class for binary arrays.
 * 
 * `BArrayDense` class objects are arbitrary dense-arrays. The data
 * is stored internally in the `el` member, which can be accessed
 * using the member function `get_data()`, by column.
 *
 * @tparam Cell_Type Type of cell (any type).
 * @tparam Data_Type Data type of the array (bool default).
 */
template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDense {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCol<Cell_Type,Data_Type>;
    friend class BArrayDenseCol_const<Cell_Type,Data_Type>;
    friend class BArrayDenseRow<Cell_Type,Data_Type>;
    friend class BArrayDenseRow_const<Cell_Type,Data_Type>;
    // friend class Support<Cell_Type,Data_Type>;
    // friend class StatsCounter<Cell_Type,Data_Type>;
private:
    size_t N;
    size_t M;
    // size_t NCells = 0u;
    std::vector< Cell_Type > el;
    std::vector< Cell_Type > el_rowsums;
    std::vector< Cell_Type > el_colsums;
    Data_Type * data = nullptr;
    bool delete_data = false;

    static Cell_Type Cell_default;
    static const bool dense = true;

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
     * @param target An size_t vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     * @param value Cell_Type defaul fill-in value (zero, by default.)
     */
    ///@{
    
    /** @brief Zero-size array */
    BArrayDense() : N(0u), M(0u), el(0u), el_rowsums(0u), el_colsums(0u) {};
    
    /** @brief Empty array */
    BArrayDense (size_t N_, size_t M_, Cell_Type value = static_cast<Cell_Type>(0)) :
        N(N_), M(M_), el(N_ * M_, value),
        el_rowsums(N_, static_cast<Cell_Type>(value * M_)), el_colsums(M_, static_cast<Cell_Type>(value * N_)) {};
    
    /** @brief Edgelist with data */
    BArrayDense (
        size_t N_,
        size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArrayDense (
        size_t N_, size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
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
    Data_Type * D_ptr();
    const Data_Type * D_ptr() const;
    Data_Type & D();
    const Data_Type & D() const;
    ///@}
    
    // Function to access the elements
    // bool check_cell
    void out_of_range(size_t i, size_t j) const;
    Cell_Type get_cell(size_t i, size_t j, bool check_bounds = true) const; 
    std::vector< Cell_Type >      get_col_vec(size_t i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(size_t i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    
    BArrayDenseRow<Cell_Type,Data_Type> & row(size_t i, bool check_bounds = true);
    const BArrayDenseRow_const<Cell_Type,Data_Type> row(size_t i, bool check_bounds = true) const;

    BArrayDenseCol<Cell_Type,Data_Type> & col(size_t j, bool check_bounds = true);
    const BArrayDenseCol_const<Cell_Type,Data_Type> col(size_t j, bool check_bounds = true) const;

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
    bool is_empty(size_t i, size_t j, bool check_bounds = true) const;
    size_t nrow() const noexcept;
    size_t ncol() const noexcept;
    size_t nnozero() const noexcept;
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
    BArrayDense<Cell_Type,Data_Type> & operator+=(const std::pair<size_t, size_t> & coords);
    BArrayDense<Cell_Type,Data_Type> & operator-=(const std::pair<size_t, size_t> & coords);
    BArrayDenseCell<Cell_Type,Data_Type> operator()(size_t i, size_t j, bool check_bounds = true);
    const Cell_Type operator()(size_t i, size_t j, bool check_bounds = true) const;
    
    void rm_cell(size_t i, size_t j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(size_t i, size_t j, const Cell< Cell_Type > & v, bool check_bounds, bool check_exists);
    // void insert_cell(size_t i, size_t j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(size_t i, size_t j, Cell_Type v, bool check_bounds, bool check_exists);
    
    void swap_cells(
        size_t i0, size_t j0, size_t i1, size_t j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(size_t i, size_t j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(size_t i, size_t j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(size_t i0, size_t i1, bool check_bounds = true);
    void swap_cols(size_t j0, size_t j1, bool check_bounds = true);
    
    void zero_row(size_t i, bool check_bounds = true);
    void zero_col(size_t j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(size_t N_, size_t M_);
    void reserve();

    // Advances operators
    // void toggle_iterator
    
    // Misc
    void print(const char * fmt = nullptr, ...) const;

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
    // operator BArrayDense<size_t,bool>() const;
    // operator BArrayDense<bool,bool>() const;
    // ///@}
    
    bool is_dense() const noexcept {return dense;};

    const std::vector< Cell_Type > & get_data() const;
    const Cell_Type rowsum(size_t i) const;
    const Cell_Type colsum(size_t i) const;
};

#endif
