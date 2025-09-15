#ifndef BARRY_COUNTERS_BONES_HPP
#define BARRY_COUNTERS_BONES_HPP 1

/**
 * @defgroup counting
 * @details `barry` includes a flexible way to generate counters based on change
 * statistics. Since most of the time we are counting many motifs in a graph,
 * change statistics make a reasonable (and efficient) way to make such counts.
 * 
 * In particular, let the motif be defined as \f$s(y)\f$, with \f$y\f$ as the
 * binary array. The change statistic when adding cell \f$y_{ij}\f$, i.e. when
 * the cell moves from being emty to have a one, is defined as
 * 
 * \f[
 * \delta(y_{ij}) = s^+_{ij}(y) - s^-_{ij}(y),
 * \f]
 * 
 * where \f$s^+_{ij}(y)\f$ and \f$s^-_{ij}(y)\f$ represent the motif statistic
 * with and without the ij-cell. For example, in the case of networks, the change
 * statistic for the number of edges is always 1. 
 * 
 * To count statistics in an array, the [Counter] class will empty the array, 
 * initialize the counters, and then start counting while adding at each step
 * a single cell, until matching the original array. 
 */

/**
  * @ingroup counting Implementation of motif counting
  * @brief A counter function based on change statistics.
  * 
  * This class is used by `CountStats` and `StatsCounter` as a way to count
  * statistics using change statistics.
  */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Counter {
public:
    
    Counter_fun_type<Array_Type,Data_Type> count_fun;
    Counter_fun_type<Array_Type,Data_Type> init_fun;
    Hasher_fun_type<Array_Type,Data_Type> hasher_fun;

    Data_Type data;
    std::string  name = "";
    std::string  desc = "";

    /**
     * @name Creator passing a counter and an initializer
     * 
     * @param count_fun_ The main counter function.
     * @param init_fun_ The initializer function can also be used to check if the
     *  `BArray` as the needed variables (see BArray::data).
     * @param data_ Data to be used with the counter.
     * @param delete_data_ When `true`, the destructor will delete the pointer
     * in the main data.
     */
    ///@{
    Counter() : count_fun(nullptr), init_fun(nullptr), hasher_fun(nullptr) {};
    
    Counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_,
        Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
        Data_Type                              data_,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
        ): count_fun(count_fun_), init_fun(init_fun_), hasher_fun(hasher_fun_), data(data_),
            name(name_), desc(desc_) {};
    
    Counter(const Counter<Array_Type,Data_Type> & counter_); ///< Copy constructor
    Counter(Counter<Array_Type,Data_Type> && counter_) noexcept; ///< Move constructor
    Counter<Array_Type,Data_Type> operator=(const Counter<Array_Type,Data_Type> & counter_); ///< Copy assignment
    Counter<Array_Type,Data_Type>& operator=(Counter<Array_Type,Data_Type> && counter_) noexcept; ///< Move assignment
    ///@}

    ~Counter() {};
    
    /***
      * ! Main functions.
      */
    double count(Array_Type & Array, size_t i, size_t j);
    double init(Array_Type & Array, size_t i, size_t j);
    std::string get_name() const;
    std::string get_description() const;
    void set_name(std::string new_name);
    void set_description(std::string new_desc);

    /**
     * @brief Get and set the hasher function
     * 
     * The hasher function is used to characterize the support of the array.
     * This way, if possible, the support enumeration is recycled.
     * 
     * @param fun 
     */
    ///@{
    void set_hasher(Hasher_fun_type<Array_Type,Data_Type> fun);
    Hasher_fun_type<Array_Type,Data_Type> get_hasher();
    ///@}
    
};

/**
  * @brief Vector of counters.
  * 
  * Various functions hold more than one counter, so this class is a helper class
  * that allows managing multiple counters efficiently. The main data is a vector
  * to pointers of counters.
  */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Counters {
    
private:
    std::vector< Counter<Array_Type,Data_Type > > data;
    Hasher_fun_type<Array_Type,Data_Type> hasher;
    
public: 
    
    // Constructors
    Counters();
    
    // Destructor needs to deal with the pointers
    ~Counters() {};

    /**
     * @brief Copy constructor
     * @param counter_ 
     */
    Counters(const Counters<Array_Type,Data_Type> & counter_);
    
    /**
     * @brief Move constructor
     * 
     * @param counters_ 
     */
    Counters(Counters<Array_Type,Data_Type> && counters_) noexcept;

    /**
     * @brief Copy assignment constructor
     * 
     * @param counter_ 
     * @return Counters<Array_Type,Data_Type> 
     */
    Counters<Array_Type,Data_Type> operator=(const Counters<Array_Type,Data_Type> & counter_);

    /**
     * @brief Move assignment constructor
     * 
     * @param counter_ 
     * @return Counters<Array_Type,Data_Type>& 
     */
    Counters<Array_Type,Data_Type> & operator=(Counters<Array_Type,Data_Type> && counter_) noexcept;
    
    /**
     * @brief Returns a pointer to a particular counter.
     * 
     * @param idx Id of the counter
     * @return Counter<Array_Type,Data_Type>* 
     */
    Counter<Array_Type,Data_Type> & operator[](size_t idx);

    /**
     * @brief Number of counters in the set.
     * 
     * @return size_t 
     */
    std::size_t size() const noexcept {
        return data.size();
        };
    
    // Functions to add counters
    void add_counter(Counter<Array_Type, Data_Type> counter);
    void add_counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_,
        Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
        Data_Type                              data_,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
    );
    
    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    /**
     * @brief Generates a hash for the given array according to the counters.
     * 
     * @param array 
     * @param add_dims When `true` (default) the dimmension of the array will
     * be added to the hash.
     * @return std::vector< double > That can be hashed later.
     */
    std::vector< double > gen_hash(
      const Array_Type & array,
      bool add_dims = true
      );

    void add_hash(
      Hasher_fun_type<Array_Type,Data_Type> fun_
    );
    
};

#endif

