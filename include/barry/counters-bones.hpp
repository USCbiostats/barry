#include "typedefs.hpp"
#include "barray-bones.hpp"

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
    Data_Type * data = nullptr;
    bool delete_data = false;
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
    Counter() : count_fun(nullptr), init_fun(nullptr) {};
    
    Counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
        Data_Type *                            data_        = nullptr,
        bool                                   delete_data_ = false,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
        ): count_fun(count_fun_), init_fun(init_fun_), data(data_),
            delete_data(delete_data_), name(name_), desc(desc_) {};
    
    Counter(const Counter<Array_Type,Data_Type> & counter_);
    ///@}

    Counter<Array_Type,Data_Type> operator=(const Counter<Array_Type,Data_Type> & counter_);

    ~Counter() {
        // delete data;
        if (delete_data)
            delete data;
    };
    
    /***
      * ! Main functions.
      */
    double count(Array_Type & Array, uint i, uint j);
    double init(Array_Type & Array, uint i, uint j);
    
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
    std::vector< Counter<Array_Type,Data_Type >* > data = {};
    std::vector< uint >                            to_be_deleted = {};
    
public: 
    
    // Constructors
    Counters() {};
    
    // Destructor needs to deal with the pointers
    ~Counters() {
        this->clear();
    }

    Counters(const Counters<Array_Type,Data_Type> & counter_);
    Counters<Array_Type,Data_Type> operator=(const Counters<Array_Type,Data_Type> & counter_);
    
    /**
     * @brief Returns a pointer to a particular counter.
     * 
     * @param idx Id of the counter
     * @return Counter<Array_Type,Data_Type>* 
     */
    Counter<Array_Type,Data_Type> & operator[](uint idx);

    /**
     * @brief Number of counters in the set.
     * 
     * @return uint 
     */
    std::size_t size() const noexcept {
        return data.size();
        };
    
    // Functions to add counters
    void add_counter(Counter<Array_Type, Data_Type> & counter);
    void add_counter(Counter<Array_Type, Data_Type> * counter);
    void add_counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
        Data_Type *                            data_        = nullptr,
        bool                                   delete_data_ = false,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
    );
    void clear();
    
};

#endif
