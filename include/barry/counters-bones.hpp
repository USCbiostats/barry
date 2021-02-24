#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRY_COUNTERS_BONES_HPP
#define BARRY_COUNTERS_BONES_HPP 1
 

/**
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

  /***
   * ! Initializers
   */
  Counter() : count_fun(nullptr), init_fun(nullptr) {};
  
  /**
   * @brief Creator passing a counter and an initializer
   * 
   * @param count_fun_ The main counter function.
   * @param init_fun_ The initializer function can also be used to check if the
   *  `BArray` as the needed variables (see BArray::data).
   * @param data_ Data to be used with the counter.
   * @param delete_data_ When `true`, the destructor will delete the pointer
   * in the main data.
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
    Data_Type *                            data_        = nullptr,
    bool                                   delete_data_ = false
    ): count_fun(count_fun_), init_fun(init_fun_), data(data_), delete_data(delete_data_) {};
  
  Counter(const Counter<Array_Type,Data_Type> & counter_);
  Counter<Array_Type,Data_Type> operator=(const Counter<Array_Type,Data_Type> & counter_);

  ~Counter() {
    // delete data;
    if (delete_data)
      delete data;
  };
  
  /***
   * ! Main functions.
   */
  double count(Array_Type * Array, uint i, uint j);
  double init(Array_Type * Array, uint i, uint j);
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
  
  Counter<Array_Type,Data_Type> * operator[](uint idx);
  uint size() const {
    return data.size();
    };
  
  // Functions to add counters
  void add_counter(Counter<Array_Type, Data_Type> & counter);
  void add_counter(Counter<Array_Type, Data_Type> * counter);
  void add_counter(
      Counter_fun_type<Array_Type,Data_Type> count_fun_,
      Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
      Data_Type *                            data_        = nullptr,
      bool                                   delete_data_ = false
    );
  void clear();
  
};

#endif
