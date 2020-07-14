#include "typedefs.hpp"
#include "barray-bones.hpp"

#ifndef BARRAY_COUNTERS_H
#define BARRAY_COUNTERS_H

  
#define COUNTER_FUNCTION(a) template <typename Array_Type, typename Data_Type> \
  inline double (a) (Array_Type * Array, uint i, uint j, Data_Type * data)\
    

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
   * @param count_fun The main counter function.
   * @param init_fun The initializer function can also be used to check if the
   *  `BArray` as the needed variables (see BArray::data).
   * @param delete_data_ When `true`, the destructor will delete the pointer
   * in the main data.
   */
  Counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_    = nullptr,
    Data_Type *                            data_        = nullptr,
    bool                                   delete_data_ = false
    ): count_fun(count_fun_), init_fun(init_fun_), data(data_), delete_data(delete_data_) {};
  
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

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::count(
    Array_Type * Array, uint i, uint j) {
  if (count_fun == nullptr)
    return 0.0;
  return count_fun(Array, i, j, data);
}

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::init(
    Array_Type * Array, uint i, uint j) {
  if (init_fun == nullptr)
    return 0.0;
  return init_fun(Array, i, j, data);
}

/**@brief Vector of counters.
 * 
 * Various functions hold more than one counter, so this class is a helper class
 * that allows managing multiple counters efficiently. The main data is a vector
 * to pointers of counters.
 */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class CounterVector {
  
private:
  std::vector< Counter<Array_Type,Data_Type >* > data = {};
  std::vector< uint >                            to_be_deleted = {};
  
public: 
  
  // Constructors
  CounterVector() {};
  
  // Destructor needs to deal with the pointers
  ~CounterVector() {
    this->clear();
  }
  
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

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type> * CounterVector<Array_Type,Data_Type>::operator[](uint idx) {
  return data[idx];
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::add_counter(
  Counter<Array_Type, Data_Type> & counter
) {
  
  to_be_deleted.push_back(data.size());
  data.push_back(new Counter<Array_Type, Data_Type>(counter));
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * counter
) {
  
  data.push_back(counter);
  return;
  
}

template <typename Array_Type, typename Data_Type>
void CounterVector<Array_Type,Data_Type>::add_counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Data_Type *                            data_,
    bool                                   delete_data_
) {
 
  to_be_deleted.push_back(data.size());
  data.push_back(new Counter<Array_Type,Data_Type>(
    count_fun_,
    init_fun_,
    data_,
    delete_data_
  ));
 
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void CounterVector<Array_Type,Data_Type>::clear() {
  
  for (auto iter = to_be_deleted.begin(); iter != to_be_deleted.end(); ++iter)
    delete data[*iter];
  
  to_be_deleted.clear();
  
  return;
  
}

#endif
