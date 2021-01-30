// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"
#include "rules-bones.hpp"

#ifndef SUPPORT_HPP 
#define SUPPORT_HPP 1

/**@brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 */ 
template <
  typename Array_Type = BArray<>,
  typename Data_Counter_Type = bool,
  typename Data_Rule_Type = bool
  >
class Support {
  
private:
  void calc_backend(
      uint pos = 0u,
      std::vector< Array_Type > * array_bank = nullptr,
      std::vector< std::vector< double > > * stats_bank = nullptr
  );
  
public:
  
  /**@brief Reference array to generate the support.
   */
  Array_Type                            EmptyArray;
  FreqTable<>                           data;
  CounterVector<Array_Type,Data_Counter_Type> * counters;
  Rules<Array_Type,Data_Rule_Type>            * rules;

  uint N, M;
  bool initialized     = false;
  bool counter_deleted = false;
  bool rules_deleted   = false;
  
  // Temp variables to reduce memory allocation
  std::vector< double >                current_stats;
  std::vector< std::pair<uint,uint> >  coordinates_free;
  std::vector< std::pair<uint,uint> >  coordinates_locked;
  std::vector< std::vector< double > > change_stats;
  
  /**@brief Constructor passing a reference Array.
   */
  Support(const Array_Type * Array_) :
    EmptyArray(*Array_),
    counters(new CounterVector<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(Array_->nrow()), M(Array_->ncol()) {
    // init_support();
    return;
    
  };
  
  /**@brief Constructor specifying the dimensions of the array (empty).
   */
  Support(uint N_, uint M_) :
    EmptyArray(N_, M_),
    counters(new CounterVector<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(N_), M(M_) {
    // init_support();
    return;
  };
  
  Support() :
    EmptyArray(0u, 0u),
    counters(new CounterVector<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(0u), M(0u) {
    // init_support();
    return;
  };
  
  ~Support() {
    if (!counter_deleted)
      delete counters;
    if (!rules_deleted)
      delete rules;
  };
  
  void init_support();
  
  
  /**@brief Resets the support calculator
   * 
   * If needed, the counters of a support object can be reused.
   * 
   * @param Array_ New array over which the support will be computed.
   */
  void reset_array();
  void reset_array(const Array_Type * Array_);
  
  void add_counter(Counter<Array_Type, Data_Counter_Type> * f_);
  void add_counter(Counter<Array_Type,Data_Counter_Type> f_);
  void set_counters(CounterVector<Array_Type,Data_Counter_Type> * counters_);
  
  void add_rule(Rule<Array_Type, Data_Rule_Type> * f_);
  void add_rule(Rule<Array_Type,Data_Rule_Type> f_);
  void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);
  
  /**@brief Computes the entire support
   * 
   * Not to be used by the user. Sets the starting point in the array
   * (column-major).
   * 
   * @param diag When `true`, includes the diagonal (i=j) in the counts.
   * 
   * @param array_bank If specified, the counter will add to the vector each 
   * possible state of the array, as it counts.
   * 
   * @param stats_bank If specified, the counter will add to the vector each
   * possible set of statistics, as it counts.
   * 
   */
  void calc(
      std::vector< Array_Type > * array_bank = nullptr,
      std::vector< std::vector< double > > * stats_bank = nullptr
    );
  
  
  Counts_type           get_counts() const;
  const MapVec_type<> * get_counts_ptr() const;
  void print() const;
  
};

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type>::init_support() {
  
  // Computing the locations
  rules->get_seq(&EmptyArray, &coordinates_free, &coordinates_locked);
  
  return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type, Data_Counter_Type, Data_Rule_Type>::reset_array() {
  
  data.clear();
  initialized = false;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type, Data_Counter_Type, Data_Rule_Type>::reset_array(const Array_Type * Array_) {
  
  data.clear();
  initialized = false;
  EmptyArray = *Array_;
  N = Array_->nrow();
  M = Array_->ncol();
  // init_support();
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type, Data_Counter_Type, Data_Rule_Type>::calc_backend(
    uint                                   pos,
    std::vector< Array_Type > *            array_bank,
    std::vector< std::vector< double > > * stats_bank
  ) {
  
  // Did we reached the end??
  if (pos >= coordinates_free.size()) {
    return;
  }
  
  // // No self ties, go to the next step and return.
  // if (!diag && (coordinates_free[pos].first == coordinates_free[pos].second))
  //   return calc(pos + 1u, diag, array_bank, stats_bank);
  
  // Initializing
  if (!initialized) {
    
    // Do we have any counter?
    if (counters->size() == 0u)
      throw std::logic_error("No counters added: Cannot compute the support without knowning what to count!");
    
    // Initializing
    initialized = true;
    EmptyArray.clear(true);
    EmptyArray.reserve();
    current_stats.resize(counters->size());
    change_stats.resize(N*M, current_stats);
    
    // Resizing support
    data.reserve(pow(2.0, N * (M - 1.0))); 
    
    for (uint n = 0u; n < counters->size(); ++n)
      current_stats[n] = counters->operator[](n)->init(
        &EmptyArray,
        coordinates_free[pos].first, 
        coordinates_free[pos].second
      );
    
    // Adding to the overall count
    data.add(current_stats);
    
    if (array_bank != nullptr) 
      array_bank->push_back(EmptyArray);
    
    if (stats_bank != nullptr)
      stats_bank->push_back(current_stats);
    
  }
  
  // We will pass it to the next step, if the iteration makes sense.
  calc_backend(pos + 1u, array_bank, stats_bank);
  
  // Once we have returned, everything will be back as it used to be, so we
  // treat the data as if nothing has changed.
  
  // Toggle the cell (we will toggle it back after calling the counter)
  EmptyArray.insert_cell(
    coordinates_free[pos].first,
    coordinates_free[pos].second,
    EmptyArray.Cell_default , false, false
    );

  // Counting
  // std::vector< double > change_stats(counters.size());
  for (uint n = 0u; n < counters->size(); ++n) {
    change_stats[pos][n] = counters->operator[](n)->count(
      &EmptyArray,
      coordinates_free[pos].first,
      coordinates_free[pos].second
      );
    current_stats[n] += change_stats[pos][n];
  }
  
  // Adding to the overall count
  data.add(current_stats);
  
  // Need to save?
  if (array_bank != nullptr)
    array_bank->push_back(EmptyArray);
  
  if (stats_bank != nullptr)
    stats_bank->push_back(current_stats);
  
  // Again, we only pass it to the next level iff the next level is not
  // passed the last step.
  calc_backend(pos + 1, array_bank, stats_bank);
  
  // We need to restore the state of the cell
  EmptyArray.rm_cell(
    coordinates_free[pos].first,
    coordinates_free[pos].second,
    false, false
    );
  
  for (uint n = 0u; n < counters->size(); ++n) 
    current_stats[n] -= change_stats[pos][n];
  
  
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type, Data_Counter_Type, Data_Rule_Type>::calc(
    std::vector< Array_Type > *            array_bank,
    std::vector< std::vector< double > > * stats_bank
) {

  this->init_support();
  calc_backend(0u, array_bank, stats_bank);
  return;
  
}
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::add_counter(
    Counter<Array_Type, Data_Counter_Type> * f_
  ) {
  
  counters->add_counter(f_);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::add_counter(
    Counter<Array_Type,Data_Counter_Type> f_
) {
  
  counters->add_counter(f_);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::set_counters(
    CounterVector<Array_Type,Data_Counter_Type> * counters_
) {
  
  // Cleaning up before replacing the memory
  if (!counter_deleted)
    delete counters;
  counter_deleted = true;
  counters = counters_;
  
  return;
  
}

/////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::add_rule(
    Rule<Array_Type, Data_Rule_Type> * f_
) {
  
  rules->add_rule(f_);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::add_rule(
    Rule<Array_Type,Data_Rule_Type> f_
) {
  
  rules->add_rule(f_);
  return;
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::set_rules(
    Rules<Array_Type,Data_Rule_Type> * rules_
) {
  
  // Cleaning up before replacing the memory
  if (!rules_deleted)
    delete rules;
  rules_deleted = true;
  rules = rules_;
  
  return;
  
}

//////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline Counts_type Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::get_counts() const {
  
  return data.as_vector(); 
  
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline const MapVec_type<> * Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::get_counts_ptr() const {
  
  return data.get_data_ptr();
   
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type>::print() const {
  data.print();
}

#endif
