// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"

#ifndef SUPPORT_HPP 
#define SUPPORT_HPP 1

/**@brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 */ 
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Support {
public:
  
  /**@brief Reference array to generate the support.
   */
  Array_Type                                       EmptyArray;
  FreqTable<>                                      data;
  std::vector< Counter<Array_Type, Data_Type>* > * counters;
  std::vector< uint >                              to_be_deleted;

  uint N, M;
  bool initialized = false;
  bool counter_deleted = false;
  
  // Temp variables to reduce memory allocation
  std::vector< double >                current_stats;
  std::vector< uint >                  pos_i;
  std::vector< uint >                  pos_j;
  std::vector< std::vector< double > > change_stats;
  
  /**@brief Constructor passing a reference Array.
   */
  Support(const Array_Type * Array_) :
    EmptyArray(*Array_), N(Array_->N), M(Array_->M) {
    
    counters = new std::vector< Counter<Array_Type, Data_Type>* >(0u);
    init_support();
    return;
    
  };
  
  /**@brief Constructor specifying the dimensions of the array (empty).
   */
  Support(uint N_, uint M_) :
    EmptyArray(N_, M_) ,N(N_), M(M_) {
    
    counters = new std::vector< Counter<Array_Type, Data_Type>* >(0u);
    init_support();
    return;
    
  };
  
  ~Support() {
    for (auto iter = to_be_deleted.begin(); iter != to_be_deleted.end(); ++iter)
      delete counters->operator[](*iter);
    
    if (!counter_deleted)
      delete counters;
  };
  
  void init_support();
  
  
  /**@brief Resets the support calculator
   * 
   * If needed, the counters of a support object can be reused.
   * 
   * @param Array_ New array over which the support will be computed.
   */
  void reset();
  void reset(const Array_Type * Array_);
  void add_counter(Counter<Array_Type, Data_Type> * f_);
  void add_counter(Counter<Array_Type,Data_Type> f_);
  void set_counters(std::vector< Counter<Array_Type,Data_Type> *> * counters_);
  
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
      uint pos = 0u,
      const bool & diag = false,
      std::vector< Array_Type > * array_bank = nullptr,
      std::vector< std::vector< double > > * stats_bank = nullptr
    );
  
  Counts_type           get_counts() const;
  const MapVec_type<> * get_counts_ptr() const;
  
};

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::init_support() {
  
  pos_i.resize(N*M);
  pos_j.resize(N*M);
  for (uint pos = 0u; pos < pos_i.size(); ++pos) {
    pos_i[pos] = (int) pos % (int) N;
    pos_j[pos] = floor((int) pos / (int) N);
  }
  
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::reset() {
  
  data.clear();
  initialized = false;
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::reset(const Array_Type * Array_) {
  
  data.clear();
  initialized = false;
  EmptyArray = *Array_;
  N = Array_->N;
  M = Array_->M;
  init_support();
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type, Data_Type>::calc(
    uint pos, const bool & diag,
    std::vector< Array_Type > * array_bank,
    std::vector< std::vector< double > > * stats_bank
  ) {
  
  if (pos >= N*M) {
    return;
  }
  
  // No self ties, go to the next step and return.
  if (!diag && (pos_i[pos] == pos_j[pos]))
    return calc(pos + 1u, diag, array_bank, stats_bank);
  
  // Initializing
  if (!initialized) {
    
    // Initializing
    initialized = true;
    EmptyArray.clear(true);
    EmptyArray.reserve();
    current_stats.resize(counters->size());
    change_stats.resize(N*M, current_stats);
    
    // Resizing support
    data.reserve(pow(2.0, N * (M - 1.0)));
    
    for (uint n = 0u; n < counters->size(); ++n) 
      current_stats[n] = counters->operator[](n)->init(&EmptyArray, pos_i[pos], pos_j[pos]);
    
    // Adding to the overall count
    data.add(current_stats);
    
    if (array_bank != nullptr)
      array_bank->push_back(EmptyArray);
    
    if (stats_bank != nullptr)
      stats_bank->push_back(current_stats);
    
  }
  
  // We will pass it to the next step, if the iteration makes sense.
  calc(pos + 1u, diag, array_bank, stats_bank);
  
  // Once we have returned, everything will be back as it used to be, so we
  // treat the data as if nothing has changed.
  
  // Toggle the cell (we will toggle it back after calling the counter)
  EmptyArray.insert_cell(pos_i[pos], pos_j[pos], true, false, false);
  
  // Counting
  // std::vector< double > change_stats(counters.size());
  for (uint n = 0u; n < counters->size(); ++n) {
    change_stats[pos][n] = counters->operator[](n)->count(&EmptyArray, pos_i[pos], pos_j[pos]);
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
  calc(pos + 1, diag, array_bank, stats_bank);
  
  // We need to restore the state of the cell
  EmptyArray.rm_cell(pos_i[pos], pos_j[pos], false, false);
  for (uint n = 0u; n < counters->size(); ++n) 
    current_stats[n] -= change_stats[pos][n];
  
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * f_
  ) {
  counters->push_back(f_);
  return;
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type,Data_Type> f_
) {
  
  to_be_deleted.push_back(counters->size());
  counters->push_back(new Counter<Array_Type,Data_Type>(f_));
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline void Support<Array_Type,Data_Type>::set_counters(
    std::vector< Counter<Array_Type,Data_Type> *> * counters_
) {
  
  for (auto iter = to_be_deleted.begin(); iter != to_be_deleted.end(); ++iter)
    delete counters->operator[](*iter);
  
  if (!counter_deleted)
    delete counters;
  
  counters = counters;
  
  return;
  
}

template <typename Array_Type, typename Data_Type>
inline Counts_type Support<Array_Type,Data_Type>::get_counts() const {
  
  return data.as_vector(); 
  
}

template <typename Array_Type, typename Data_Type>
inline const MapVec_type<> * Support<Array_Type,Data_Type>::get_counts_ptr() const {
  
  return data.get_data_ptr();
   
}


#endif
